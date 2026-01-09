#!/usr/bin/env python3

import argparse
import csv
import json
import math
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np


C_KMS = 299792.458
G_MPC_KMS2_PER_MSUN = 4.30091e-9  # (km/s)^2 * Mpc / Msun


def _print(msg: str) -> None:
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}")


@dataclass
class Profile:
    r_mpc: np.ndarray
    y: np.ndarray
    yerr: Optional[np.ndarray] = None


def _read_csv_profile(path: Path, r_col: str, y_col: str, yerr_col: str = "") -> Profile:
    with open(path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        rows = list(reader)

    if not rows:
        raise ValueError(f"Empty CSV: {path}")

    def _get(name: str, row: Dict[str, str]) -> float:
        v = row.get(name, "")
        if v is None:
            return float("nan")
        try:
            return float(v)
        except Exception:
            return float("nan")

    r_list: List[float] = []
    y_list: List[float] = []
    e_list: List[float] = []

    for row in rows:
        r = _get(r_col, row)
        y = _get(y_col, row)
        if not (math.isfinite(r) and math.isfinite(y)):
            continue
        r_list.append(r)
        y_list.append(y)
        if yerr_col:
            e = _get(yerr_col, row)
            e_list.append(e)

    if len(r_list) < 3:
        raise ValueError(f"Too few valid rows in {path}")

    r = np.asarray(r_list, dtype=float)
    y = np.asarray(y_list, dtype=float)
    order = np.argsort(r)
    r = r[order]
    y = y[order]

    yerr = None
    if yerr_col:
        e = np.asarray(e_list, dtype=float)
        e = e[order]
        yerr = e

    return Profile(r_mpc=r, y=y, yerr=yerr)


def _potential_from_enclosed_mass(r_mpc: np.ndarray, m_msun: np.ndarray) -> np.ndarray:
    r = np.asarray(r_mpc, dtype=float)
    m = np.asarray(m_msun, dtype=float)
    if r.size < 3:
        raise ValueError("need >=3 points")
    if not np.all(np.diff(r) > 0):
        raise ValueError("r must be strictly increasing")

    # Integrate inward with boundary condition Phi(r_max)=0
    phi = np.zeros_like(r)
    for i in range(r.size - 2, -1, -1):
        r_hi = r[i + 1]
        r_lo = r[i]
        dr = r_hi - r_lo

        # Use midpoint radius/mass for the integrand GM/r^2
        r_mid = 0.5 * (r_hi + r_lo)
        m_mid = np.interp(r_mid, r, m)

        dphi = -G_MPC_KMS2_PER_MSUN * m_mid * dr / (r_mid * r_mid)
        phi[i] = phi[i + 1] + dphi

    return phi


def _dv_from_phi(phi: np.ndarray) -> np.ndarray:
    # Convention: Δv is the relative shift of galaxies at radius r with respect to
    # the central galaxy (BCG). Since Phi(0) is deeper (more negative), galaxies at
    # larger r appear blueshifted relative to the BCG, giving a negative Δv.
    return -(phi - float(phi[0])) / C_KMS


def _weighted_fit_scale_offset(x: np.ndarray, y: np.ndarray, yerr: np.ndarray) -> Tuple[float, float, float, int]:
    good = np.isfinite(x) & np.isfinite(y) & np.isfinite(yerr) & (yerr > 0)
    if np.count_nonzero(good) < 3:
        return float("nan"), float("nan"), float("nan"), 0

    xx = x[good]
    yy = y[good]
    ww = 1.0 / (yerr[good] ** 2)

    X = np.column_stack([xx, np.ones_like(xx)])
    W = np.diag(ww)

    try:
        beta = np.linalg.solve(X.T @ W @ X, X.T @ W @ yy)
    except np.linalg.LinAlgError:
        return float("nan"), float("nan"), float("nan"), int(np.count_nonzero(good))

    s = float(beta[0])
    o = float(beta[1])

    yhat = s * xx + o
    chi2 = float(np.sum(ww * (yy - yhat) ** 2))
    dof = int(xx.size - 2)

    return s, o, chi2, dof


def _weighted_fit_offset(y_model: np.ndarray, y_obs: np.ndarray, yerr: np.ndarray) -> Tuple[float, float, int]:
    good = np.isfinite(y_model) & np.isfinite(y_obs) & np.isfinite(yerr) & (yerr > 0)
    if np.count_nonzero(good) < 2:
        return float("nan"), float("nan"), 0
    ww = 1.0 / (yerr[good] ** 2)
    offset = float(np.sum(ww * (y_obs[good] - y_model[good])) / np.sum(ww))
    chi2 = float(np.sum(ww * (y_obs[good] - (y_model[good] + offset)) ** 2))
    dof = int(np.count_nonzero(good) - 1)
    return offset, chi2, dof


def _chi2(y_model: np.ndarray, y_obs: np.ndarray, yerr: np.ndarray) -> Tuple[float, int]:
    good = np.isfinite(y_model) & np.isfinite(y_obs) & np.isfinite(yerr) & (yerr > 0)
    if np.count_nonzero(good) < 1:
        return float("nan"), 0
    ww = 1.0 / (yerr[good] ** 2)
    chi2 = float(np.sum(ww * (y_obs[good] - y_model[good]) ** 2))
    return chi2, int(np.count_nonzero(good))


def _interp_to_target(r_src: np.ndarray, y_src: np.ndarray, r_tgt: np.ndarray) -> np.ndarray:
    return np.interp(r_tgt, r_src, y_src)


def _nfw_enclosed_mass_msun(r_mpc: np.ndarray, m500_msun: float, r500_mpc: float, c500: float) -> np.ndarray:
    r = np.asarray(r_mpc, dtype=float)
    x = r / float(r500_mpc)
    c = float(c500)
    if not (math.isfinite(c) and c > 0):
        raise ValueError("c500 must be > 0")

    def f(u: np.ndarray) -> np.ndarray:
        return np.log(1.0 + u) - (u / (1.0 + u))

    norm = f(np.array([c]))[0]
    if not (math.isfinite(norm) and norm > 0):
        raise ValueError("invalid NFW normalization")

    m = float(m500_msun) * f(c * x) / norm
    m = np.where(np.isfinite(m), m, np.nan)
    return m


def main() -> None:
    p = argparse.ArgumentParser(description="Step 3.2: Cluster three-leg closure (lensing + dynamics + gravitational redshift)")

    p.add_argument("--lensing-mass-csv", default="")
    p.add_argument("--dynamics-mass-csv", default="")
    p.add_argument("--grz-csv", default="")

    p.add_argument("--mass-r-col", default="r_mpc")
    p.add_argument("--mass-m-col", default="m_enclosed_msun")

    p.add_argument("--grz-r-col", default="r_mpc")
    p.add_argument("--grz-dv-col", default="dv_kms")
    p.add_argument("--grz-dv-err-col", default="dv_err_kms")

    p.add_argument("--nfw-m500-msun", type=float, default=0.0)
    p.add_argument("--nfw-r500-mpc", type=float, default=0.0)
    p.add_argument("--nfw-c500", type=float, default=3.0)
    p.add_argument("--nfw-rmax-r500", type=float, default=4.0)
    p.add_argument("--nfw-npts", type=int, default=400)

    p.add_argument("--fit-mode", choices=["fixed", "offset_only", "scale_offset"], default="fixed")
    p.add_argument("--corr-csv", default="")
    p.add_argument("--corr-r-col", default="r_mpc")
    p.add_argument("--corr-dv-col", default="dv_kms")

    p.add_argument("--output", default="results/outputs/step_3_2_cluster_closure.json")

    args = p.parse_args()

    if not args.grz_csv:
        raise ValueError("--grz-csv is required")

    grz = _read_csv_profile(Path(args.grz_csv), args.grz_r_col, args.grz_dv_col, args.grz_dv_err_col)
    if grz.yerr is None:
        raise ValueError("grz CSV must include an uncertainty column")

    corr_profile: Optional[Profile] = None
    if args.corr_csv:
        corr_profile = _read_csv_profile(Path(args.corr_csv), args.corr_r_col, args.corr_dv_col)

    results: Dict[str, object] = {
        "analysis_date": datetime.now().isoformat(),
        "inputs": {
            "lensing_mass_csv": args.lensing_mass_csv,
            "dynamics_mass_csv": args.dynamics_mass_csv,
            "grz_csv": args.grz_csv,
            "corr_csv": args.corr_csv,
        },
        "fit_mode": str(args.fit_mode),
        "fit": {},
        "profiles": {
            "grz": {
                "r_mpc": grz.r_mpc.tolist(),
                "dv_kms": grz.y.tolist(),
                "dv_err_kms": grz.yerr.tolist(),
            },
        },
    }

    if corr_profile is not None:
        results["profiles"]["corr"] = {
            "r_mpc": corr_profile.r_mpc.tolist(),
            "dv_kms": corr_profile.y.tolist(),
        }

    def _build_nfw_mass_profile() -> Optional[Profile]:
        m500 = float(args.nfw_m500_msun)
        r500 = float(args.nfw_r500_mpc)
        if not (math.isfinite(m500) and math.isfinite(r500) and m500 > 0 and r500 > 0):
            return None
        c500 = float(args.nfw_c500)
        rmax = float(args.nfw_rmax_r500) * r500
        npts = int(args.nfw_npts)
        if npts < 50:
            npts = 50
        r = np.geomspace(max(1e-3 * r500, 1e-4), rmax, npts)
        m = _nfw_enclosed_mass_msun(r, m500_msun=m500, r500_mpc=r500, c500=c500)
        return Profile(r_mpc=r, y=m, yerr=None)

    def _process_mass(label: str, path_str: str) -> None:
        prof: Optional[Profile] = None
        if path_str:
            prof = _read_csv_profile(Path(path_str), args.mass_r_col, args.mass_m_col)
        else:
            prof = _build_nfw_mass_profile()
        if prof is None:
            return
        phi = _potential_from_enclosed_mass(prof.r_mpc, prof.y)
        dv_pred = _dv_from_phi(phi)

        dv_corr = np.zeros_like(dv_pred)
        if corr_profile is not None:
            dv_corr = np.interp(prof.r_mpc, corr_profile.r_mpc, corr_profile.y)

        dv_model = dv_pred + dv_corr

        dv_on_grz = _interp_to_target(prof.r_mpc, dv_model, grz.r_mpc)

        fixed_chi2, fixed_n = _chi2(dv_on_grz, grz.y, grz.yerr)
        off, off_chi2, off_dof = _weighted_fit_offset(dv_on_grz, grz.y, grz.yerr)
        s, o, so_chi2, so_dof = _weighted_fit_scale_offset(dv_on_grz, grz.y, grz.yerr)

        if args.fit_mode == "fixed":
            s_use = 1.0
            o_use = 0.0
            chi2_use = fixed_chi2
            dof_use = fixed_n
        elif args.fit_mode == "offset_only":
            s_use = 1.0
            o_use = off
            chi2_use = off_chi2
            dof_use = off_dof
        else:
            s_use = s
            o_use = o
            chi2_use = so_chi2
            dof_use = so_dof

        dv_fit = s_use * dv_on_grz + o_use
        resid = grz.y - dv_fit

        results["profiles"][label] = {
            "r_mpc": prof.r_mpc.tolist(),
            "m_enclosed_msun": prof.y.tolist(),
            "phi_kms2": phi.tolist(),
            "dv_pred_kms": dv_pred.tolist(),
            "dv_corr_kms": dv_corr.tolist(),
            "dv_model_kms": dv_model.tolist(),
        }
        results["fit"][label] = {
            "fixed": {
                "scale_s": 1.0,
                "offset_kms": 0.0,
                "chi2": fixed_chi2,
                "dof": fixed_n,
                "chi2_red": (fixed_chi2 / fixed_n) if (fixed_n > 0 and math.isfinite(fixed_chi2)) else float("nan"),
            },
            "offset_only": {
                "scale_s": 1.0,
                "offset_kms": off,
                "chi2": off_chi2,
                "dof": off_dof,
                "chi2_red": (off_chi2 / off_dof) if (off_dof > 0 and math.isfinite(off_chi2)) else float("nan"),
            },
            "scale_offset": {
                "scale_s": s,
                "offset_kms": o,
                "chi2": so_chi2,
                "dof": so_dof,
                "chi2_red": (so_chi2 / so_dof) if (so_dof > 0 and math.isfinite(so_chi2)) else float("nan"),
            },
            "selected": {
                "scale_s": s_use,
                "offset_kms": o_use,
                "chi2": chi2_use,
                "dof": dof_use,
                "chi2_red": (chi2_use / dof_use) if (dof_use > 0 and math.isfinite(chi2_use)) else float("nan"),
            },
            "residuals_at_grz": {
                "r_mpc": grz.r_mpc.tolist(),
                "dv_obs_kms": grz.y.tolist(),
                "dv_model_kms": dv_fit.tolist(),
                "dv_resid_kms": resid.tolist(),
            },
        }

    _process_mass("lensing", args.lensing_mass_csv)
    _process_mass("dynamics", args.dynamics_mass_csv)

    if "lensing" in results["profiles"] and "dynamics" in results["profiles"]:
        dv_l = np.asarray(results["profiles"]["lensing"]["dv_pred_kms"], dtype=float)
        r_l = np.asarray(results["profiles"]["lensing"]["r_mpc"], dtype=float)
        dv_d = np.asarray(results["profiles"]["dynamics"]["dv_pred_kms"], dtype=float)
        r_d = np.asarray(results["profiles"]["dynamics"]["r_mpc"], dtype=float)

        r_common = grz.r_mpc
        dv_lc = np.interp(r_common, r_l, dv_l)
        dv_dc = np.interp(r_common, r_d, dv_d)

        diff = dv_lc - dv_dc
        results["closure"] = {
            "r_mpc": r_common.tolist(),
            "dv_lensing_minus_dynamics_kms": diff.tolist(),
            "rms_kms": float(np.sqrt(np.mean(diff * diff))),
        }

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(results, indent=2), encoding="utf-8")

    _print(f"Wrote: {out_path}")


if __name__ == "__main__":
    main()
