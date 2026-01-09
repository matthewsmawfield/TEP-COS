#!/usr/bin/env python3

import argparse
import json
import math
from pathlib import Path
from typing import Dict, Optional, Tuple

import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.io import fits
import astropy.units as u
from sklearn.linear_model import HuberRegressor, LinearRegression

CMB_RA_DEG = 168.0
CMB_DEC_DEG = -7.0

DEFAULT_DAPTYPE = "HYB10-MILESHC-MASTARSSP"
DEFAULT_RRE_MIN = 0.8
DEFAULT_RRE_MAX = 1.2
DEFAULT_RRE_SYS_MAX = 0.1


def _get_hdu_data(hdul: fits.HDUList, name: str) -> np.ndarray:
    if name not in hdul:
        raise KeyError(f"Missing HDU '{name}'")
    data = hdul[name].data
    if data is None:
        raise RuntimeError(f"Empty HDU '{name}'")
    return np.asarray(data)


def _get_channel_index(header: fits.Header, desired_name: str) -> Optional[int]:
    desired = desired_name.strip()
    for k, v in header.items():
        if not k.startswith("C"):
            continue
        try:
            idx_1 = int(k[1:])
        except Exception:
            continue
        if str(v).strip() == desired:
            return idx_1 - 1
    return None


def _find_channel_index_by_substring(header: fits.Header, substrings: list[str]) -> Optional[int]:
    subs = [s.lower() for s in substrings]
    for k, v in header.items():
        if not k.startswith("C"):
            continue
        try:
            idx_1 = int(k[1:])
        except Exception:
            continue
        name = str(v).strip().lower()
        if any(s in name for s in subs):
            return idx_1 - 1
    return None


def compute_x_cmb(ra_deg: float, dec_deg: float) -> float:
    cmb = SkyCoord(ra=CMB_RA_DEG * u.deg, dec=CMB_DEC_DEG * u.deg)
    gal = SkyCoord(ra=ra_deg * u.deg, dec=dec_deg * u.deg)
    return float(np.cos(gal.separation(cmb).rad))


def robust_slope(x: np.ndarray, y: np.ndarray, w: Optional[np.ndarray] = None) -> float:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if w is None:
        w = np.ones_like(x)
    w = np.asarray(w, dtype=float)
    ok = np.isfinite(x) & np.isfinite(y) & np.isfinite(w) & (w > 0)
    x = x[ok]
    y = y[ok]
    w = w[ok]
    if x.size < 20:
        return float("nan")

    X = x.reshape(-1, 1)
    try:
        model = HuberRegressor(epsilon=1.35, max_iter=200)
        model.fit(X, y, sample_weight=w)
        return float(model.coef_[0])
    except Exception:
        model = LinearRegression()
        model.fit(X, y, sample_weight=w)
        return float(model.coef_[0])


def bootstrap_slope(
    x: np.ndarray,
    y: np.ndarray,
    w: Optional[np.ndarray] = None,
    *,
    n_boot: int = 2000,
    seed: int = 42,
) -> Tuple[float, float, Tuple[float, float]]:
    rng = np.random.default_rng(seed)
    x = np.asarray(x)
    y = np.asarray(y)
    if w is None:
        w = np.ones_like(x)
    w = np.asarray(w)

    ok = np.isfinite(x) & np.isfinite(y) & np.isfinite(w) & (w > 0)
    x = x[ok]
    y = y[ok]
    w = w[ok]
    n = x.size
    if n < 30:
        return float("nan"), float("nan"), (float("nan"), float("nan"))

    slopes = []
    for _ in range(int(n_boot)):
        idx = rng.integers(0, n, size=n)
        s = robust_slope(x[idx], y[idx], w[idx])
        if np.isfinite(s):
            slopes.append(s)

    if len(slopes) < 50:
        return float("nan"), float("nan"), (float("nan"), float("nan"))

    arr = np.array(slopes, dtype=float)
    return float(np.mean(arr)), float(np.std(arr, ddof=1)), (float(np.percentile(arr, 2.5)), float(np.percentile(arr, 97.5)))


def permutation_p_value(
    x: np.ndarray,
    y: np.ndarray,
    w: Optional[np.ndarray] = None,
    *,
    n_perm: int = 5000,
    seed: int = 42,
) -> float:
    rng = np.random.default_rng(seed)
    x = np.asarray(x)
    y = np.asarray(y)
    if w is None:
        w = np.ones_like(x)
    w = np.asarray(w)

    ok = np.isfinite(x) & np.isfinite(y) & np.isfinite(w) & (w > 0)
    x = x[ok]
    y = y[ok]
    w = w[ok]

    a_obs = robust_slope(x, y, w)
    if not np.isfinite(a_obs):
        return float("nan")

    count = 0
    for _ in range(int(n_perm)):
        y_perm = rng.permutation(y)
        a_perm = robust_slope(x, y_perm, w)
        if np.isfinite(a_perm) and abs(a_perm) >= abs(a_obs):
            count += 1
    return float(count / float(n_perm))


def axis_randomization_p_value(
    nvecs: np.ndarray,
    y: np.ndarray,
    w: Optional[np.ndarray],
    a_obs: float,
    *,
    n_axes: int = 2000,
    seed: int = 43,
) -> float:
    rng = np.random.default_rng(seed)
    if w is None:
        w = np.ones_like(y)
    y = np.asarray(y, dtype=float)
    w = np.asarray(w, dtype=float)

    ok = np.isfinite(y) & np.isfinite(w) & (w > 0) & np.all(np.isfinite(nvecs), axis=1)
    nvecs = nvecs[ok]
    y = y[ok]
    w = w[ok]

    if nvecs.shape[0] < 30 or not np.isfinite(a_obs):
        return float("nan")

    hits = 0
    for _ in range(int(n_axes)):
        v = rng.normal(size=3)
        v /= np.linalg.norm(v)
        x_rand = nvecs @ v
        a_rand = robust_slope(x_rand, y, w)
        if np.isfinite(a_rand) and abs(a_rand) >= abs(a_obs):
            hits += 1
    return float(hits / float(n_axes))


def _extract_theta_map(hdul: fits.HDUList) -> np.ndarray:
    ell = _get_hdu_data(hdul, "SPX_ELLCOO")
    if ell.ndim == 2:
        raise RuntimeError("SPX_ELLCOO expected to have channel dimension")

    header = hdul["SPX_ELLCOO"].header
    idx = (
        _get_channel_index(header, "Elliptical azimuth")
        or _get_channel_index(header, "Azimuth")
        or _find_channel_index_by_substring(header, ["azimuth", "az"])
    )
    if idx is None:
        idx = 3
    if idx < 0 or idx >= ell.shape[0]:
        idx = 3
    return np.deg2rad(ell[idx, :, :])


def _extract_rre_map(hdul: fits.HDUList) -> np.ndarray:
    ell = _get_hdu_data(hdul, "SPX_ELLCOO")
    if ell.ndim == 2:
        raise RuntimeError("SPX_ELLCOO expected to have channel dimension")

    header = hdul["SPX_ELLCOO"].header
    idx = (
        _get_channel_index(header, "R/Re")
        or _get_channel_index(header, "R/Reff")
        or _find_channel_index_by_substring(header, ["r/re", "r/reff"])
    )
    if idx is None:
        idx = 1
    if idx < 0 or idx >= ell.shape[0]:
        idx = 1
    return ell[idx, :, :]


def compute_m1_lopsidedness_from_maps(
    maps_path: Path,
    *,
    daptype: str,
    rre_min: float,
    rre_max: float,
    rre_sys_max: float,
    min_spaxels: int = 80,
) -> Optional[Dict[str, float]]:
    try:
        with fits.open(maps_path) as hdul:
            v_map = _get_hdu_data(hdul, "STELLAR_VEL")
            m_map = _get_hdu_data(hdul, "STELLAR_VEL_MASK")
            if v_map.ndim == 3:
                v_map = v_map[0]
            if m_map.ndim == 3:
                m_map = m_map[0]

            rre = _extract_rre_map(hdul)
            theta = _extract_theta_map(hdul)

        good_vr = (m_map == 0) & np.isfinite(v_map) & np.isfinite(rre)
        sys_sel = good_vr & (rre <= rre_sys_max)
        if np.count_nonzero(sys_sel) < 5:
            return None
        v0 = float(np.nanmedian(v_map[sys_sel]))
        v = v_map - v0

        ann = good_vr & (rre >= rre_min) & (rre <= rre_max) & np.isfinite(theta)
        if np.count_nonzero(ann) < min_spaxels:
            return None

        th = theta[ann].astype(float)
        vv = v[ann].astype(float)

        X = np.column_stack([np.sin(th), np.cos(th), np.ones_like(th)])
        try:
            model = HuberRegressor(epsilon=1.35, max_iter=200)
            model.fit(X[:, :2], vv)
            pred = model.predict(X[:, :2])
            resid = vv - pred
        except Exception:
            model = LinearRegression()
            model.fit(X, vv)
            resid = vv - model.predict(X)

        resid = resid - float(np.mean(resid))

        cos_th = np.cos(th)
        sin_th = np.sin(th)
        a1 = 2.0 * float(np.mean(resid * cos_th))
        b1 = 2.0 * float(np.mean(resid * sin_th))
        amp = float(math.sqrt(a1 * a1 + b1 * b1))
        phase_rad = float(math.atan2(b1, a1))

        return {
            "m1_amp": amp,
            "m1_phase_deg": float(np.degrees(phase_rad)),
            "n_spaxels": float(np.count_nonzero(ann)),
            "resid_rms": float(np.std(resid)),
        }
    except Exception:
        return None


def main() -> None:
    parser = argparse.ArgumentParser(description="Step 2.11: Residual-field m=1 lopsidedness test")
    parser.add_argument("--daptype", default=DEFAULT_DAPTYPE)
    parser.add_argument("--maps-dir", default=None)
    parser.add_argument("--dapall", default="data/dapall/dapall-v3_1_1-3.1.0.fits")
    parser.add_argument("--plateifu-list", default="results/outputs/step_1_0_plateifu_selection.txt")
    parser.add_argument("--output-dir", default="results/outputs/m1_lopsidedness")
    parser.add_argument("--max-galaxies", type=int, default=2000)
    parser.add_argument("--rre-min", type=float, default=DEFAULT_RRE_MIN)
    parser.add_argument("--rre-max", type=float, default=DEFAULT_RRE_MAX)
    parser.add_argument("--rre-sys-max", type=float, default=DEFAULT_RRE_SYS_MAX)
    parser.add_argument("--min-spaxels", type=int, default=80)
    parser.add_argument("--n-perm", type=int, default=5000)
    parser.add_argument("--n-axis-rand", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    daptype = str(args.daptype)
    maps_dir = Path(args.maps_dir) if args.maps_dir else Path(f"data/maps/{daptype}")

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    meta: Dict[str, Dict[str, float]] = {}
    with fits.open(args.dapall) as hdul:
        data = hdul[1].data
        for row in data:
            plateifu = row["PLATEIFU"].strip()
            ra = float(row["OBJRA"])
            dec = float(row["OBJDEC"])
            meta[plateifu] = {
                "ra": ra,
                "dec": dec,
                "x_cmb": compute_x_cmb(ra, dec),
            }

    with open(args.plateifu_list) as f:
        plateifus = [ln.strip() for ln in f if ln.strip()]

    rows = []

    for plateifu in plateifus[: int(args.max_galaxies)]:
        plate, ifu = plateifu.split("-")
        maps_path = maps_dir / plate / ifu / f"manga-{plateifu}-MAPS-{daptype}.fits.gz"
        if not maps_path.exists():
            continue
        res = compute_m1_lopsidedness_from_maps(
            maps_path,
            daptype=daptype,
            rre_min=float(args.rre_min),
            rre_max=float(args.rre_max),
            rre_sys_max=float(args.rre_sys_max),
            min_spaxels=int(args.min_spaxels),
        )
        if res is None:
            continue
        m = meta.get(plateifu)
        if m is None:
            continue
        rows.append({
            "plateifu": plateifu,
            "ra_deg": m["ra"],
            "dec_deg": m["dec"],
            "x_cmb": m["x_cmb"],
            **res,
        })

    df = pd.DataFrame(rows)
    df.to_csv(out_dir / "step_2_11_m1_per_galaxy.csv", index=False)

    summary: Dict[str, object] = {
        "n_galaxies": int(df.shape[0]),
        "rre_min": float(args.rre_min),
        "rre_max": float(args.rre_max),
        "rre_sys_max": float(args.rre_sys_max),
        "min_spaxels": int(args.min_spaxels),
    }

    if df.shape[0] >= 30:
        x = df["x_cmb"].to_numpy(float)
        y = df["m1_amp"].to_numpy(float)
        w = np.ones_like(x)

        a_obs = robust_slope(x, y, w)
        mean_a, std_a, ci = bootstrap_slope(x, y, w, n_boot=2000, seed=int(args.seed))
        p_perm = permutation_p_value(x, y, w, n_perm=int(args.n_perm), seed=int(args.seed))

        nvecs = np.stack([
            np.cos(np.deg2rad(df["dec_deg"].to_numpy(float))) * np.cos(np.deg2rad(df["ra_deg"].to_numpy(float))),
            np.cos(np.deg2rad(df["dec_deg"].to_numpy(float))) * np.sin(np.deg2rad(df["ra_deg"].to_numpy(float))),
            np.sin(np.deg2rad(df["dec_deg"].to_numpy(float))),
        ], axis=1)

        p_axis = axis_randomization_p_value(
            nvecs,
            y,
            w,
            a_obs,
            n_axes=int(args.n_axis_rand),
            seed=int(args.seed) + 101,
        )

        summary["m1_amp_slope"] = {
            "slope": float(a_obs),
            "boot_mean": float(mean_a),
            "boot_err": float(std_a),
            "boot_ci_95": [float(ci[0]), float(ci[1])],
            "p_perm": float(p_perm),
            "p_axis_rand": float(p_axis),
        }

    with open(out_dir / "step_2_11_m1_summary.json", "w") as f:
        json.dump(summary, f, indent=2)


if __name__ == "__main__":
    main()
