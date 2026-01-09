#!/usr/bin/env python3

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import argparse
import json
import math
from dataclasses import dataclass
from typing import Dict, Optional, Tuple

import numpy as np
from astropy.io import fits

from scripts.utils.logger import TEPLogger, print_status, set_step_logger


DEFAULT_DAPTYPE = "HYB10-MILESHC-MASTARSSP"
DEFAULT_DRPVER = "v3_1_1"
DEFAULT_CMB_RA_DEG = 168.0
DEFAULT_CMB_DEC_DEG = -7.0
DEFAULT_RRE_MIN = 0.8
DEFAULT_RRE_MAX = 1.2
DEFAULT_RRE_SYS_MAX = 0.1


@dataclass
class ResidualResult:
    plateifu: str
    ra_deg: float
    dec_deg: float
    z: float
    x_cmb: float
    mag_r: float
    re_arcsec: float
    ba: float
    sigma_1re: float
    tf_logvrot: float
    tf_mabs: float
    tf_resid: float
    fp_logre: float
    fp_mu_e: float
    fp_resid: float


def _vec_from_radec_deg(ra_deg: float, dec_deg: float) -> np.ndarray:
    ra = np.deg2rad(ra_deg)
    dec = np.deg2rad(dec_deg)
    return np.array([
        math.cos(dec) * math.cos(ra),
        math.cos(dec) * math.sin(ra),
        math.sin(dec),
    ])


def _plateifu_from_path(path: Path) -> Optional[str]:
    name = path.name
    if not name.startswith("manga-"):
        return None
    parts = name.split("-")
    if len(parts) < 3:
        return None
    try:
        plate = int(parts[1])
        ifu = int(parts[2])
    except Exception:
        return None
    return f"{plate}-{ifu}"


def find_maps_files(data_maps_root: Path, daptype: str) -> list[Path]:
    base = data_maps_root / daptype
    if not base.exists():
        return []
    return sorted(base.rglob(f"manga-*-MAPS-{daptype}.fits.gz"))


def _get_hdu_data(hdul: fits.HDUList, name: str) -> np.ndarray:
    if name not in hdul:
        raise KeyError(f"Missing HDU '{name}'")
    data = hdul[name].data
    if data is None:
        raise RuntimeError(f"Empty HDU '{name}'")
    return np.asarray(data)


def _get_channel_index(header: fits.Header, target: str) -> Optional[int]:
    t = target.strip().lower()
    for k, v in header.items():
        if not str(k).startswith("C"):
            continue
        if str(k).endswith("NAME") and str(v).strip().lower() == t:
            try:
                return int(str(k)[1:-4]) - 1
            except Exception:
                return None
    return None


def _find_channel_index_by_substring(header: fits.Header, substrings: list[str]) -> Optional[int]:
    subs = [s.strip().lower() for s in substrings]
    for k, v in header.items():
        if not str(k).startswith("C"):
            continue
        if not str(k).endswith("NAME"):
            continue
        vv = str(v).strip().lower()
        if any(s in vv for s in subs):
            try:
                return int(str(k)[1:-4]) - 1
            except Exception:
                return None
    return None


def _extract_velocity_map(hdul: fits.HDUList, velocity_source: str) -> Tuple[np.ndarray, np.ndarray]:
    if velocity_source == "stellar":
        v = _get_hdu_data(hdul, "STELLAR_VEL")
        m = _get_hdu_data(hdul, "STELLAR_VEL_MASK")
        return v, m

    raise ValueError("velocity_source must be 'stellar' for Step 2.1 TF")


def _extract_rre_map(hdul: fits.HDUList) -> np.ndarray:
    ell = _get_hdu_data(hdul, "SPX_ELLCOO")
    if ell.ndim == 2:
        raise RuntimeError("SPX_ELLCOO expected to have channel dimension")
    if ell.shape[0] < 2:
        raise RuntimeError("SPX_ELLCOO has fewer than 2 channels")

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


def _estimate_center_xy(good: np.ndarray, rre_map: np.ndarray, rre_sys_max: float) -> Optional[Tuple[float, float]]:
    sel = good & (rre_map <= float(rre_sys_max))
    if not np.any(sel):
        return None
    idx = np.argwhere(sel)
    if idx.size == 0:
        return None
    y0 = float(np.median(idx[:, 0]))
    x0 = float(np.median(idx[:, 1]))
    if not np.isfinite(x0) or not np.isfinite(y0):
        return None
    return x0, y0


def _vrot_from_wedge(
    v_map: np.ndarray,
    mask_map: np.ndarray,
    rre_map: np.ndarray,
    rre_min: float,
    rre_max: float,
    rre_sys_max: float,
    delta_phi_deg: float,
) -> Optional[float]:
    good = (mask_map == 0) & np.isfinite(v_map) & np.isfinite(rre_map)

    sys_sel = good & (rre_map <= float(rre_sys_max))
    if not np.any(sys_sel):
        return None
    v0 = float(np.nanmedian(v_map[sys_sel]))
    v = v_map - v0

    center = _estimate_center_xy(good, rre_map, rre_sys_max)
    if center is None:
        return None
    x0, y0 = center

    ann = good & (rre_map >= float(rre_min)) & (rre_map <= float(rre_max))
    if not np.any(ann):
        return None

    idx = np.argwhere(ann)
    if idx.size == 0:
        return None
    yy = idx[:, 0].astype(float)
    xx = idx[:, 1].astype(float)
    X = xx - x0
    Y = yy - y0
    vv = v[ann].astype(float)

    A = np.column_stack([X, Y, np.ones_like(X)])
    try:
        coef, _, _, _ = np.linalg.lstsq(A, vv, rcond=None)
    except Exception:
        return None

    a = float(coef[0])
    b = float(coef[1])
    gnorm = math.hypot(a, b)
    if not np.isfinite(gnorm) or gnorm <= 0:
        return None

    ux = a / gnorm
    uy = b / gnorm
    along = ux * X + uy * Y
    perp = -uy * X + ux * Y
    phi = np.arctan2(perp, along)

    dphi = math.radians(float(delta_phi_deg))
    w_a = np.abs(phi) <= dphi
    w_b = np.abs(phi) >= (math.pi - dphi)

    n_a = int(np.count_nonzero(w_a))
    n_b = int(np.count_nonzero(w_b))
    if n_a < 20 or n_b < 20:
        return None

    va = vv[w_a]
    vb = vv[w_b]
    m_a = float(np.mean(va))
    m_b = float(np.mean(vb))
    if m_a * m_b >= 0:
        return None

    v_rec = m_a if m_a > 0 else m_b
    v_app = m_b if m_a > 0 else m_a
    vrot = 0.5 * (abs(v_rec) + abs(v_app))
    if not np.isfinite(vrot) or vrot <= 0:
        return None
    return float(vrot)


def _sin_inclination_from_ba(ba: float, q0: float) -> float:
    if not np.isfinite(ba):
        return float("nan")
    b = float(ba)
    q = float(q0)
    if b <= 0 or b > 1 or q <= 0 or q >= 1:
        return float("nan")
    num = b * b - q * q
    den = 1.0 - q * q
    if num <= 0:
        return 1.0
    c2 = num / den
    c2 = min(max(c2, 0.0), 1.0)
    s2 = 1.0 - c2
    return float(math.sqrt(max(s2, 0.0)))


def _find_col(names: list[str], desired: list[str]) -> Optional[str]:
    u = {n.upper(): n for n in names}
    for d in desired:
        if d.upper() in u:
            return u[d.upper()]
    return None


def _norm_plateifu(v: object) -> str:
    if isinstance(v, (bytes, np.bytes_)):
        try:
            return v.decode("utf-8").strip()
        except Exception:
            return str(v).strip()
    return str(v).strip()


def _first_col_containing(names: list[str], substrings: list[str]) -> Optional[str]:
    subs = [s.upper() for s in substrings]
    for n in names:
        nu = n.upper()
        if any(s in nu for s in subs):
            return n
    return None


def _sanitize_output_tag(tag: Optional[str]) -> str:
    if not tag:
        return ""
    t = "".join(ch for ch in str(tag).strip() if (ch.isalnum() or ch in "-_"))
    if not t:
        return ""
    return f"_{t}"


def _load_selection(selection_path: Path) -> list[str]:
    return [
        s.strip()
        for s in selection_path.read_text(encoding="utf-8").splitlines()
        if s.strip() and not s.strip().startswith("#")
    ]


def _robust_huber_fit(x: np.ndarray, y: np.ndarray, w: np.ndarray, huber_c: float = 1.345, max_iter: int = 30) -> Tuple[float, float]:
    wsum = float(np.sum(w))
    xbar = float(np.sum(w * x) / wsum)
    ybar = float(np.sum(w * y) / wsum)
    xx = float(np.sum(w * (x - xbar) ** 2))
    if xx == 0:
        raise RuntimeError("Zero variance in x")
    a = float(np.sum(w * (x - xbar) * (y - ybar)) / xx)
    b = float(ybar - a * xbar)

    for _ in range(int(max_iter)):
        r = y - (a * x + b)
        mad = float(np.median(np.abs(r)))
        s = max(1.4826 * mad, 1e-12)
        u = r / (s * float(huber_c))
        wr = np.ones_like(u)
        m = np.abs(u) > 1.0
        wr[m] = 1.0 / np.abs(u[m])
        w_eff = w * wr

        wsum = float(np.sum(w_eff))
        xbar = float(np.sum(w_eff * x) / wsum)
        ybar = float(np.sum(w_eff * y) / wsum)
        xx = float(np.sum(w_eff * (x - xbar) ** 2))
        if xx == 0:
            break
        a_new = float(np.sum(w_eff * (x - xbar) * (y - ybar)) / xx)
        b_new = float(ybar - a_new * xbar)

        if abs(a_new - a) <= 1e-10 * max(1.0, abs(a)) and abs(b_new - b) <= 1e-10 * max(1.0, abs(b)):
            a, b = a_new, b_new
            break
        a, b = a_new, b_new

    return float(a), float(b)


def _wls_slope_only(x: np.ndarray, y: np.ndarray, w: np.ndarray) -> float:
    wsum = float(np.sum(w))
    xbar = float(np.sum(w * x) / wsum)
    ybar = float(np.sum(w * y) / wsum)
    xx = float(np.sum(w * (x - xbar) ** 2))
    if xx == 0:
        raise RuntimeError("Zero variance in x")
    return float(np.sum(w * (x - xbar) * (y - ybar)) / xx)


def _axis_randomization_test(
    nvecs: np.ndarray,
    y: np.ndarray,
    w: np.ndarray,
    a_obs: float,
    n_axes: int,
    seed: int,
) -> Dict[str, object]:
    rng = np.random.default_rng(int(seed))
    if int(n_axes) <= 0:
        return {"status": "skipped", "reason": "n_axes<=0"}

    n = int(y.size)
    if n < 3:
        return {"status": "skipped", "reason": "insufficient_n"}

    count = 0
    best_abs_a = -1.0
    best_axis: Optional[np.ndarray] = None

    a_abs_obs = abs(float(a_obs))
    for _ in range(int(n_axes)):
        u = float(rng.uniform(0.0, 1.0))
        v = float(rng.uniform(0.0, 1.0))
        phi = 2.0 * math.pi * u
        cos_t = 2.0 * v - 1.0
        sin_t = math.sqrt(max(0.0, 1.0 - cos_t * cos_t))
        axis = np.array([sin_t * math.cos(phi), sin_t * math.sin(phi), cos_t], dtype=float)

        x_rand = nvecs @ axis
        try:
            a_rand = _wls_slope_only(x_rand, y, w)
        except Exception:
            continue

        a_abs = abs(float(a_rand))
        if a_abs >= a_abs_obs:
            count += 1
        if a_abs > best_abs_a:
            best_abs_a = a_abs
            best_axis = axis

    p = (count + 1) / (int(n_axes) + 1)
    out: Dict[str, object] = {
        "n_axes": int(n_axes),
        "seed": int(seed),
        "a_obs": float(a_obs),
        "p_value": float(p),
        "best_abs_a": float(best_abs_a),
    }
    if best_axis is not None:
        ra = float((math.degrees(math.atan2(float(best_axis[1]), float(best_axis[0]))) + 360.0) % 360.0)
        dec = float(math.degrees(math.asin(float(best_axis[2]))))
        out["best_axis_radec_deg"] = {"ra": ra, "dec": dec}
    return out


def _perm_p_value(x: np.ndarray, y: np.ndarray, w: np.ndarray, n_perm: int, seed: int, robust: bool) -> Tuple[float, float, float]:
    rng = np.random.default_rng(seed)

    if robust:
        a_obs, _ = _robust_huber_fit(x, y, w)
    else:
        wsum = float(np.sum(w))
        xbar = float(np.sum(w * x) / wsum)
        ybar = float(np.sum(w * y) / wsum)
        xx = float(np.sum(w * (x - xbar) ** 2))
        if xx == 0:
            raise RuntimeError("Zero variance in x")
        a_obs = float(np.sum(w * (x - xbar) * (y - ybar)) / xx)

    n = int(y.size)
    count_y = 0
    count_pair = 0
    for _ in range(int(n_perm)):
        perm = rng.permutation(n)
        yp = y[perm]
        wp = w[perm]

        if robust:
            ap_y, _ = _robust_huber_fit(x, yp, w)
            ap_p, _ = _robust_huber_fit(x, yp, wp)
        else:
            wsum = float(np.sum(w))
            xbar = float(np.sum(w * x) / wsum)
            xx = float(np.sum(w * (x - xbar) ** 2))
            ybar = float(np.sum(w * yp) / wsum)
            ap_y = float(np.sum(w * (x - xbar) * (yp - ybar)) / xx) if xx != 0 else 0.0

            wsum_p = float(np.sum(wp))
            xbar_p = float(np.sum(wp * x) / wsum_p)
            xx_p = float(np.sum(wp * (x - xbar_p) ** 2))
            ybar_p = float(np.sum(wp * yp) / wsum_p)
            ap_p = float(np.sum(wp * (x - xbar_p) * (yp - ybar_p)) / xx_p) if xx_p != 0 else 0.0

        if abs(ap_y) >= abs(a_obs):
            count_y += 1
        if abs(ap_p) >= abs(a_obs):
            count_pair += 1

    p_y = (count_y + 1) / (int(n_perm) + 1)
    p_pair = (count_pair + 1) / (int(n_perm) + 1)
    return float(a_obs), float(p_pair), float(p_y)


def _re_kpc_from_re_arcsec(z: np.ndarray, re_arcsec: np.ndarray) -> np.ndarray:
    from astropy.cosmology import Planck18

    out = np.full_like(np.asarray(z, dtype=float), np.nan, dtype=float)
    zf = np.asarray(z, dtype=float)
    rf = np.asarray(re_arcsec, dtype=float)
    ok = np.isfinite(zf) & np.isfinite(rf) & (zf > 0.0) & (zf < 1.0) & (rf > 0.0)
    if not np.any(ok):
        return out

    da_kpc = Planck18.angular_diameter_distance(zf[ok]).to_value("kpc")
    arcsec_rad = np.deg2rad(1.0 / 3600.0)
    out[ok] = rf[ok] * arcsec_rad * da_kpc
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description="TEP-COS Step 2.1 - Distance-residual dipole")

    parser.add_argument("--daptype", default=DEFAULT_DAPTYPE)
    parser.add_argument("--drpver", default=DEFAULT_DRPVER)

    parser.add_argument("--velocity-source", choices=["stellar"], default="stellar")
    parser.add_argument("--rre-min", type=float, default=DEFAULT_RRE_MIN)
    parser.add_argument("--rre-max", type=float, default=DEFAULT_RRE_MAX)
    parser.add_argument("--rre-sys-max", type=float, default=DEFAULT_RRE_SYS_MAX)
    parser.add_argument("--delta-phi-deg", type=float, default=20.0)
    parser.add_argument("--tf-q0", type=float, default=0.2)
    parser.add_argument("--tf-ba-min", type=float, default=0.2)
    parser.add_argument("--tf-ba-max", type=float, default=0.85)

    parser.add_argument("--cmb-ra-deg", type=float, default=DEFAULT_CMB_RA_DEG)
    parser.add_argument("--cmb-dec-deg", type=float, default=DEFAULT_CMB_DEC_DEG)

    parser.add_argument("--n-perm", type=int, default=10000)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--min-galaxies", type=int, default=20)
    parser.add_argument("--n-axis-rand", type=int, default=2000, help="Number of random axes for look-elsewhere control")

    parser.add_argument("--plateifu-list", type=str, default=None)
    parser.add_argument("--ignore-step1-selection", action="store_true")

    parser.add_argument(
        "--output-tag",
        type=str,
        default=None,
        help="Optional tag appended to output filenames (CSV/JSON/MD and figures) to prevent overwrites across reruns.",
    )

    args = parser.parse_args()

    logger = TEPLogger(
        "step_2_1_distance_residual_dipole",
        log_file_path=PROJECT_ROOT / "logs" / "step_2_1_distance_residual_dipole.log",
    )
    set_step_logger(logger)

    out_dir = PROJECT_ROOT / "results" / "outputs"
    fig_dir = PROJECT_ROOT / "results" / "figures"
    out_dir.mkdir(parents=True, exist_ok=True)
    fig_dir.mkdir(parents=True, exist_ok=True)

    tag = _sanitize_output_tag(args.output_tag)

    selection_path = PROJECT_ROOT / "results" / "outputs" / "step_1_0_plateifu_selection.txt"
    if args.plateifu_list:
        plateifus = _load_selection(Path(args.plateifu_list))
    elif selection_path.exists() and not args.ignore_step1_selection:
        plateifus = _load_selection(selection_path)
    else:
        raise FileNotFoundError("No PLATEIFU selection list found. Provide --plateifu-list or run Step 1.")

    drpall_dir = PROJECT_ROOT / "data" / "drpall"
    drpall_path = drpall_dir / f"drpall-{args.drpver}.fits"
    if not drpall_path.exists():
        raise FileNotFoundError(f"Missing drpall: {drpall_path}. Re-run Step 1 with --download-drpall.")

    dapall_dir = PROJECT_ROOT / "data" / "dapall"
    dapall_candidates = sorted(dapall_dir.glob("dapall-*.fits*"))
    dapall_path = dapall_candidates[0] if dapall_candidates else None

    with fits.open(drpall_path, memmap=True) as hdul:
        data = hdul[1].data
        if data is None:
            raise RuntimeError("drpall table missing")
        cols = list(data.names)

        c_plateifu = _find_col(cols, ["plateifu"])
        c_ra = _find_col(cols, ["objra", "ifura", "ra"])
        c_dec = _find_col(cols, ["objdec", "ifudec", "dec"])
        c_z = _find_col(cols, ["nsa_z", "z"])
        c_re = _find_col(cols, ["nsa_elpetro_th50_r", "nsa_sersic_th50", "nsa_petro_th50"])
        c_ba = _find_col(cols, ["nsa_elpetro_ba", "nsa_sersic_ba"])
        c_absmag = _find_col(cols, ["nsa_elpetro_absmag", "nsa_sersic_absmag"])
        c_mass = _find_col(cols, ["nsa_elpetro_mass", "nsa_sersic_mass"])

        if c_plateifu is None or c_ra is None or c_dec is None or c_z is None or c_re is None or c_absmag is None:
            raise RuntimeError(
                f"drpall missing required columns: plateifu/ra/dec/z/re/absmag. Found cols={len(cols)}"
            )

        idx_map: Dict[str, int] = {}
        for i, v in enumerate(data[c_plateifu]):
            idx_map[_norm_plateifu(v)] = int(i)

        cmb_vec = _vec_from_radec_deg(args.cmb_ra_deg, args.cmb_dec_deg)

        rows: list[ResidualResult] = []
        for pl in plateifus:
            if pl not in idx_map:
                continue
            i = idx_map[pl]
            ra = float(data[c_ra][i])
            dec = float(data[c_dec][i])
            z = float(data[c_z][i])
            if not np.isfinite(ra) or not np.isfinite(dec) or not np.isfinite(z):
                continue
            if z <= 0.0 or z >= 1.0:
                continue
            re = float(data[c_re][i])
            if not np.isfinite(re) or re <= 0.0:
                continue
            ba = float(data[c_ba][i]) if c_ba is not None else float("nan")
            sigma = float("nan")

            absmag = data[c_absmag][i]
            if isinstance(absmag, np.ndarray):
                mag_r = float(absmag[4]) if absmag.size > 4 else float(absmag[0])
            else:
                mag_r = float(absmag)

            x = float(np.dot(_vec_from_radec_deg(ra, dec), cmb_vec))

            rows.append(
                ResidualResult(
                    plateifu=str(pl),
                    ra_deg=ra,
                    dec_deg=dec,
                    z=z,
                    x_cmb=x,
                    mag_r=mag_r,
                    re_arcsec=re,
                    ba=ba,
                    sigma_1re=sigma,
                    tf_logvrot=float("nan"),
                    tf_mabs=float("nan"),
                    tf_resid=float("nan"),
                    fp_logre=float("nan"),
                    fp_mu_e=float("nan"),
                    fp_resid=float("nan"),
                )
            )

    if len(rows) < max(2, int(args.min_galaxies)):
        raise RuntimeError(f"Insufficient matched drpall rows: N={len(rows)}")

    mag_r = np.array([r.mag_r for r in rows], dtype=float)
    re = np.array([r.re_arcsec for r in rows], dtype=float)
    ba = np.array([r.ba for r in rows], dtype=float)

    sigma = np.array([r.sigma_1re for r in rows], dtype=float)
    if dapall_path is not None and (not np.all(np.isfinite(sigma))):
        with fits.open(dapall_path, memmap=True) as hdul:
            ext = args.daptype
            if ext not in hdul:
                raise RuntimeError(f"dapall missing extension {ext}")
            t = hdul[ext].data
            if t is None:
                raise RuntimeError("dapall extension table missing")
            cols_d = list(t.names)
            c_pl = _find_col(cols_d, ["PLATEIFU"]) or _find_col(cols_d, ["plateifu"])
            c_sig = _find_col(cols_d, ["STELLAR_SIGMA_1RE"]) or _find_col(cols_d, ["stellar_sigma_1re"])
            if c_pl is None or c_sig is None:
                raise RuntimeError("dapall missing PLATEIFU or STELLAR_SIGMA_1RE")
            idx = {_norm_plateifu(v): i for i, v in enumerate(t[c_pl])}
            for j, rr in enumerate(rows):
                k = idx.get(rr.plateifu)
                if k is None:
                    continue
                try:
                    sigma[j] = float(t[c_sig][k])
                except Exception:
                    continue

    z = np.array([r.z for r in rows], dtype=float)
    re_kpc = _re_kpc_from_re_arcsec(z, re)

    with np.errstate(invalid="ignore", divide="ignore"):
        mu_e = mag_r + 2.5 * np.log10(2.0 * math.pi * np.maximum(re_kpc, 1e-6) ** 2)

    ok_fp = np.isfinite(re_kpc) & (re_kpc > 0) & np.isfinite(mu_e) & np.isfinite(sigma) & (sigma > 0)
    fp_logre = np.log10(np.maximum(re_kpc, 1e-12))

    x = np.array([r.x_cmb for r in rows], dtype=float)

    vrot = np.full_like(mag_r, np.nan, dtype=float)
    maps_root = PROJECT_ROOT / "data" / "maps"
    maps_files = find_maps_files(maps_root, args.daptype)
    if maps_files:
        sel = set(r.plateifu for r in rows)
        mp: Dict[str, Path] = {}
        for p in maps_files:
            pl = _plateifu_from_path(p)
            if pl is None:
                continue
            if pl in sel:
                mp[pl] = p

        for i, rr in enumerate(rows):
            path = mp.get(rr.plateifu)
            if path is None:
                continue
            try:
                with fits.open(path, memmap=True) as hdul:
                    v_map, v_mask = _extract_velocity_map(hdul, args.velocity_source)
                    rre_map = _extract_rre_map(hdul)
                v0 = _vrot_from_wedge(
                    v_map=v_map,
                    mask_map=v_mask,
                    rre_map=rre_map,
                    rre_min=args.rre_min,
                    rre_max=args.rre_max,
                    rre_sys_max=args.rre_sys_max,
                    delta_phi_deg=args.delta_phi_deg,
                )
                if v0 is None:
                    continue
                sin_i = _sin_inclination_from_ba(rr.ba, args.tf_q0)
                if not np.isfinite(sin_i) or sin_i <= 0:
                    continue
                if not (args.tf_ba_min <= rr.ba <= args.tf_ba_max):
                    continue
                vrot[i] = float(v0 / max(sin_i, 1e-3))
            except Exception:
                continue

    out_csv = out_dir / f"step_2_1_distance_residual_per_galaxy{tag}.csv"
    out_json = out_dir / f"step_2_1_distance_residual_summary{tag}.json"
    out_md = out_dir / f"step_2_1_distance_residual_report{tag}.md"

    tf_logv = np.full_like(vrot, np.nan)
    with np.errstate(invalid="ignore", divide="ignore"):
        tf_logv = np.log10(vrot)
    tf_mabs = mag_r
    tf_resid = np.full_like(tf_mabs, np.nan)
    tf_summary: Optional[Dict[str, object]] = None
    tf_dipole: Optional[Dict[str, float]] = None
    tf_dipole_robust: Optional[Dict[str, float]] = None
    tf_axis_rand: Optional[Dict[str, object]] = None

    ok_tf = np.isfinite(tf_mabs) & np.isfinite(tf_logv)
    min_tf_n = max(20, int(args.min_galaxies) // 2)
    if np.count_nonzero(ok_tf) >= min_tf_n:
        Xtf = np.column_stack([tf_logv[ok_tf], np.ones(int(np.count_nonzero(ok_tf)))])
        ytf = tf_mabs[ok_tf]
        coef, _, _, _ = np.linalg.lstsq(Xtf, ytf, rcond=None)
        alpha = float(coef[0])
        beta = float(coef[1])
        yhat = Xtf @ coef
        resid = ytf - yhat
        tf_resid[ok_tf] = resid

        xf = x[ok_tf]
        w = np.ones_like(resid)
        a_dip, p_pair, p_y = _perm_p_value(xf, resid, w, n_perm=args.n_perm, seed=args.seed + 7, robust=False)
        a_dip_r, p_pair_r, p_y_r = _perm_p_value(xf, resid, w, n_perm=args.n_perm, seed=args.seed + 7, robust=True)

        tf_summary = {"n": int(xf.size), "tf": {"alpha": alpha, "beta": beta}}
        tf_dipole = {"a_obs": float(a_dip), "p_value_pair": float(p_pair), "p_value_y_only": float(p_y)}
        tf_dipole_robust = {"a_obs": float(a_dip_r), "p_value_pair": float(p_pair_r), "p_value_y_only": float(p_y_r)}

        nvecs_tf = np.array([_vec_from_radec_deg(rows[j].ra_deg, rows[j].dec_deg) for j in np.where(ok_tf)[0]], dtype=float)
        tf_axis_rand = _axis_randomization_test(nvecs_tf, resid, w, a_dip, n_axes=args.n_axis_rand, seed=args.seed + 201)

    lines = [
        "plateifu,ra_deg,dec_deg,z,x_cmb,mag_r,re_arcsec,ba,sigma_1re,tf_logvrot,tf_mabs,tf_resid,fp_logre,fp_mu_e,fp_resid"
    ]

    fp_resid = np.full_like(fp_logre, np.nan)
    fp_summary: Optional[Dict[str, object]] = None
    fp_dipole: Optional[Dict[str, float]] = None
    fp_dipole_robust: Optional[Dict[str, float]] = None
    fp_axis_rand: Optional[Dict[str, object]] = None

    min_fp_n = max(20, int(args.min_galaxies) // 2)
    if np.count_nonzero(ok_fp) >= min_fp_n:
        xf = x[ok_fp]
        yf = fp_logre[ok_fp]
        mf = mu_e[ok_fp]
        sf = np.log10(sigma[ok_fp])

        X = np.column_stack([sf, mf, np.ones_like(sf)])
        coef, _, _, _ = np.linalg.lstsq(X, yf, rcond=None)
        a_fp = float(coef[0])
        b_fp = float(coef[1])
        c_fp = float(coef[2])
        yhat = X @ coef
        resid = yf - yhat
        fp_resid[ok_fp] = resid

        w = np.ones_like(resid)
        a_dip, p_pair, p_y = _perm_p_value(xf, resid, w, n_perm=args.n_perm, seed=args.seed, robust=False)
        a_dip_r, p_pair_r, p_y_r = _perm_p_value(xf, resid, w, n_perm=args.n_perm, seed=args.seed, robust=True)

        fp_summary = {
            "n": int(xf.size),
            "fp": {"a": a_fp, "b": b_fp, "c": c_fp},
        }
        fp_dipole = {"a_obs": float(a_dip), "p_value_pair": float(p_pair), "p_value_y_only": float(p_y)}
        fp_dipole_robust = {"a_obs": float(a_dip_r), "p_value_pair": float(p_pair_r), "p_value_y_only": float(p_y_r)}

        nvecs_fp = np.array([_vec_from_radec_deg(rows[j].ra_deg, rows[j].dec_deg) for j in np.where(ok_fp)[0]], dtype=float)
        fp_axis_rand = _axis_randomization_test(nvecs_fp, resid, w, a_dip, n_axes=args.n_axis_rand, seed=args.seed + 200)

    for i, r in enumerate(rows):
        lines.append(
            ",".join(
                [
                    r.plateifu,
                    f"{r.ra_deg:.8f}",
                    f"{r.dec_deg:.8f}",
                    f"{r.z:.8f}",
                    f"{x[i]:.8f}",
                    f"{mag_r[i]:.6f}",
                    f"{re[i]:.6f}",
                    f"{ba[i]:.6f}",
                    f"{sigma[i]:.6f}",
                    f"{tf_logv[i]:.8g}" if np.isfinite(tf_logv[i]) else "nan",
                    f"{tf_mabs[i]:.6f}",
                    f"{tf_resid[i]:.8g}" if np.isfinite(tf_resid[i]) else "nan",
                    f"{fp_logre[i]:.6f}",
                    f"{mu_e[i]:.6f}",
                    f"{fp_resid[i]:.8g}" if np.isfinite(fp_resid[i]) else "nan",
                ]
            )
        )

    out_csv.write_text("\n".join(lines) + "\n", encoding="utf-8")

    summary: Dict[str, object] = {
        "n_galaxies": int(len(rows)),
        "drpver": str(args.drpver),
        "daptype": str(args.daptype),
        "cmb_ra_deg": float(args.cmb_ra_deg),
        "cmb_dec_deg": float(args.cmb_dec_deg),
        "tf_q0": float(args.tf_q0),
        "tf_ba_min": float(args.tf_ba_min),
        "tf_ba_max": float(args.tf_ba_max),
        "rre_min": float(args.rre_min),
        "rre_max": float(args.rre_max),
        "rre_sys_max": float(args.rre_sys_max),
        "delta_phi_deg": float(args.delta_phi_deg),
        "n_perm": int(args.n_perm),
        "seed": int(args.seed),
        "fp": fp_summary,
        "fp_dipole": fp_dipole,
        "fp_dipole_robust": fp_dipole_robust,
        "fp_axis_randomization": fp_axis_rand,
        "tf": tf_summary,
        "tf_dipole": tf_dipole,
        "tf_dipole_robust": tf_dipole_robust,
        "tf_axis_randomization": tf_axis_rand,
        "note": "TF-like residuals use a wedge-based rotation proxy from MAPS (inclination-corrected via b/a), and absolute magnitudes from drpall. FP-like residuals use dapall sigma_1re and drpall photometry/sizes.",
    }

    out_json.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    out_md.write_text(
        "\n".join(
            [
                "# Distance-Residual Dipole Report (Step 2.1)",
                "",
                f"- Galaxies matched in drpall: N = {int(len(rows))}",
                f"- drpver: {args.drpver}",
                f"- CMB axis (RA, Dec): ({args.cmb_ra_deg:.3f} deg, {args.cmb_dec_deg:.3f} deg)",
                "",
                "## Fundamental Plane (FP-like)",
                f"- status: {'computed' if fp_summary else 'insufficient data'}",
                (f"- dipole p_pair = {fp_dipole['p_value_pair']:.6g}" if fp_dipole else ""),
                (f"- dipole robust p_pair = {fp_dipole_robust['p_value_pair']:.6g}" if fp_dipole_robust else ""),
                (f"- axis-rand p = {fp_axis_rand['p_value']:.6g}" if fp_axis_rand and 'p_value' in fp_axis_rand else ""),
                "",
                "## Tully-Fisher (TF-like)",
                f"- status: {'computed' if tf_summary else 'insufficient data or missing MAPS'}",
                (f"- dipole p_pair = {tf_dipole['p_value_pair']:.6g}" if tf_dipole else ""),
                (f"- dipole robust p_pair = {tf_dipole_robust['p_value_pair']:.6g}" if tf_dipole_robust else ""),
                (f"- axis-rand p = {tf_axis_rand['p_value']:.6g}" if tf_axis_rand and 'p_value' in tf_axis_rand else ""),
                "",
                "## Deliverables",
                f"- Per-galaxy table: `{out_csv.name}`",
                f"- Summary JSON: `{out_json.name}`",
                "",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    print_status(f"Saved per-galaxy table: {out_csv}", "SUCCESS")
    print_status(f"Saved summary: {out_json}", "SUCCESS")
    print_status(f"Saved report: {out_md}", "SUCCESS")
    print_status("Step 2.1 complete", "SUCCESS")

    try:
        if fp_dipole and fp_dipole_robust and np.count_nonzero(ok_fp) >= min_fp_n:
            import matplotlib

            matplotlib.use("Agg")
            import matplotlib.pyplot as plt

            xf = x[ok_fp]
            rf = fp_resid[ok_fp]
            xs = np.linspace(-1, 1, 200)

            fig = plt.figure(figsize=(6.5, 4.0))
            ax = fig.add_subplot(1, 1, 1)
            ax.scatter(xf, rf, s=18, alpha=0.75)

            edges = np.linspace(-1.0, 1.0, 11)
            bx = []
            by = []
            be = []
            for i in range(10):
                lo = edges[i]
                hi = edges[i + 1]
                m = (xf >= lo) & (xf < hi) if i < 9 else (xf >= lo) & (xf <= hi)
                if not np.any(m):
                    continue
                yy = rf[m]
                bx.append(float(np.mean(xf[m])))
                by.append(float(np.mean(yy)))
                be.append(float(np.std(yy, ddof=1) / math.sqrt(max(int(yy.size), 1))) if yy.size > 1 else 0.0)
            if bx:
                ax.errorbar(bx, by, yerr=be, fmt="o", color="#555555", markersize=4, capsize=2, alpha=0.9)

            a = float(fp_dipole["a_obs"])
            ar = float(fp_dipole_robust["a_obs"])
            ax.plot(xs, a * xs + 0.0, color="black", linewidth=2)
            ax.plot(xs, ar * xs + 0.0, color="#C0392B", linewidth=2)

            ax.set_xlabel("x = n · n_CMB")
            ax.set_ylabel("FP residual (log Re - fit)")
            ax.set_title(f"Distance-residual dipole (FP-like; N={int(np.count_nonzero(ok_fp))})")
            ax.grid(True, alpha=0.25)
            ax.text(
                0.02,
                0.98,
                f"OLS a={a:.3g}  P_pair={fp_dipole['p_value_pair']:.2g}  P_y={fp_dipole['p_value_y_only']:.2g}\n"
                f"Robust a={ar:.3g}  P_pair={fp_dipole_robust['p_value_pair']:.2g}  P_y={fp_dipole_robust['p_value_y_only']:.2g}",
                transform=ax.transAxes,
                ha="left",
                va="top",
                fontsize=10,
                bbox={"facecolor": "white", "alpha": 0.85, "edgecolor": "none"},
            )

            fig_path = fig_dir / f"step_2_1_fp_residual_dipole_fit{tag}.png"
            fig_path_pdf = fig_dir / f"step_2_1_fp_residual_dipole_fit{tag}.pdf"
            fig.tight_layout()
            fig.savefig(fig_path, dpi=200)
            fig.savefig(fig_path_pdf)
            plt.close(fig)
            print_status(f"Saved figure: {fig_path}", "SUCCESS")
            print_status(f"Saved figure: {fig_path_pdf}", "SUCCESS")
    except Exception:
        pass

    try:
        if tf_dipole and tf_dipole_robust and np.count_nonzero(ok_tf) >= min_tf_n:
            import matplotlib

            matplotlib.use("Agg")
            import matplotlib.pyplot as plt

            xf = x[ok_tf]
            rf = tf_resid[ok_tf]
            xs = np.linspace(-1, 1, 200)

            fig = plt.figure(figsize=(6.5, 4.0))
            ax = fig.add_subplot(1, 1, 1)
            ax.scatter(xf, rf, s=18, alpha=0.75)

            edges = np.linspace(-1.0, 1.0, 11)
            bx = []
            by = []
            be = []
            for i in range(10):
                lo = edges[i]
                hi = edges[i + 1]
                m = (xf >= lo) & (xf < hi) if i < 9 else (xf >= lo) & (xf <= hi)
                if not np.any(m):
                    continue
                yy = rf[m]
                bx.append(float(np.mean(xf[m])))
                by.append(float(np.mean(yy)))
                be.append(float(np.std(yy, ddof=1) / math.sqrt(max(int(yy.size), 1))) if yy.size > 1 else 0.0)
            if bx:
                ax.errorbar(bx, by, yerr=be, fmt="o", color="#555555", markersize=4, capsize=2, alpha=0.9)

            a = float(tf_dipole["a_obs"])
            ar = float(tf_dipole_robust["a_obs"])
            ax.plot(xs, a * xs + 0.0, color="black", linewidth=2)
            ax.plot(xs, ar * xs + 0.0, color="#C0392B", linewidth=2)

            ax.set_xlabel("x = n · n_CMB")
            ax.set_ylabel("TF residual (M_r - fit)")
            ax.set_title(f"Distance-residual dipole (TF-like; N={int(np.count_nonzero(ok_tf))})")
            ax.grid(True, alpha=0.25)
            ax.text(
                0.02,
                0.98,
                f"OLS a={a:.3g}  P_pair={tf_dipole['p_value_pair']:.2g}  P_y={tf_dipole['p_value_y_only']:.2g}\n"
                f"Robust a={ar:.3g}  P_pair={tf_dipole_robust['p_value_pair']:.2g}  P_y={tf_dipole_robust['p_value_y_only']:.2g}",
                transform=ax.transAxes,
                ha="left",
                va="top",
                fontsize=10,
                bbox={"facecolor": "white", "alpha": 0.85, "edgecolor": "none"},
            )

            fig_path = fig_dir / f"step_2_1_tf_residual_dipole_fit{tag}.png"
            fig_path_pdf = fig_dir / f"step_2_1_tf_residual_dipole_fit{tag}.pdf"
            fig.tight_layout()
            fig.savefig(fig_path, dpi=200)
            fig.savefig(fig_path_pdf)
            plt.close(fig)
            print_status(f"Saved figure: {fig_path}", "SUCCESS")
            print_status(f"Saved figure: {fig_path_pdf}", "SUCCESS")
    except Exception:
        pass


if __name__ == "__main__":
    main()
