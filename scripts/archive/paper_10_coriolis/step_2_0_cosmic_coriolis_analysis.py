#!/usr/bin/env python3

import sys
from pathlib import Path

# Ensure repo root is importable when executing this script directly
PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import argparse
import json
import math
from dataclasses import asdict, dataclass
from typing import Dict, Optional, Tuple

import numpy as np
from astropy.io import fits

from scripts.utils.logger import TEPLogger, print_status, set_step_logger

DEFAULT_DAPTYPE = "HYB10-MILESHC-MASTARSSP"
DEFAULT_RRE_MIN = 0.8
DEFAULT_RRE_MAX = 1.2
DEFAULT_RRE_SYS_MAX = 0.1

DEFAULT_CMB_RA_DEG = 168.0
DEFAULT_CMB_DEC_DEG = -7.0


@dataclass
class GalaxyResult:
    plateifu: str
    ra_deg: float
    dec_deg: float
    x_cmb: float
    n_spaxels: int
    n_pos: int
    n_neg: int
    vpos_mean_abs: float
    vneg_mean_abs: float
    delta_v: float
    delta_v_norm: float
    delta_v_sigma: float
    delta_v_axis: float
    delta_v_axis_norm: float
    delta_v_axis_sigma: float
    n_side_pos: int
    n_side_neg: int
    kpa_deg: float
    delta_v_wedge: float
    delta_v_wedge_norm: float
    delta_v_wedge_sigma: float
    n_wedge_a: int
    n_wedge_b: int


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


def _vec_from_radec_deg(ra_deg: float, dec_deg: float) -> np.ndarray:
    ra = np.deg2rad(ra_deg)
    dec = np.deg2rad(dec_deg)
    return np.array([
        math.cos(dec) * math.cos(ra),
        math.cos(dec) * math.sin(ra),
        math.sin(dec),
    ])


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


def _extract_velocity_map(
    hdul: fits.HDUList,
    velocity_source: str,
    gas_line: str,
) -> Tuple[np.ndarray, np.ndarray]:
    if velocity_source == "stellar":
        v = _get_hdu_data(hdul, "STELLAR_VEL")
        m = _get_hdu_data(hdul, "STELLAR_VEL_MASK")
        return v, m

    if velocity_source != "gas":
        raise ValueError("velocity_source must be 'stellar' or 'gas'")

    v_all = _get_hdu_data(hdul, "EMLINE_GVEL")
    m_all = _get_hdu_data(hdul, "EMLINE_GVEL_MASK")

    header = hdul["EMLINE_GVEL"].header
    idx = _get_channel_index(header, gas_line)
    if idx is None:
        idx = 0

    if v_all.ndim != 3:
        raise RuntimeError("Unexpected EMLINE_GVEL dimensions")

    v = v_all[idx, :, :]
    m = m_all[idx, :, :]
    return v, m


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


def compute_delta_v(
    v_map: np.ndarray,
    mask_map: np.ndarray,
    rre_map: np.ndarray,
    rre_min: float,
    rre_max: float,
    rre_sys_max: float,
) -> Optional[Tuple[float, float, float, int, int, int, float, float]]:
    good = (mask_map == 0) & np.isfinite(v_map) & np.isfinite(rre_map)

    sys_sel = good & (rre_map <= rre_sys_max)
    if not np.any(sys_sel):
        return None

    v0 = float(np.nanmedian(v_map[sys_sel]))
    v = v_map - v0

    ann = good & (rre_map >= rre_min) & (rre_map <= rre_max)
    if not np.any(ann):
        return None

    v_ann = v[ann]

    pos = v_ann > 0
    neg = v_ann < 0

    n_total = int(v_ann.size)
    n_pos = int(np.count_nonzero(pos))
    n_neg = int(np.count_nonzero(neg))

    if n_pos < 20 or n_neg < 20:
        return None

    vpos = np.abs(v_ann[pos])
    vneg = np.abs(v_ann[neg])

    vpos_mean = float(np.mean(vpos))
    vneg_mean = float(np.mean(vneg))

    delta_v = vpos_mean - vneg_mean

    denom = vpos_mean + vneg_mean
    delta_v_norm = float(delta_v / denom) if denom != 0 else 0.0

    vpos_var = float(np.var(vpos, ddof=1)) if vpos.size > 1 else 0.0
    vneg_var = float(np.var(vneg, ddof=1)) if vneg.size > 1 else 0.0

    sigma = math.sqrt(vpos_var / max(n_pos, 1) + vneg_var / max(n_neg, 1))

    return delta_v, delta_v_norm, sigma, n_total, n_pos, n_neg, vpos_mean, vneg_mean


def _estimate_center_xy(good_mask: np.ndarray, rre_map: np.ndarray, rre_sys_max: float) -> Optional[Tuple[float, float]]:
    sys_sel = good_mask & (rre_map <= rre_sys_max)
    if not np.any(sys_sel):
        return None
    idx = np.argwhere(sys_sel)
    if idx.size == 0:
        return None
    y0 = float(np.median(idx[:, 0]))
    x0 = float(np.median(idx[:, 1]))
    return x0, y0


def compute_delta_v_axis(
    v_map: np.ndarray,
    mask_map: np.ndarray,
    rre_map: np.ndarray,
    rre_min: float,
    rre_max: float,
    rre_sys_max: float,
) -> Optional[Tuple[float, float, float, int, int, int, float]]:
    good = (mask_map == 0) & np.isfinite(v_map) & np.isfinite(rre_map)

    sys_sel = good & (rre_map <= rre_sys_max)
    if not np.any(sys_sel):
        return None
    v0 = float(np.nanmedian(v_map[sys_sel]))
    v = v_map - v0

    center = _estimate_center_xy(good, rre_map, rre_sys_max)
    if center is None:
        return None
    x0, y0 = center

    ann = good & (rre_map >= rre_min) & (rre_map <= rre_max)
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
    kpa_deg = float((np.degrees(np.arctan2(uy, ux)) + 360.0) % 180.0)

    s = ux * X + uy * Y
    side_pos = s > 0
    side_neg = s < 0

    n_total = int(vv.size)
    n_pos = int(np.count_nonzero(side_pos))
    n_neg = int(np.count_nonzero(side_neg))
    if n_pos < 20 or n_neg < 20:
        return None

    vpos = np.abs(vv[side_pos])
    vneg = np.abs(vv[side_neg])

    vpos_mean = float(np.mean(vpos))
    vneg_mean = float(np.mean(vneg))
    delta_v = vpos_mean - vneg_mean
    denom = vpos_mean + vneg_mean
    delta_v_norm = float(delta_v / denom) if denom != 0 else 0.0

    vpos_var = float(np.var(vpos, ddof=1)) if vpos.size > 1 else 0.0
    vneg_var = float(np.var(vneg, ddof=1)) if vneg.size > 1 else 0.0
    sigma = math.sqrt(vpos_var / max(n_pos, 1) + vneg_var / max(n_neg, 1))

    return delta_v, delta_v_norm, sigma, n_total, n_pos, n_neg, kpa_deg


def compute_delta_v_wedge(
    v_map: np.ndarray,
    mask_map: np.ndarray,
    rre_map: np.ndarray,
    rre_min: float,
    rre_max: float,
    rre_sys_max: float,
    delta_phi_deg: float,
) -> Optional[Tuple[float, float, float, int, int]]:
    good = (mask_map == 0) & np.isfinite(v_map) & np.isfinite(rre_map)

    sys_sel = good & (rre_map <= rre_sys_max)
    if not np.any(sys_sel):
        return None
    v0 = float(np.nanmedian(v_map[sys_sel]))
    v = v_map - v0

    center = _estimate_center_xy(good, rre_map, rre_sys_max)
    if center is None:
        return None
    x0, y0 = center

    ann = good & (rre_map >= rre_min) & (rre_map <= rre_max)
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

    if m_a > 0:
        v_rec = m_a
        v_app = m_b
        v_rec_s = va
        v_app_s = vb
    else:
        v_rec = m_b
        v_app = m_a
        v_rec_s = vb
        v_app_s = va

    delta_v = v_rec + v_app
    denom = v_rec - v_app
    delta_v_norm = float(delta_v / denom) if denom != 0 else 0.0

    v_rec_var = float(np.var(v_rec_s, ddof=1)) if v_rec_s.size > 1 else 0.0
    v_app_var = float(np.var(v_app_s, ddof=1)) if v_app_s.size > 1 else 0.0
    sigma = math.sqrt(v_rec_var / max(int(v_rec_s.size), 1) + v_app_var / max(int(v_app_s.size), 1))

    return float(delta_v), float(delta_v_norm), float(sigma), int(n_a), int(n_b)

def weighted_linear_fit(x: np.ndarray, y: np.ndarray, w: np.ndarray) -> Tuple[float, float, float, float]:
    if x.size != y.size or x.size != w.size:
        raise ValueError("x, y, w must have same size")

    wsum = np.sum(w)
    xbar = np.sum(w * x) / wsum
    ybar = np.sum(w * y) / wsum

    xx = np.sum(w * (x - xbar) ** 2)
    xy = np.sum(w * (x - xbar) * (y - ybar))

    if xx == 0:
        raise RuntimeError("Zero variance in x")

    a = xy / xx
    b = ybar - a * xbar

    yhat = a * x + b
    resid = y - yhat

    dof = max(int(x.size) - 2, 1)
    s2 = float(np.sum(w * resid**2) / dof)

    a_se = math.sqrt(s2 / xx)

    b_se = math.sqrt(s2 * (1.0 / wsum + (xbar**2) / xx))

    return float(a), float(b), float(a_se), float(b_se)


def robust_huber_fit(
    x: np.ndarray,
    y: np.ndarray,
    w: np.ndarray,
    *,
    huber_c: float = 1.345,
    max_iter: int = 30,
    tol: float = 1e-10,
) -> Tuple[float, float, float, int]:
    if x.size != y.size or x.size != w.size:
        raise ValueError("x, y, w must have same size")

    a, b, _, _ = weighted_linear_fit(x, y, w)

    s = 1.0
    for it in range(int(max_iter)):
        r = y - (a * x + b)
        mad = float(np.median(np.abs(r)))
        s = max(1.4826 * mad, 1e-12)
        u = r / (s * float(huber_c))
        wr = np.ones_like(u)
        m = np.abs(u) > 1.0
        wr[m] = 1.0 / np.abs(u[m])
        w_eff = w * wr

        a_new, b_new, _, _ = weighted_linear_fit(x, y, w_eff)
        if abs(a_new - a) <= tol * max(1.0, abs(a)) and abs(b_new - b) <= tol * max(1.0, abs(b)):
            a, b = a_new, b_new
            return float(a), float(b), float(s), int(it + 1)
        a, b = a_new, b_new

    return float(a), float(b), float(s), int(max_iter)


def robust_permutation_p_value(
    x: np.ndarray,
    y: np.ndarray,
    w: np.ndarray,
    n_perm: int,
    seed: int,
) -> Tuple[float, float, float]:
    rng = np.random.default_rng(seed)
    a_obs, _, _, _ = robust_huber_fit(x, y, w)

    count_y_only = 0
    count_pair = 0
    n = int(y.size)
    for _ in range(int(n_perm)):
        perm = rng.permutation(n)
        yp = y[perm]
        wp = w[perm]

        ap_y, _, _, _ = robust_huber_fit(x, yp, w)
        if abs(ap_y) >= abs(a_obs):
            count_y_only += 1

        ap_p, _, _, _ = robust_huber_fit(x, yp, wp)
        if abs(ap_p) >= abs(a_obs):
            count_pair += 1

    p_y_only = (count_y_only + 1) / (int(n_perm) + 1)
    p_pair = (count_pair + 1) / (int(n_perm) + 1)
    return float(a_obs), float(p_pair), float(p_y_only)


def binned_means(x: np.ndarray, y: np.ndarray, n_bins: int = 10) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    if x.size != y.size:
        raise ValueError("x and y must have same size")

    edges = np.linspace(-1.0, 1.0, int(n_bins) + 1)
    xc: list[float] = []
    ym: list[float] = []
    ys: list[float] = []
    for i in range(int(n_bins)):
        lo = edges[i]
        hi = edges[i + 1]
        if i == int(n_bins) - 1:
            m = (x >= lo) & (x <= hi)
        else:
            m = (x >= lo) & (x < hi)
        if not np.any(m):
            continue
        yy = y[m]
        xc.append(float(np.mean(x[m])))
        ym.append(float(np.mean(yy)))
        ys.append(float(np.std(yy, ddof=1) / math.sqrt(max(int(yy.size), 1))) if yy.size > 1 else 0.0)
    return np.array(xc, dtype=float), np.array(ym, dtype=float), np.array(ys, dtype=float)


def permutation_p_value(
    x: np.ndarray,
    y: np.ndarray,
    w: np.ndarray,
    n_perm: int,
    seed: int,
) -> Tuple[float, float, float]:
    rng = np.random.default_rng(seed)

    wsum = float(np.sum(w))
    xbar = float(np.sum(w * x) / wsum)
    xx = float(np.sum(w * (x - xbar) ** 2))
    if xx == 0:
        raise RuntimeError("Zero variance in x")
    ybar = float(np.sum(w * y) / wsum)
    a_obs = float(np.sum(w * (x - xbar) * (y - ybar)) / xx)

    count_y_only = 0
    count_pair = 0
    n = int(y.size)
    for _ in range(int(n_perm)):
        perm = rng.permutation(n)
        yp = y[perm]
        wp = w[perm]

        ybar_p_y_only = float(np.sum(w * yp) / wsum)
        ap_y_only = float(np.sum(w * (x - xbar) * (yp - ybar_p_y_only)) / xx)
        if abs(ap_y_only) >= abs(a_obs):
            count_y_only += 1

        wsum_p = float(np.sum(wp))
        xbar_p = float(np.sum(wp * x) / wsum_p)
        xx_p = float(np.sum(wp * (x - xbar_p) ** 2))
        if xx_p == 0:
            continue
        ybar_p = float(np.sum(wp * yp) / wsum_p)
        ap = float(np.sum(wp * (x - xbar_p) * (yp - ybar_p)) / xx_p)
        if abs(ap) >= abs(a_obs):
            count_pair += 1

    p_y_only = (count_y_only + 1) / (int(n_perm) + 1)
    p_pair = (count_pair + 1) / (int(n_perm) + 1)
    return a_obs, float(p_pair), float(p_y_only)


def ordinary_linear_fit(x: np.ndarray, y: np.ndarray) -> Tuple[float, float]:
    xbar = float(np.mean(x))
    ybar = float(np.mean(y))
    xx = float(np.sum((x - xbar) ** 2))
    if xx == 0:
        raise RuntimeError("Zero variance in x")
    a = float(np.sum((x - xbar) * (y - ybar)) / xx)
    b = float(ybar - a * xbar)
    return a, b


def _wls_slope_only(x: np.ndarray, y: np.ndarray, w: np.ndarray) -> float:
    wsum = float(np.sum(w))
    xbar = float(np.sum(w * x) / wsum)
    ybar = float(np.sum(w * y) / wsum)
    xx = float(np.sum(w * (x - xbar) ** 2))
    if xx == 0:
        raise RuntimeError("Zero variance in x")
    return float(np.sum(w * (x - xbar) * (y - ybar)) / xx)


def axis_randomization_test(
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


def robustness_diagnostics(x: np.ndarray, y: np.ndarray, w: np.ndarray) -> Dict[str, object]:
    out: Dict[str, object] = {}

    try:
        a_ols, b_ols = ordinary_linear_fit(x, y)
        out["ols"] = {"a": a_ols, "b": b_ols}
    except Exception:
        out["ols"] = None

    n = int(y.size)
    if n <= 250:
        loo: list[float] = []
        for i in range(n):
            mask = np.ones(n, dtype=bool)
            mask[i] = False
            try:
                loo.append(_wls_slope_only(x[mask], y[mask], w[mask]))
            except Exception:
                continue
        if loo:
            loo_arr = np.array(loo, dtype=float)
            out["loo_wls"] = {
                "n": int(loo_arr.size),
                "min": float(np.min(loo_arr)),
                "max": float(np.max(loo_arr)),
                "mean": float(np.mean(loo_arr)),
                "std": float(np.std(loo_arr, ddof=1)) if loo_arr.size > 1 else 0.0,
            }
        else:
            out["loo_wls"] = None
    else:
        out["loo_wls"] = {"skipped": True, "reason": "n>250"}

    abs_y = np.abs(y)
    order = np.argsort(abs_y)
    for k in [1, 2]:
        if n <= k:
            continue
        keep = order[:-k]
        try:
            out[f"trim_abs_y_{k}"] = {
                "a_wls": _wls_slope_only(x[keep], y[keep], w[keep]),
                "a_ols": ordinary_linear_fit(x[keep], y[keep])[0],
            }
        except Exception:
            out[f"trim_abs_y_{k}"] = None

    return out


def find_maps_files(data_maps_root: Path, daptype: str) -> list[Path]:
    base = data_maps_root / daptype
    if not base.exists():
        return []
    return sorted(base.rglob(f"manga-*-MAPS-{daptype}.fits.gz"))


def load_radec_from_header(hdul: fits.HDUList) -> Tuple[float, float]:
    hdr = hdul[0].header
    ra = hdr.get("OBJRA")
    dec = hdr.get("OBJDEC")
    if ra is None or dec is None:
        raise KeyError("OBJRA/OBJDEC not found in PRIMARY header")
    return float(ra), float(dec)


def main() -> None:
    parser = argparse.ArgumentParser(description="TEP-COS Step 2.0 - Cosmic Coriolis feasibility analysis")

    parser.add_argument("--daptype", default=DEFAULT_DAPTYPE)
    parser.add_argument("--velocity-source", choices=["stellar", "gas"], default="stellar")
    parser.add_argument("--gas-line", default="Ha-6564")

    parser.add_argument("--rre-min", type=float, default=DEFAULT_RRE_MIN)
    parser.add_argument("--rre-max", type=float, default=DEFAULT_RRE_MAX)
    parser.add_argument("--rre-sys-max", type=float, default=DEFAULT_RRE_SYS_MAX)
    parser.add_argument("--delta-phi-deg", type=float, default=20.0)

    parser.add_argument("--cmb-ra-deg", type=float, default=DEFAULT_CMB_RA_DEG)
    parser.add_argument("--cmb-dec-deg", type=float, default=DEFAULT_CMB_DEC_DEG)

    parser.add_argument("--min-galaxies", type=int, default=20)

    parser.add_argument("--n-perm", type=int, default=2000)
    parser.add_argument("--n-axis-rand", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=42)

    parser.add_argument(
        "--output-tag",
        type=str,
        default=None,
        help="Optional tag appended to output filenames (CSV/JSON/MD and figures) to prevent overwrites across reruns.",
    )

    parser.add_argument(
        "--plateifu-list",
        type=str,
        default=None,
        help="Optional newline-delimited PLATEIFU list. If provided, only MAPS files matching this list are analyzed.",
    )
    parser.add_argument(
        "--ignore-step1-selection",
        action="store_true",
        help="Ignore results/outputs/step_1_0_plateifu_selection.txt even if it exists.",
    )

    parser.add_argument("--use-norm", action="store_true", help="Fit using delta_v_norm instead of delta_v")

    args = parser.parse_args()

    logger = TEPLogger(
        "step_2_0_cosmic_coriolis_analysis",
        log_file_path=PROJECT_ROOT / "logs" / "step_2_0_cosmic_coriolis_analysis.log",
    )
    set_step_logger(logger)

    maps_root = PROJECT_ROOT / "data" / "maps"
    out_dir = PROJECT_ROOT / "results" / "outputs"
    fig_dir = PROJECT_ROOT / "results" / "figures"

    out_dir.mkdir(parents=True, exist_ok=True)
    fig_dir.mkdir(parents=True, exist_ok=True)

    maps_files = find_maps_files(maps_root, args.daptype)
    if not maps_files:
        raise FileNotFoundError(
            f"No MAPS files found under {maps_root / args.daptype}. Run step_1_0_data_acquisition.py with --download-maps first."
        )

    selection_path = PROJECT_ROOT / "results" / "outputs" / "step_1_0_plateifu_selection.txt"
    plateifu_set: Optional[set[str]] = None
    if args.plateifu_list:
        plateifu_set = set(
            s.strip()
            for s in Path(args.plateifu_list).read_text(encoding="utf-8").splitlines()
            if s.strip() and not s.strip().startswith("#")
        )
    elif selection_path.exists() and not args.ignore_step1_selection:
        plateifu_set = set(
            s.strip()
            for s in selection_path.read_text(encoding="utf-8").splitlines()
            if s.strip() and not s.strip().startswith("#")
        )

    if plateifu_set is not None:
        filtered: list[Path] = []
        missing = 0
        for p in maps_files:
            pl = _plateifu_from_path(p)
            if pl is None:
                continue
            if pl in plateifu_set:
                filtered.append(p)
        missing = max(int(len(plateifu_set)) - int(len(filtered)), 0)
        maps_files = filtered
        print_status(
            f"Found MAPS files (filtered by selection): N={len(maps_files)} | missing_from_list={missing}",
            "SUCCESS",
        )
    else:
        print_status(f"Found MAPS files: N={len(maps_files)}", "SUCCESS")

    cmb_vec = _vec_from_radec_deg(args.cmb_ra_deg, args.cmb_dec_deg)

    results: list[GalaxyResult] = []

    for i, path in enumerate(maps_files, start=1):
        plateifu = _plateifu_from_path(path) or path.stem
        print_status(f"[{i}/{len(maps_files)}] Processing {plateifu}", "PROCESS")

        try:
            with fits.open(path, memmap=True) as hdul:
                ra_deg, dec_deg = load_radec_from_header(hdul)

                v_map, v_mask = _extract_velocity_map(
                    hdul,
                    velocity_source=args.velocity_source,
                    gas_line=args.gas_line,
                )

                rre_map = _extract_rre_map(hdul)

            d = compute_delta_v(
                v_map=v_map,
                mask_map=v_mask,
                rre_map=rre_map,
                rre_min=args.rre_min,
                rre_max=args.rre_max,
                rre_sys_max=args.rre_sys_max,
            )
            if d is None:
                continue

            d_axis = compute_delta_v_axis(
                v_map=v_map,
                mask_map=v_mask,
                rre_map=rre_map,
                rre_min=args.rre_min,
                rre_max=args.rre_max,
                rre_sys_max=args.rre_sys_max,
            )
            if d_axis is None:
                continue

            d_wedge = compute_delta_v_wedge(
                v_map=v_map,
                mask_map=v_mask,
                rre_map=rre_map,
                rre_min=args.rre_min,
                rre_max=args.rre_max,
                rre_sys_max=args.rre_sys_max,
                delta_phi_deg=args.delta_phi_deg,
            )

            delta_v, delta_v_norm, sigma, n_total, n_pos, n_neg, vpos_mean, vneg_mean = d
            delta_v_axis, delta_v_axis_norm, sigma_axis, n_total_axis, n_side_pos, n_side_neg, kpa_deg = d_axis
            if n_total_axis != n_total:
                n_total = int(min(n_total, n_total_axis))

            if d_wedge is None:
                delta_v_wedge = float("nan")
                delta_v_wedge_norm = float("nan")
                sigma_wedge = float("nan")
                n_wedge_a = 0
                n_wedge_b = 0
            else:
                delta_v_wedge, delta_v_wedge_norm, sigma_wedge, n_wedge_a, n_wedge_b = d_wedge

            x_cmb = float(np.dot(_vec_from_radec_deg(ra_deg, dec_deg), cmb_vec))

            results.append(
                GalaxyResult(
                    plateifu=str(plateifu),
                    ra_deg=float(ra_deg),
                    dec_deg=float(dec_deg),
                    x_cmb=x_cmb,
                    n_spaxels=int(n_total),
                    n_pos=int(n_pos),
                    n_neg=int(n_neg),
                    vpos_mean_abs=float(vpos_mean),
                    vneg_mean_abs=float(vneg_mean),
                    delta_v=float(delta_v),
                    delta_v_norm=float(delta_v_norm),
                    delta_v_sigma=float(sigma if sigma > 0 else 1.0),
                    delta_v_axis=float(delta_v_axis),
                    delta_v_axis_norm=float(delta_v_axis_norm),
                    delta_v_axis_sigma=float(sigma_axis if sigma_axis > 0 else 1.0),
                    n_side_pos=int(n_side_pos),
                    n_side_neg=int(n_side_neg),
                    kpa_deg=float(kpa_deg),
                    delta_v_wedge=float(delta_v_wedge),
                    delta_v_wedge_norm=float(delta_v_wedge_norm),
                    delta_v_wedge_sigma=float(sigma_wedge if np.isfinite(sigma_wedge) and sigma_wedge > 0 else 1.0),
                    n_wedge_a=int(n_wedge_a),
                    n_wedge_b=int(n_wedge_b),
                )
            )
        except Exception as e:
            print_status(f"Failed {plateifu}: {e}", "WARNING")
            continue

    tag = "" if not args.output_tag else f"_{str(args.output_tag).strip()}"
    out_csv = out_dir / f"step_2_0_per_galaxy{tag}.csv"
    out_json = out_dir / f"step_2_0_cosmic_coriolis_summary{tag}.json"
    out_md = out_dir / f"step_2_0_cosmic_coriolis_report{tag}.md"

    header = [
        "plateifu",
        "ra_deg",
        "dec_deg",
        "x_cmb",
        "n_spaxels",
        "n_pos",
        "n_neg",
        "vpos_mean_abs",
        "vneg_mean_abs",
        "delta_v",
        "delta_v_norm",
        "delta_v_sigma",
        "delta_v_axis",
        "delta_v_axis_norm",
        "delta_v_axis_sigma",
        "n_side_pos",
        "n_side_neg",
        "kpa_deg",
        "delta_v_wedge",
        "delta_v_wedge_norm",
        "delta_v_wedge_sigma",
        "n_wedge_a",
        "n_wedge_b",
    ]

    lines = [",".join(header)]
    for r in results:
        row = [
            r.plateifu,
            f"{r.ra_deg:.8f}",
            f"{r.dec_deg:.8f}",
            f"{r.x_cmb:.8f}",
            str(r.n_spaxels),
            str(r.n_pos),
            str(r.n_neg),
            f"{r.vpos_mean_abs:.6f}",
            f"{r.vneg_mean_abs:.6f}",
            f"{r.delta_v:.6f}",
            f"{r.delta_v_norm:.8f}",
            f"{r.delta_v_sigma:.6f}",
            f"{r.delta_v_axis:.6f}",
            f"{r.delta_v_axis_norm:.8f}",
            f"{r.delta_v_axis_sigma:.6f}",
            str(r.n_side_pos),
            str(r.n_side_neg),
            f"{r.kpa_deg:.3f}",
            f"{r.delta_v_wedge:.6f}",
            f"{r.delta_v_wedge_norm:.8f}",
            f"{r.delta_v_wedge_sigma:.6f}",
            str(r.n_wedge_a),
            str(r.n_wedge_b),
        ]
        lines.append(",".join(row))

    out_csv.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print_status(f"Saved per-galaxy table: {out_csv}", "SUCCESS")

    if len(results) < max(2, int(args.min_galaxies)):
        summary: Dict[str, object] = {
            "status": "fit_skipped_insufficient_n",
            "n_galaxies": len(results),
            "min_galaxies": int(args.min_galaxies),
            "daptype": args.daptype,
            "velocity_source": args.velocity_source,
            "gas_line": args.gas_line,
            "rre_min": args.rre_min,
            "rre_max": args.rre_max,
            "rre_sys_max": args.rre_sys_max,
            "cmb_ra_deg": args.cmb_ra_deg,
            "cmb_dec_deg": args.cmb_dec_deg,
            "use_norm": bool(args.use_norm),
            "note": "Insufficient galaxies for dipole fit/permutation test. Increase --sample-n and/or lower --min-galaxies for exploratory runs.",
        }
        out_json.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
        out_md.write_text(
            "\n".join(
                [
                    "# Cosmic Coriolis Feasibility Report (Step 2.0)",
                    "",
                    f"- Dataset: MaNGA DAP MAPS ({args.daptype})",
                    f"- Velocity source: {args.velocity_source}" + (f" ({args.gas_line})" if args.velocity_source == "gas" else ""),
                    f"- Radial window: R/Re in [{args.rre_min}, {args.rre_max}] (systematics cutoff {args.rre_sys_max})",
                    f"- CMB axis (RA, Dec): ({args.cmb_ra_deg:.3f} deg, {args.cmb_dec_deg:.3f} deg)",
                    f"- Galaxies passing QC: N = {len(results)} (min required: {int(args.min_galaxies)})",
                    "",
                    "## Outcome", 
                    "Insufficient sample size for dipole fitting and permutation testing.",
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
        print_status(
            f"Fit skipped (insufficient N): N={len(results)} < min={int(args.min_galaxies)}. Wrote summary/report and exiting cleanly.",
            "INFO",
        )
        print_status("Step 2.0 complete", "SUCCESS")
        return

    x = np.array([r.x_cmb for r in results], dtype=float)

    y_legacy = np.array([r.delta_v_norm if args.use_norm else r.delta_v for r in results], dtype=float)
    sigma_legacy = np.array([r.delta_v_sigma for r in results], dtype=float)
    # BUG FIX: Use minimum sigma of 0.1 km/s to prevent infinite weights from sigma=0 galaxies
    sigma_legacy_clipped = np.maximum(sigma_legacy, 0.1)
    w_legacy = 1.0 / sigma_legacy_clipped ** 2

    y_axis = np.array([r.delta_v_axis_norm if args.use_norm else r.delta_v_axis for r in results], dtype=float)
    sigma_axis = np.array([r.delta_v_axis_sigma for r in results], dtype=float)
    # BUG FIX: Use minimum sigma of 0.1 km/s
    sigma_axis_clipped = np.maximum(sigma_axis, 0.1)
    w_axis = 1.0 / sigma_axis_clipped ** 2

    y_wedge_all = np.array([r.delta_v_wedge_norm if args.use_norm else r.delta_v_wedge for r in results], dtype=float)
    sigma_wedge_all = np.array([r.delta_v_wedge_sigma for r in results], dtype=float)
    # BUG FIX: Use minimum sigma of 0.1 km/s
    sigma_wedge_clipped = np.maximum(sigma_wedge_all, 0.1)
    w_wedge_all = 1.0 / sigma_wedge_clipped ** 2
    ok_wedge = np.isfinite(y_wedge_all) & np.isfinite(sigma_wedge_all)
    x_wedge = x[ok_wedge]
    y_wedge = y_wedge_all[ok_wedge]
    w_wedge = w_wedge_all[ok_wedge]

    a, b, a_se, b_se = weighted_linear_fit(x, y_legacy, w_legacy)
    a_obs, p_perm_pair, p_perm_y_only = permutation_p_value(x, y_legacy, w_legacy, n_perm=args.n_perm, seed=args.seed)
    a_r, b_r, s_r, it_r = robust_huber_fit(x, y_legacy, w_legacy)
    a_r_obs, p_r_pair, p_r_y_only = robust_permutation_p_value(x, y_legacy, w_legacy, n_perm=args.n_perm, seed=args.seed)

    a_ax, b_ax, a_ax_se, b_ax_se = weighted_linear_fit(x, y_axis, w_axis)
    a_ax_obs, p_ax_pair, p_ax_y_only = permutation_p_value(x, y_axis, w_axis, n_perm=args.n_perm, seed=args.seed)
    a_ax_r, b_ax_r, s_ax_r, it_ax_r = robust_huber_fit(x, y_axis, w_axis)
    a_ax_r_obs, p_ax_r_pair, p_ax_r_y_only = robust_permutation_p_value(x, y_axis, w_axis, n_perm=args.n_perm, seed=args.seed)

    wedge_fit: Optional[Dict[str, float]] = None
    wedge_perm: Optional[Dict[str, float]] = None
    wedge_robust_fit: Optional[Dict[str, float]] = None
    wedge_robust_perm: Optional[Dict[str, float]] = None
    min_wedge_n = max(20, int(args.min_galaxies) // 2)
    if x_wedge.size >= min_wedge_n:
        a_w, b_w, a_w_se, b_w_se = weighted_linear_fit(x_wedge, y_wedge, w_wedge)
        a_w_obs, p_w_pair, p_w_y = permutation_p_value(x_wedge, y_wedge, w_wedge, n_perm=args.n_perm, seed=args.seed)
        a_w_r, b_w_r, s_w_r, it_w_r = robust_huber_fit(x_wedge, y_wedge, w_wedge)
        a_w_r_obs, p_w_r_pair, p_w_r_y = robust_permutation_p_value(x_wedge, y_wedge, w_wedge, n_perm=args.n_perm, seed=args.seed)
        wedge_fit = {"a": float(a_w), "b": float(b_w), "a_se": float(a_w_se), "b_se": float(b_w_se), "n": int(x_wedge.size)}
        wedge_perm = {"a_obs": float(a_w_obs), "p_value_pair": float(p_w_pair), "p_value_y_only": float(p_w_y), "n": int(x_wedge.size)}
        wedge_robust_fit = {"a": float(a_w_r), "b": float(b_w_r), "scale": float(s_w_r), "n_iter": int(it_w_r), "n": int(x_wedge.size)}
        wedge_robust_perm = {"a_obs": float(a_w_r_obs), "p_value_pair": float(p_w_r_pair), "p_value_y_only": float(p_w_r_y), "n": int(x_wedge.size)}

    nvecs = np.array([_vec_from_radec_deg(r.ra_deg, r.dec_deg) for r in results], dtype=float)
    axis_rand_legacy = axis_randomization_test(nvecs, y_legacy, w_legacy, a_obs, n_axes=args.n_axis_rand, seed=args.seed + 101)
    axis_rand_axis = axis_randomization_test(nvecs, y_axis, w_axis, a_ax_obs, n_axes=args.n_axis_rand, seed=args.seed + 102)
    axis_rand_wedge = (
        axis_randomization_test(
            nvecs[ok_wedge],
            y_wedge,
            w_wedge,
            float(wedge_perm["a_obs"]),
            n_axes=args.n_axis_rand,
            seed=args.seed + 103,
        )
        if (x_wedge.size and wedge_perm)
        else {"status": "skipped", "reason": "wedge_fit_not_available"}
    )

    summary: Dict[str, object] = {
        "n_galaxies": len(results),
        "daptype": args.daptype,
        "velocity_source": args.velocity_source,
        "gas_line": args.gas_line,
        "rre_min": args.rre_min,
        "rre_max": args.rre_max,
        "rre_sys_max": args.rre_sys_max,
        "delta_phi_deg": float(args.delta_phi_deg),
        "cmb_ra_deg": args.cmb_ra_deg,
        "cmb_dec_deg": args.cmb_dec_deg,
        "use_norm": bool(args.use_norm),
        "n_axis_rand": int(args.n_axis_rand),
        "fit": {
            "a": a,
            "b": b,
            "a_se": a_se,
            "b_se": b_se,
        },
        "permutation": {
            "n_perm": args.n_perm,
            "seed": args.seed,
            "a_obs": a_obs,
            "p_value_pair": p_perm_pair,
            "p_value_y_only": p_perm_y_only,
        },
        "robust_fit": {
            "a": a_r,
            "b": b_r,
            "scale": s_r,
            "n_iter": int(it_r),
        },
        "robust_permutation": {
            "n_perm": args.n_perm,
            "seed": args.seed,
            "a_obs": a_r_obs,
            "p_value_pair": p_r_pair,
            "p_value_y_only": p_r_y_only,
        },
        "fit_axis": {
            "a": a_ax,
            "b": b_ax,
            "a_se": a_ax_se,
            "b_se": b_ax_se,
        },
        "permutation_axis": {
            "n_perm": args.n_perm,
            "seed": args.seed,
            "a_obs": a_ax_obs,
            "p_value_pair": p_ax_pair,
            "p_value_y_only": p_ax_y_only,
        },
        "robust_fit_axis": {
            "a": a_ax_r,
            "b": b_ax_r,
            "scale": s_ax_r,
            "n_iter": int(it_ax_r),
        },
        "robust_permutation_axis": {
            "n_perm": args.n_perm,
            "seed": args.seed,
            "a_obs": a_ax_r_obs,
            "p_value_pair": p_ax_r_pair,
            "p_value_y_only": p_ax_r_y_only,
        },
        "fit_wedge": wedge_fit,
        "permutation_wedge": wedge_perm,
        "robust_fit_wedge": wedge_robust_fit,
        "robust_permutation_wedge": wedge_robust_perm,
        "robustness": {
            "legacy": robustness_diagnostics(x, y_legacy, w_legacy),
            "axis": robustness_diagnostics(x, y_axis, w_axis),
            "wedge": robustness_diagnostics(x_wedge, y_wedge, w_wedge) if x_wedge.size else None,
        },
        "axis_randomization": {
            "legacy": axis_rand_legacy,
            "axis": axis_rand_axis,
            "wedge": axis_rand_wedge,
        },
    }

    out_json.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print_status(f"Saved summary: {out_json}", "SUCCESS")

    out_md.write_text(
        "\n".join(
            [
                "# Cosmic Coriolis Feasibility Report (Step 2.0)",
                "",
                f"- Dataset: MaNGA DAP MAPS ({args.daptype})",
                f"- Velocity source: {args.velocity_source}" + (f" ({args.gas_line})" if args.velocity_source == "gas" else ""),
                f"- Radial window: R/Re in [{args.rre_min}, {args.rre_max}] (systematics cutoff {args.rre_sys_max})",
                f"- CMB axis (RA, Dec): ({args.cmb_ra_deg:.3f} deg, {args.cmb_dec_deg:.3f} deg)",
                f"- Galaxies passing QC: N = {len(results)}",
                "",
                "## Dipole-fit model", 
                "We fit weighted linear models of the form: y = a x + b, where x = n · n_CMB.",
                "",
                "### Legacy observable (sign-split deltaV)",
                f"- a = {a:.6g} ± {a_se:.3g}",
                f"- b = {b:.6g} ± {b_se:.3g}",
                f"- robust a = {a_r:.6g} (P_pair={p_r_pair:.6g}, P_y={p_r_y_only:.6g})",
                "",
                "### Kinematic-axis observable (hemisphere deltaV)",
                f"- a = {a_ax:.6g} ± {a_ax_se:.3g}",
                f"- b = {b_ax:.6g} ± {b_ax_se:.3g}",
                f"- robust a = {a_ax_r:.6g} (P_pair={p_ax_r_pair:.6g}, P_y={p_ax_r_y_only:.6g})",
                "",
                "### Wedge observable (major-axis wedges)",
                f"- N_wedge = {int(x_wedge.size)}",
                f"- a = {wedge_fit['a']:.6g} ± {wedge_fit['a_se']:.3g}" if wedge_fit else "- a = n/a",
                f"- b = {wedge_fit['b']:.6g} ± {wedge_fit['b_se']:.3g}" if wedge_fit else "- b = n/a",
                f"- robust a = {wedge_robust_fit['a']:.6g} (P_pair={wedge_robust_perm['p_value_pair']:.6g}, P_y={wedge_robust_perm['p_value_y_only']:.6g})" if wedge_robust_fit and wedge_robust_perm else "- robust a = n/a",
                "",
                "## Permutation test", 
                f"- n_perm = {int(args.n_perm)} (seed={int(args.seed)})",
                "### Legacy observable (sign-split)",
                f"- p_value_pair(a) = {p_perm_pair:.6g}",
                f"- p_value_y_only(a) = {p_perm_y_only:.6g}",
                "",
                "### Kinematic-axis observable (hemisphere)",
                f"- p_value_pair(a) = {p_ax_pair:.6g}",
                f"- p_value_y_only(a) = {p_ax_y_only:.6g}",
                "",
                "## Robustness diagnostics", 
                "The summary JSON includes OLS fits, leave-one-out WLS slope stability (if N≤250), and trimmed |y| sensitivity (k=1,2).",
                "",
                "## Look-elsewhere control (axis randomization)",
                f"- n_axis_rand = {int(args.n_axis_rand)}",
                f"- P_axis_rand(|a|>=|a_CMB|), legacy = {axis_rand_legacy.get('p_value', float('nan')):.6g}" if isinstance(axis_rand_legacy, dict) else "- legacy = n/a",
                f"- P_axis_rand(|a|>=|a_CMB|), kinematic-axis = {axis_rand_axis.get('p_value', float('nan')):.6g}" if isinstance(axis_rand_axis, dict) else "- axis = n/a",
                f"- P_axis_rand(|a|>=|a_CMB|), wedge = {axis_rand_wedge.get('p_value', float('nan')):.6g}" if isinstance(axis_rand_wedge, dict) else "- wedge = n/a",
                "",
                "## Deliverables", 
                f"- Per-galaxy table: `{out_csv.name}`",
                f"- Summary JSON: `{out_json.name}`",
                f"- Figures: `step_2_0_cmb_dipole_fit{tag}.(png|pdf)`, `step_2_0_cmb_dipole_fit_axis{tag}.(png|pdf)`, and `step_2_0_cmb_dipole_fit_wedge{tag}.(png|pdf)`",
                "",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    print_status(f"Saved report: {out_md}", "SUCCESS")

    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        try:
            plt.rcParams.update(
                {
                    "font.family": "serif",
                    "font.serif": ["Times New Roman", "STIXGeneral", "DejaVu Serif"],
                    "mathtext.fontset": "stix",
                    "axes.labelsize": 11,
                    "axes.titlesize": 12,
                    "xtick.labelsize": 10,
                    "ytick.labelsize": 10,
                }
            )
        except Exception:
            pass

        def _save_one(
            xv: np.ndarray,
            yv: np.ndarray,
            aa: float,
            bb: float,
            aa_se: float,
            pp: float,
            py: float,
            aa_r: float,
            bb_r: float,
            pp_r: float,
            py_r: float,
            stem: str,
            title: str,
        ) -> None:
            fig = plt.figure(figsize=(6.5, 4.0))
            ax = fig.add_subplot(1, 1, 1)
            ax.scatter(xv, yv, s=18, alpha=0.75)

            bx, by, be = binned_means(xv, yv, n_bins=10)
            if bx.size:
                ax.errorbar(bx, by, yerr=be, fmt="o", color="#555555", markersize=4, capsize=2, alpha=0.9)

            xs = np.linspace(-1, 1, 200)
            ax.plot(xs, aa * xs + bb, color="black", linewidth=2)
            ax.plot(xs, aa_r * xs + bb_r, color="#C0392B", linewidth=2)

            ax.set_xlabel("x = n · n_CMB")
            ax.set_ylabel("deltaV_norm" if args.use_norm else "deltaV")
            ax.set_title(title)
            ax.grid(True, alpha=0.25)

            ax.text(
                0.02,
                0.98,
                f"WLS a={aa:.3g}±{aa_se:.2g}  P_pair={pp:.2g}  P_y={py:.2g}\nRobust a={aa_r:.3g}  P_pair={pp_r:.2g}  P_y={py_r:.2g}",
                transform=ax.transAxes,
                ha="left",
                va="top",
                fontsize=10,
                bbox={"facecolor": "white", "alpha": 0.85, "edgecolor": "none"},
            )

            fig_path = fig_dir / f"{stem}.png"
            fig_path_pdf = fig_dir / f"{stem}.pdf"
            fig.tight_layout()
            fig.savefig(fig_path, dpi=200)
            fig.savefig(fig_path_pdf)
            plt.close(fig)

            print_status(f"Saved figure: {fig_path}", "SUCCESS")
            print_status(f"Saved figure: {fig_path_pdf}", "SUCCESS")

        _save_one(
            x,
            y_legacy,
            a,
            b,
            a_se,
            p_perm_pair,
            p_perm_y_only,
            a_r,
            b_r,
            p_r_pair,
            p_r_y_only,
            f"step_2_0_cmb_dipole_fit{tag}",
            f"Cosmic Coriolis (legacy sign-split; N={len(results)})",
        )
        _save_one(
            x,
            y_axis,
            a_ax,
            b_ax,
            a_ax_se,
            p_ax_pair,
            p_ax_y_only,
            a_ax_r,
            b_ax_r,
            p_ax_r_pair,
            p_ax_r_y_only,
            f"step_2_0_cmb_dipole_fit_axis{tag}",
            f"Cosmic Coriolis (kinematic-axis hemispheres; N={len(results)})",
        )

        if wedge_fit and wedge_perm and wedge_robust_fit and wedge_robust_perm:
            _save_one(
                x_wedge,
                y_wedge,
                float(wedge_fit["a"]),
                float(wedge_fit["b"]),
                float(wedge_fit["a_se"]),
                float(wedge_perm["p_value_pair"]),
                float(wedge_perm["p_value_y_only"]),
                float(wedge_robust_fit["a"]),
                float(wedge_robust_fit["b"]),
                float(wedge_robust_perm["p_value_pair"]),
                float(wedge_robust_perm["p_value_y_only"]),
                f"step_2_0_cmb_dipole_fit_wedge{tag}",
                f"Cosmic Coriolis (major-axis wedges; N={int(x_wedge.size)})",
            )
    except Exception as e:
        print_status(f"Figure generation failed: {e}", "WARNING")

    print_status("Step 2.0 complete", "SUCCESS")


if __name__ == "__main__":
    main()
