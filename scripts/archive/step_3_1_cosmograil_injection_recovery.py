#!/usr/bin/env python3

import argparse
import json
import math
import sys
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
from scipy import interpolate, signal
from scipy.ndimage import gaussian_filter1d


@dataclass
class LightCurve:
    label: str
    t: np.ndarray
    y: np.ndarray
    yerr: np.ndarray


def _print(msg: str) -> None:
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}")


def _read_rdb(path: Path) -> Tuple[List[str], Dict[str, LightCurve]]:
    lines = path.read_text(encoding="utf-8").splitlines()

    header: Optional[List[str]] = None
    data_start: Optional[int] = None
    for i, line in enumerate(lines):
        s = line.strip()
        if not s:
            continue
        if header is None and ("mhjd" in s.lower() or "mjd" in s.lower()):
            header = s.split()
            continue
        if header is not None and ("====" in s or "----" in s):
            continue
        if header is not None:
            data_start = i
            break

    if header is None or data_start is None:
        raise ValueError(f"Could not parse header/data in {path}")

    labels: List[str] = []
    mag_col: Dict[str, int] = {}
    err_col: Dict[str, int] = {}

    for j, c in enumerate(header):
        cl = c.lower()
        if cl.startswith("mag_") and "err" not in cl:
            lab = c.split("_", 1)[1]
            labels.append(lab)
            mag_col[lab] = j
        if cl.startswith("magerr_"):
            lab = c.split("_", 1)[1]
            err_col[lab] = j

    data: Dict[str, Dict[str, List[float]]] = {lab: {"t": [], "y": [], "e": []} for lab in labels}

    for line in lines[data_start:]:
        s = line.strip()
        if not s or "====" in s or s.startswith("%"):
            continue
        parts = s.split()
        if len(parts) < 2:
            continue
        try:
            t = float(parts[0])
        except Exception:
            continue

        for lab in labels:
            mi = mag_col.get(lab)
            ei = err_col.get(lab)
            if mi is None or mi >= len(parts):
                continue
            try:
                y = float(parts[mi])
                e = float(parts[ei]) if ei is not None and ei < len(parts) else 0.01
            except Exception:
                continue
            if not (math.isfinite(t) and math.isfinite(y) and math.isfinite(e)):
                continue
            data[lab]["t"].append(t)
            data[lab]["y"].append(y)
            data[lab]["e"].append(e)

    lcs: Dict[str, LightCurve] = {}
    for lab in labels:
        t = np.asarray(data[lab]["t"], dtype=float)
        y = np.asarray(data[lab]["y"], dtype=float)
        e = np.asarray(data[lab]["e"], dtype=float)
        if t.size < 10:
            continue
        order = np.argsort(t)
        lcs[lab] = LightCurve(label=lab, t=t[order], y=y[order], yerr=e[order])

    return labels, lcs


def _uniform_grid(t1: np.ndarray, t2: np.ndarray, dt: float) -> np.ndarray:
    t0 = max(float(np.min(t1)), float(np.min(t2)))
    t1m = min(float(np.max(t1)), float(np.max(t2)))
    if not (t1m > t0):
        return np.array([], dtype=float)
    n = int(math.floor((t1m - t0) / dt)) + 1
    if n < 50:
        return np.array([], dtype=float)
    return t0 + dt * np.arange(n, dtype=float)


def _interp_to_grid(t: np.ndarray, y: np.ndarray, grid: np.ndarray) -> np.ndarray:
    f = interpolate.interp1d(t, y, kind="linear", bounds_error=False, fill_value=np.nan)
    return f(grid)


def _fill_nans(y: np.ndarray) -> np.ndarray:
    out = y.copy()
    good = np.isfinite(out)
    if not np.any(good):
        return np.zeros_like(out)
    med = float(np.nanmedian(out[good]))
    out[~good] = med
    return out


def _lowpass_detrend(y: np.ndarray, sigma_samples: float) -> Tuple[np.ndarray, np.ndarray]:
    if sigma_samples <= 0:
        return y.copy(), np.full_like(y, np.nan)

    good = np.isfinite(y)
    if not np.any(good):
        return y.copy(), np.full_like(y, np.nan)

    y_filled = _fill_nans(y)
    low = gaussian_filter1d(y_filled, sigma=sigma_samples, mode="nearest")
    out = y - low
    out[~good] = np.nan
    low[~good] = np.nan
    return out, low


def _bandpass_dog(y: np.ndarray, sigma1: float, sigma2: float) -> np.ndarray:
    good = np.isfinite(y)
    if not np.any(good):
        return y.copy()

    y_filled = _fill_nans(y)
    g1 = gaussian_filter1d(y_filled, sigma=sigma1, mode="nearest")
    g2 = gaussian_filter1d(y_filled, sigma=sigma2, mode="nearest")
    out = g1 - g2
    out[~good] = np.nan
    return out


def _corrcoef_nan(a: np.ndarray, b: np.ndarray) -> Tuple[float, int]:
    good = np.isfinite(a) & np.isfinite(b)
    n = int(np.count_nonzero(good))
    if n < 10:
        return float("nan"), n
    aa = a[good]
    bb = b[good]
    sa = float(np.std(aa))
    sb = float(np.std(bb))
    if sa <= 0 or sb <= 0:
        return float("nan"), n
    r = float(np.corrcoef(aa, bb)[0, 1])
    return r, n


def _estimate_delay_corr(
    a: np.ndarray,
    b: np.ndarray,
    dt: float,
    lag_max_days: float,
    min_points: int,
) -> Tuple[float, float, int]:
    nlag = int(round(lag_max_days / dt))
    if nlag < 1:
        return 0.0, float("nan"), 0

    best_r = -1.0
    best_k = 0
    best_n = 0
    r_by_k: Dict[int, float] = {}

    for k in range(-nlag, nlag + 1):
        if k >= 0:
            aa = a[k:]
            bb = b[: a.size - k]
        else:
            kk = -k
            aa = a[: a.size - kk]
            bb = b[kk:]

        r, n = _corrcoef_nan(aa, bb)
        r_by_k[k] = r
        if n < min_points or not math.isfinite(r):
            continue
        if r > best_r:
            best_r = r
            best_k = k
            best_n = n

    if not math.isfinite(best_r):
        return best_k * dt, best_r, best_n

    r0 = r_by_k.get(best_k, float("nan"))
    r1 = r_by_k.get(best_k - 1, float("nan"))
    r2 = r_by_k.get(best_k + 1, float("nan"))

    delta = 0.0
    if math.isfinite(r0) and math.isfinite(r1) and math.isfinite(r2):
        denom = (r1 - 2.0 * r0 + r2)
        if denom != 0.0:
            delta = 0.5 * (r1 - r2) / denom
            delta = float(np.clip(delta, -1.0, 1.0))

    return (best_k + delta) * dt, float(best_r), int(best_n)


def _shift_series(y: np.ndarray, grid: np.ndarray, shift_days: float) -> np.ndarray:
    f = interpolate.interp1d(grid, y, kind="linear", bounds_error=False, fill_value=np.nan)
    return f(grid - shift_days)


def _inject_two_scale_warp(
    src: np.ndarray,
    grid: np.ndarray,
    base_delay_days: float,
    extra_slow_delay_days: float,
    slow_sigma_days: float,
    dt: float,
) -> np.ndarray:
    sigma = slow_sigma_days / dt
    if sigma < 1.0:
        sigma = 1.0

    src_filled = _fill_nans(src)
    slow = gaussian_filter1d(src_filled, sigma=sigma, mode="nearest")
    fast = src_filled - slow

    slow_shift = _shift_series(slow, grid, base_delay_days + extra_slow_delay_days)
    fast_shift = _shift_series(fast, grid, base_delay_days)

    out = slow_shift + fast_shift
    out[~np.isfinite(out)] = np.nan
    return out


def _fit_delta_delay(
    tau: List[float],
    delays: List[float],
    corrs: List[float],
    tau_fast_max: float,
    tau_slow_min: float,
    min_corr: float,
) -> Tuple[float, int, int]:
    fast = [delays[i] for i in range(len(tau)) if tau[i] <= tau_fast_max and math.isfinite(delays[i]) and corrs[i] >= min_corr]
    slow = [delays[i] for i in range(len(tau)) if tau[i] >= tau_slow_min and math.isfinite(delays[i]) and corrs[i] >= min_corr]
    if len(fast) == 0 or len(slow) == 0:
        return float("nan"), len(fast), len(slow)
    return float(np.nanmedian(slow) - np.nanmedian(fast)), len(fast), len(slow)


def _spectral_delay_delta(
    a: np.ndarray,
    b: np.ndarray,
    dt_days: float,
    min_coh: float,
    f_fast_min: float,
    f_slow_max: float,
) -> Tuple[float, int, int]:
    a0 = _fill_nans(a)
    b0 = _fill_nans(b)

    a0 = a0 - float(np.mean(a0))
    b0 = b0 - float(np.mean(b0))

    n = int(a0.size)
    # We care about very low frequencies (long timescales), so we prefer larger
    # segments to improve frequency resolution.
    if n < 128:
        return float("nan"), 0, 0
    nperseg = min(1024, n)
    if nperseg < 256:
        nperseg = max(128, (nperseg // 4) * 4)

    fs = 1.0 / dt_days
    f, pxy = signal.csd(a0, b0, fs=fs, nperseg=nperseg)
    _, pxx = signal.welch(a0, fs=fs, nperseg=nperseg)
    _, pyy = signal.welch(b0, fs=fs, nperseg=nperseg)

    with np.errstate(divide="ignore", invalid="ignore"):
        coh = (np.abs(pxy) ** 2) / (pxx * pyy)

    good = np.isfinite(f) & np.isfinite(coh) & np.isfinite(pxy) & (f > 0)
    if not np.any(good):
        return float("nan"), 0, 0

    f = f[good]
    pxy = pxy[good]
    coh = coh[good]

    # IMPORTANT: unwrap phase on the full, frequency-ordered grid.
    # Unwrapping only on discontiguous subsets (e.g. coherence-selected bins)
    # produces non-physical jumps.
    phase_full = np.unwrap(np.angle(pxy))
    dt_full = phase_full / (2.0 * math.pi * f)

    def _select_band(mask: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        # Apply coherence cut within the band; if too few points pass, fall back
        # to the highest-coherence bins *within the band*.
        if np.count_nonzero(mask) < 8:
            return np.array([], dtype=float), np.array([], dtype=float)
        band_idx = np.where(mask)[0]
        coh_band = coh[band_idx]
        keep = coh_band >= min_coh
        if np.count_nonzero(keep) < 8:
            k = min(32, coh_band.size)
            top = np.argsort(coh_band)[-k:]
            keep = np.zeros_like(coh_band, dtype=bool)
            keep[top] = True
        use_idx = band_idx[keep]
        return dt_full[use_idx], coh[use_idx]

    slow_mask = f <= float(f_slow_max)
    fast_mask = f >= float(f_fast_min)

    slow_dt, slow_w = _select_band(slow_mask)
    fast_dt, fast_w = _select_band(fast_mask)

    slow_dt = slow_dt[np.isfinite(slow_dt)]
    fast_dt = fast_dt[np.isfinite(fast_dt)]

    if slow_dt.size < 5 or fast_dt.size < 5:
        return float("nan"), int(fast_dt.size), int(slow_dt.size)

    # Robust: medians (less sensitive to outliers / residual phase issues)
    return float(np.median(slow_dt) - np.median(fast_dt)), int(fast_dt.size), int(slow_dt.size)


def _analyze_pair_multiscale(
    grid: np.ndarray,
    a: np.ndarray,
    b: np.ndarray,
    dt: float,
    lag_max: float,
    detrend_days: float,
    taus: List[float],
    dog_width: float,
    min_points: int,
    min_corr: float,
) -> Dict[str, object]:
    sigma_det = detrend_days / dt
    a_det, _ = _lowpass_detrend(a, sigma_det)
    b_det, _ = _lowpass_detrend(b, sigma_det)

    delays: List[float] = []
    corrs: List[float] = []
    npts: List[int] = []

    for tau in taus:
        s1 = max(1.0, (tau * (1.0 - dog_width)) / dt)
        s2 = max(s1 + 1.0, (tau * (1.0 + dog_width)) / dt)
        aa = _bandpass_dog(a_det, s1, s2)
        bb = _bandpass_dog(b_det, s1, s2)
        d, r, n = _estimate_delay_corr(aa, bb, dt, lag_max, min_points)
        if not math.isfinite(r) or r < min_corr:
            delays.append(float("nan"))
        else:
            delays.append(float(d))
        corrs.append(float(r) if math.isfinite(r) else float("nan"))
        npts.append(int(n))

    delta, n_fast, n_slow = _fit_delta_delay(
        taus,
        delays,
        corrs,
        tau_fast_max=40.0,
        tau_slow_min=80.0,
        min_corr=min_corr,
    )

    return {
        "detrend_days": float(detrend_days),
        "taus": [float(x) for x in taus],
        "delays": delays,
        "corrs": corrs,
        "npts": npts,
        "delta_slow_minus_fast": float(delta),
        "n_fast": int(n_fast),
        "n_slow": int(n_slow),
    }


def _simulate_once(
    grid: np.ndarray,
    a_grid: np.ndarray,
    b_grid: np.ndarray,
    a_err: np.ndarray,
    b_err: np.ndarray,
    dt: float,
    lag_max: float,
    base_delay: float,
    extra_slow_delay: float,
    slow_sigma_days: float,
    microlens_sigma_days: float,
    noise_scale: float,
    detrend_days_list: List[float],
    taus: List[float],
    dog_width: float,
    min_points: int,
    min_corr: float,
    min_coh: float,
    spec_fast_days: float,
    spec_slow_days: float,
) -> Dict[str, object]:
    sigma_ml = microlens_sigma_days / dt
    if sigma_ml < 1.0:
        sigma_ml = 1.0

    a0 = a_grid.copy()
    b0 = b_grid.copy()

    a_fill = _fill_nans(a0)
    b_fill = _fill_nans(b0)

    b_shift_back = _shift_series(b_fill, grid, -base_delay)
    ml_raw = b_fill - _shift_series(a_fill, grid, base_delay)
    ml = gaussian_filter1d(ml_raw, sigma=sigma_ml, mode="nearest")

    src = 0.5 * (a_fill + (b_shift_back - _shift_series(ml, grid, -base_delay)))

    b_inj = _inject_two_scale_warp(src, grid, base_delay, extra_slow_delay, slow_sigma_days, dt)
    b_inj = b_inj + ml

    mask_a = np.isfinite(a0)
    mask_b = np.isfinite(b0)
    a_sim = a_fill
    b_sim = b_inj
    a_sim[~mask_a] = np.nan
    b_sim[~mask_b] = np.nan

    rng = np.random.default_rng()

    ea = a_err.copy()
    eb = b_err.copy()
    ea = np.where(np.isfinite(ea), ea, np.nanmedian(ea[np.isfinite(ea)]) if np.any(np.isfinite(ea)) else 0.02)
    eb = np.where(np.isfinite(eb), eb, np.nanmedian(eb[np.isfinite(eb)]) if np.any(np.isfinite(eb)) else 0.02)

    a_sim = a_sim + rng.standard_normal(a_sim.size) * (noise_scale * ea)
    b_sim = b_sim + rng.standard_normal(b_sim.size) * (noise_scale * eb)

    out: Dict[str, object] = {
        "base_delay": float(base_delay),
        "extra_slow_delay": float(extra_slow_delay),
    }

    per_detrend: Dict[str, object] = {}
    for dd in detrend_days_list:
        res_ms = _analyze_pair_multiscale(
            grid,
            a_sim,
            b_sim,
            dt,
            lag_max,
            dd,
            taus,
            dog_width,
            min_points,
            min_corr,
        )

        dt_delta_spec, nf_fast, nf_slow = _spectral_delay_delta(
            _lowpass_detrend(a_sim, dd / dt)[0],
            _lowpass_detrend(b_sim, dd / dt)[0],
            dt,
            min_coh=float(min_coh),
            f_fast_min=1.0 / float(spec_fast_days),
            f_slow_max=1.0 / float(spec_slow_days),
        )

        res_ms["delta_spec_slow_minus_fast"] = float(dt_delta_spec)
        res_ms["n_fast_freq"] = int(nf_fast)
        res_ms["n_slow_freq"] = int(nf_slow)

        per_detrend[str(int(dd))] = res_ms

    out["by_detrend"] = per_detrend
    return out


def _summarize_runs(
    runs: List[Dict[str, object]],
    detrend_days_list: List[float],
) -> Dict[str, object]:
    summary: Dict[str, object] = {}
    for dd in detrend_days_list:
        key = str(int(dd))
        deltas = []
        deltas_spec = []
        for r in runs:
            by = r.get("by_detrend", {})
            if key not in by:
                continue
            v = by[key].get("delta_slow_minus_fast")
            vs = by[key].get("delta_spec_slow_minus_fast")
            if isinstance(v, (int, float)) and math.isfinite(float(v)):
                deltas.append(float(v))
            if isinstance(vs, (int, float)) and math.isfinite(float(vs)):
                deltas_spec.append(float(vs))

        def _stat(x: List[float]) -> Dict[str, float]:
            if len(x) == 0:
                return {"n": 0, "mean": float("nan"), "std": float("nan"), "median": float("nan")}
            arr = np.asarray(x, dtype=float)
            return {
                "n": int(arr.size),
                "mean": float(np.mean(arr)),
                "std": float(np.std(arr)),
                "median": float(np.median(arr)),
            }

        summary[key] = {
            "delta_time_domain": _stat(deltas),
            "delta_spectral": _stat(deltas_spec),
        }

    return summary


def _parse_float_list(s: str) -> List[float]:
    out: List[float] = []
    for part in s.split(","):
        part = part.strip()
        if not part:
            continue
        out.append(float(part))
    return out


def main() -> None:
    p = argparse.ArgumentParser(description="Step 3.1: COSMOGRAIL injection-recovery time-warp power study")
    p.add_argument("--rdb", default="data/cosmograil/HE0435_Bonvin2016.rdb")
    p.add_argument("--pair", default="A-B")
    p.add_argument("--dt", type=float, default=1.0)
    p.add_argument("--lag-max", type=float, default=50.0)
    p.add_argument("--detrend-days", default="100,200,400")
    p.add_argument("--taus", default="10,20,40,80,160")
    p.add_argument("--dog-width", type=float, default=0.5)
    p.add_argument("--min-points", type=int, default=80)
    p.add_argument("--min-corr", type=float, default=0.3)

    p.add_argument("--min-coh", type=float, default=0.2)
    p.add_argument("--spec-fast-days", type=float, default=40.0)
    p.add_argument("--spec-slow-days", type=float, default=120.0)

    p.add_argument("--slow-sigma-days", type=float, default=120.0)
    p.add_argument("--microlens-sigma-days", type=float, default=400.0)

    p.add_argument("--extra-slow-delays", default="0,0.5,1,2,5")
    p.add_argument("--n-mc", type=int, default=30)
    p.add_argument("--noise-scale", type=float, default=1.0)

    p.add_argument("--output", default="results/outputs/step_3_1_injection_recovery.json")

    args = p.parse_args()

    rdb = Path(args.rdb)
    if not rdb.exists():
        raise FileNotFoundError(rdb)

    labels, lcs = _read_rdb(rdb)

    pair = args.pair.strip()
    if "-" not in pair:
        raise ValueError("--pair must be like A-B")
    a_lab, b_lab = pair.split("-", 1)
    a_lab = a_lab.strip()
    b_lab = b_lab.strip()

    if a_lab not in lcs or b_lab not in lcs:
        raise ValueError(f"pair {pair} not in file; have labels {sorted(lcs.keys())}")

    lc_a = lcs[a_lab]
    lc_b = lcs[b_lab]

    dt = float(args.dt)
    grid = _uniform_grid(lc_a.t, lc_b.t, dt)
    if grid.size == 0:
        raise RuntimeError("empty overlap grid")

    a_grid = _interp_to_grid(lc_a.t, lc_a.y, grid)
    b_grid = _interp_to_grid(lc_b.t, lc_b.y, grid)
    a_err = _interp_to_grid(lc_a.t, lc_a.yerr, grid)
    b_err = _interp_to_grid(lc_b.t, lc_b.yerr, grid)

    detrend_days_list = _parse_float_list(args.detrend_days)
    taus = _parse_float_list(args.taus)
    extra_list = _parse_float_list(args.extra_slow_delays)

    a_det0, _ = _lowpass_detrend(a_grid, detrend_days_list[0] / dt)
    b_det0, _ = _lowpass_detrend(b_grid, detrend_days_list[0] / dt)
    base_delay, base_r, base_n = _estimate_delay_corr(a_det0, b_det0, dt, float(args.lag_max), int(args.min_points))

    _print(f"RDB: {rdb}")
    _print(f"Pair: {pair} | overlap grid N={grid.size} | dt={dt} d")
    _print(f"Estimated base delay (from detrend {detrend_days_list[0]} d): {base_delay:.2f} d (r={base_r:.3f}, n={base_n})")

    results: Dict[str, object] = {
        "rdb": str(rdb),
        "pair": pair,
        "dt": dt,
        "lag_max": float(args.lag_max),
        "taus": taus,
        "detrend_days": detrend_days_list,
        "slow_sigma_days": float(args.slow_sigma_days),
        "microlens_sigma_days": float(args.microlens_sigma_days),
        "min_points": int(args.min_points),
        "min_corr": float(args.min_corr),
        "dog_width": float(args.dog_width),
        "base_delay_est": float(base_delay),
        "n_mc": int(args.n_mc),
        "noise_scale": float(args.noise_scale),
        "runs": {},
        "analysis_date": datetime.now().isoformat(),
    }

    for extra in extra_list:
        _print(f"Running injection extra_slow_delay={extra} d ...")
        runs: List[Dict[str, object]] = []
        for _ in range(int(args.n_mc)):
            runs.append(
                _simulate_once(
                    grid,
                    a_grid,
                    b_grid,
                    a_err,
                    b_err,
                    dt,
                    float(args.lag_max),
                    float(base_delay),
                    float(extra),
                    float(args.slow_sigma_days),
                    float(args.microlens_sigma_days),
                    float(args.noise_scale),
                    detrend_days_list,
                    taus,
                    float(args.dog_width),
                    int(args.min_points),
                    float(args.min_corr),
                    float(args.min_coh),
                    float(args.spec_fast_days),
                    float(args.spec_slow_days),
                )
            )

        summary = _summarize_runs(runs, detrend_days_list)
        results["runs"][str(extra)] = {
            "extra_slow_delay": float(extra),
            "summary": summary,
        }

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(results, indent=2), encoding="utf-8")
    _print(f"Wrote: {out_path}")


if __name__ == "__main__":
    main()
