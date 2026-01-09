#!/usr/bin/env python3

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import argparse
import csv
import json
import math
from dataclasses import asdict, dataclass
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np
from scipy.interpolate import interp1d
from scipy.ndimage import gaussian_filter1d

from scripts.utils.logger import TEPLogger, print_status, set_step_logger


@dataclass
class LightCurve:
    t: np.ndarray
    y: np.ndarray
    yerr: np.ndarray


def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Step 2.1: COSMOGRAIL multi-band time-delay consistency analysis")

    p.add_argument("--data-dir", default=str(PROJECT_ROOT / "data" / "cosmograil"))
    p.add_argument("--manifest", default="")
    p.add_argument("--output-dir", default=str(PROJECT_ROOT / "results" / "outputs"))

    p.add_argument("--grid-dt", type=float, default=1.0)
    p.add_argument("--lag-max", type=float, default=250.0)
    p.add_argument("--detrend-days", type=float, default=200.0)
    p.add_argument("--taus", default="5,10,20,40,80")

    p.add_argument("--n-mc", type=int, default=200)
    p.add_argument("--min-points", type=int, default=50)

    return p.parse_args()


def _read_light_curve_csv(path: Path) -> LightCurve:
    with open(path, "r", encoding="utf-8") as f:
        reader = csv.reader(f)
        header: Optional[List[str]] = None
        rows: List[List[str]] = []
        for row in reader:
            if not row:
                continue
            if header is None:
                header = [c.strip() for c in row]
                continue
            rows.append(row)

    if header is None:
        raise ValueError(f"Empty light curve: {path}")

    cols = {name: i for i, name in enumerate(header)}

    def _idx(name: str) -> int:
        if name in cols:
            return cols[name]
        lower = {k.lower(): v for k, v in cols.items()}
        if name.lower() in lower:
            return lower[name.lower()]
        raise KeyError(f"Missing column '{name}' in {path.name} (have: {header})")

    ti, yi, ei = _idx("t"), _idx("y"), _idx("yerr")

    t_list: List[float] = []
    y_list: List[float] = []
    e_list: List[float] = []

    for r in rows:
        try:
            t = float(r[ti])
            y = float(r[yi])
            e = float(r[ei])
        except Exception:
            continue
        if not (math.isfinite(t) and math.isfinite(y) and math.isfinite(e)):
            continue
        t_list.append(t)
        y_list.append(y)
        e_list.append(e)

    if len(t_list) < 3:
        raise ValueError(f"Too few valid rows in {path}")

    t = np.asarray(t_list, dtype=float)
    y = np.asarray(y_list, dtype=float)
    yerr = np.asarray(e_list, dtype=float)

    order = np.argsort(t)
    return LightCurve(t=t[order], y=y[order], yerr=yerr[order])


def _build_uniform_grid(lc_a: LightCurve, lc_b: LightCurve, dt: float) -> np.ndarray:
    t0 = max(float(np.min(lc_a.t)), float(np.min(lc_b.t)))
    t1 = min(float(np.max(lc_a.t)), float(np.max(lc_b.t)))
    if not (t1 > t0):
        return np.array([], dtype=float)
    n = int(math.floor((t1 - t0) / dt)) + 1
    if n <= 2:
        return np.array([], dtype=float)
    return t0 + dt * np.arange(n, dtype=float)


def _interp_to_grid(lc: LightCurve, grid_t: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    f = interp1d(lc.t, lc.y, kind="linear", bounds_error=False, fill_value=np.nan)
    y = f(grid_t)
    fe = interp1d(lc.t, lc.yerr, kind="linear", bounds_error=False, fill_value=np.nan)
    yerr = fe(grid_t)
    return y, yerr


def _detrend_lowpass(y: np.ndarray, sigma_samples: float) -> np.ndarray:
    if sigma_samples <= 0:
        return y
    good = np.isfinite(y)
    if not np.any(good):
        return y

    y_filled = y.copy()
    med = float(np.nanmedian(y_filled[good]))
    y_filled[~good] = med

    low = gaussian_filter1d(y_filled, sigma=sigma_samples, mode="nearest")
    out = y - low
    out[~good] = np.nan
    return out


def _bandpass_dog(y: np.ndarray, sigma1: float, sigma2: float) -> np.ndarray:
    if sigma1 <= 0 or sigma2 <= 0:
        return y
    good = np.isfinite(y)
    if not np.any(good):
        return y
    y_filled = y.copy()
    med = float(np.nanmedian(y_filled[good]))
    y_filled[~good] = med

    g1 = gaussian_filter1d(y_filled, sigma=sigma1, mode="nearest")
    g2 = gaussian_filter1d(y_filled, sigma=sigma2, mode="nearest")
    out = g1 - g2
    out[~good] = np.nan
    return out


def _corrcoef_nan(a: np.ndarray, b: np.ndarray) -> Tuple[float, int]:
    good = np.isfinite(a) & np.isfinite(b)
    n = int(np.count_nonzero(good))
    if n < 3:
        return float("nan"), n
    aa = a[good]
    bb = b[good]
    sa = float(np.std(aa))
    sb = float(np.std(bb))
    if sa == 0.0 or sb == 0.0:
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
    n_by_k: Dict[int, int] = {}

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
        n_by_k[k] = n
        if n < min_points or not math.isfinite(r):
            continue
        if r > best_r:
            best_r = r
            best_k = k
            best_n = n

    if best_r < -0.5:
        return best_k * dt, best_r, best_n

    k0 = best_k
    r0 = r_by_k.get(k0, float("nan"))
    r1 = r_by_k.get(k0 - 1, float("nan"))
    r2 = r_by_k.get(k0 + 1, float("nan"))

    delta = 0.0
    if math.isfinite(r0) and math.isfinite(r1) and math.isfinite(r2):
        denom = (r1 - 2.0 * r0 + r2)
        if denom != 0.0:
            delta = 0.5 * (r1 - r2) / denom
            delta = float(np.clip(delta, -1.0, 1.0))

    return (k0 + delta) * dt, best_r, best_n


def _parse_manifest(path: Path) -> Dict[str, Dict[str, Dict[str, Path]]]:
    obj = json.loads(path.read_text(encoding="utf-8"))
    systems = obj.get("systems", [])
    out: Dict[str, Dict[str, Dict[str, Path]]] = {}

    for s in systems:
        sid = str(s.get("system_id", "")).strip()
        if not sid:
            continue
        bands = s.get("bands", {})
        if not isinstance(bands, dict):
            continue
        out[sid] = {}
        for band, images in bands.items():
            if not isinstance(images, dict):
                continue
            out[sid][str(band)] = {}
            for image_label, rel in images.items():
                p = Path(rel)
                out[sid][str(band)][str(image_label)] = p

    return out


def _auto_discover(data_dir: Path) -> Dict[str, Dict[str, Dict[str, Path]]]:
    out: Dict[str, Dict[str, Dict[str, Path]]] = {}
    for path in sorted(data_dir.rglob("*.csv")):
        name = path.name
        if "__" not in name:
            continue
        stem = name[:-4]
        parts = stem.split("__")
        if len(parts) != 3:
            continue
        system_id, image_label, band = parts
        out.setdefault(system_id, {}).setdefault(band, {})[image_label] = path
    return out


def _iter_image_pairs(images: List[str]) -> Iterable[Tuple[str, str]]:
    for i in range(len(images)):
        for j in range(i + 1, len(images)):
            yield images[i], images[j]


def main() -> None:
    args = _parse_args()

    log = TEPLogger("step_2_1_cosmograil", log_file_path=PROJECT_ROOT / "logs" / "step_2_1_cosmograil_time_delay.log")
    set_step_logger(log)

    data_dir = Path(args.data_dir)
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    taus = [float(x.strip()) for x in args.taus.split(",") if x.strip()]
    if not taus:
        raise ValueError("No taus provided")

    if args.manifest:
        manifest_path = Path(args.manifest)
    else:
        manifest_path = data_dir / "manifest.json"

    mapping: Dict[str, Dict[str, Dict[str, Path]]]
    if manifest_path.exists():
        print_status(f"Using manifest: {manifest_path}", "PROCESS")
        mapping = _parse_manifest(manifest_path)
        for sid in list(mapping.keys()):
            for band in list(mapping[sid].keys()):
                for img in list(mapping[sid][band].keys()):
                    p = mapping[sid][band][img]
                    if not p.is_absolute():
                        mapping[sid][band][img] = data_dir / p
    else:
        print_status(f"Manifest not found, auto-discovering CSVs under: {data_dir}", "PROCESS")
        if not data_dir.exists():
            raise FileNotFoundError(f"Missing data directory: {data_dir}")
        mapping = _auto_discover(data_dir)

    if not mapping:
        raise RuntimeError("No COSMOGRAIL light curves discovered")

    per_band_rows: List[Dict[str, object]] = []
    multiscale_rows: List[Dict[str, object]] = []
    gamma_rows: List[Dict[str, object]] = []

    cfg = {
        "grid_dt": args.grid_dt,
        "lag_max": args.lag_max,
        "detrend_days": args.detrend_days,
        "taus": taus,
        "n_mc": args.n_mc,
        "min_points": args.min_points,
    }

    detrend_sigma = args.detrend_days / args.grid_dt

    for system_id, bands in mapping.items():
        print_status(f"System: {system_id}", "PROCESS")
        for band, images in bands.items():
            image_labels = sorted(images.keys())
            if len(image_labels) < 2:
                continue

            lc_by_img: Dict[str, LightCurve] = {}
            for img in image_labels:
                p = images[img]
                if not p.exists():
                    print_status(f"Missing light curve: {p}", "WARNING")
                    continue
                try:
                    lc_by_img[img] = _read_light_curve_csv(p)
                except Exception as e:
                    print_status(f"Failed to read {p.name}: {e}", "WARNING")

            labels = sorted(lc_by_img.keys())
            if len(labels) < 2:
                continue

            for img_a, img_b in _iter_image_pairs(labels):
                lc_a = lc_by_img[img_a]
                lc_b = lc_by_img[img_b]

                grid_t = _build_uniform_grid(lc_a, lc_b, args.grid_dt)
                if grid_t.size == 0:
                    continue

                y_a, _ = _interp_to_grid(lc_a, grid_t)
                y_b, _ = _interp_to_grid(lc_b, grid_t)

                y_a_d = _detrend_lowpass(y_a, detrend_sigma)
                y_b_d = _detrend_lowpass(y_b, detrend_sigma)

                delay, peak_r, n_eff = _estimate_delay_corr(
                    y_a_d,
                    y_b_d,
                    dt=args.grid_dt,
                    lag_max_days=args.lag_max,
                    min_points=args.min_points,
                )

                delays_mc: List[float] = []
                if args.n_mc > 0:
                    fa = interp1d(lc_a.t, lc_a.y, kind="linear", bounds_error=False, fill_value=np.nan)
                    fb = interp1d(lc_b.t, lc_b.y, kind="linear", bounds_error=False, fill_value=np.nan)
                    fea = interp1d(lc_a.t, lc_a.yerr, kind="linear", bounds_error=False, fill_value=np.nan)
                    feb = interp1d(lc_b.t, lc_b.yerr, kind="linear", bounds_error=False, fill_value=np.nan)

                    base_a = fa(grid_t)
                    base_b = fb(grid_t)
                    err_a = fea(grid_t)
                    err_b = feb(grid_t)

                    rng = np.random.default_rng(42)
                    for _ in range(args.n_mc):
                        a_s = base_a + rng.normal(0.0, err_a)
                        b_s = base_b + rng.normal(0.0, err_b)
                        a_s = _detrend_lowpass(a_s, detrend_sigma)
                        b_s = _detrend_lowpass(b_s, detrend_sigma)
                        d_s, _, n_s = _estimate_delay_corr(
                            a_s,
                            b_s,
                            dt=args.grid_dt,
                            lag_max_days=args.lag_max,
                            min_points=args.min_points,
                        )
                        if n_s >= args.min_points and math.isfinite(d_s):
                            delays_mc.append(float(d_s))

                delay_sigma = float(np.std(np.asarray(delays_mc), ddof=1)) if len(delays_mc) > 2 else float("nan")

                per_band_rows.append(
                    {
                        "system_id": system_id,
                        "band": band,
                        "img_a": img_a,
                        "img_b": img_b,
                        "delay_days": delay,
                        "delay_sigma_days": delay_sigma,
                        "peak_r": peak_r,
                        "n_eff": n_eff,
                        "grid_dt": args.grid_dt,
                        "detrend_days": args.detrend_days,
                    }
                )

                tau_delays: List[Tuple[float, float]] = []
                for tau in taus:
                    sigma1 = max(0.5, 0.5 * tau / args.grid_dt)
                    sigma2 = max(sigma1 + 0.5, 2.0 * sigma1)
                    a_bp = _bandpass_dog(y_a_d, sigma1=sigma1, sigma2=sigma2)
                    b_bp = _bandpass_dog(y_b_d, sigma1=sigma1, sigma2=sigma2)
                    d_tau, r_tau, n_tau = _estimate_delay_corr(
                        a_bp,
                        b_bp,
                        dt=args.grid_dt,
                        lag_max_days=args.lag_max,
                        min_points=args.min_points,
                    )
                    if n_tau >= args.min_points and math.isfinite(d_tau):
                        tau_delays.append((tau, float(d_tau)))
                    multiscale_rows.append(
                        {
                            "system_id": system_id,
                            "band": band,
                            "img_a": img_a,
                            "img_b": img_b,
                            "tau_days": tau,
                            "delay_days": d_tau,
                            "peak_r": r_tau,
                            "n_eff": n_tau,
                        }
                    )

                if len(tau_delays) >= 2:
                    xs = np.log(np.asarray([t for t, _ in tau_delays], dtype=float))
                    ys = np.asarray([d for _, d in tau_delays], dtype=float)
                    x0 = float(np.mean(xs))
                    y0 = float(np.mean(ys))
                    num = float(np.sum((xs - x0) * (ys - y0)))
                    den = float(np.sum((xs - x0) ** 2))
                    gamma = float(num / den) if den != 0.0 else float("nan")
                else:
                    gamma = float("nan")

                gamma_rows.append(
                    {
                        "system_id": system_id,
                        "band": band,
                        "img_a": img_a,
                        "img_b": img_b,
                        "gamma": gamma,
                        "n_tau": len(tau_delays),
                    }
                )

    out_per_band = out_dir / "step_2_1_delay_per_band.csv"
    out_multiscale = out_dir / "step_2_1_delay_multiscale.csv"
    out_gamma = out_dir / "step_2_1_delay_gamma.csv"
    out_summary = out_dir / "step_2_1_summary.json"

    def _write_csv(path: Path, rows: List[Dict[str, object]]) -> None:
        if not rows:
            path.write_text("", encoding="utf-8")
            return
        keys = list(rows[0].keys())
        with open(path, "w", encoding="utf-8", newline="") as f:
            w = csv.DictWriter(f, fieldnames=keys)
            w.writeheader()
            w.writerows(rows)

    _write_csv(out_per_band, per_band_rows)
    _write_csv(out_multiscale, multiscale_rows)
    _write_csv(out_gamma, gamma_rows)

    summary = {
        "config": cfg,
        "n_systems": len(mapping),
        "n_pairs_per_band": len(per_band_rows),
        "outputs": {
            "delay_per_band_csv": str(out_per_band),
            "delay_multiscale_csv": str(out_multiscale),
            "delay_gamma_csv": str(out_gamma),
        },
    }
    out_summary.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print_status(f"Wrote: {out_per_band}", "SUCCESS")
    print_status(f"Wrote: {out_multiscale}", "SUCCESS")
    print_status(f"Wrote: {out_gamma}", "SUCCESS")
    print_status(f"Wrote: {out_summary}", "SUCCESS")


if __name__ == "__main__":
    main()
