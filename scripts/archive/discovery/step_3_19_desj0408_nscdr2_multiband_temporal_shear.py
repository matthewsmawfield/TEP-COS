#!/usr/bin/env python3
"""Step 3.19: DESJ0408 NSC DR2 multi-band temporal shear + chromaticity

This step consumes the per-band light curves produced by:
  scripts/steps/step_3_18_desj0408_nscdr2_selfextract.py

It then:
- builds LensSystem/LightCurve objects compatible with step_3_0_cosmograil_temporal_shear
- runs analyze_system per band
- computes chromaticity diagnostics: ΔΓ between bands for each image pair

Outputs
- results/outputs/step_3_19_desj0408_nscdr2_multiband_temporal_shear.json

Notes
- NSC DR2 provides calibrated magnitudes and magnitude errors; we use those directly.
- We only keep rows with finite mjd, mag, and mag_err>0.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Dict, List

import numpy as np


# Allow running as a standalone script by ensuring repo root is on sys.path.
_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from scripts.steps.step_3_0_cosmograil_temporal_shear import LightCurve, LensSystem, analyze_system


def load_band_csv(csv_path: Path) -> Dict[str, Dict[str, np.ndarray]]:
    data = np.genfromtxt(
        csv_path,
        delimiter=",",
        names=True,
        dtype=None,
        encoding="utf-8",
    )

    if data.size == 0:
        raise RuntimeError(f"Empty CSV: {csv_path}")

    if data.shape == ():
        data = np.array([data])

    cols = set(data.dtype.names or [])
    required = {"mjd", "image", "mag", "mag_err"}
    if not required.issubset(cols):
        raise RuntimeError(f"CSV missing required columns {required} in {csv_path}")

    out: Dict[str, Dict[str, List[float]]] = {}
    for row in data:
        img = str(row["image"]).strip()
        mjd = float(row["mjd"]) if row["mjd"] is not None else np.nan
        mag = float(row["mag"]) if row["mag"] is not None else np.nan
        magerr = float(row["mag_err"]) if row["mag_err"] is not None else np.nan

        if not (np.isfinite(mjd) and np.isfinite(mag) and np.isfinite(magerr) and magerr > 0):
            continue

        out.setdefault(img, {"t": [], "mag": [], "magerr": []})
        out[img]["t"].append(mjd)
        out[img]["mag"].append(mag)
        out[img]["magerr"].append(magerr)

    return {k: {kk: np.array(vv, dtype=float) for kk, vv in d.items()} for k, d in out.items()}


def build_system(system_id: str, band: str, per_image: Dict[str, Dict[str, np.ndarray]], min_epochs: int) -> LensSystem:
    lcs: Dict[str, LightCurve] = {}
    for label, arrs in per_image.items():
        if arrs["t"].size < min_epochs:
            continue
        lcs[label] = LightCurve(label=label, t=arrs["t"], mag=arrs["mag"], magerr=arrs["magerr"])

    if len(lcs) < 2:
        raise RuntimeError(f"Not enough images with >={min_epochs} usable epochs for band {band} (got {list(lcs.keys())})")

    return LensSystem(system_id=system_id, light_curves=lcs, band=band)


def compute_delta_gamma(per_band_results: Dict[str, Dict]) -> Dict[str, Dict[str, float]]:
    gamma_by_band: Dict[str, Dict[str, float]] = {}
    for band, res in per_band_results.items():
        gamma_by_band[band] = {}
        for pair_key, pair_data in res["pairs"].items():
            gamma_by_band[band][pair_key] = float(pair_data["gamma"]["value"])

    bands = sorted(gamma_by_band.keys())
    if len(bands) < 2:
        return {}

    pairs = set()
    for b in bands:
        pairs |= set(gamma_by_band[b].keys())

    out: Dict[str, Dict[str, float]] = {}
    for pair in sorted(pairs):
        out[pair] = {}
        for i in range(len(bands)):
            for j in range(i + 1, len(bands)):
                b1, b2 = bands[i], bands[j]
                g1 = gamma_by_band[b1].get(pair, np.nan)
                g2 = gamma_by_band[b2].get(pair, np.nan)
                key = f"{b1}-{b2}"
                out[pair][key] = float(g1 - g2) if (np.isfinite(g1) and np.isfinite(g2)) else float("nan")

    return out


def main():
    parser = argparse.ArgumentParser(description="DESJ0408 NSC DR2 multi-band temporal shear + chromaticity")
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=Path("data/selfextract/desj0408_nscdr2"),
        help="Directory containing lightcurves_desj0408_nscdr2_<band>.csv",
    )
    parser.add_argument(
        "--bands",
        type=str,
        default="g,r,i,z",
        help="Comma-separated bands (filenames use these exactly)",
    )
    parser.add_argument(
        "--min-epochs",
        type=int,
        default=10,
        help="Minimum usable epochs per image to include that image in a band analysis.",
    )
    parser.add_argument(
        "--tau-values",
        type=str,
        default="",
        help="Comma-separated tau values in days (overrides estimator default if provided)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("results/outputs/step_3_19_desj0408_nscdr2_multiband_temporal_shear.json"),
        help="Output JSON path",
    )
    parser.add_argument("--detrend-window", type=float, default=200.0)
    parser.add_argument("--lag-min", type=float, default=-150.0, help="Minimum lag (days) for correlation scan")
    parser.add_argument("--lag-max", type=float, default=150.0, help="Maximum lag (days) for correlation scan")
    parser.add_argument("--mode-lock-window", type=float, default=50.0, help="Half-width (days) around broadband delay")
    parser.add_argument(
        "--min-variance-fraction",
        type=float,
        default=0.02,
        help="Minimum fraction of original variance preserved by bandpass filter",
    )
    args = parser.parse_args()

    bands = [b.strip() for b in args.bands.split(",") if b.strip()]
    tau_values = None
    if args.tau_values.strip():
        tau_values = [float(x.strip()) for x in args.tau_values.split(",") if x.strip()]

    per_band_results: Dict[str, Dict] = {}
    per_band_qc: Dict[str, Dict] = {}
    per_band_skipped: Dict[str, str] = {}

    for band in bands:
        csv_path = args.input_dir / f"lightcurves_desj0408_nscdr2_{band}.csv"
        if not csv_path.exists():
            per_band_skipped[band] = f"Missing expected input CSV: {csv_path}"
            continue

        per_image = load_band_csv(csv_path)
        try:
            system = build_system("DESJ0408_NSCDR2", band, per_image, min_epochs=args.min_epochs)
        except Exception as e:
            per_band_skipped[band] = str(e)
            continue

        per_band_qc[band] = {
            "n_images": system.n_images,
            "image_labels": system.image_labels,
            "n_epochs": {k: int(v.n_epochs) for k, v in system.light_curves.items()},
            "baseline_days": {k: float(v.baseline_days) for k, v in system.light_curves.items()},
        }

        res = analyze_system(
            system,
            detrend_window=args.detrend_window,
            tau_values=tau_values,
            lag_range=(args.lag_min, args.lag_max),
            mode_lock_window=args.mode_lock_window,
            min_variance_fraction=args.min_variance_fraction,
        )
        per_band_results[band] = res

    if not per_band_results:
        raise RuntimeError(f"No bands had sufficient usable data to analyze. Skipped: {per_band_skipped}")

    delta_gamma = compute_delta_gamma(per_band_results)

    out = {
        "system_id": "DESJ0408_NSCDR2",
        "bands": bands,
        "bands_analyzed": sorted(per_band_results.keys()),
        "bands_skipped": per_band_skipped,
        "lag_range": [args.lag_min, args.lag_max],
        "tau_values": tau_values,
        "mode_lock_window_days": args.mode_lock_window,
        "min_variance_fraction": args.min_variance_fraction,
        "qc": per_band_qc,
        "per_band": per_band_results,
        "delta_gamma": delta_gamma,
    }

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w") as f:
        json.dump(out, f, indent=2, default=lambda x: None if (isinstance(x, float) and not np.isfinite(x)) else x)


if __name__ == "__main__":
    main()
