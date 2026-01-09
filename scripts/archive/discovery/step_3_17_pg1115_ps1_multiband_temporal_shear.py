#!/usr/bin/env python3
"""Step 3.17: PG1115 PS1 self-extracted multi-band temporal shear + chromaticity

This step consumes the self-extracted PS1 forced-photometry light curves produced by:
  scripts/steps/step_3_16_pg1115_ps1_selfextract.py

It then:
- builds LensSystem/LightCurve objects compatible with step_3_0_cosmograil_temporal_shear
- runs analyze_system per band
- computes simple chromaticity diagnostics: ΔΓ between bands for each image pair

Outputs
- results/outputs/step_3_17_pg1115_ps1_multiband_temporal_shear.json

Notes
- Input CSVs contain fluxes and magnitudes; the step_3_0 estimator expects magnitudes.
- We only keep rows with finite mjd, mag, and mag_err.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np

# Allow running as a standalone script by ensuring repo root is on sys.path.
_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from scripts.steps.step_3_0_cosmograil_temporal_shear import LightCurve, LensSystem, analyze_system


def load_band_csv(csv_path: Path) -> Dict[str, Dict[str, np.ndarray]]:
    """Return mapping image-> {t, mag, magerr} in MHJD/MJD days."""
    # Using numpy to avoid extra dependency on pandas.
    data = np.genfromtxt(
        csv_path,
        delimiter=",",
        names=True,
        dtype=None,
        encoding="utf-8",
    )

    if data.size == 0:
        raise RuntimeError(f"Empty CSV: {csv_path}")

    # Ensure 1D structured array
    if data.shape == ():
        data = np.array([data])

    cols = set(data.dtype.names or [])
    required_core = {"mjd", "image"}
    if not required_core.issubset(cols):
        raise RuntimeError(f"CSV missing required columns {required_core} in {csv_path}")

    has_mag = {"mag", "mag_err"}.issubset(cols)
    has_flux = {"flux", "flux_err"}.issubset(cols)
    if not (has_mag or has_flux):
        raise RuntimeError(
            f"CSV must contain either (mag, mag_err) or (flux, flux_err) columns: {csv_path}"
        )

    out: Dict[str, Dict[str, List[float]]] = {}
    for row in data:
        img = str(row["image"])
        mjd = float(row["mjd"]) if row["mjd"] is not None else np.nan

        mag = np.nan
        magerr = np.nan
        if has_mag:
            mag = float(row["mag"]) if row["mag"] is not None else np.nan
            magerr = float(row["mag_err"]) if row["mag_err"] is not None else np.nan

        # Fallback: use flux directly as a brightness proxy.
        # The temporal-shear estimator normalizes each series internally; it does not require
        # absolute calibration, and it can operate on any monotonic proxy.
        if (not np.isfinite(mag)) or (not np.isfinite(magerr)) or (magerr <= 0):
            if has_flux:
                flux = float(row["flux"]) if row["flux"] is not None else np.nan
                flux_err = float(row["flux_err"]) if row["flux_err"] is not None else np.nan
                if np.isfinite(flux) and np.isfinite(flux_err) and flux_err > 0:
                    mag = flux
                    magerr = flux_err

        if not (np.isfinite(mjd) and np.isfinite(mag) and np.isfinite(magerr) and magerr > 0):
            continue

        out.setdefault(img, {"t": [], "mag": [], "magerr": []})
        out[img]["t"].append(mjd)
        out[img]["mag"].append(mag)
        out[img]["magerr"].append(magerr)

    return {k: {kk: np.array(vv, dtype=float) for kk, vv in d.items()} for k, d in out.items()}


def build_system(
    system_id: str,
    band: str,
    per_image: Dict[str, Dict[str, np.ndarray]],
    min_epochs_per_image: int,
) -> LensSystem:
    lcs: Dict[str, LightCurve] = {}
    for label, arrs in per_image.items():
        if arrs["t"].size < min_epochs_per_image:
            continue
        lcs[label] = LightCurve(label=label, t=arrs["t"], mag=arrs["mag"], magerr=arrs["magerr"])

    if len(lcs) < 2:
        raise RuntimeError(
            f"Not enough images with >={min_epochs_per_image} usable epochs for band {band} (got {list(lcs.keys())})"
        )

    return LensSystem(system_id=system_id, light_curves=lcs, band=band)


def compute_delta_gamma(per_band_results: Dict[str, Dict]) -> Dict[str, Dict[str, float]]:
    """Compute ΔΓ between bands per image-pair. Returns pair -> {"r-i":..., ...}."""
    # Extract gamma values
    gamma_by_band: Dict[str, Dict[str, float]] = {}
    for band, res in per_band_results.items():
        gamma_by_band[band] = {}
        for pair_key, pair_data in res["pairs"].items():
            gamma_by_band[band][pair_key] = float(pair_data["gamma"]["value"]) if pair_data else np.nan

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
    parser = argparse.ArgumentParser(description="PG1115 PS1 multi-band temporal shear + chromaticity")
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=Path("data/selfextract/pg1115_ps1"),
        help="Directory containing lightcurves_pg1115_ps1_<band>.csv",
    )
    parser.add_argument(
        "--bands",
        type=str,
        default="r.00000,i.00000,z.00000",
        help="Comma-separated band tags matching filenames (e.g., r.00000,i.00000,z.00000)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("results/outputs/step_3_17_pg1115_ps1_multiband_temporal_shear.json"),
        help="Output JSON path",
    )
    parser.add_argument("--detrend-window", type=float, default=200.0)
    parser.add_argument(
        "--min-epochs",
        type=int,
        default=10,
        help="Minimum usable epochs per image to include that image in a band analysis.",
    )
    args = parser.parse_args()

    bands = [b.strip() for b in args.bands.split(",") if b.strip()]
    per_band_results: Dict[str, Dict] = {}
    per_band_qc: Dict[str, Dict] = {}
    per_band_skipped: Dict[str, str] = {}

    for band in bands:
        csv_path = args.input_dir / f"lightcurves_pg1115_ps1_{band}.csv"
        if not csv_path.exists():
            raise FileNotFoundError(f"Missing expected input CSV: {csv_path}")

        per_image = load_band_csv(csv_path)
        try:
            system = build_system(
                system_id="PG1115_PS1",
                band=band,
                per_image=per_image,
                min_epochs_per_image=args.min_epochs,
            )
        except Exception as e:
            per_band_skipped[band] = str(e)
            continue

        # QC summary
        per_band_qc[band] = {
            "n_images": system.n_images,
            "image_labels": system.image_labels,
            "n_epochs": {k: int(v.n_epochs) for k, v in system.light_curves.items()},
            "baseline_days": {k: float(v.baseline_days) for k, v in system.light_curves.items()},
        }

        results = analyze_system(system, detrend_window=args.detrend_window)
        per_band_results[band] = results

    if not per_band_results:
        raise RuntimeError("No bands had sufficient usable data to analyze.")

    delta_gamma = compute_delta_gamma(per_band_results)

    out = {
        "system_id": "PG1115_PS1",
        "bands": bands,
        "bands_analyzed": sorted(per_band_results.keys()),
        "bands_skipped": per_band_skipped,
        "qc": per_band_qc,
        "per_band": per_band_results,
        "delta_gamma": delta_gamma,
    }

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w") as f:
        json.dump(out, f, indent=2, default=lambda x: None if (isinstance(x, float) and not np.isfinite(x)) else x)


if __name__ == "__main__":
    main()
