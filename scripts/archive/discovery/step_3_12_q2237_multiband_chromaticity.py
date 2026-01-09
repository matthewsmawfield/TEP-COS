#!/usr/bin/env python3
"""Step 3.12: Multi-band chromaticity test (Q2237+0305)

Goal
- Execute an explicit multi-band test of temporal shear chromaticity using a lens with
  public multi-band monitoring.

Dataset
- VizieR: J/A+A/637/A89 (Goicoechea+ 2020): QSO 2237+0305 (Einstein Cross) light curves
  in multiple bands.

Strategy
- Download the light-curve tables for g and r bands (optionally V and I).
- For each band, run the canonical Step 3.0 temporal-shear estimator.
- Compute ΔΓ between bands for the same image-pair (default: A-B, also report all pairs).

Outputs
- data/cosmograil/q2237_JAA637A89_{band}.csv
- results/outputs/step_3_12_q2237_multiband_chromaticity.json
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, Tuple

import numpy as np

from astroquery.vizier import Vizier

from step_3_0_cosmograil_temporal_shear import LightCurve, LensSystem, analyze_system


DATA_DIR = Path("data/cosmograil")
OUT_DIR = Path("results/outputs")
OUT_JSON = OUT_DIR / "step_3_12_q2237_multiband_chromaticity.json"

CAT_ID = "J/A+A/637/A89"

# Table mapping discovered via astroquery:
# table10: g band light curves
# table12: r band light curves
# table11: V band light curves
# table14: I band light curves
TABLES = {
    "g": "table10",
    "r": "table12",
    "V": "table11",
    "I": "table14",
}


def fetch_band_table(band: str):
    Vizier.ROW_LIMIT = -1
    tables = Vizier.get_catalogs(CAT_ID)

    for t in tables:
        name = t.meta.get("name", "")
        if name.endswith(TABLES[band]):
            return t
    raise RuntimeError(f"Could not find {TABLES[band]} in {CAT_ID}")


def table_to_system(table, system_id: str, band: str) -> LensSystem:
    # Required columns: MJD, mA, e_mA, mB, e_mB, mC, e_mC, mD, e_mD
    t = np.array(table["MJD"], dtype=float)

    lcs = {}
    for img in ["A", "B", "C", "D"]:
        mcol = f"m{img}"
        ecol = f"e_m{img}"
        if mcol not in table.colnames or ecol not in table.colnames:
            continue
        mag = np.array(table[mcol], dtype=float)
        err = np.array(table[ecol], dtype=float)
        lcs[img] = LightCurve(label=img, t=t, mag=mag, magerr=err)

    return LensSystem(system_id=system_id, light_curves=lcs, band=band)


def compare_bands(res1: Dict, res2: Dict) -> Dict:
    pairs = set(res1["pairs"].keys()) & set(res2["pairs"].keys())
    out = {}

    for pair in sorted(pairs):
        g1 = res1["pairs"][pair]["gamma"]
        g2 = res2["pairs"][pair]["gamma"]

        v1 = float(g1["value"]) if np.isfinite(g1["value"]) else np.nan
        e1 = float(g1["uncertainty"]) if np.isfinite(g1["uncertainty"]) else np.nan
        v2 = float(g2["value"]) if np.isfinite(g2["value"]) else np.nan
        e2 = float(g2["uncertainty"]) if np.isfinite(g2["uncertainty"]) else np.nan

        if not (np.isfinite(v1) and np.isfinite(v2) and np.isfinite(e1) and np.isfinite(e2) and e1 > 0 and e2 > 0):
            continue

        delta = v1 - v2
        sigma = float(np.sqrt(e1 * e1 + e2 * e2))
        sig = float(abs(delta) / sigma) if sigma > 0 else np.nan

        out[pair] = {
            "delta_gamma": float(delta),
            "delta_uncertainty": float(sigma),
            "delta_sigma": float(sig),
            "band1_gamma": float(v1),
            "band1_uncertainty": float(e1),
            "band2_gamma": float(v2),
            "band2_uncertainty": float(e2),
        }

    return out


def main() -> None:
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # Download and persist band tables
    band_tables = {}
    for band in ["g", "r", "V", "I"]:
        t = fetch_band_table(band)
        band_tables[band] = t
        out_csv = DATA_DIR / f"q2237_JAA637A89_{band}.csv"
        t.write(out_csv, format="csv", overwrite=True)

    # Run canonical estimator per band
    results = {}
    for band in ["g", "r", "V", "I"]:
        sys = table_to_system(band_tables[band], system_id=f"Q2237_{band}", band=band)
        # Q2237 has short physical delays; keep lag window modest.
        res = analyze_system(sys, detrend_window=200.0, tau_values=[5, 10, 20, 40, 80, 160, 320], lag_range=(-50, 50))
        results[band] = res

    # Primary chromaticity comparisons
    comparisons = {
        "g_minus_r": compare_bands(results["g"], results["r"]),
        "V_minus_I": compare_bands(results["V"], results["I"]),
    }

    out = {
        "catalog": CAT_ID,
        "system": "Q2237+0305",
        "bands": {
            b: {
                "n_images": results[b]["n_images"],
                "pairs": results[b]["pairs"],
                "tau_values": results[b]["tau_values"],
                "detrend_window_days": results[b]["detrend_window_days"],
            }
            for b in results
        },
        "comparisons": comparisons,
    }

    with open(OUT_JSON, "w") as f:
        json.dump(out, f, indent=2)

    # Console summary: A-B if present
    for key in ["g_minus_r", "V_minus_I"]:
        comp = comparisons[key]
        if "A-B" in comp:
            d = comp["A-B"]
            print(f"{key} A-B: ΔΓ = {d['delta_gamma']:.2f} ± {d['delta_uncertainty']:.2f} (σ={d['delta_sigma']:.2f})")

    print(f"Saved: {OUT_JSON}")


if __name__ == "__main__":
    main()
