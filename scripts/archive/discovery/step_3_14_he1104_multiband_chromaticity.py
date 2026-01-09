#!/usr/bin/env python3
"""Step 3.14: Multi-band chromaticity test (HE1104-1805)

Goal
- Extend multi-band chromaticity validation beyond Q2237+0305 and HE0435-1223.
- Use a public, image-resolved multi-band light curve for HE1104-1805.

Dataset
- VizieR: J/ApJ/798/95 (Blackburne+ 2015): HE1104-1805 BVRIJ light curves.
  Columns: HJD, A, e_A, B, e_B, Filt.

Method
- Split into per-filter light curves.
- Run canonical temporal-shear estimator (Step 3.0) for the A-B pair in each band.
- Compute ΔΓ between bands for A-B with propagated uncertainty.

Outputs
- data/cosmograil/he1104_JApJ798_95_{band}.csv
- results/outputs/step_3_14_he1104_multiband_chromaticity.json
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, List

import numpy as np

from astroquery.vizier import Vizier

from step_3_0_cosmograil_temporal_shear import LightCurve, LensSystem, analyze_system


DATA_DIR = Path("data/cosmograil")
OUT_DIR = Path("results/outputs")
OUT_JSON = OUT_DIR / "step_3_14_he1104_multiband_chromaticity.json"

CAT_ID = "J/ApJ/798/95"


def fetch_table():
    Vizier.ROW_LIMIT = -1
    tables = Vizier.get_catalogs(CAT_ID)
    if not tables:
        raise RuntimeError(f"No tables returned for {CAT_ID}")
    return tables[0]


def to_system(rows, band: str) -> LensSystem:
    # rows is an astropy Table already filtered
    t = np.array(rows["HJD"], dtype=float)
    lcA = LightCurve("A", t=t, mag=np.array(rows["A"], dtype=float), magerr=np.array(rows["e_A"], dtype=float))
    lcB = LightCurve("B", t=t, mag=np.array(rows["B"], dtype=float), magerr=np.array(rows["e_B"], dtype=float))
    return LensSystem(system_id=f"HE1104_{band}", light_curves={"A": lcA, "B": lcB}, band=band)


def compare_pair_gamma(res1: Dict, res2: Dict, pair: str = "A-B") -> Dict:
    g1 = res1["pairs"][pair]["gamma"]
    g2 = res2["pairs"][pair]["gamma"]

    v1 = float(g1["value"]) if np.isfinite(g1["value"]) else np.nan
    e1 = float(g1["uncertainty"]) if np.isfinite(g1["uncertainty"]) else np.nan
    v2 = float(g2["value"]) if np.isfinite(g2["value"]) else np.nan
    e2 = float(g2["uncertainty"]) if np.isfinite(g2["uncertainty"]) else np.nan

    if not (np.isfinite(v1) and np.isfinite(v2) and np.isfinite(e1) and np.isfinite(e2) and e1 > 0 and e2 > 0):
        return {}

    delta = v1 - v2
    sigma = float(np.sqrt(e1 * e1 + e2 * e2))
    sig = float(abs(delta) / sigma) if sigma > 0 else np.nan

    return {
        "delta_gamma": float(delta),
        "delta_uncertainty": float(sigma),
        "delta_sigma": float(sig),
        "band1_gamma": float(v1),
        "band1_uncertainty": float(e1),
        "band2_gamma": float(v2),
        "band2_uncertainty": float(e2),
    }


def main() -> None:
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    t = fetch_table()

    # Filters in this catalog: B, I, J, R, V (V has very few points)
    bands = ["R", "B", "I", "J"]

    band_tables = {}
    for b in bands:
        mask = np.array([str(x).strip() == b for x in t["Filt"]])
        tb = t[mask]
        band_tables[b] = tb
        tb.write(DATA_DIR / f"he1104_JApJ798_95_{b}.csv", format="csv", overwrite=True)

    # Canonical estimator settings
    # HE1104 has O(100d) delay, so allow a wide lag window.
    lag_range = (-300, 300)
    tau_values: List[float] = [5, 10, 20, 40, 80, 160, 320]

    results: Dict[str, Dict] = {}
    for b in bands:
        sys = to_system(band_tables[b], b)
        results[b] = analyze_system(sys, detrend_window=200.0, tau_values=tau_values, lag_range=lag_range)

    # Comparisons vs R (most sampled)
    comparisons: Dict[str, Dict] = {}
    base = "R"
    for b in [x for x in bands if x != base]:
        comp = compare_pair_gamma(results[b], results[base])
        if comp:
            comparisons[f"{b}_minus_{base}"] = comp

    out = {
        "catalog": CAT_ID,
        "system": "HE1104-1805",
        "bands_analyzed": bands,
        "tau_values": tau_values,
        "lag_range": lag_range,
        "results": results,
        "comparisons": comparisons,
    }

    with open(OUT_JSON, "w") as f:
        json.dump(out, f, indent=2)

    # Console summary
    for k, v in comparisons.items():
        print(f"{k} (A-B): ΔΓ = {v['delta_gamma']:.2f} ± {v['delta_uncertainty']:.2f} (σ={v['delta_sigma']:.2f})")

    gR = results["R"]["pairs"]["A-B"]["gamma"]
    print(f"R-band A-B: Gamma = {gR['value']:.2f} ± {gR['uncertainty']:.2f} (sigma={gR['sigma']:.2f})")
    print(f"Saved: {OUT_JSON}")


if __name__ == "__main__":
    main()
