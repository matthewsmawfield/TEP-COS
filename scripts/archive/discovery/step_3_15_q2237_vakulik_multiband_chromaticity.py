#!/usr/bin/env python3
"""Step 3.15: Multi-band chromaticity test (Q2237+0305 Vakulik VRI)

Goal
- Provide a fourth independent multi-band chromaticity test using a different Q2237 dataset.

Dataset
- VizieR: J/A+A/420/447 (Vakulik+ 2004): Q2237+0305 VRI photometry 1995-2000.
  Table6: VRI photometry of Q2237+0305 A,B,C,D.

Method
- Split into per-band light curves.
- Run canonical temporal-shear estimator (Step 3.0) for each image pair in each band.
- Compute ΔΓ between bands for A-B with propagated uncertainty.

Outputs
- data/cosmograil/q2237_vakulik_JAA420_447_{band}.csv
- results/outputs/step_3_15_q2237_vakulik_multiband_chromaticity.json
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
OUT_JSON = OUT_DIR / "step_3_15_q2237_vakulik_multiband_chromaticity.json"

CAT_ID = "J/A+A/420/447"


def fetch_table6():
    Vizier.ROW_LIMIT = -1
    tables = Vizier.get_catalogs(CAT_ID)
    for t in tables:
        name = t.meta.get("name", "")
        if "table6" in name:
            return t
    raise RuntimeError(f"Could not find table6 in {CAT_ID}")


def to_system(rows, band: str) -> LensSystem:
    # Convert date string to MJD (approximate: YYYY-MM-DD -> JD)
    from astropy.time import Time
    dates = [str(d) for d in rows["Obs.Date"]]
    t = np.array([Time(d, format="iso").mjd for d in dates])

    lcs = {}
    for img in ["A", "B", "C", "D"]:
        mcol = f"mag{img}"
        ecol = f"e_mag{img}"
        mag = np.array(rows[mcol], dtype=float)
        err = np.array(rows[ecol], dtype=float)
        lcs[img] = LightCurve(label=img, t=t, mag=mag, magerr=err)

    return LensSystem(system_id=f"Q2237_Vakulik_{band}", light_curves=lcs, band=band)


def compare_pair_gamma(res1: Dict, res2: Dict, pair: str = "A-B") -> Dict:
    if pair not in res1["pairs"] or pair not in res2["pairs"]:
        return {}

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

    t = fetch_table6()

    bands = ["V", "R", "I"]
    band_tables = {}
    for b in bands:
        mask = np.array([str(x).strip() == b for x in t["Band"]])
        tb = t[mask]
        band_tables[b] = tb
        tb.write(DATA_DIR / f"q2237_vakulik_JAA420_447_{b}.csv", format="csv", overwrite=True)

    # Q2237 has very short physical delays (~hours), so use a narrow lag window
    lag_range = (-10, 10)
    tau_values: List[float] = [5, 10, 20, 40, 80, 160]

    results: Dict[str, Dict] = {}
    for b in bands:
        sys = to_system(band_tables[b], b)
        results[b] = analyze_system(sys, detrend_window=200.0, tau_values=tau_values, lag_range=lag_range)

    # Comparisons vs R
    comparisons: Dict[str, Dict] = {}
    base = "R"
    for b in [x for x in bands if x != base]:
        comp = compare_pair_gamma(results[b], results[base])
        if comp:
            comparisons[f"{b}_minus_{base}"] = comp

    out = {
        "catalog": CAT_ID,
        "system": "Q2237+0305 (Vakulik)",
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

    if "R" in results and "A-B" in results["R"]["pairs"]:
        gR = results["R"]["pairs"]["A-B"]["gamma"]
        print(f"R-band A-B: Gamma = {gR['value']:.2f} ± {gR['uncertainty']:.2f} (sigma={gR['sigma']:.2f})")

    print(f"Saved: {OUT_JSON}")


if __name__ == "__main__":
    main()
