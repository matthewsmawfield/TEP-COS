#!/usr/bin/env python3
"""Step 3.13: Multi-band chromaticity test (HE0435-1223)

Goal
- Provide a second multi-band chromaticity test using a lens with public V and R light curves.

Dataset
- VizieR: J/A+A/703/A250 (Sorgenfrei+ 2025): HE0435-1223 V and R light curves.

Method
- Download the V and R tables.
- Run canonical temporal-shear estimator (Step 3.0) independently on each band.
- Compute ΔΓ between V and R per image pair.

Outputs
- data/cosmograil/he0435_JAA703A250_V.csv
- data/cosmograil/he0435_JAA703A250_R.csv
- results/outputs/step_3_13_he0435_multiband_chromaticity.json
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Dict

import numpy as np

from astroquery.vizier import Vizier

from step_3_0_cosmograil_temporal_shear import LightCurve, LensSystem, analyze_system


DATA_DIR = Path("data/cosmograil")
OUT_DIR = Path("results/outputs")
OUT_JSON = OUT_DIR / "step_3_13_he0435_multiband_chromaticity.json"

CAT_ID = "J/A+A/703/A250"
TABLES = {
    "V": "he0435v",
    "R": "he0435r",
}


def fetch_table(suffix: str):
    Vizier.ROW_LIMIT = -1
    tables = Vizier.get_catalogs(CAT_ID)
    for t in tables:
        name = t.meta.get("name", "")
        if name.endswith(suffix):
            return t
    raise RuntimeError(f"Could not find {suffix} in {CAT_ID}")


def table_to_system(table, band: str) -> LensSystem:
    t = np.array(table["MJD"], dtype=float)

    if band == "V":
        mapping = {
            "A": ("VmagA", "e_VmagA"),
            "B": ("VmagB", "e_VmagB"),
            "C": ("VmagC", "e_VmagC"),
            "D": ("VmagD", "e_VmagD"),
        }
    elif band == "R":
        mapping = {
            "A": ("RmagA", "e_RmagA"),
            "B": ("RmagB", "e_RmagB"),
            "C": ("RmagC", "e_RmagC"),
            "D": ("RmagD", "e_RmagD"),
        }
    else:
        raise ValueError(band)

    lcs = {}
    for img, (mcol, ecol) in mapping.items():
        mag = np.array(table[mcol], dtype=float)
        err = np.array(table[ecol], dtype=float)
        lcs[img] = LightCurve(label=img, t=t, mag=mag, magerr=err)

    return LensSystem(system_id=f"HE0435_{band}", light_curves=lcs, band=band)


def compare_bands(res_v: Dict, res_r: Dict) -> Dict:
    pairs = set(res_v["pairs"].keys()) & set(res_r["pairs"].keys())
    out: Dict = {}

    for pair in sorted(pairs):
        gv = res_v["pairs"][pair]["gamma"]
        gr = res_r["pairs"][pair]["gamma"]

        v1 = float(gv["value"]) if np.isfinite(gv["value"]) else np.nan
        e1 = float(gv["uncertainty"]) if np.isfinite(gv["uncertainty"]) else np.nan
        v2 = float(gr["value"]) if np.isfinite(gr["value"]) else np.nan
        e2 = float(gr["uncertainty"]) if np.isfinite(gr["uncertainty"]) else np.nan

        if not (np.isfinite(v1) and np.isfinite(v2) and np.isfinite(e1) and np.isfinite(e2) and e1 > 0 and e2 > 0):
            continue

        delta = v1 - v2
        sigma = float(np.sqrt(e1 * e1 + e2 * e2))
        sig = float(abs(delta) / sigma) if sigma > 0 else np.nan

        out[pair] = {
            "delta_gamma": float(delta),
            "delta_uncertainty": float(sigma),
            "delta_sigma": float(sig),
            "V_gamma": float(v1),
            "V_uncertainty": float(e1),
            "R_gamma": float(v2),
            "R_uncertainty": float(e2),
        }

    return out


def main() -> None:
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    t_v = fetch_table(TABLES["V"])
    t_r = fetch_table(TABLES["R"])

    # Persist inputs
    t_v.write(DATA_DIR / "he0435_JAA703A250_V.csv", format="csv", overwrite=True)
    t_r.write(DATA_DIR / "he0435_JAA703A250_R.csv", format="csv", overwrite=True)

    sys_v = table_to_system(t_v, "V")
    sys_r = table_to_system(t_r, "R")

    # HE0435 delays are O(days) to O(weeks), so a modest lag range is sufficient.
    tau_values = [5, 10, 20, 40, 80, 160]

    res_v = analyze_system(sys_v, detrend_window=200.0, tau_values=tau_values, lag_range=(-50, 50))
    res_r = analyze_system(sys_r, detrend_window=200.0, tau_values=tau_values, lag_range=(-50, 50))

    comp = compare_bands(res_v, res_r)

    out = {
        "catalog": CAT_ID,
        "system": "HE0435-1223",
        "bands": {
            "V": res_v,
            "R": res_r,
        },
        "comparisons": {
            "V_minus_R": comp,
        },
    }

    with open(OUT_JSON, "w") as f:
        json.dump(out, f, indent=2)

    if "A-B" in comp:
        d = comp["A-B"]
        print(f"V_minus_R A-B: ΔΓ = {d['delta_gamma']:.2f} ± {d['delta_uncertainty']:.2f} (σ={d['delta_sigma']:.2f})")
    print(f"Saved: {OUT_JSON}")


if __name__ == "__main__":
    main()
