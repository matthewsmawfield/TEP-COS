#!/usr/bin/env python3
"""Step 3.11: Q0957+561 Temporal Shear (GLENDAMA)

Purpose
- Use the canonical COSMOGRAIL temporal-shear estimator (Step 3.0) on GLENDAMA QSO B0957+561.
- Compute multiscale delays and fit Γ = d(Δt)/d(log τ) for the A-B pair.

Input
- data/cosmograil/glendama_J_ApA_616_A118_table6.csv
  Columns: MJD, mA, e_mA, mB, e_mB

Notes on time axis
- The GLENDAMA VizieR tables use a compact time coordinate. For table6, values span 117..7528.
  This matches the common convention JD-2450000 (1996..2016). The estimator only requires
  relative time units in days, so we keep the coordinate as-is (days).

Outputs
- results/outputs/step_3_11_q0957_glendama_temporal_shear.json
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from step_3_0_cosmograil_temporal_shear import LightCurve, LensSystem, analyze_system


DATA_FILE = Path("data/cosmograil/glendama_J_ApA_616_A118_table6.csv")
OUT_FILE = Path("results/outputs/step_3_11_q0957_glendama_temporal_shear.json")


def load_q0957_glendama() -> LensSystem:
    if not DATA_FILE.exists():
        raise FileNotFoundError(str(DATA_FILE))

    df = pd.read_csv(DATA_FILE)

    required = {"MJD", "mA", "e_mA", "mB", "e_mB"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing columns in {DATA_FILE.name}: {sorted(missing)}")

    # Keep finite points only
    t = df["MJD"].to_numpy(dtype=float)

    mA = df["mA"].to_numpy(dtype=float)
    eA = df["e_mA"].to_numpy(dtype=float)

    mB = df["mB"].to_numpy(dtype=float)
    eB = df["e_mB"].to_numpy(dtype=float)

    lcA = LightCurve(label="A", t=t, mag=mA, magerr=eA)
    lcB = LightCurve(label="B", t=t, mag=mB, magerr=eB)

    return LensSystem(system_id="Q0957_GLENDAMA", light_curves={"A": lcA, "B": lcB}, band="R")


def main() -> None:
    system = load_q0957_glendama()

    # Q0957 delay is ~417 days, so widen lag range well beyond the default.
    # We also extend tau_values toward longer timescales, since Q0957 has a long baseline.
    results = analyze_system(
        system,
        detrend_window=200.0,
        tau_values=[5, 10, 20, 40, 80, 160, 320],
        lag_range=(-700, 700),
    )

    OUT_FILE.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT_FILE, "w") as f:
        json.dump(results, f, indent=2)

    # Minimal console summary for the single pair
    pair_key = "A-B"
    pair = results["pairs"].get(pair_key)
    if pair:
        bb = pair["broadband"]
        g = pair["gamma"]
        print("Q0957+561 (GLENDAMA) summary")
        print(f"Broadband delay: {bb['delay_days']:.2f} d (r={bb['correlation']:.3f}, σ≈{bb['uncertainty_days']:.2f} d)")
        print(f"Gamma: {g['value']:.2f} ± {g['uncertainty']:.2f} days/decade (sigma={g['sigma']:.2f})")
        print(f"Saved: {OUT_FILE}")


if __name__ == "__main__":
    main()
