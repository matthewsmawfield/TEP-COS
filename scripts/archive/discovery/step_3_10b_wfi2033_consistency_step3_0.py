#!/usr/bin/env python3
"""Step 3.10b: Instrumental Consistency (WFI2033) using canonical Step 3.0 estimator.

This recomputes Γ separately on EulerCAM and SMARTS WFI2033 datasets using
`scripts/steps/step_3_0_cosmograil_temporal_shear.py`.

Outputs:
- results/outputs/step_3_10b_wfi2033_consistency_step3_0.json

Notes:
- WFI2033_ecam.dat and WFI2033_smarts.dat format: MHJD A errA B errB C errC
- We report Γ for each pair and also a simple ΔΓ for pair A-B (if available).
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, Tuple

import numpy as np

from step_3_0_cosmograil_temporal_shear import LightCurve, LensSystem, analyze_system


DATA_DIR = Path("data/cosmograil")
OUT_FILE = Path("results/outputs/step_3_10b_wfi2033_consistency_step3_0.json")


def load_wfi_dat(path: Path, system_id: str) -> LensSystem:
    raw = np.loadtxt(path)
    if raw.ndim != 2 or raw.shape[1] < 7:
        raise ValueError(f"Unexpected format for {path.name}: shape={raw.shape}")

    t = raw[:, 0]
    lcs = {
        "A": LightCurve("A", t=t, mag=raw[:, 1], magerr=raw[:, 2]),
        "B": LightCurve("B", t=t, mag=raw[:, 3], magerr=raw[:, 4]),
        "C": LightCurve("C", t=t, mag=raw[:, 5], magerr=raw[:, 6]),
    }

    return LensSystem(system_id=system_id, light_curves=lcs, band="R")


def extract_pair_gamma(results: Dict, pair: str) -> Tuple[float, float, float]:
    """Return (gamma, gamma_err, sigma)."""
    g = results["pairs"][pair]["gamma"]
    return float(g["value"]), float(g["uncertainty"]), float(g["sigma"])


def main() -> None:
    ecam_path = DATA_DIR / "WFI2033_ecam.dat"
    smarts_path = DATA_DIR / "WFI2033_smarts.dat"

    if not ecam_path.exists() or not smarts_path.exists():
        raise FileNotFoundError("Missing WFI2033_ecam.dat or WFI2033_smarts.dat")

    sys_ecam = load_wfi_dat(ecam_path, "WFI2033_EulerCAM")
    sys_smarts = load_wfi_dat(smarts_path, "WFI2033_SMARTS")

    # Use the same tau grid as Step 3.0 defaults.
    # Lag range stays at the canonical default (±200 d), since WFI2033 delays are small.
    res_ecam = analyze_system(sys_ecam, detrend_window=200.0, tau_values=[5, 10, 20, 40, 80, 160], lag_range=(-200, 200))
    res_smarts = analyze_system(sys_smarts, detrend_window=200.0, tau_values=[5, 10, 20, 40, 80, 160], lag_range=(-200, 200))

    out: Dict = {
        "ecam": res_ecam,
        "smarts": res_smarts,
        "summary": {},
    }

    # Prefer the A-B pair for a single ΔΓ comparison (if present)
    pair_key = "A-B"
    if pair_key in res_ecam["pairs"] and pair_key in res_smarts["pairs"]:
        g1, e1, s1 = extract_pair_gamma(res_ecam, pair_key)
        g2, e2, s2 = extract_pair_gamma(res_smarts, pair_key)
        out["summary"]["pair"] = pair_key
        out["summary"]["ecam_gamma"] = g1
        out["summary"]["ecam_gamma_err"] = e1
        out["summary"]["ecam_sigma"] = s1
        out["summary"]["smarts_gamma"] = g2
        out["summary"]["smarts_gamma_err"] = e2
        out["summary"]["smarts_sigma"] = s2
        out["summary"]["delta_gamma"] = abs(g1 - g2)

    OUT_FILE.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT_FILE, "w") as f:
        json.dump(out, f, indent=2)

    print("WFI2033 instrumental consistency (Step 3.0 estimator)")
    if "pair" in out["summary"]:
        print(f"Pair {out['summary']['pair']}")
        print(f"EulerCAM: {out['summary']['ecam_gamma']:.2f} ± {out['summary']['ecam_gamma_err']:.2f} (σ={out['summary']['ecam_sigma']:.2f})")
        print(f"SMARTS:   {out['summary']['smarts_gamma']:.2f} ± {out['summary']['smarts_gamma_err']:.2f} (σ={out['summary']['smarts_sigma']:.2f})")
        print(f"ΔΓ = {out['summary']['delta_gamma']:.2f} days/decade")
    print(f"Saved: {OUT_FILE}")


if __name__ == "__main__":
    main()
