#!/usr/bin/env python3
"""
Step 24: Equal Cluster Weighting Analysis for Density Scaling

This script performs a robustness test of the suppressed density scaling result
by comparing weighting schemes on the ACTUAL pipeline data (step_07 controlled
residuals and step_02 raw cluster means), not hardcoded representative values.

Two dependent variables are analyzed:
1. Raw logPdot_abs (cluster mean observed spin-down) — comparable to step_12
2. Controlled residual (GC − matched field, from step_07) — field-subtracted

Purpose: Address the concern that extreme clusters (Terzan 5, NGC 6517) with
large pulsar populations may dominate the statistics.

Author: TEP Collaboration
Date: 2026-03-30
"""

import json
import numpy as np
import pandas as pd
from scipy import stats
from pathlib import Path
from typing import Dict, Tuple

REPO = Path(__file__).resolve().parent.parent.parent

# Cluster central densities (log10 rho_c in L_sun/pc^3)
CLUSTER_DENSITIES = {
    "Terzan 5": 5.50, "47 Tuc (NGC 104)": 4.88, "NGC 6517": 5.80,
    "M28 (NGC 6626)": 4.52, "M62 (NGC 6266)": 5.16, "M13 (NGC 6205)": 3.79,
    "M15 (NGC 7078)": 5.05, "M5 (NGC 5904)": 3.53, "Terzan 1": 5.00,
    "NGC 6752": 4.30, "M2 (NGC 7089)": 4.15, "Omega Centauri (NGC 5139)": 3.12,
    "M53 (NGC 5024)": 2.96, "M3 (NGC 5272)": 3.68, "M71 (NGC 6838)": 2.29,
    "NGC 6397": 5.68, "NGC 1851": 5.09, "NGC 6522": 5.50,
    "NGC 6544": 5.20, "NGC 6624": 5.60, "NGC 6760": 3.80,
    "M22 (NGC 6656)": 2.97, "M80 (NGC 6093)": 4.79, "M92 (NGC 6341)": 4.30,
    "NGC 6712": 3.70, "NGC 6652": 4.50, "M14 (NGC 6402)": 3.44,
    "NGC 6539": 3.30, "M4 (NGC 6121)": 2.85, "NGC 6440": 5.1,
    "NGC 6441": 5.0,
}


def ols_with_stderr(x: np.ndarray, y: np.ndarray) -> Tuple[float, float, float, float]:
    """OLS regression with standard errors."""
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
    return float(slope), float(intercept), float(std_err), float(p_value)


def wls_by_n(x: np.ndarray, y: np.ndarray, n: np.ndarray) -> Tuple[float, float, float, float]:
    """Weighted least squares by sample size (repeat observations)."""
    xw = np.repeat(x, n)
    yw = np.repeat(y, n)
    slope, intercept, r_value, p_value, std_err = stats.linregress(xw, yw)
    return float(slope), float(intercept), float(std_err), float(p_value)


def analyze_density_scaling(x, y, weights, label):
    """Run OLS and WLS and return results dict."""
    x = np.array(x)
    y = np.array(y)
    w = np.array(weights)

    slope_ols, int_ols, se_ols, p_ols = ols_with_stderr(x, y)
    slope_wls, int_wls, se_wls, p_wls = wls_by_n(x, y, w)

    return {
        "label": label,
        "n_clusters": len(x),
        "ols_slope": slope_ols,
        "ols_stderr": se_ols,
        "ols_p": p_ols,
        "wls_slope": slope_wls,
        "wls_stderr": se_wls,
        "wls_p": p_wls,
    }


def main():
    print("--- Step 24: Equal Cluster Weighting Analysis (Real Data) ---")

    # Load actual pipeline outputs
    with open(REPO / "results/outputs/step_07_per_cluster_controlled_residuals.json") as f:
        s31 = json.load(f)

    df_pulsars = pd.read_csv(REPO / "results/outputs/step_02_pulsar_population_controls.csv")
    gc_df = df_pulsars[df_pulsars['environment'] == 'globular_cluster']

    # Build cluster-level datasets
    clusters_ctrl = []   # controlled residuals
    clusters_raw = []    # raw logPdot means
    for name, data in s31["clusters"].items():
        if name not in CLUSTER_DENSITIES:
            continue
        clusters_ctrl.append({
            "cluster": name,
            "log_rho": CLUSTER_DENSITIES[name],
            "y": data["controlled_residual"],
            "n": data["n_pulsars"],
        })
        gc_cluster = gc_df[gc_df['cluster'] == name]
        if len(gc_cluster) > 0:
            clusters_raw.append({
                "cluster": name,
                "log_rho": CLUSTER_DENSITIES[name],
                "y": gc_cluster['logPdot_abs'].mean(),
                "n": len(gc_cluster),
            })

    # Analysis A: Raw logPdot (comparable to step_12 primary)
    raw_results = analyze_density_scaling(
        [c["log_rho"] for c in clusters_raw],
        [c["y"] for c in clusters_raw],
        [c["n"] for c in clusters_raw],
        "Raw cluster-mean log|Ṗ|"
    )

    # Analysis B: Controlled residuals (field-subtracted)
    ctrl_results = analyze_density_scaling(
        [c["log_rho"] for c in clusters_ctrl],
        [c["y"] for c in clusters_ctrl],
        [c["n"] for c in clusters_ctrl],
        "Controlled residual (GC − matched field)"
    )

    # Load Newtonian prediction from literature consensus
    s48_path = Path('results/outputs/step_14_cmc_literature.json')
    if s48_path.exists():
        with open(s48_path) as f:
            s48 = json.load(f)
        comp = s48.get('comparison', {})
        newt_slope = comp.get('cmc_predicted_slope', 0.748)
        newt_err = comp.get('cmc_predicted_error', 0.039)
    else:
        newt_slope = 0.748
        newt_err = 0.039
        print(f"  Warning: {s48_path} not found; using literature consensus fallback.")

    def tension(slope, stderr):
        diff = slope - newt_slope
        combined = np.sqrt(stderr**2 + newt_err**2)
        return abs(diff) / combined

    output = {
        "analysis_raw_logpdot": {
            **raw_results,
            "ols_tension_sigma": tension(raw_results["ols_slope"], raw_results["ols_stderr"]),
            "wls_tension_sigma": tension(raw_results["wls_slope"], raw_results["wls_stderr"]),
        },
        "analysis_controlled_residual": {
            **ctrl_results,
            "ols_tension_sigma": tension(ctrl_results["ols_slope"], ctrl_results["ols_stderr"]),
            "wls_tension_sigma": tension(ctrl_results["wls_slope"], ctrl_results["wls_stderr"]),
        },
        "note": (
            "The raw-logPdot and controlled-residual slopes differ because the "
            "matched field mean itself correlates with cluster density (period/B-proxy "
            "distributions differ across clusters). The primary discriminant uses raw "
            "logPdot with mixed-effects (step_12)."
        ),
    }

    out_path = REPO / "results/outputs/step_24_equal_cluster_weighting.json"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(output, f, indent=2)

    print(f"\nResults saved to: {out_path}")
    print("\n" + "=" * 70)
    print("DENSITY SLOPE DEFINITIONS SUMMARY")
    print("=" * 70)
    print(f"\nA) Raw log|Ṗ| (comparable to step_12 primary):")
    print(f"   OLS (equal weight):  {raw_results['ols_slope']:.3f} ± {raw_results['ols_stderr']:.3f}")
    print(f"   WLS (by N):          {raw_results['wls_slope']:.3f} ± {raw_results['wls_stderr']:.3f}")
    print(f"   Tension vs Newtonian: {tension(raw_results['ols_slope'], raw_results['ols_stderr']):.2f}σ")

    print(f"\nB) Controlled residual (GC − matched field):")
    print(f"   OLS (equal weight):  {ctrl_results['ols_slope']:.3f} ± {ctrl_results['ols_stderr']:.3f}")
    print(f"   WLS (by N):          {ctrl_results['wls_slope']:.3f} ± {ctrl_results['wls_stderr']:.3f}")
    print(f"   Tension vs Newtonian: {tension(ctrl_results['ols_slope'], ctrl_results['ols_stderr']):.2f}σ")
    print("\n" + "=" * 70)


if __name__ == "__main__":
    main()
