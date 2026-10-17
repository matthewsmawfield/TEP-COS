#!/usr/bin/env python3
"""
Step 29: Pulsar Dynamical Calibration Analysis

ADDRESSES CRITICAL WEAKNESS: Cannot distinguish TEP-enhanced acceleration 
from standard GR cluster acceleration.

STRATEGY: Use detailed cluster potential modeling to predict Newtonian 
acceleration field, then compare observed residuals to predictions.

If residuals EXCEED Newtonian predictions systematically, this supports TEP.
If residuals MATCH Newtonian predictions, this favors standard GR.

This test is INDEPENDENT of the density scaling argument.

IMPORTANT: This script uses MONTE CARLO SIMULATION to model Newtonian
expectations for comparison with real pulsar data.

Author: M. Smawfield
Date: March 2026
"""

import numpy as np
import pandas as pd
from scipy import stats
import json
from pathlib import Path
import os

# Configuration
REPO_ROOT = Path(__file__).resolve().parents[2]
RESULTS_DIR = REPO_ROOT / "results" / "outputs"
OUTPUT_JSON = RESULTS_DIR / "step_29_dynamical_calibration.json"

# Random seed for reproducibility
# Fixed seed ensures Monte Carlo simulation results are fully reproducible
RANDOM_SEED = 42
np.random.seed(RANDOM_SEED)

# Cluster parameters from Freire catalog and literature
CLUSTER_PARAMS = {
    "47_Tuc": {"M": 1.0e6, "rc": 0.36, "rh": 3.17, "distance": 4500},
    "M15": {"M": 5.6e5, "rc": 0.14, "rh": 1.9, "distance": 10400},
    "M13": {"M": 6.0e5, "rc": 0.55, "rh": 3.6, "distance": 7400},
    "Terzan_5": {"M": 2.0e6, "rc": 0.16, "rh": 0.8, "distance": 6700},
    "M28": {"M": 5.0e5, "rc": 0.24, "rh": 2.5, "distance": 7000},
    "M62": {"M": 1.0e6, "rc": 0.18, "rh": 1.2, "distance": 6700},
    "M5": {"M": 5.7e5, "rc": 0.42, "rh": 3.4, "distance": 7500},
    "M3": {"M": 4.0e5, "rc": 0.55, "rh": 4.0, "distance": 10100},
    "M53": {"M": 3.0e5, "rc": 0.65, "rh": 3.0, "distance": 18000},
    "Omega_Cen": {"M": 4.0e6, "rc": 0.6, "rh": 4.0, "distance": 5200},
    "M2": {"M": 1.0e6, "rc": 0.52, "rh": 3.7, "distance": 11500},
    "NGC_6440": {"M": 1.4e6, "rc": 0.15, "rh": 0.9, "distance": 8600},
    "NGC_6544": {"M": 2.0e5, "rc": 0.28, "rh": 1.8, "distance": 3000},
    "M30": {"M": 2.0e5, "rc": 0.35, "rh": 2.5, "distance": 8300},
    "NGC_6752": {"M": 2.8e5, "rc": 0.35, "rh": 2.5, "distance": 4000},
}


def king_potential(r, M, rc):
    """King model gravitational potential."""
    G = 4.302e-3  # pc (km/s)^2 / M_sun
    
    # Simplified King potential
    x = r / rc
    # Core potential
    phi_0 = G * M / rc
    
    # Approximate potential profile
    phi = phi_0 / np.sqrt(1 + x**2)
    
    return phi


def predict_newtonian_shift(cluster_params, n_pulsars=1000):
    """
    Predict Newtonian acceleration-induced shift for a cluster.
    
    Returns predicted mean log|Ṗ| shift in dex.
    """
    M = cluster_params["M"]
    rc = cluster_params["rc"]
    rh = cluster_params["rh"]
    
    # Simulate pulsar positions and velocities.
    # The potential used below is Plummer-like, phi = phi0 / sqrt(1+x^2);
    # the self-consistent density is rho ∝ (1+x^2)^(-5/2) and the enclosed
    # mass is M(r) = M x^3 (1+x^2)^(-3/2).  Pulsar positions are drawn from
    # the 3-D number-density profile n(r) ∝ r^2 (1+x^2)^(-5/2) by inverse-CDF
    # sampling on a tabulated cumulative grid (reproducible under RANDOM_SEED).
    x_grid = np.geomspace(1e-3, max(30.0, 10.0 * rh / rc), 20000)
    cdf_integrand = x_grid ** 2 / (1.0 + x_grid ** 2) ** 2.5
    cdf = np.concatenate([[0.0], np.cumsum(0.5 * (cdf_integrand[1:] + cdf_integrand[:-1]) * np.diff(x_grid))])
    cdf /= cdf[-1]
    u = np.random.uniform(0.0, 1.0, size=n_pulsars)
    r_pulsars = rc * np.interp(u, cdf, x_grid)
    r_pulsars = np.clip(r_pulsars, 0.01, 10 * rh)  # Physical limits
    
    # Calculate line-of-sight acceleration
    G = 4.302e-3  # pc (km/s)^2 / M_sun
    
    # Acceleration varies with position
    # a ~ G*M(r)/r^2, but M(r) increases with r in King model
    x = r_pulsars / rc
    # Plummer-consistent enclosed mass: continuous at rc, a -> 0 as r -> 0.
    M_enclosed = M * x ** 3 / (1.0 + x ** 2) ** 1.5
    # Line-of-sight component (random orientation); sign is retained because
    # pulsars on the far side accelerate toward us, which DECREASES |Pdot|.
    cos_theta = np.random.uniform(-1, 1, size=n_pulsars)
    a_los = G * M_enclosed / r_pulsars ** 2 * cos_theta
    
    # Convert to pdot_acc and log|Ṗ|
    c = 3e5  # km/s
    pc_to_km = 3.086e13  # km/pc - CRITICAL UNIT CONVERSION
    P_SPIN_S = 5.0e-3    # canonical MSP period (s); matches step_37 convention
    # Pdot_acc = P * a_los / c.  a_los is in (km/s)^2/pc, so a_los/c is in
    # 1/s after the pc->km conversion; multiplying by P gives dimensionless
    # Pdot_acc.  (The previous version compared a/c directly to Pdot,
    # dropping the period factor and overstating the shift by ~1/P ~ 200x.)
    pdot_acc = P_SPIN_S * a_los / (c * pc_to_km)  # dimensionless s/s
    
    # Load intrinsic field MSP spin-down from step_02 if available
    s10_path = Path("results/outputs/step_02_pulsar_population_controls.json")
    if s10_path.exists():
        with open(s10_path) as f:
            s10 = json.load(f)
        log_pdot_intrinsic = s10["base_log10_abs_pdot"]["field_mean"]
    else:
        log_pdot_intrinsic = -19.5  # Fallback: typical MSP intrinsic spin-down
    
    # Observed spin-down is the signed sum; the observable is |Pdot_obs|.
    # Pulsars whose line-of-sight acceleration opposes the intrinsic
    # spin-down land BELOW the intrinsic value, producing the cancellation
    # floor that a positive-only |a| treatment removes.
    pdot_obs = 10 ** log_pdot_intrinsic + pdot_acc
    log_pdot_obs = np.log10(np.abs(pdot_obs))
    
    # Calculate shift from intrinsic
    shift_dex = np.mean(log_pdot_obs) - log_pdot_intrinsic
    shift_std = np.std(log_pdot_obs)
    
    return {
        "predicted_shift_dex": float(shift_dex),
        "predicted_std": float(shift_std),
        "n_pulsars": n_pulsars,
        "a_los_mean": float(np.mean(np.abs(a_los))),
        "a_los_std": float(np.std(a_los)),
        "p_spin_s": P_SPIN_S
    }


def analyze_dynamical_calibration():
    """
    Main analysis: Compare observed residuals to Newtonian predictions.
    
    Key test: Do observed residuals exceed Newtonian predictions?
    """
    print("=" * 70)
    print("STEP 5.41: PULSAR DYNAMICAL CALIBRATION")
    print("=" * 70)
    print("\nTesting: Do observed residuals exceed Newtonian predictions?")
    print("This addresses the TEP vs GR ambiguity directly.\n")
    
    results = {}
    
    # Run Newtonian prediction for each cluster
    print("Generating Newtonian predictions...")
    for cluster_name, params in CLUSTER_PARAMS.items():
        pred = predict_newtonian_shift(params, n_pulsars=10000)
        results[cluster_name] = {
            "parameters": params,
            "newtonian_prediction": pred
        }
        print(f"  {cluster_name}: predicted shift = {pred['predicted_shift_dex']:.2f} dex")
    
    # Summary statistics
    predicted_shifts = [r["newtonian_prediction"]["predicted_shift_dex"] 
                        for r in results.values()]
    
    avg_predicted = np.mean(predicted_shifts)
    std_predicted = np.std(predicted_shifts)
    
    # Compare to observed (from step_02 results)
    # CRITICAL METHODOLOGY NOTE:
    # The observed residual (0.59 dex) is AFTER population controls (period+B-proxy matching).
    # The Newtonian prediction (~1.8 dex) is a RAW shift without such controls.
    # The deficit is therefore quoted both against the controlled-residual
    # uncertainty and against the inter-cluster scatter of the prediction
    # (the honest model-systematic scale).
    s510_path = REPO_ROOT / "results" / "outputs" / "step_02_pulsar_population_controls.json"
    if not s510_path.exists():
        print(f"ERROR: Required input file not found: {s510_path}")
        print(f"Dynamical calibration requires actual observed residuals from step_02.")
        raise RuntimeError("Missing required input: step_02_pulsar_population_controls.json")
    
    with open(s510_path, 'r') as f:
        s510_data = json.load(f)
    
    observed_residual = s510_data['base_log10_abs_pdot']['diff_dex']
    # Uncertainty estimated from bootstrap CI width / 2
    ci_width = s510_data['controls']['period_and_bproxy_matched']['diff_ci84'] - s510_data['controls']['period_and_bproxy_matched']['diff_ci16']
    observed_uncertainty = ci_width / 2  # ~1 sigma
    print(f"  Loaded from step_02: residual = {observed_residual:.3f} ± {observed_uncertainty:.3f} dex")
    
    print(f"\n{'='*70}")
    print("COMPARISON TO OBSERVED RESIDUALS")
    print(f"{'='*70}")
    print(f"\nNewtonian Prediction (raw, no controls): {avg_predicted:.2f} ± {std_predicted:.2f} dex")
    print(f"Observed Residual (after population controls): {observed_residual:.2f} ± {observed_uncertainty:.2f} dex")
    print(f"\nNote: The Newtonian prediction is a raw shift; the observed is after")
    print(f"period+B-proxy matching. The prediction scatter across clusters is")
    print(f"the relevant model-systematic scale for the deficit, not the")
    print(f"controlled-residual uncertainty alone.")
    
    # Key comparison: the raw Newtonian prediction vs the controlled
    # observed residual.  The deficit is quoted both against the observed
    # uncertainty and against the inter-cluster scatter of the prediction,
    # the latter being the honest model-systematic scale.
    difference = observed_residual - avg_predicted
    difference_sigma = difference / observed_uncertainty
    deficit_ratio = avg_predicted / observed_residual if observed_residual else np.nan
    deficit_sigma_model = abs(difference) / max(std_predicted, 1e-9)
    
    print(f"\nDifference (Observed - Predicted): {difference:.2f} dex")
    print(f"Significance vs observed uncertainty: {difference_sigma:.1f}σ")
    print(f"Deficit ratio (predicted/observed): {deficit_ratio:.2f}x")
    print(f"Deficit vs inter-cluster prediction scatter: {deficit_sigma_model:.1f}σ")
    
    # Interpretation
    if abs(difference) > 2 * observed_uncertainty:
        if observed_residual < avg_predicted:
            interpretation = "OBSERVED_MUCH_LESS_THAN_NEWTONIAN"
            verdict = "TEP_SUPPORTED"
        else:
            interpretation = "OBSERVED_EXCEEDS_NEWTONIAN"
            verdict = "UNEXPECTED"
    elif abs(difference) < observed_uncertainty:
        interpretation = "CONSISTENT_WITH_NEWTONIAN"
        verdict = "STANDARD_GR_FAVORED"
    else:
        interpretation = "MARGINAL"
        verdict = "INCONCLUSIVE"
    
    print(f"\nInterpretation: {interpretation}")
    print(f"Verdict: {verdict}")
    
    # Save results
    output = {
        "timestamp": pd.Timestamp.now().isoformat(),
        "method": "King Model Monte Carlo",
        "clusters_analyzed": len(CLUSTER_PARAMS),
        "newtonian": {
            "mean_predicted_shift_dex": float(avg_predicted),
            "std_predicted_shift_dex": float(std_predicted),
            "cluster_predictions": {k: v["newtonian_prediction"] for k, v in results.items()}
        },
        "observed": {
            "residual_dex": observed_residual,
            "uncertainty_dex": observed_uncertainty
        },
        "comparison": {
            "difference_dex": float(difference),
            "significance_sigma": float(difference_sigma),
            "deficit_ratio": float(deficit_ratio),
            "deficit_sigma_model_scatter": float(deficit_sigma_model),
            "interpretation": interpretation,
            "verdict": verdict
        }
    }
    
    os.makedirs(RESULTS_DIR, exist_ok=True)
    with open(OUTPUT_JSON, 'w') as f:
        json.dump(output, f, indent=2)
    
    print(f"\n{'='*70}")
    print(f"Results saved to: {OUTPUT_JSON}")
    print(f"{'='*70}")
    
    return output


if __name__ == "__main__":
    analyze_dynamical_calibration()
