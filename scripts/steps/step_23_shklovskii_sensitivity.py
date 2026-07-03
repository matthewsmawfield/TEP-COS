#!/usr/bin/env python3
"""
Step 23: Shklovskii Cancellation Sensitivity Analysis

Tests whether Shklovskii effect amplification could explain the observed
suppressed density scaling by simulating cancellation scenarios.

Author: TEP-COS Analysis Pipeline
Date: 2026
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import json
import os
from pathlib import Path
from scipy import stats


def load_upstream_values():
    """Load dynamically computed values from upstream pipeline outputs."""
    results_dir = Path("results/outputs")

    # Observed slope from hierarchical density scaling (step_12)
    hier_file = results_dir / "step_12_hierarchical_density_results.json"
    if hier_file.exists():
        with open(hier_file, 'r') as f:
            hier_data = json.load(f)
        observed_slope = hier_data.get('model_b_mixed_slope', 0.393)
    else:
        raise FileNotFoundError(f"Required hierarchical input missing: {hier_file}")

    # Newtonian slope from CMC literature consensus (step_14)
    cmc_file = results_dir / "step_14_cmc_literature.json"
    if cmc_file.exists():
        with open(cmc_file, 'r') as f:
            cmc_data = json.load(f)
        newtonian_slope = cmc_data.get('cmc_consensus', {}).get('weighted_mean', 0.748)
    else:
        raise FileNotFoundError(f"Required CMC input missing: {cmc_file}")

    return observed_slope, newtonian_slope


def run_cancellation_analysis():
    print("--- Step 23: Shklovskii Cancellation Analysis ---")

    # Load upstream values dynamically
    observed_slope, newtonian_slope = load_upstream_values()

    # Hypothesis: Shklovskii effect (positive term) cancels Cluster Acceleration (negative term)
    # causing the observed suppressed density scaling.
    # Question: By what factor must Shklovskii be amplified to reduce the density slope
    # from Newtonian prediction to observed?

    # CORRECTED APPROACH: Combine quantities in LINEAR space, then measure log slope
    # When P_net = P_acc + K * P_shk, the log slope is NOT slope_acc + K * slope_shk
    # Instead, we must: (1) sum in linear space, (2) take log, (3) regress

    # Simulation Parameters
    n_clusters = 50
    # Log Density range: 2.3 to 5.8
    log_densities = np.linspace(2.3, 5.8, n_clusters)
    densities = 10**log_densities  # Convert to linear density

    # 1. Newtonian Acceleration (in linear space)
    # P_acc scales with cluster potential ~ rho^newtonian_slope in dex space
    # In linear space: P_acc ~ rho^newtonian_slope, but we need absolute values
    # Normalize: at log_rho = 4 (rho = 10^4), set P_acc = 1.0 as reference
    ref_log_density = 4.0
    P_acc_ref = 1.0
    slope_acc_dex = newtonian_slope
    # P_acc = P_acc_ref * (rho / rho_ref)^slope_acc_dex
    P_acc = P_acc_ref * (densities / 10**ref_log_density)**slope_acc_dex

    # 2. Shklovskii Effect (in linear space)
    # Shklovskii: P_shk ~ v^2/D where v^2 ~ sigma_v^2 ~ M/R (virialized)
    # For homologous systems: M/R ~ rho^(2/3), giving sigma_v^2 ~ rho^(2/3) ~ rho^0.67
    # However, observed velocity dispersion profiles show weaker scaling in cluster cores.
    # Empirically, Shklovskii has shallower density dependence than acceleration.
    # We model this as: P_shk ~ rho^0.50 (flatter than P_acc)

    slope_shk_dex = 0.50  # Flatter density scaling than acceleration
    base_shk_fraction = 0.15  # Shklovskii is 15% of acceleration at reference point
    # P_shk has different density scaling than P_acc
    P_shk_base = base_shk_fraction * P_acc_ref * (densities / 10**ref_log_density)**slope_shk_dex

    print(f"Base Acceleration Slope (input):       {slope_acc_dex:.3f} (from CMC consensus)")
    print(f"Base Shklovskii Slope (input):         {slope_shk_dex:.2f}")
    print(f"Base Shklovskii Fraction at rho_ref:   {base_shk_fraction*100:.1f}%")
    print(f"Target Net Slope:                      {observed_slope:.3f} (Mixed-Effects Observed)")

    # 3. Find amplification factor K such that combined slope = observed_slope
    # P_net = P_acc + K * P_shk
    # Measure slope of log(P_net) vs log(density)

    target_slope = observed_slope
    
    def measure_combined_slope(K):
        """Combine in linear space, then measure log slope via OLS."""
        P_net = P_acc + K * P_shk_base
        # Take logarithms for slope measurement
        log_P_net = np.log10(P_net)
        # OLS regression: log_P vs log_density
        slope, intercept, r_value, p_value, std_err = stats.linregress(log_densities, log_P_net)
        return slope, r_value**2
    
    # Search for K
    amplification_factors = np.linspace(0, 20, 200)
    measured_slopes = []
    
    for K in amplification_factors:
        slope, r2 = measure_combined_slope(K)
        measured_slopes.append(slope)
    
    measured_slopes = np.array(measured_slopes)
    
    # Find crossing with target slope
    # Slope decreases as K increases (Shklovskii adds flat component)
    # Find where measured slope = target_slope
    idx = np.argmin(np.abs(measured_slopes - target_slope))
    required_K = amplification_factors[idx]
    final_slope = measured_slopes[idx]

    # If the closest slope within the search range is still far from target,
    # the required amplification exceeds our search boundary.
    exceeds_range = abs(final_slope - target_slope) > 0.05

    print(f"Measured Acc-only slope (verification): {measured_slopes[0]:.3f}")
    if exceeds_range:
        print(f"Required Shklovskii Amplification K:    >{amplification_factors.max():.1f} (exceeds search range)")
        print(f"Closest net slope within range:           {final_slope:.3f} (target: {target_slope:.3f})")
    else:
        print(f"Required Shklovskii Amplification K:    {required_K:.2f}")
        print(f"Resulting net slope at K={required_K:.2f}:     {final_slope:.3f}")

    # 4. Physical plausibility assessment
    # Shklovskii = 2.43e-21 * P * mu^2 * D
    # To get amplification K, we need:
    # - Distance error factor: D_true = K * D_catalog (since Shklovskii ~ D)
    # - Proper motion error factor: mu_true = sqrt(K) * mu_catalog (since Shklovskii ~ mu^2)

    distance_error_factor = required_K
    pm_error_factor = np.sqrt(required_K)

    # Physical constraints (literature references for comparison)
    # Distance: Baumgardt et al. (2021) reports GC distance uncertainties ~5-15%
    typical_distance_uncertainty = 0.10  # Conservative 10% representative value
    # Proper motion: Gaia EDR3 astrometric precision ~0.02 mas/yr; typical GC pulsar PM ~5 mas/yr
    typical_pm_uncertainty = 0.02 / 5.0  # ~0.4% relative PM uncertainty

    conclusion = ""
    if exceeds_range:
        conclusion = f"Physically excluded. Required amplification exceeds {amplification_factors.max():.0f}x. Even at K={amplification_factors.max():.0f}x the net slope is {final_slope:.2f}, still far from observed {target_slope:.2f}. Shklovskii cannot explain the suppressed scaling."
    elif required_K > 5.0:
        conclusion = f"Physically excluded. Requires distances off by {distance_error_factor:.1f}x (typical error ~10%) or proper motions off by {pm_error_factor:.1f}x (Gaia precision <1%)."
    elif required_K > 2.0:
        conclusion = f"Highly implausible. Requires distances off by {distance_error_factor:.1f}x or proper motions off by {pm_error_factor:.1f}x."
    else:
        conclusion = "Plausible systematic."

    print(conclusion)

    # Save results
    results = {
        "acc_slope": float(slope_acc_dex),
        "shk_slope": float(slope_shk_dex),
        "observed_slope": float(observed_slope),
        "newtonian_slope": float(newtonian_slope),
        "base_shk_fraction": float(base_shk_fraction),
        "required_amplification_factor": float(required_K),
        "amplification_exceeds_range": bool(exceeds_range),
        "required_distance_error_factor": float(required_K),
        "required_pm_error_factor": float(np.sqrt(required_K)),
        "conclusion": conclusion
    }

    with open("results/outputs/step_23_sensitivity_results.json", "w") as f:
        json.dump(results, f, indent=4)

    # Plot
    plt.figure(figsize=(8,5))
    plt.plot(amplification_factors, measured_slopes, label='Net Density Slope')
    plt.axhline(observed_slope, color='r', linestyle='--', label=f'Observed Slope ({observed_slope:.3f})')
    plt.axhline(newtonian_slope, color='g', linestyle='--', label=f'Newtonian Prediction ({newtonian_slope:.3f})')
    plt.axvline(required_K, color='k', linestyle=':', label=f'Required Amp ({required_K:.1f}x)')
    plt.xlabel('Shklovskii Amplification Factor')
    plt.ylabel('Resulting Density Slope')
    plt.title('Cancellation Analysis: Can Shklovskii Explain Suppressed Scaling?')
    plt.legend()
    plt.grid(True, alpha=0.3)

    os.makedirs('site/figures', exist_ok=True)
    plt.savefig('site/figures/step_23_shklovskii_cancellation.png')

if __name__ == "__main__":
    run_cancellation_analysis()
