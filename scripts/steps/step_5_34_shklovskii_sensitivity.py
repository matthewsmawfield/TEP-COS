import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import json
import os
from scipy import stats

def run_cancellation_analysis():
    print("--- Step 5.34: Shklovskii Cancellation Analysis ---")
    
    # Hypothesis: Shklovskii effect (positive term) cancels Cluster Acceleration (negative term)
    # causing the observed suppressed density scaling.
    # Question: By what factor must Shklovskii be amplified to reduce the density slope
    # from 0.82 (Newtonian) to 0.39 (Observed)?
    
    # CORRECTED APPROACH: Combine quantities in LINEAR space, then measure log slope
    # When P_net = P_acc + K * P_shk, the log slope is NOT slope_acc + K * slope_shk
    # Instead, we must: (1) sum in linear space, (2) take log, (3) regress
    
    # Simulation Parameters
    n_clusters = 50
    # Log Density range: 2.3 to 5.8
    log_densities = np.linspace(2.3, 5.8, n_clusters)
    densities = 10**log_densities  # Convert to linear density
    
    # 1. Newtonian Acceleration (in linear space)
    # P_acc scales with cluster potential ~ rho^(0.82) in dex space
    # In linear space: P_acc ~ rho^0.82, but we need absolute values
    # Normalize: at log_rho = 4 (rho = 10^4), set P_acc = 1.0 as reference
    ref_log_density = 4.0
    P_acc_ref = 1.0
    slope_acc_dex = 0.82
    # P_acc = P_acc_ref * (rho / rho_ref)^slope_acc_dex
    P_acc = P_acc_ref * (densities / 10**ref_log_density)**slope_acc_dex
    
    # 2. Shklovskii Effect (in linear space)
    # Shklovskii: P_shk ~ v^2/D where v^2 ~ sigma_v^2 ~ M/R (virialized)
    # For homologous systems: M/R ~ rho^(2/3), giving sigma_v^2 ~ rho^(2/3) ~ rho^0.67
    # However, observed velocity dispersion profiles show weaker scaling in cluster cores.
    # Empirically, Shklovskii has shallower density dependence than acceleration.
    # We model this as: P_shk ~ rho^0.50 (flatter than P_acc ~ rho^0.82)
    
    slope_shk_dex = 0.50  # Flatter density scaling than acceleration
    base_shk_fraction = 0.15  # Shklovskii is 15% of acceleration at reference point
    # P_shk has different density scaling than P_acc
    P_shk_base = base_shk_fraction * P_acc_ref * (densities / 10**ref_log_density)**slope_shk_dex
    
    print(f"Base Acceleration Slope (input):       {slope_acc_dex:.2f}")
    print(f"Base Shklovskii Slope (input):         {slope_shk_dex:.2f}")
    print(f"Base Shklovskii Fraction at rho_ref:   {base_shk_fraction*100:.1f}%")
    print(f"Target Net Slope:                      0.39 (Mixed-Effects Observed)")
    
    # 3. Find amplification factor K such that combined slope = 0.39
    # P_net = P_acc + K * P_shk
    # Measure slope of log(P_net) vs log(density)
    
    target_slope = 0.39
    
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
    
    print(f"Measured Acc-only slope (verification): {measured_slopes[0]:.3f}")
    print(f"Required Shklovskii Amplification K:    {required_K:.2f}")
    print(f"Resulting net slope at K={required_K:.2f}:     {final_slope:.3f}")
    
    # 4. Physical plausibility assessment
    # Shklovskii = 2.43e-21 * P * mu^2 * D
    # To get amplification K, we need:
    # - Distance error factor: D_true = K * D_catalog (since Shklovskii ~ D)
    # - Proper motion error factor: mu_true = sqrt(K) * mu_catalog (since Shklovskii ~ mu^2)
    
    distance_error_factor = required_K
    pm_error_factor = np.sqrt(required_K)
    
    # Physical constraints
    typical_distance_uncertainty = 0.10  # 10% for GCs (Baumgardt et al.)
    typical_pm_uncertainty = 0.02 / 5.0  # Gaia EDR3: ~0.02 mas/yr vs typical 5 mas/yr
    
    conclusion = ""
    if required_K > 5.0:
        conclusion = f"Physically excluded. Requires distances off by {distance_error_factor:.1f}x (typical error ~10%) or proper motions off by {pm_error_factor:.1f}x (Gaia precision <1%)."
    elif required_K > 2.0:
        conclusion = f"Highly implausible. Requires distances off by {distance_error_factor:.1f}x or proper motions off by {pm_error_factor:.1f}x."
    else:
        conclusion = "Plausible systematic."
        
    print(conclusion)
    
    # Save results
    results = {
        "acc_slope": slope_acc_dex,
        "shk_slope": slope_shk_dex,
        "observed_slope": 0.39,
        "base_shk_fraction": base_shk_fraction,
        "required_amplification_factor": float(required_K),
        "required_distance_error_factor": float(required_K),
        "required_pm_error_factor": float(np.sqrt(required_K)),
        "conclusion": conclusion
    }
    
    with open("results/outputs/step_5_34_sensitivity_results.json", "w") as f:
        json.dump(results, f, indent=4)
        
    # Plot
    plt.figure(figsize=(8,5))
    plt.plot(amplification_factors, measured_slopes, label='Net Density Slope')
    plt.axhline(0.39, color='r', linestyle='--', label='Observed Slope (0.39)')
    plt.axhline(0.82, color='g', linestyle='--', label='Newtonian Prediction (0.82)')
    plt.axvline(required_K, color='k', linestyle=':', label=f'Required Amp ({required_K:.1f}x)')
    plt.xlabel('Shklovskii Amplification Factor')
    plt.ylabel('Resulting Density Slope')
    plt.title('Cancellation Analysis: Can Shklovskii Explain Suppressed Scaling?')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    os.makedirs('site/figures', exist_ok=True)
    plt.savefig('site/figures/shklovskii_cancellation.png')

if __name__ == "__main__":
    run_cancellation_analysis()
