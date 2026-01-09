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
    # Question: By what factor must Shklovskii be amplified (e.g. underestimated distance/proper motion)
    # to reduce the density slope from 0.82 (Newtonian) to 0.35 (Observed)?
    
    # Simulation Parameters
    n_clusters = 50
    # Log Density range: 2.3 to 5.8
    densities = np.linspace(2.3, 5.8, n_clusters)
    
    # 1. Newtonian Acceleration Scaling (Cluster Potential)
    # Scales as rho (approx 0.82 slope in dex space)
    # This increases |Pdot| (shifts it positive in dex)
    slope_acc = 0.82
    intercept_acc = -0.5 # Normalized so dense clusters have high shift
    shift_acc = slope_acc * densities + intercept_acc
    
    # 2. Shklovskii Scaling
    # Shklovskii term P_shk ~ v^2/D
    # v ~ sigma_v ~ sqrt(M/R) ~ rho^...
    # Empirically, dense clusters have higher velocity dispersion.
    # If sigma_v^2 scales with potential depth, and potential depth scales with density...
    # Let's assume Shklovskii also scales with density, but acts to REDUCE |Pdot|.
    # (Since Shklovskii is positive and Pdot_int is negative).
    
    # In dex space, adding a positive term to a negative value reduces its absolute magnitude.
    # So Shklovskii contributes a NEGATIVE slope to log|Pdot|.
    
    # Estimate Shklovskii slope relative to Acc slope
    # In standard models, Shklovskii is sub-dominant (~10% of Acc in cores).
    # So initial slope_shk ~ -0.08 (approx 10% of 0.82)
    base_shk_fraction = 0.15 
    slope_shk_base = -1 * base_shk_fraction * slope_acc
    
    print(f"Base Acceleration Slope: {slope_acc:.2f}")
    print(f"Base Shklovskii Slope:   {slope_shk_base:.2f} (Assumed {base_shk_fraction*100}% of Acc)")
    print(f"Target Net Slope:        0.35 (Mixed-Effects Observed)")
    
    # 3. Solve for Amplification Factor K
    # Net Slope = Slope_Acc + K * Slope_Shk_Base
    # 0.35 = 0.82 + K * (-0.12)
    # K * 0.12 = 0.82 - 0.35 = 0.47
    # K = 0.47 / 0.12 ~ 3.9
    
    amplification_factors = np.linspace(0, 10, 100)
    net_slopes = []
    
    for K in amplification_factors:
        # Net shift = Acc Shift + K * Shk Shift
        # Note: We simulate this by combining slopes directly for robustness
        net_slope = slope_acc + K * slope_shk_base
        net_slopes.append(net_slope)
        
    # Find crossing
    net_slopes = np.array(net_slopes)
    idx = np.argmin(np.abs(net_slopes - 0.35))
    required_K = amplification_factors[idx]
    
    print(f"Required Shklovskii Amplification K: {required_K:.2f}")
    
    # Validate with physics
    # Shklovskii = 2.43e-21 * P * mu^2 * D
    # To get K=3.9, we need D * mu^2 to be 3.9x larger.
    # If D is error source: D_true = 3.9 * D_catalog? Or D_true = D_cat / 3.9?
    # Shklovskii proportional to D. So D would need to be 3.9x LARGER.
    # OR mu would need to be sqrt(3.9) ~ 2.0x LARGER.
    
    # Implication:
    # Are distances to Globular Clusters underestimated by factor 3.9?
    # Typical uncertainty is 5-10%. 390% is impossible.
    # Are proper motions underestimated by factor 2.0?
    # Gaia EDR3 precision is ~0.02 mas/yr. PMs are ~5 mas/yr. Error is <1%.
    
    conclusion = ""
    if required_K > 2.0:
        conclusion = f"Physically excluded. Requires distances to be off by factor {required_K:.1f}x or proper motions by {np.sqrt(required_K):.1f}x."
    else:
        conclusion = "Plausible systematic."
        
    print(conclusion)
    
    # Save results
    results = {
        "acc_slope": slope_acc,
        "observed_slope": 0.35,
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
    plt.plot(amplification_factors, net_slopes, label='Net Density Slope')
    plt.axhline(0.35, color='r', linestyle='--', label='Observed Slope (0.35)')
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
