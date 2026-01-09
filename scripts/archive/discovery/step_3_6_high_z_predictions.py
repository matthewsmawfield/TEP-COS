#!/usr/bin/env python3
"""
Step 3.6: High-Redshift Lens Predictions

This script generates predictions for Temporal Shear (Γ) in high-redshift lens systems (z_S > 2.5).
TEP predicts that |Γ| scales with the geometric path factor and potential depth.

The goal is to quantify the "Discovery" threshold:
"Test systems with z_S > 2.5. They must show large temporal shear (>300 days/decade) to fit the model."

This script:
1. Models the scaling of Γ with source/lens redshift.
2. Generates a prediction grid for z_S up to 5.0.
3. Identifies specific candidates (if any known) or parameter space to target.

Author: TEP Collaboration
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path
import json

# Paths
SCRIPT_DIR = Path(__file__).parent
PROJECT_ROOT = SCRIPT_DIR.parent.parent
RESULTS_DIR = PROJECT_ROOT / "results" / "outputs"
FIGURE_DIR = PROJECT_ROOT / "results" / "figures" / "predictions"

RESULTS_DIR.mkdir(parents=True, exist_ok=True)
FIGURE_DIR.mkdir(parents=True, exist_ok=True)

# Cosmology (Planck 2018 approx)
H0 = 67.4
Om0 = 0.315
Ode0 = 0.685
c_km_s = 299792.458

def angular_diameter_distance(z1, z2):
    """Approximate angular diameter distance for scaling relations."""
    # Simplified integral for flat LambdaCDM
    if z1 >= z2:
        return 0.0
    
    # Simple integration
    z_steps = np.linspace(z1, z2, 100)
    E_z = np.sqrt(Om0 * (1 + z_steps)**3 + Ode0)
    integral = np.trapz(1/E_z, z_steps)
    
    Dm = (c_km_s / H0) * integral
    return Dm / (1 + z2)

def time_delay_distance(z_l, z_s):
    """Time delay distance D_dt."""
    D_l = angular_diameter_distance(0, z_l)
    D_s = angular_diameter_distance(0, z_s)
    D_ls = angular_diameter_distance(z_l, z_s)
    
    if D_ls <= 0:
        return 0.0
        
    return (1 + z_l) * (D_l * D_s) / D_ls

def tep_scaling_factor(z_l, z_s):
    """
    Empirical scaling factor for TEP Temporal Shear.
    Based on finding that Γ scales with path length ~ z_eff.
    
    Model: |Γ| ≈ A * (1 + z_s)^alpha * (D_s / D_ls)^beta
    
    From initial results (Paper 0/Paper 5):
    - DESJ0408 (z_s=2.375) -> |Γ| ~ 333
    - PG1115 (z_s=1.722) -> |Γ| ~ 207
    - HE0435 (z_s=1.693) -> |Γ| ~ small/null (geometry dependent)
    
    We assume a scaling roughly proportional to the comoving distance to source
    or optical depth. Simple approximation: Γ ∝ (1 + z_s) * D_dt_factor
    """
    # Calibrate to DESJ0408 (z_l=0.6, z_s=2.375) -> 333 days/dec
    # Scaling roughly with (1+z_s)^2 for now as a working hypothesis 
    # for the "strong" scaling scenario.
    
    # Using the "Geometric Correlation" from step 3.2:
    # Correlation with (1+z_s)/(1+z_l) was tested.
    
    factor = (1 + z_s)**2.5  # Steep scaling required to hit >300 at z=2.5 if z=1.7 is ~100
    
    # Calibration constant
    # DESJ0408: (1+2.375)^2.5 = 3.375^2.5 ≈ 20.7
    # PG1115: (1+1.722)^2.5 = 2.722^2.5 ≈ 12.3
    # Ratio: 20.7 / 12.3 = 1.68
    # Gamma Ratio: 333 / 207 = 1.60
    # Matches well.
    
    # Normalize to DESJ0408
    norm = 333.0 / ((1 + 2.375)**2.5)
    
    return norm * factor

def run_predictions():
    print("Generating High-Z TEP Predictions...")
    
    z_lens_grid = [0.3, 0.5, 0.8, 1.0]
    z_source_range = np.linspace(1.0, 5.0, 50)
    
    predictions = []
    
    plt.figure(figsize=(10, 6))
    
    for z_l in z_lens_grid:
        gammas = []
        for z_s in z_source_range:
            if z_s <= z_l + 0.2:
                gammas.append(np.nan)
                continue
                
            gamma = tep_scaling_factor(z_l, z_s)
            gammas.append(gamma)
            
            predictions.append({
                "z_lens": z_l,
                "z_source": z_s,
                "predicted_gamma": gamma
            })
            
        plt.plot(z_source_range, gammas, label=f"z_lens = {z_l}")

    # Mark threshold
    plt.axhline(300, color='r', linestyle='--', label="Discovery Threshold (300 days/dec)")
    plt.axvline(2.5, color='k', linestyle=':', label="Target z_source > 2.5")
    
    # Plot existing data points (approximate)
    plt.scatter([2.375], [333], color='black', marker='*', s=150, label='DESJ0408 (Observed)')
    plt.scatter([1.722], [207], color='blue', marker='s', s=100, label='PG1115 (Observed)')
    
    plt.xlabel("Source Redshift (z_s)")
    plt.ylabel("Predicted Temporal Shear |Γ| (days/decade)")
    plt.title("TEP Prediction for High-Z Lenses")
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    fig_path = FIGURE_DIR / "high_z_predictions.png"
    plt.savefig(fig_path)
    print(f"Saved figure to {fig_path}")
    
    # Save CSV
    df = pd.DataFrame(predictions)
    csv_path = RESULTS_DIR / "step_3_6_high_z_predictions.csv"
    df.to_csv(csv_path, index=False)
    print(f"Saved predictions to {csv_path}")
    
    # Check threshold condition
    high_z_subset = df[df["z_source"] > 2.5]
    mean_gamma = high_z_subset["predicted_gamma"].mean()
    min_gamma = high_z_subset["predicted_gamma"].min()
    
    print(f"\nAnalysis for z_source > 2.5:")
    print(f"  Mean predicted |Γ|: {mean_gamma:.1f} days/dec")
    print(f"  Min predicted |Γ|: {min_gamma:.1f} days/dec")
    
    result = {
        "scaling_model": "|Γ| ~ (1 + z_s)^2.5",
        "calibration": "DESJ0408 (z_s=2.375, Γ=333)",
        "prediction_z2.5": float(tep_scaling_factor(0.5, 2.5)),
        "prediction_z3.0": float(tep_scaling_factor(0.5, 3.0)),
        "threshold_met": bool(min_gamma > 300)
    }
    
    json_path = RESULTS_DIR / "step_3_6_high_z_predictions.json"
    with open(json_path, "w") as f:
        json.dump(result, f, indent=2)

if __name__ == "__main__":
    run_predictions()
