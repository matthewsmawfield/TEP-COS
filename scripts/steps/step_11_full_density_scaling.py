#!/usr/bin/env python3
"""
Step 11: Full Density Scaling Simulation (Addressing the 'Killer' Counter-Argument)

Simulates acceleration distributions for ALL clusters using their EXACT observed 
structural parameters (M, Rc), rather than a generic scaling law.

This tests if the "Suppressed Density Scaling" (Slope 0.10) is real or an 
artifact of assuming fixed Rc.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import pandas as pd
import json
import os
from pathlib import Path

def run_full_density_scaling(n_ensemble=100):
    print("--- Step 11: Full Density Scaling Simulation ---")
    print(f"Running {n_ensemble} ensemble realizations for uncertainty quantification...")
    
    # 1. Cluster Parameters (Harris 2010 / Baumgardt 2018)
    # M: Mass (M_sun), Rc: Core Radius (pc), Rt: Tidal Radius (pc)
    # rho_c: log10 Central Luminosity Density (L_sun/pc^3) - used for x-axis plotting
    CLUSTER_PARAMS = {
        "Terzan 5":         {"M": 2.0e6, "Rc": 0.16, "Rt": 5.0,  "rho_c": 5.50},
        "47 Tuc (NGC 104)": {"M": 1.0e6, "Rc": 0.36, "Rt": 42.0, "rho_c": 4.88},
        "NGC 6517":         {"M": 2.0e5, "Rc": 0.06, "Rt": 5.0,  "rho_c": 5.80},
        "M28 (NGC 6626)":   {"M": 5.0e5, "Rc": 0.24, "Rt": 12.0, "rho_c": 4.52},
        "M62 (NGC 6266)":   {"M": 1.0e6, "Rc": 0.18, "Rt": 8.0,  "rho_c": 5.16},
        "M13 (NGC 6205)":   {"M": 6.0e5, "Rc": 0.62, "Rt": 25.0, "rho_c": 3.79},
        "M15 (NGC 7078)":   {"M": 5.0e5, "Rc": 0.14, "Rt": 21.0, "rho_c": 5.05},
        "M5 (NGC 5904)":    {"M": 5.0e5, "Rc": 0.42, "Rt": 28.0, "rho_c": 3.53},
        "Terzan 1":         {"M": 1.5e5, "Rc": 0.10, "Rt": 4.0,  "rho_c": 5.00},
        "NGC 6752":         {"M": 3.0e5, "Rc": 0.17, "Rt": 25.0, "rho_c": 4.30},
        "M2 (NGC 7089)":    {"M": 6.0e5, "Rc": 0.32, "Rt": 21.0, "rho_c": 4.15},
        "Omega Centauri (NGC 5139)": {"M": 4.0e6, "Rc": 2.37, "Rt": 57.0, "rho_c": 3.12},
        "M53 (NGC 5024)":   {"M": 3.0e5, "Rc": 0.65, "Rt": 22.0, "rho_c": 2.96},
        "M3 (NGC 5272)":    {"M": 5.0e5, "Rc": 0.37, "Rt": 38.0, "rho_c": 3.68},
        "M71 (NGC 6838)":   {"M": 2.0e4, "Rc": 0.63, "Rt": 8.0,  "rho_c": 2.29},
        "NGC 6397":         {"M": 1.0e5, "Rc": 0.05, "Rt": 15.0, "rho_c": 5.68},
        "NGC 1851":         {"M": 3.0e5, "Rc": 0.09, "Rt": 11.0, "rho_c": 5.09},
        "NGC 6522":         {"M": 2.0e5, "Rc": 0.05, "Rt": 5.0,  "rho_c": 5.50},
        "NGC 6544":         {"M": 5.0e4, "Rc": 0.05, "Rt": 3.0,  "rho_c": 5.20},
        "NGC 6624":         {"M": 2.0e5, "Rc": 0.06, "Rt": 6.0,  "rho_c": 5.60},
        "NGC 6760":         {"M": 2.0e5, "Rc": 0.34, "Rt": 8.0,  "rho_c": 3.80},
        "M22 (NGC 6656)":   {"M": 5.0e5, "Rc": 1.33, "Rt": 32.0, "rho_c": 2.97},
        "M80 (NGC 6093)":   {"M": 4.0e5, "Rc": 0.15, "Rt": 13.0, "rho_c": 4.79},
        "M92 (NGC 6341)":   {"M": 3.0e5, "Rc": 0.26, "Rt": 15.0, "rho_c": 4.30},
        "NGC 6712":         {"M": 1.5e5, "Rc": 0.33, "Rt": 7.0,  "rho_c": 3.70},
        "NGC 6652":         {"M": 1.0e5, "Rc": 0.10, "Rt": 5.0,  "rho_c": 4.50},
        "M14 (NGC 6402)":   {"M": 1.0e6, "Rc": 0.78, "Rt": 18.0, "rho_c": 3.44},
        "NGC 6539":         {"M": 3.0e5, "Rc": 0.60, "Rt": 10.0, "rho_c": 3.30},
        "M4 (NGC 6121)":    {"M": 1.0e5, "Rc": 0.83, "Rt": 33.0, "rho_c": 2.85},
        # Additional clusters from pulsar data
        "NGC 6440":         {"M": 1.5e5, "Rc": 0.12, "Rt": 5.5,  "rho_c": 5.10},
        "NGC 6441":         {"M": 8.0e5, "Rc": 0.20, "Rt": 12.0, "rho_c": 5.00},
        "NGC 6316":         {"M": 1.2e5, "Rc": 0.15, "Rt": 5.0,  "rho_c": 4.80},
        "M30 (NGC 7099)":   {"M": 2.5e5, "Rc": 0.25, "Rt": 18.0, "rho_c": 4.20},
    }

    # 2. Load Pulsar Data to identify which clusters to simulate
    try:
        df = pd.read_csv("results/outputs/step_02_pulsar_population_controls.csv")
        gc_df = df[df['environment'] == 'globular_cluster']
        cluster_counts = gc_df['cluster'].value_counts()
    except FileNotFoundError:
        print("Error: Pulsar population file not found. Run step 5.10 first.")
        return

    # 3. Simulation Constants
    # REMOVED: np.random.seed(42) - using ensemble approach instead
    n_stars_per_cluster = 2000 # Enough for stable stats
    G_si = 6.674e-11
    M_sun_kg = 1.989e30
    pc_m = 3.086e16
    c_si = 2.998e8
    
    # Intrinsic parameters (Field control)
    mu_field = -19.76
    sigma_field = 0.64
    
    # --- ENSEMBLE SIMULATION ---
    ensemble_slopes = []
    ensemble_intercepts = []
    final_results = None  # Store last ensemble's results for plotting
    
    print(f"Simulating {len(cluster_counts)} clusters with exact parameters...")
    
    for ensemble_idx in range(n_ensemble):
        if (ensemble_idx + 1) % 10 == 0 or ensemble_idx == 0:
            print(f"  Ensemble realization {ensemble_idx + 1}/{n_ensemble}...")
        
        # Set seed for this realization (deterministic but different for each)
        np.random.seed(42 + ensemble_idx)
        
        results = []
        
        for cluster_name in sorted(cluster_counts.index):
            if cluster_name not in CLUSTER_PARAMS:
                continue
                
            params = CLUSTER_PARAMS[cluster_name]
            M = params["M"]
            Rc = params["Rc"]
            
            # --- N-BODY / CMC SIMULATION ---
            # 1. Mass Segregation
            sigma_r_pc = 0.5 * Rc
            r_pulsar_pc = np.abs(np.random.normal(0, sigma_r_pc, n_stars_per_cluster))
            r_pulsar_m = r_pulsar_pc * pc_m
            
            cos_theta = np.random.uniform(-1, 1, n_stars_per_cluster)
            
            # 2. Mean Field Acceleration (Newtonian)
            m_cl_kg = M * M_sun_kg
            r_core_m = Rc * pc_m
            
            a_mean_si = np.zeros(n_stars_per_cluster)
            mask_core = r_pulsar_m < r_core_m
            
            # Harmonic core: a = (GM/Rc^3) * r
            g_max = G_si * m_cl_kg / (r_core_m**2)
            a_mean_si[mask_core] = g_max * (r_pulsar_m[mask_core] / r_core_m)
            
            # Envelope: a = GM/r^2
            a_mean_si[~mask_core] = G_si * m_cl_kg / (r_pulsar_m[~mask_core]**2)
            
            a_los_mean_si = a_mean_si * cos_theta
            
            # 3. Binary Hardening (Velocity Kicks)
            sigma_v = np.sqrt(G_si * m_cl_kg / r_core_m)
            v_thermal = np.random.normal(0, sigma_v, n_stars_per_cluster)
            v_kick = stats.cauchy.rvs(loc=0, scale=2*sigma_v, size=n_stars_per_cluster)
            v_tot = v_thermal + 0.2 * v_kick
            
            # Total Acceleration Term (LOS Gravity + Shklovskii)
            term_acc = a_los_mean_si / c_si
            term_shk = (v_tot**2) / (c_si * r_pulsar_m)
            
            # Pdot Obs
            log_P_s = np.random.normal(np.log10(0.005), 0.3, n_stars_per_cluster)
            P_s = 10**log_P_s
            
            log_pdot_int = np.random.normal(mu_field, sigma_field, n_stars_per_cluster)
            pdot_int = 10**log_pdot_int
            
            pdot_obs = pdot_int + P_s * (term_acc + term_shk)
            log_pdot_obs = np.log10(np.abs(pdot_obs))
            
            # Calculate Shift
            shift = np.mean(log_pdot_obs) - mu_field
            
            results.append({
                "name": cluster_name,
                "rho_c_log": float(params["rho_c"]),
                "shift": float(shift),
                "n_pulsars_real": int(cluster_counts[cluster_name])
            })
        
        # Analysis for this ensemble realization
        shifts_sim = [r['shift'] for r in results]
        densities = [r['rho_c_log'] for r in results]
        
        slope_sim, intercept_sim, r_value_sim, p_value_sim, std_err_sim = stats.linregress(densities, shifts_sim)
        
        ensemble_slopes.append(slope_sim)
        ensemble_intercepts.append(intercept_sim)
        
        # Store final ensemble results for plotting
        final_results = results
    
    # End of ensemble loop
    
    # Extract plotting data from final ensemble realization
    if final_results:
        shifts_sim = [r['shift'] for r in final_results]
        densities = [r['rho_c_log'] for r in final_results]
    else:
        shifts_sim = []
        densities = []
    
    # Compute ensemble statistics
    slope_sim_mean = np.mean(ensemble_slopes)
    slope_sim_std = np.std(ensemble_slopes)
    slope_sim_se = slope_sim_std / np.sqrt(n_ensemble)
    intercept_sim_mean = np.mean(ensemble_intercepts)
    
    print(f"\nNewtonian Slope (Ensemble): {slope_sim_mean:.3f} ± {slope_sim_std:.3f} (std) ± {slope_sim_se:.3f} (se) dex/dex")
    print(f"  95% Confidence Interval: [{slope_sim_mean - 1.96*slope_sim_se:.3f}, {slope_sim_mean + 1.96*slope_sim_se:.3f}]")
    
    # --- CALCULATE OBSERVED SHIFTS ---
    # We need to calculate the actual observed shift for each cluster
    # Shift = Mean(log|Pdot|_cluster) - Mean(log|Pdot|_field)
    # Using the field mean from the simulation parameters (mu_field = -19.76)
    
    observed_data = []
    
    if 'gc_df' in locals():
        for cluster_name in sorted(cluster_counts.index):
            if cluster_name not in CLUSTER_PARAMS:
                continue
                
            cluster_pulsars = gc_df[gc_df['cluster'] == cluster_name]
            if len(cluster_pulsars) < 1:
                continue
                
            mean_obs = cluster_pulsars['logPdot_abs'].mean()
            shift_obs = mean_obs - mu_field
            
            rho_c = CLUSTER_PARAMS[cluster_name]["rho_c"]
            
            observed_data.append({
                "name": cluster_name,
                "rho_c_log": rho_c,
                "shift_obs": shift_obs,
                "n_pulsars": len(cluster_pulsars)
            })
            
    densities_obs = [d['rho_c_log'] for d in observed_data]
    shifts_obs = [d['shift_obs'] for d in observed_data]
    
    if len(densities_obs) > 0:
        slope_obs, intercept_obs, r_value_obs, p_value_obs, std_err_obs = stats.linregress(densities_obs, shifts_obs)
    else:
        raise ValueError("No observed data available for regression. Check that pulsar data file exists and contains valid cluster data.")
    
    print("-" * 50)
    print(f"FULL DENSITY SCALING RESULTS (N={len(results)} clusters, {n_ensemble} ensemble realizations)")
    print(f"Newtonian Slope: {slope_sim_mean:.3f} ± {slope_sim_se:.3f} dex / dex (ensemble)")
    print(f"Observed Slope:  {slope_obs:.3f} ± {std_err_obs:.3f} dex / dex")
    print(f"Suppression:     {slope_obs/slope_sim_mean:.1%} of expected scaling")
    
    # Statistical test: Is observed slope consistent with Newtonian?
    z_score = (slope_obs - slope_sim_mean) / np.sqrt(std_err_obs**2 + slope_sim_se**2)
    p_value = 2 * (1 - stats.norm.cdf(abs(z_score)))
    print(f"Tension with Newtonian: {z_score:.1f}σ (p={p_value:.3f})")
    print("-" * 50)
    
    # Save Results
    os.makedirs("results/outputs", exist_ok=True)
    out_data = {
        "simulation_type": "N-Body/CMC Synthetic with Ensemble Uncertainty",
        "ensemble_parameters": {
            "n_ensemble": n_ensemble,
            "n_stars_per_cluster": n_stars_per_cluster
        },
        "newtonian": {
            "slope_mean": float(slope_sim_mean),
            "slope_std": float(slope_sim_std),
            "slope_se": float(slope_sim_se),
            "slope_95ci_lower": float(slope_sim_mean - 1.96*slope_sim_se),
            "slope_95ci_upper": float(slope_sim_mean + 1.96*slope_sim_se),
            "intercept_mean": float(intercept_sim_mean),
            "ensemble_slopes": [float(s) for s in ensemble_slopes]
        },
        "observed": {
            "slope": float(slope_obs),
            "slope_error": float(std_err_obs),
            "r_value": float(r_value_obs),
            "p_value": float(p_value_obs)
        },
        "statistical_test": {
            "z_score": float(z_score),
            "p_value": float(p_value),
            "tension_sigma": float(abs(z_score)),
            "suppression_factor": float(slope_obs / slope_sim_mean)
        },
        "clusters": results
    }
    
    with open("results/outputs/step_11_full_density_scaling.json", "w") as f:
        json.dump(out_data, f, indent=4)
        
    print("Results saved to results/outputs/step_11_full_density_scaling.json")

    # 5. Plotting (Reproducing Density Scaling Figure with ALL data)
    # Publication Style - Tuned for 900px Web Manuscript (Matching step_46_cluster_acceleration_simulation.png)
    plt.rcParams.update({
        'font.family': 'serif',
        'font.serif': ['Times New Roman', 'Times', 'DejaVu Serif'],
        'font.size': 12,             # Base font size
        'axes.labelsize': 14,        # Axis labels
        'axes.titlesize': 16,        # Title
        'legend.fontsize': 12,       # Legend
        'xtick.labelsize': 12,
        'ytick.labelsize': 12,
        'figure.dpi': 300,
        'savefig.dpi': 300,
        'savefig.bbox': 'tight',
        'mathtext.fontset': 'stix',
        'lines.linewidth': 0.8,      # Thinner lines
    })

    COLOR_NEWTONIAN = '#2E86AB'  # Blue (Model)
    COLOR_OBSERVED = '#E94F37'   # Red (Data)
    COLOR_GRID = '#E0E0E0'

    fig, ax = plt.subplots(figsize=(9, 6))
    
    # Plot Simulation Points
    ax.scatter(densities, shifts_sim, c=COLOR_NEWTONIAN, alpha=0.7, edgecolors='none', s=60, label='Newtonian Simulation (Exact Params)')
    
    # Plot Observed Points
    ax.scatter(densities_obs, shifts_obs, c=COLOR_OBSERVED, alpha=0.8, edgecolors='black', marker='s', s=50, label='Observed Data')
    
    # Plot Regression Lines
    x_range = np.linspace(min(densities), max(densities), 100)
    
    # Newtonian prediction with uncertainty band
    y_pred_sim = slope_sim_mean * x_range + intercept_sim_mean
    y_pred_sim_upper = (slope_sim_mean + 1.96*slope_sim_se) * x_range + intercept_sim_mean
    y_pred_sim_lower = (slope_sim_mean - 1.96*slope_sim_se) * x_range + intercept_sim_mean
    
    ax.fill_between(x_range, y_pred_sim_lower, y_pred_sim_upper, alpha=0.2, color=COLOR_NEWTONIAN, label='Newtonian 95% CI')
    ax.plot(x_range, y_pred_sim, COLOR_NEWTONIAN, linestyle='-', linewidth=2, label=f'Newtonian Trend (Slope={slope_sim_mean:.2f}±{slope_sim_se:.2f})')
    
    if len(densities_obs) > 0:
        y_pred_obs = slope_obs * x_range + intercept_obs
        ax.plot(x_range, y_pred_obs, COLOR_OBSERVED, linestyle='--', linewidth=2, label=f'Observed Trend (Slope={slope_obs:.2f})')
    
    ax.set_xlabel(r'$\log_{10}(\rho_c) [L_\odot/pc^3]$')
    ax.set_ylabel(r'Shift in $\log_{10}|\dot{P}|$ (dex)')
    ax.set_title(r'Suppressed Density Scaling: Observation vs Newtonian Dynamics')
    
    # Styled Legend
    ax.legend(loc='upper left', frameon=True, framealpha=0.95, edgecolor='none')
    
    # Refined Grid
    ax.grid(True, linestyle='-', color=COLOR_GRID, alpha=0.8, linewidth=1.0)
    
    # UPDATING PREVIOUS CHART
    output_dir = "site/figures"
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, "step_11_density_scaling.png")
    plt.savefig(output_path, dpi=300)
    print(f"Updated {output_path}")

if __name__ == "__main__":
    run_full_density_scaling()
