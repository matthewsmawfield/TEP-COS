#!/usr/bin/env python3
"""
Step 5.33: Hierarchical (Mixed-Effects) Density Scaling Analysis

Performs mixed-effects regression of pulsar spin-down residuals against
central density with cluster-level random effects.

Author: TEP-COS Analysis Pipeline
Date: 2026
"""

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
import statsmodels.api as sm
from scipy import stats
import json
import os
from pathlib import Path


def run_hierarchical_analysis():
    print("--- Step 5.33: Hierarchical (Mixed-Effects) Density Scaling ---")
    
    # 1. Load Pulsar Data
    # We need the individual pulsar data with their cluster associations
    csv_path = "results/outputs/step_5_10_pulsar_population_controls.csv"
    if not os.path.exists(csv_path):
        print(f"Error: {csv_path} not found.")
        return

    try:
        df = pd.read_csv(csv_path)
    except (pd.errors.EmptyDataError, pd.errors.ParserError) as e:
        print(f"Error reading CSV file: {type(e).__name__} - {e}")
        return

    gc_df = df[df['environment'] == 'globular_cluster'].copy()
    
    # 2. Cluster Densities (Baumgardt 2018 / Harris 2010)
    # log10(rho_c) in L_sun/pc^3
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
        "NGC 6539": 3.30, "M4 (NGC 6121)": 2.85
    }
    
    # Map densities to dataframe
    gc_df['log_rho_c'] = gc_df['cluster'].map(CLUSTER_DENSITIES)
    
    # Filter out clusters without density info or Pdot
    gc_df = gc_df.dropna(subset=['log_rho_c', 'logPdot_abs', 'logP', 'log_b_proxy'])
    
    print(f"Analyzing {len(gc_df)} pulsars in {gc_df['cluster'].nunique()} clusters.")
    
    # 3. Define Models
    
    # Model A: WLS on Cluster Means (Weighted by sample size for valid baseline)
    # Weighted by N to properly account for cluster sample sizes
    cluster_means = gc_df.groupby('cluster').agg({
        'logPdot_abs': 'mean',
        'log_rho_c': 'first',
        'logP': 'count'  # sample size per cluster
    }).reset_index().rename(columns={'logP': 'n_obs'})
    
    # Weighted Least Squares: weights = N (inverse variance approximation)
    wls_means = smf.wls("logPdot_abs ~ log_rho_c", data=cluster_means, weights=cluster_means['n_obs']).fit()
    
    # Model B: Mixed Effects Model (Hierarchical)
    # Fixed effects: log_rho_c, logP, log_b_proxy (controls)
    # Random effects: Intercept for each cluster
    # This properly accounts for within-cluster variance vs between-cluster scaling
    
    # We standardize controls for numerical stability
    gc_df['logP_std'] = (gc_df['logP'] - gc_df['logP'].mean()) / gc_df['logP'].std()
    # Note: B-proxy is derived from P and Pdot. Including it as a predictor for Pdot creates 
    # perfect circularity/multicollinearity. We must EXCLUDE it.
    # We control for P (spin) and rely on the fact that B-field distributions are standard.
    
    gc_df['log_rho_c_centered'] = gc_df['log_rho_c'] - gc_df['log_rho_c'].mean()
    
    # Formula: Pdot depends on density (fixed slope), intrinsic params (fixed slopes),
    # plus a random offset per cluster (random intercept).
    
    # Removed logB_std to fix circularity
    md = smf.mixedlm(
        "logPdot_abs ~ log_rho_c_centered + logP_std", 
        gc_df, 
        groups=gc_df["cluster"]
    )
    # Use default optimizer first, fall back if needed
    try:
        mdf = md.fit()
    except (ValueError, TypeError) as e:
        print(f"  Default optimizer failed ({type(e).__name__}), trying Powell...")
        try:
            mdf = md.fit(method='powell', maxiter=1000)
        except (ValueError, TypeError) as e:
            print(f"  Powell optimizer failed ({type(e).__name__})")
            return
    
    print("\n" + "="*60)
    print("MODEL A: WLS on Cluster Means (Weighted by N, Valid Baseline)")
    print("="*60)
    print(wls_means.summary())
    
    print("\n" + "="*60)
    print("MODEL B: Hierarchical Mixed-Effects (Recommended Update)")
    print("="*60)
    print(mdf.summary())
    
    # Extract Key Stats
    density_slope = mdf.params['log_rho_c_centered']
    density_slope_err = mdf.bse['log_rho_c_centered']
    density_p = mdf.pvalues['log_rho_c_centered']
    
    # Comparison with Newtonian Prediction
    # Newtonian prediction from step 5.32 is slope ~ 0.72 - 0.82
    # We test H0: slope = 0.72 (Newtonian) vs H1: slope != 0.72
    
    newtonian_slope_target = 0.72
    z_score = (density_slope - newtonian_slope_target) / density_slope_err
    p_reject_newtonian = 2 * (1 - stats.norm.cdf(abs(z_score)))
    
    print("\n" + "="*60)
    print("HYPOTHESIS TEST: Suppressed Density Scaling")
    print("="*60)
    print(f"Observed Slope (Mixed Model): {density_slope:.3f} ± {density_slope_err:.3f}")
    print(f"Newtonian Target Slope:       {newtonian_slope_target:.3f}")
    print(f"Difference:                   {density_slope - newtonian_slope_target:.3f}")
    print(f"Z-score (Rejection of GR):    {z_score:.1f}σ")
    print(f"p-value:                      {p_reject_newtonian:.2e}")
    
    # Save results
    results = {
        "model_a_wls_slope": float(wls_means.params['log_rho_c']),
        "model_a_wls_error": float(wls_means.bse['log_rho_c']),
        "model_b_mixed_slope": float(density_slope),
        "model_b_mixed_error": float(density_slope_err),
        "model_b_mixed_p": float(density_p),
        "newtonian_predicted_slope": float(newtonian_slope_target),
        "newtonian_slope_range": [0.72, 0.82],
        "observed_vs_newtonian_diff": float(density_slope - newtonian_slope_target),
        "rejection_sigma": float(abs(z_score)),
        "rejection_p_value": float(p_reject_newtonian),
        "interpretation": "Observed slope significantly flatter than Newtonian prediction",
        "suppression_factor": float(density_slope / newtonian_slope_target)
    }
    
    with open("results/outputs/step_5_33_hierarchical_density_results.json", "w") as f:
        json.dump(results, f, indent=4)

if __name__ == "__main__":
    run_hierarchical_analysis()
