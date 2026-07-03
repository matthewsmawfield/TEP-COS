#!/usr/bin/env python3
"""
Step 6.10: Test E - Light-Weighted vs Mass-Weighted Age Discrepancy

TEP HYPOTHESIS:
MaNGA Firefly provides both light-weighted (LW) and mass-weighted (MW) stellar ages. 
These probe different stellar populations. Under TEP, time dilation is a field effect 
that should affect the underlying clock rates. If dilation depends on potential depth, 
we might expect systematic differences in how these two age measures diverge at high sigma, 
potentially due to the different effective formation epochs they probe.

TEP PREDICTION:
  Define: ΔAge = log(Age_LW) - log(Age_MW)
  At fixed metallicity:
    r(ΔAge, σ) ≠ 0 (Systematic divergence tracking potential)

DATA:
  - mangaFirefly_miles (Ages, Z, Mass)
  - mangaDAPall (Sigma, Redshift)

Author: M. Smawfield
Date: January 2026
"""

import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
import json
import os
import requests
import time

# Configuration
SDSS_URL = "https://skyserver.sdss.org/dr18/SkyServerWS/SearchTools/SqlSearch"
DATA_DIR = os.path.join(os.path.dirname(__file__), '..', '..', 'data', 'sdss')
RESULTS_DIR = os.path.join(os.path.dirname(__file__), '..', '..', 'results', 'outputs')
FIGURES_DIR = os.path.join(os.path.dirname(__file__), '..', '..', 'results', 'figures')

os.makedirs(DATA_DIR, exist_ok=True)
os.makedirs(RESULTS_DIR, exist_ok=True)
os.makedirs(FIGURES_DIR, exist_ok=True)

def query_sdss(sql, max_retries=3):
    for attempt in range(max_retries):
        try:
            response = requests.get(
                SDSS_URL,
                params={"cmd": sql, "format": "json"},
                timeout=300
            )
            if response.status_code == 200:
                data = response.json()
                if data and len(data) > 0 and "Rows" in data[0]:
                    return pd.DataFrame(data[0]["Rows"])
            else:
                print(f"  HTTP {response.status_code}")
                print(f"  Response: {response.text[:500]}")
        except Exception as e:
            print(f"  Error: {e}")
        if attempt < max_retries - 1:
            time.sleep(2)
    return None

def download_data(limit=10000):
    print(f"Querying SDSS for Test E (Limit: {limit})...")
    
    # Note: Using mangaFirefly_miles and mangaDAPall
    # Joining on PLATEIFU
    
    sql = f"""
    SELECT TOP {limit}
        f.PLATEIFU,
        d.z AS redshift,
        
        -- Ages (Linear Gyr)
        f.LW_AGE_1RE,
        f.LW_AGE_1RE_ERROR,
        f.MW_AGE_1RE,
        f.MW_AGE_1RE_ERROR,
        
        -- Metallicity
        f.LW_Z_1RE,
        f.MW_Z_1RE,
        
        -- Mass
        f.PHOTOMETRIC_MASS,
        
        -- Kinematics from DAPall
        d.stellar_sigma_1re
        
    FROM mangaFirefly_miles f
    JOIN mangaDAPall d ON f.PLATEIFU = d.plateifu
    
    WHERE 
        d.z BETWEEN 0.01 AND 0.15
        AND d.stellar_sigma_1re > 50 AND d.stellar_sigma_1re < 400
        AND f.LW_AGE_1RE > 0 AND f.MW_AGE_1RE > 0
        AND f.LW_AGE_1RE_ERROR > 0
        AND f.PHOTOMETRIC_MASS > 0
    """
    
    return query_sdss(sql)

def analyze_age_discrepancy(df):
    print("Analyzing Age Discrepancy...")
    
    # 1. Compute Log Ages and Delta
    # Ages are in Gyr.
    df['log_Age_LW'] = np.log10(df['LW_AGE_1RE'])
    df['log_Age_MW'] = np.log10(df['MW_AGE_1RE'])
    df['delta_log_age'] = df['log_Age_LW'] - df['log_Age_MW']
    
    df['log_sigma'] = np.log10(df['stellar_sigma_1re'])
    
    # 2. Simple Correlation
    r_simple, p_simple = stats.pearsonr(df['log_sigma'], df['delta_log_age'])
    print(f"Simple r(ΔAge, σ): {r_simple:.4f} (p={p_simple:.2e})")
    
    # 3. Controlled Correlation (Metallicity & Mass)
    # MW age is naturally older than LW age. The question is if the GAP depends on sigma.
    # Control for Metallicity (LW_Z_1RE) and Mass (PHOTOMETRIC_MASS)

    # Prepare regressors
    # Handle NaNs if any (though SQL should filter)
    df_clean = df.dropna(subset=['delta_log_age', 'log_sigma', 'LW_Z_1RE', 'PHOTOMETRIC_MASS']).copy()

    # Standardize regressors to prevent numerical ill-conditioning
    X = df_clean[['LW_Z_1RE', 'PHOTOMETRIC_MASS']].values
    y = df_clean['delta_log_age'].values
    X_mean = X.mean(axis=0)
    X_std = X.std(axis=0)
    X_std[X_std == 0] = 1.0  # Guard against zero-variance columns
    X_scaled = (X - X_mean) / X_std

    # Use numpy.linalg.lstsq instead of sklearn LinearRegression to avoid
    # a spurious matmul warning triggered by sklearn/numpy interaction on
    # some platforms (divide-by-zero in matmul despite clean data).
    X_design = np.column_stack([np.ones(len(X_scaled)), X_scaled])
    coeffs, *_ = np.linalg.lstsq(X_design, y, rcond=None)
    intercept = coeffs[0]
    beta = coeffs[1:]
    # Suppress spurious BLAS divide-by-zero/overflow warnings on macOS/ARM
    with np.errstate(divide='ignore', invalid='ignore', over='ignore'):
        pred = X_design @ coeffs
    df_clean['delta_age_resid'] = y - pred

    r_controlled, p_controlled = stats.pearsonr(df_clean['log_sigma'], df_clean['delta_age_resid'])
    print(f"Controlled r(ΔAge_resid, σ): {r_controlled:.4f} (p={p_controlled:.2e})")
    
    return {
        'r_simple': float(r_simple),
        'p_simple': float(p_simple),
        'r_controlled': float(r_controlled),
        'p_controlled': float(p_controlled),
        'n_sample': int(len(df_clean))
    }

def create_figure(df, results):
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    # Raw Delta vs Sigma
    ax = axes[0]
    ax.hexbin(df['log_sigma'], df['delta_log_age'], gridsize=30, cmap='magma_r', mincnt=1)
    ax.set_xlabel(r'$\log(\sigma_{1Re})$')
    ax.set_ylabel(r'$\log(\text{Age}_{LW}) - \log(\text{Age}_{MW})$')
    ax.set_title(fr'Age Discrepancy vs $\sigma$ (r={results["r_simple"]:.3f})')
    
    # Controlled Delta vs Sigma
    # Need to re-compute residual for plotting if I want to plot the residual
    # (Just reusing the logic from analyze for the plot data if available, 
    # but 'df' here might not have the 'delta_age_resid' col if I used df_clean inside the function)
    # Let's quickly re-calc for the plot or just pass the clean df.
    # For simplicity, I will assume df has the col or just skip if not present (it won't be in original df)
    
    # Re-calculate residuals for plotting
    df_clean = df.dropna(subset=['delta_log_age', 'log_sigma', 'LW_Z_1RE', 'PHOTOMETRIC_MASS']).copy()
    X = df_clean[['LW_Z_1RE', 'PHOTOMETRIC_MASS']].values
    y = df_clean['delta_log_age'].values
    X_mean = X.mean(axis=0)
    X_std = X.std(axis=0)
    X_std[X_std == 0] = 1.0
    X_scaled = (X - X_mean) / X_std
    X_design = np.column_stack([np.ones(len(X_scaled)), X_scaled])
    coeffs, *_ = np.linalg.lstsq(X_design, y, rcond=None)
    with np.errstate(divide='ignore', invalid='ignore', over='ignore'):
        pred = X_design @ coeffs
    df_clean['delta_age_resid'] = y - pred
    
    ax = axes[1]
    ax.hexbin(df_clean['log_sigma'], df_clean['delta_age_resid'], gridsize=30, cmap='magma_r', mincnt=1)
    ax.axhline(0, color='k', ls='--')
    ax.set_xlabel(r'$\log(\sigma_{1Re})$')
    ax.set_ylabel(r'Residual $\Delta$Age (after Z, M* control)')
    ax.set_title(rf'Controlled Residual vs $\sigma$ (r={results["r_controlled"]:.3f})')
    
    plt.tight_layout()
    out_file = os.path.join(FIGURES_DIR, 'step_60_sdss_test_e_age_discrepancy.png')
    plt.savefig(out_file, dpi=150)
    print(f"Figure saved to {out_file}")

def main():
    print("="*60)
    print("STEP 6.10: TEST E - LW vs MW AGE DISCREPANCY")
    print("="*60)
    
    cache_path = os.path.join(DATA_DIR, 'manga_age_data.csv')
    if os.path.exists(cache_path):
        print("Loading cached data...")
        df = pd.read_csv(cache_path)
    else:
        df = download_data()
        if df is not None:
            df.to_csv(cache_path, index=False)
        else:
            print("Failed to download data.")
            return
            
    print(f"Data size: {len(df)}")
    
    results = analyze_age_discrepancy(df)
    create_figure(df, results)
    
    out_path = os.path.join(RESULTS_DIR, 'step_60_manga_test_e_age_discrepancy.json')
    with open(out_path, 'w') as f:
        json.dump(results, f, indent=2)
        
    print("\nSUMMARY:")
    print(f"Discrepancy Detected: {results['p_controlled'] < 0.05}")
    print(f"Correlation Strength: {results['r_controlled']:.4f}")

if __name__ == "__main__":
    main()
