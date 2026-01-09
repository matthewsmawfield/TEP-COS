#!/usr/bin/env python3
"""
Step 6.9: Test D - Photometric vs Spectroscopic Redshift Residual

TEP HYPOTHESIS:
Photometric redshift estimation relies on SED fitting, which assumes standard stellar 
evolution timescales. If TEP time dilation is real, high-σ galaxies should have 
systematically biased photo-z relative to spec-z because their SEDs will correspond 
to "older" populations (or different evolutionary stages) than standard templates 
predict for their redshift.

TEP PREDICTION:
  Define: Δz = z_photo - z_spec
  At fixed color and stellar mass:
    r(Δz, σ) ≠ 0 (Systematic bias)

DATA:
  - Photoz
  - SpecPhotoAll

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
from datetime import datetime

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
                print(f"  Response: {response.text}")
        except Exception as e:
            print(f"  Error: {e}")
        if attempt < max_retries - 1:
            time.sleep(2)
    return None

def download_data(limit=10000):
    print(f"Querying SDSS for Test D (Limit: {limit})...")
    sql = f"""
    SELECT TOP {limit}
        sp.specObjID,
        sp.z AS z_spec,
        pz.z AS z_photo,
        pz.zErr AS z_photo_err,
        e.sigmaStars AS sigma_stars,
        s.logMass,
        
        -- Colors
        sp.dered_g - sp.dered_r AS g_minus_r,
        sp.dered_r - sp.dered_i AS r_minus_i
        
    FROM SpecPhotoAll sp
    JOIN Photoz pz ON sp.objID = pz.objID
    JOIN emissionLinesPort e ON sp.specObjID = e.specObjID
    JOIN stellarMassFSPSGranWideDust s ON sp.specObjID = s.specObjID
    
    WHERE 
        sp.z BETWEEN 0.02 AND 0.30
        AND sp.zWarning = 0
        AND pz.zErr < 0.1
        AND e.sigmaStars > 50 AND e.sigmaStars < 400
        AND pz.photoErrorClass = 1
    """
    return query_sdss(sql)

def analyze_photoz_residual(df):
    print("Analyzing Photo-z residuals...")
    
    # 1. Compute Residual
    df['delta_z'] = df['z_photo'] - df['z_spec']
    df['delta_z_norm'] = df['delta_z'] / (1 + df['z_spec'])
    df['log_sigma'] = np.log10(df['sigma_stars'])
    
    # 2. Simple Correlation
    r_resid, p_resid = stats.pearsonr(df['log_sigma'], df['delta_z_norm'])
    print(f"Simple r(Δz_norm, σ): {r_resid:.4f} (p={p_resid:.2e})")
    
    # 3. Controlled Correlation (Partial)
    # Control for color and mass, as photo-z depends heavily on these
    # Regress out color and mass from delta_z_norm
    
    from sklearn.linear_model import LinearRegression
    X_control = df[['g_minus_r', 'r_minus_i', 'logMass', 'z_spec']].values
    y = df['delta_z_norm'].values
    
    reg = LinearRegression().fit(X_control, y)
    df['delta_z_resid'] = y - reg.predict(X_control)
    
    r_controlled, p_controlled = stats.pearsonr(df['log_sigma'], df['delta_z_resid'])
    print(f"Controlled r(Δz_resid, σ): {r_controlled:.4f} (p={p_controlled:.2e})")
    
    return {
        'r_simple': float(r_resid),
        'p_simple': float(p_resid),
        'r_controlled': float(r_controlled),
        'p_controlled': float(p_controlled)
    }

def create_figure(df, results):
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    # Delta z vs Sigma (Simple)
    ax = axes[0]
    ax.hexbin(df['log_sigma'], df['delta_z_norm'], gridsize=40, cmap='RdBu_r', mincnt=1)
    ax.axhline(0, color='k', ls='--')
    ax.set_xlabel(r'$\log(\sigma)$')
    ax.set_ylabel(r'$\Delta z / (1+z_{spec})$')
    ax.set_title(f'Photo-z Residual vs $\sigma$ (r={results["r_simple"]:.3f})')
    
    # Delta z vs Sigma (Controlled)
    ax = axes[1]
    ax.hexbin(df['log_sigma'], df['delta_z_resid'], gridsize=40, cmap='RdBu_r', mincnt=1)
    ax.axhline(0, color='k', ls='--')
    ax.set_xlabel(r'$\log(\sigma)$')
    ax.set_ylabel(r'Residual $\Delta z$ (after M*, color correction)')
    ax.set_title(f'Controlled Residual vs $\sigma$ (r={results["r_controlled"]:.3f})')
    
    plt.tight_layout()
    out_file = os.path.join(FIGURES_DIR, 'sdss_test_d_photoz_residual.png')
    plt.savefig(out_file, dpi=150)
    print(f"Figure saved to {out_file}")

def main():
    print("="*60)
    print("STEP 6.9: TEST D - PHOTO-Z RESIDUAL")
    print("="*60)
    
    cache_path = os.path.join(DATA_DIR, 'sdss_photoz_data.csv')
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
    
    results = analyze_photoz_residual(df)
    create_figure(df, results)
    
    out_path = os.path.join(RESULTS_DIR, 'sdss_test_d_results.json')
    with open(out_path, 'w') as f:
        json.dump(results, f, indent=2)
        
    print("\nSUMMARY:")
    print(f"Bias Detected: {results['p_controlled'] < 0.05}")
    print(f"Magnitude: {results['r_controlled']:.4f}")

if __name__ == "__main__":
    main()
