#!/usr/bin/env python3
"""
Step 6.7: Test B - Mass Estimation Method Discrepancy

TEP HYPOTHESIS:
Different stellar mass estimation methods have different dependencies on age-sensitive 
spectral features. If TEP time dilation is real, methods relying on age-sensitive features 
(SED colors, D4000) should underestimate masses for high-σ galaxies relative to methods 
using age-insensitive features (spectral PCA) or those less sensitive to the specific 
star formation history timescale.

TEP PREDICTION:
  Define: ΔM* = log(M*_SED) - log(M*_PCA)
  At fixed color and redshift:
    r(ΔM*, σ) < 0  (SED-based masses relatively lower at high σ)

DATA:
  Comparison of:
  - stellarMassFSPSGranWideDust (SED-based)
  - stellarMassPCAWiscBC03 (Spectral PCA)
  - stellarMassPortsmouth (Template fitting)

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
    """Execute SQL query against SDSS SkyServer."""
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
    """Download mass comparison data."""
    print(f"Querying SDSS for Test B (Limit: {limit})...")
    
    # Note: Using a subset of tables to ensure query success and relevance
    # stellarMassFSPSGranWideDust: SED based
    # stellarMassPCAWiscBC03: PCA based
    # stellarMassStarformingPort: Portsmouth
    
    sql = f"""
    SELECT TOP {limit}
        s1.specObjID,
        z.z AS redshift,
        p.sigmaStars AS sigma_stars,
        p.sigmaStarsErr AS sigma_stars_err,
        
        -- Method 1: FSPS Granada Wide Dust (SED-based)
        s1.logMass AS logM_FSPS,
        
        -- Method 2: Wisconsin PCA BC03 (Spectral PCA)
        s5.mstellar_median AS logM_PCA,
        
        -- Method 3: Portsmouth Star-forming
        s4.logMass AS logM_Port,
        
        -- Color (from SpecPhotoAll)
        z.dered_g - z.dered_r AS g_minus_r
        
    FROM stellarMassFSPSGranWideDust s1
    JOIN stellarMassPCAWiscBC03 s5 ON s1.specObjID = s5.specObjID
    JOIN stellarMassStarformingPort s4 ON s1.specObjID = s4.specObjID
    JOIN emissionLinesPort p ON s1.specObjID = p.specObjID
    JOIN SpecPhotoAll z ON s1.specObjID = z.specObjID
    
    WHERE 
        z.z BETWEEN 0.02 AND 0.25
        AND p.sigmaStars > 50 AND p.sigmaStars < 400
        AND p.sigmaStarsErr < 30
        AND s1.logMass > 9.0
        AND s5.mstellar_median > 9.0
        AND s4.logMass > 9.0
        AND z.dered_g > -100 AND z.dered_r > -100
    """
    
    return query_sdss(sql)

def analyze_mass_discrepancy(df):
    """Analyze the discrepancy between mass methods vs sigma."""
    print("Analyzing mass discrepancies...")
    
    # Calculate deltas
    # Prediction: SED methods (FSPS) underestimate mass at high sigma relative to PCA
    # Delta = FSPS - PCA. Expect negative correlation with sigma.
    df['delta_FSPS_PCA'] = df['logM_FSPS'] - df['logM_PCA']
    df['delta_Port_PCA'] = df['logM_Port'] - df['logM_PCA']
    df['log_sigma'] = np.log10(df['sigma_stars'])
    
    # 1. Simple Correlation
    r_fsps_pca, p_fsps_pca = stats.pearsonr(df['log_sigma'], df['delta_FSPS_PCA'])
    r_port_pca, p_port_pca = stats.pearsonr(df['log_sigma'], df['delta_Port_PCA'])
    
    print(f"Correlation r(ΔM_FSPS-PCA, σ): {r_fsps_pca:.4f} (p={p_fsps_pca:.2e})")
    print(f"Correlation r(ΔM_Port-PCA, σ): {r_port_pca:.4f} (p={p_port_pca:.2e})")
    
    # 2. Color-Controlled Analysis
    # Bin by color
    df['color_bin'] = pd.qcut(df['g_minus_r'], 5, labels=['Blue', 'Green-Blue', 'Green', 'Green-Red', 'Red'])
    
    binned_results = []
    for cbin in df['color_bin'].cat.categories:
        sub = df[df['color_bin'] == cbin]
        if len(sub) < 100: continue
        
        r, p = stats.pearsonr(sub['log_sigma'], sub['delta_FSPS_PCA'])
        binned_results.append({
            'bin': cbin,
            'n': len(sub),
            'r': r,
            'p': p,
            'mean_color': sub['g_minus_r'].mean()
        })
        print(f"  Bin {cbin:<10}: r={r:.4f}")
        
    return {
        'r_fsps_pca': float(r_fsps_pca),
        'p_fsps_pca': float(p_fsps_pca),
        'r_port_pca': float(r_port_pca),
        'binned': binned_results,
        'tep_consistent': bool(r_fsps_pca < 0)
    }

def create_figure(df, results):
    """Create summary figure."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    # Scatter plot
    ax = axes[0]
    ax.scatter(df['log_sigma'], df['delta_FSPS_PCA'], alpha=0.1, s=2, c='k')
    
    # Trend line
    slope, intercept = np.polyfit(df['log_sigma'], df['delta_FSPS_PCA'], 1)
    x = np.linspace(df['log_sigma'].min(), df['log_sigma'].max(), 100)
    ax.plot(x, slope*x + intercept, 'r-', lw=2, label=f'Slope={slope:.3f}')
    
    ax.set_xlabel(r'$\log(\sigma_*)$')
    ax.set_ylabel(r'$\Delta \log M_*$ (FSPS - PCA)')
    ax.set_title(f'Mass Discrepancy vs $\sigma$ (r={results["r_fsps_pca"]:.3f})')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Binned results
    ax = axes[1]
    bins = [b['bin'] for b in results['binned']]
    rs = [b['r'] for b in results['binned']]
    ax.bar(bins, rs, color='steelblue')
    ax.axhline(0, color='k')
    ax.set_ylim(-0.5, 0.5)
    ax.set_ylabel('Correlation Coefficient r')
    ax.set_title('Correlation by Color Bin')
    
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, 'sdss_test_b_mass_discrepancy.png'), dpi=150)
    print(f"Figure saved to {os.path.join(FIGURES_DIR, 'sdss_test_b_mass_discrepancy.png')}")

def main():
    print("="*60)
    print("STEP 6.7: TEST B - MASS METHOD DISCREPANCY")
    print("="*60)
    
    cache_path = os.path.join(DATA_DIR, 'sdss_mass_comparison.csv')
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
    
    results = analyze_mass_discrepancy(df)
    create_figure(df, results)
    
    # Save results
    out_path = os.path.join(RESULTS_DIR, 'sdss_test_b_results.json')
    with open(out_path, 'w') as f:
        json.dump(results, f, indent=2)
        
    print("\nSUMMARY:")
    print(f"TEP Prediction (r < 0): {'CONFIRMED' if results['tep_consistent'] else 'REJECTED'}")
    print(f"Global Correlation: {results['r_fsps_pca']:.4f}")

if __name__ == "__main__":
    main()
