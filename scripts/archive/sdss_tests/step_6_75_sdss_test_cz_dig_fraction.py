#!/usr/bin/env python3
"""
Step 6.75: SDSS Test CZ - Diffuse Ionized Gas (LIER) Fraction

Hypothesis:
DIG is powered largely by hot evolved stars (pAGB, WD). In deep potentials (high sigma), 
if time dilation makes populations effectively "older", the population of these stars 
might be enhanced, increasing the hardness of the ionizing field and the prevalence of LIERs.

Prediction:
Fraction of LIER-like galaxies (or LIER spectral features) correlates with Sigma.

Data:
- mangaDAPall: stellar_sigma_1re, emline_gew_1re_ha_6564, emline_gflux_1re_nii_6583, emline_gflux_1re_ha_6564
- mangaTarget: nsa_elpetro_mass

Method:
1. Query global spectral properties (1Re aperture).
2. Classify using WHAN diagram (Cid Fernandes et al. 2011):
   - SF: log([NII]/Ha) < -0.4
   - Seyfert: log([NII]/Ha) > -0.4 & EW(Ha) > 6
   - LIER: log([NII]/Ha) > -0.4 & 3 < EW(Ha) < 6
   - Retired: EW(Ha) < 3
3. Calculate the fraction of (LIER + Retired) galaxies in bins of Sigma.
4. Control for Mass? LIERs are known to be massive. We need to see if Sigma adds predictive power.
"""

import pandas as pd
import numpy as np
from scipy import stats
import os
import json
import matplotlib.pyplot as plt
import requests
import time

# Paths
DATA_DIR = os.path.join(os.path.dirname(__file__), '..', '..', 'data', 'sdss')
RESULTS_DIR = os.path.join(os.path.dirname(__file__), '..', '..', 'results', 'outputs')
FIGURES_DIR = os.path.join(os.path.dirname(__file__), '..', '..', 'results', 'figures')

SDSS_URL = "https://skyserver.sdss.org/dr18/SkyServerWS/SearchTools/SqlSearch"

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
        except Exception as e:
            print(f"  Error: {e}")
        if attempt < max_retries - 1:
            time.sleep(2)
    return None

def download_data(limit=100):
    print(f"Querying SDSS for Test CZ (Limit: {limit})...")
    
    # Note: EW in mangaDAP is often positive for emission? Or negative?
    # Standard convention: Emission is positive in flux, but EW sign varies.
    # In MaNGA DAP, emission lines usually have positive EW (need to verify, but assuming positive for emission).
    
    sql = f"""
    SELECT TOP {limit}
        d.mangaid,
        d.stellar_sigma_1re as sigma,
        d.emline_gew_1re_ha_6564 as ew_ha,
        d.emline_gflux_1re_ha_6564 as flux_ha,
        d.emline_gflux_1re_nii_6583 as flux_nii,
        s.nsa_elpetro_mass as logmass
    FROM mangaDAPall d
    JOIN mangaTarget s ON d.mangaid = s.mangaid
    WHERE d.drp3qual = 0
      AND d.stellar_sigma_1re > 0
      AND d.emline_gflux_1re_ha_6564 > 0
    """
    return query_sdss(sql)

def analyze_lier_fraction(df):
    print("Analyzing LIER Fraction...")
    
    if df is None or len(df) < 100:
        print("  Insufficient data.")
        return None
        
    df = df.copy()
    
    # Calculate Ratios
    # Protect against zero division
    df = df[df['flux_ha'] > 0]
    df['n2ha'] = np.log10(df['flux_nii'] / df['flux_ha'])
    
    # WHAN Classification
    # Note: MaNGA DAP EW is generally positive for emission.
    # Check bounds:
    # SF: n2ha < -0.4
    # AGN (Seyfert): n2ha > -0.4 & ew_ha > 6
    # LIER: n2ha > -0.4 & 3 < ew_ha < 6
    # Retired: ew_ha < 3
    
    def classify(row):
        n2ha = row['n2ha']
        ew = row['ew_ha']
        
        if pd.isna(n2ha) or pd.isna(ew): return 'Unknown'
        
        if ew < 3:
            return 'Retired' # Passive / Weak LIER
        
        if n2ha < -0.4:
            return 'SF'
        
        if ew > 6:
            return 'Seyfert'
            
        return 'LIER' # Strong LIER
        
    df['class'] = df.apply(classify, axis=1)
    
    print("  Class counts:")
    print(df['class'].value_counts())
    
    # Group LIER and Retired as "DIG-dominated" / "Old-Star Ionized"
    df['is_lier_like'] = df['class'].isin(['LIER', 'Retired']).astype(int)
    
    # Bin by Sigma
    df['sigma_bin'] = pd.qcut(df['sigma'], q=10, labels=False)
    
    grouped = df.groupby('sigma_bin').agg({
        'sigma': 'mean',
        'is_lier_like': ['mean', 'sem', 'count']
    }).reset_index()
    grouped.columns = ['bin', 'sigma', 'frac', 'sem', 'count']
    
    # Regression
    slope, intercept, r_val, p_val, std_err = stats.linregress(grouped['sigma'], grouped['frac'])
    
    print(f"  Correlation (Sigma vs LIER Fraction): r={r_val:.3f}, p={p_val:.2e}, slope={slope:.5f}")
    
    # Plot
    fig, ax = plt.subplots(1, 2, figsize=(12, 5))
    
    # WHAN Diagram
    # Subsample for speed
    plot_df = df.sample(min(len(df), 2000))
    ax[0].scatter(plot_df['n2ha'], np.log10(plot_df['ew_ha'].clip(lower=0.1)), c='gray', alpha=0.3, s=5)
    ax[0].axvline(-0.4, color='k', linestyle='--')
    ax[0].axhline(np.log10(3), color='k', linestyle='--')
    ax[0].axhline(np.log10(6), color='k', linestyle='--')
    ax[0].set_xlabel('log([NII]/Ha)')
    ax[0].set_ylabel('log(EW Ha)')
    ax[0].set_title('WHAN Diagram')
    
    # Fraction vs Sigma
    ax[1].errorbar(grouped['sigma'], grouped['frac'], yerr=grouped['sem'], fmt='o-', color='purple')
    
    # Fit line
    x_range = np.linspace(grouped['sigma'].min(), grouped['sigma'].max(), 100)
    ax[1].plot(x_range, intercept + slope*x_range, 'k--', label=f'r={r_val:.2f}')
    
    ax[1].set_xlabel('Velocity Dispersion (km/s)')
    ax[1].set_ylabel('Fraction of LIER/Retired')
    ax[1].set_title('Test CZ: LIER Fraction vs Sigma')
    ax[1].legend()
    ax[1].grid(True, alpha=0.3)
    
    out_path = os.path.join(FIGURES_DIR, 'sdss_test_cz_lier_fraction.png')
    plt.savefig(out_path, dpi=150)
    print(f"Figure saved to {out_path}")
    
    return {
        'correlation_r': r_val,
        'p_value': p_val,
        'slope': slope,
        'n_galaxies': int(len(df))
    }

def main():
    cache_path = os.path.join(DATA_DIR, 'sdss_manga_dap.csv')
    
    if os.path.exists(cache_path):
        print("Loading cached data...")
        df = pd.read_csv(cache_path)
    else:
        df = download_data()
        if df is not None:
            df.to_csv(cache_path, index=False)
        else:
            print("Download failed.")
            return

    results = analyze_lier_fraction(df)
    
    out_path = os.path.join(RESULTS_DIR, 'sdss_test_cz_results.json')
    if results:
        with open(out_path, 'w') as f:
            json.dump(results, f, indent=2)
            
        print("\nSUMMARY TEST CZ:")
        print(f"Correlation (Sigma vs LIER Frac): {results['correlation_r']:.3f} (p={results['p_value']:.2e})")
        
        if results['p_value'] < 0.05 and results['correlation_r'] > 0.1:
             print("RESULT: SIGNAL (Positive correlation)")
        elif results['p_value'] < 0.05 and results['correlation_r'] < -0.1:
             print("RESULT: CONTRADICTED (Negative correlation)")
        else:
             print("RESULT: NULL (No significant correlation)")

if __name__ == "__main__":
    main()
