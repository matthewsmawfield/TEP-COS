#!/usr/bin/env python3
"""
Step 6.89: SDSS Test DR - Brown Dwarf Desert (Formation Rates)

Hypothesis:
The paucity of brown dwarf companions ("desert") is linked to formation mechanisms. 
If TEP modifies the Jeans instability or fragmentation rates in high-metallicity/
deep-potential environments, the boundaries of this desert might shift.

Prediction:
Fraction of Substellar Companions varies with [Fe/H] or Galactic Location.

Data:
- marvelsStar: starname, teff, feh, ra, dec
- marvelsVelocityCurveUF1D: minMass, period

Method:
1. Select stars with velocity curves.
2. Identify companions in the Brown Dwarf range (13 < M < 80 M_Jup).
   Note: minMass is likely in Jupiter masses.
3. Compare the distribution of host metallicity for BD-hosts vs non-BD-hosts.
   Or correlate companion mass with host metallicity.
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
    print(f"Querying SDSS for Test DR (Limit: {limit})...")
    
    # We need all stars with velocity curves to calculate a fraction, 
    # or just the distribution of BDs?
    # To check "Fraction varies with Fe/H", we need the parent population.
    # Join marvelsStar with marvelsVelocityCurveUF1D.
    # marvelsVelocityCurveUF1D might have multiple entries per star (different planets?)
    # or one best fit? 
    # Usually UF1D is 1-Dimensional (single planet) fit.
    
    sql = f"""
    SELECT TOP {limit}
        s.starname, 
        s.feh, 
        s.teff,
        v.minMass, 
        v.period
    FROM marvelsStar s
    JOIN marvelsVelocityCurveUF1D v ON s.starname = v.starname
    WHERE s.feh > -999
      AND v.minMass > 0
    """
    return query_sdss(sql)

def analyze_bd_desert(df):
    print("Analyzing Brown Dwarf Desert...")
    
    if df is None or len(df) < 50:
        print("  Insufficient data.")
        return None
        
    df = df.dropna().copy()
    
    # Define Brown Dwarfs: 13 < M < 80 M_Jup
    # Define Planets: M < 13 M_Jup
    # Define Stars: M > 80 M_Jup
    
    df['type'] = 'Star'
    df.loc[df['minMass'] < 80, 'type'] = 'BrownDwarf'
    df.loc[df['minMass'] < 13, 'type'] = 'Planet'
    
    print("  Companion Counts:")
    print(df['type'].value_counts())
    
    bds = df[df['type'] == 'BrownDwarf']
    planets = df[df['type'] == 'Planet']
    
    if len(bds) < 5:
        print("  Not enough Brown Dwarfs for statistical analysis.")
        return None
        
    # Test 1: Metallicity Distribution of BD hosts vs Planet hosts
    # T-test
    t_stat, p_val = stats.ttest_ind(bds['feh'], planets['feh'], equal_var=False)
    
    print(f"  Metallicity Comparison (BD vs Planet hosts): t={t_stat:.3f}, p={p_val:.3f}")
    print(f"  Mean [Fe/H]: BD={bds['feh'].mean():.3f}, Planet={planets['feh'].mean():.3f}")
    
    # Test 2: Correlation of Mass with Metallicity within BD range
    slope, intercept, r_val, p_reg, std_err = stats.linregress(bds['feh'], bds['minMass'])
    print(f"  Correlation (Fe/H vs BD Mass): r={r_val:.3f}, p={p_reg:.3f}, slope={slope:.3f}")

    # Plot
    fig, ax = plt.subplots(1, 2, figsize=(12, 5))
    
    # Histograms
    ax[0].hist(planets['feh'], bins=20, alpha=0.5, label='Planets', density=True, color='blue')
    ax[0].hist(bds['feh'], bins=10, alpha=0.5, label='Brown Dwarfs', density=True, color='orange')
    ax[0].set_xlabel('[Fe/H]')
    ax[0].set_ylabel('Density')
    ax[0].set_title('Host Metallicity Distribution')
    ax[0].legend()
    
    # Scatter
    ax[1].scatter(df['feh'], df['minMass'], alpha=0.6, s=10, c='gray', label='All')
    ax[1].scatter(bds['feh'], bds['minMass'], alpha=0.8, s=20, c='orange', label='Brown Dwarfs')
    ax[1].set_yscale('log')
    ax[1].set_xlabel('[Fe/H]')
    ax[1].set_ylabel('Companion Mass (M_Jup)')
    ax[1].set_title('Mass vs Metallicity')
    ax[1].axhline(13, color='k', linestyle='--', alpha=0.3)
    ax[1].axhline(80, color='k', linestyle='--', alpha=0.3)
    ax[1].legend()
    
    out_path = os.path.join(FIGURES_DIR, 'sdss_test_dr_bd_desert.png')
    plt.savefig(out_path, dpi=150)
    print(f"Figure saved to {out_path}")
    
    return {
        't_stat': t_stat,
        'p_value_ttest': p_val,
        'slope_mass_feh': slope,
        'n_bd': int(len(bds)),
        'n_planet': int(len(planets))
    }

def main():
    cache_path = os.path.join(DATA_DIR, 'sdss_bd_desert.csv')
    
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

    results = analyze_bd_desert(df)
    
    out_path = os.path.join(RESULTS_DIR, 'sdss_test_dr_results.json')
    if results:
        with open(out_path, 'w') as f:
            json.dump(results, f, indent=2)
            
        print("\nSUMMARY TEST DR:")
        
        # If p < 0.05, distributions are different.
        if results['p_value_ttest'] < 0.05:
             print("RESULT: SIGNAL (BD hosts distinct from Planet hosts)")
        else:
             print("RESULT: NULL (No significant difference)")

if __name__ == "__main__":
    main()
