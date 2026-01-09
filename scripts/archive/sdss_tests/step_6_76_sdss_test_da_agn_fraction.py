#!/usr/bin/env python3
"""
Step 6.76: SDSS Test DA - AGN Type 1/2 Fraction (Torus Geometry)

Hypothesis:
The unified model explains Type 1 vs Type 2 AGN via viewing angle and an obscuring torus. 
The torus scale height depends on vertical velocity support. If dynamical heating/velocities 
are affected by TEP (phantom mass), the opening angle of the torus might vary with 
potential depth (sigma), changing the observed ratio of Type 1 to Type 2 AGN.

Prediction:
Fraction of Type 1 AGN (Broad Line) varies with Host Sigma.

Data:
- SpecObjAll: subClass (contains 'Broad' for Type 1), z
- emissionLinesPort: sigma_stars

Method:
1. Select objects classified as AGN (Broad vs Narrow).
2. Bin by Host Sigma.
3. Calculate Fraction(Type 1) = N(Broad) / N(Total).
4. Analyze trend.
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
    print(f"Querying SDSS for Test DA (Limit: {limit})...")
    
    sql = f"""
    SELECT TOP {limit}
        e.sigma_stars,
        CASE WHEN s.subClass LIKE '%Broad%' THEN 1 ELSE 0 END as is_broad,
        s.z
    FROM SpecObjAll s
    JOIN emissionLinesPort e ON s.specObjID = e.specObjID
    WHERE (s.class = 'QSO' OR s.class = 'GALAXY')
      AND s.subClass LIKE '%AGN%'
      AND e.sigma_stars > 0
      AND e.sigma_stars < 400 -- Exclude bad fits
    """
    return query_sdss(sql)

def analyze_agn_fraction(df):
    print("Analyzing AGN Type 1/2 Fraction...")
    
    if df is None or len(df) < 100:
        print("  Insufficient data.")
        return None
        
    df = df.dropna().copy()
    
    # Bin by Sigma
    df['sigma_bin'] = pd.qcut(df['sigma_stars'], q=10, labels=False)
    
    grouped = df.groupby('sigma_bin').agg({
        'sigma_stars': 'mean',
        'is_broad': ['mean', 'sem', 'count']
    }).reset_index()
    grouped.columns = ['bin', 'sigma', 'frac_type1', 'sem', 'count']
    
    # Regression
    slope, intercept, r_val, p_val, std_err = stats.linregress(grouped['sigma'], grouped['frac_type1'])
    
    print(f"  Correlation (Sigma vs Type 1 Frac): r={r_val:.3f}, p={p_val:.2e}, slope={slope:.5f}")
    
    # Plot
    plt.figure(figsize=(8, 6))
    plt.errorbar(grouped['sigma'], grouped['frac_type1'], yerr=grouped['sem'], fmt='o-', color='crimson')
    
    # Fit line
    x_range = np.linspace(grouped['sigma'].min(), grouped['sigma'].max(), 100)
    plt.plot(x_range, intercept + slope*x_range, 'k--', label=f'r={r_val:.2f}, p={p_val:.3f}')
    
    plt.xlabel('Host Velocity Dispersion $\sigma$ (km/s)')
    plt.ylabel('Fraction of Broad-Line (Type 1) AGN')
    plt.title('Test DA: AGN Type Fraction vs Potential Depth')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    out_path = os.path.join(FIGURES_DIR, 'sdss_test_da_agn_fraction.png')
    plt.savefig(out_path, dpi=150)
    print(f"Figure saved to {out_path}")
    
    return {
        'correlation_r': r_val,
        'p_value': p_val,
        'slope': slope,
        'n_galaxies': int(len(df))
    }

def main():
    cache_path = os.path.join(DATA_DIR, 'sdss_agn_types.csv')
    
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

    results = analyze_agn_fraction(df)
    
    out_path = os.path.join(RESULTS_DIR, 'sdss_test_da_results.json')
    if results:
        with open(out_path, 'w') as f:
            json.dump(results, f, indent=2)
            
        print("\nSUMMARY TEST DA:")
        print(f"Correlation (Sigma vs Type 1 Frac): {results['correlation_r']:.3f}")
        
        if results['p_value'] < 0.05 and abs(results['correlation_r']) > 0.1:
             print("RESULT: SIGNAL (AGN Type Fraction depends on Sigma)")
        else:
             print("RESULT: NULL (No significant dependence)")

if __name__ == "__main__":
    main()
