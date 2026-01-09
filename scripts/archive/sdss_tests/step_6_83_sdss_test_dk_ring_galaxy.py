#!/usr/bin/env python3
"""
Step 6.83: SDSS Test DK - Ring Galaxy Fraction (Collisional Clock)

Hypothesis:
Collisional Ring Galaxies form via the passage of a companion through the disk. 
The ring expands and fades on a dynamical timescale. In deep potentials, dynamical time 
is dilated. Rings should persist longer (in cosmic time) in high-sigma galaxies 
compared to low-sigma ones. We expect a higher fraction of observable rings 
in massive/high-sigma hosts.

Prediction:
Fraction of Ring Galaxies increases with Velocity Dispersion.

Data:
- zoo2MainSpecz: ring (0/1 or vote fraction), total_votes
- emissionLinesPort: sigma_stars
- stellarMassFSPSGranWideDust: logMass

Method:
1. Join Zoo 2 table with physical properties.
2. Define Ring Galaxy: ring flag or high vote fraction.
3. Calculate Fraction(Ring) in bins of Sigma.
4. Control for Mass? (Higher mass -> larger disks -> more rings?).
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
    print(f"Querying SDSS for Test DK (Limit: {limit})...")
    
    # Zoo 2: 't08_odd_feature_a24_ring_debiased' > threshold? Or 'ring' column?
    # Query Plan uses 'ring' but checking columns is safer. 
    # Let's try simplified query first.
    
    sql = f"""
    SELECT TOP {limit}
        z.dr7objid as specObjID,
        z.t08_odd_feature_a24_ring_debiased as p_ring,
        z.total_votes,
        e.sigma_stars,
        s.logMass
    FROM zoo2MainSpecz z
    JOIN emissionLinesPort e ON z.dr7objid = e.specObjID
    JOIN stellarMassFSPSGranWideDust s ON z.dr7objid = s.specObjID
    WHERE z.total_votes > 20
      AND e.sigma_stars > 0
      AND s.logMass > 9
    """
    return query_sdss(sql)

def analyze_ring_fraction(df):
    print("Analyzing Ring Galaxy Fraction...")
    
    if df is None or len(df) < 100:
        print("  Insufficient data.")
        return None
        
    df = df.dropna().copy()
    
    print(f"  Data points: {len(df)}")
    
    # Define Ring Galaxy
    # p_ring is probability. Use threshold.
    threshold = 0.5
    df['is_ring'] = (df['p_ring'] > threshold).astype(int)
    
    print("  Ring Galaxy Counts:")
    print(df['is_ring'].value_counts())
    
    n_rings = df['is_ring'].sum()
    if n_rings < 10:
        print("  Not enough ring galaxies for analysis.")
        return None
        
    # Bin by Sigma
    df['sigma_bin'] = pd.qcut(df['sigma_stars'], q=10, labels=False)
    
    grouped = df.groupby('sigma_bin').agg({
        'sigma_stars': 'mean',
        'is_ring': ['mean', 'sem', 'count']
    }).reset_index()
    grouped.columns = ['bin', 'sigma', 'frac', 'sem', 'count']
    
    # Regression
    slope, intercept, r_val, p_val, std_err = stats.linregress(grouped['sigma'], grouped['frac'])
    
    print(f"  Correlation (Sigma vs Ring Fraction): r={r_val:.3f}, p={p_val:.2e}, slope={slope:.5f}")
    
    # Plot
    plt.figure(figsize=(8, 6))
    plt.errorbar(grouped['sigma'], grouped['frac'], yerr=grouped['sem'], fmt='o-', color='orange')
    
    x_range = np.linspace(grouped['sigma'].min(), grouped['sigma'].max(), 100)
    plt.plot(x_range, intercept + slope*x_range, 'k--', label=f'r={r_val:.2f}')
    
    plt.xlabel('Velocity Dispersion (km/s)')
    plt.ylabel('Fraction of Ring Galaxies (p_ring > 0.5)')
    plt.title('Test DK: Ring Galaxy Lifetime')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    out_path = os.path.join(FIGURES_DIR, 'sdss_test_dk_ring_fraction.png')
    plt.savefig(out_path, dpi=150)
    print(f"Figure saved to {out_path}")
    
    return {
        'correlation_r': r_val,
        'p_value': p_val,
        'slope': slope,
        'n_gal': int(len(df)),
        'n_ring': int(n_rings)
    }

def main():
    cache_path = os.path.join(DATA_DIR, 'sdss_ring_fraction.csv')
    
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

    results = analyze_ring_fraction(df)
    
    out_path = os.path.join(RESULTS_DIR, 'sdss_test_dk_results.json')
    if results:
        with open(out_path, 'w') as f:
            json.dump(results, f, indent=2)
            
        print("\nSUMMARY TEST DK:")
        print(f"Slope (Sigma vs Ring Frac): {results['slope']:.5f}")
        
        if results['p_value'] < 0.05 and results['slope'] > 0:
             print("RESULT: SIGNAL (Higher Ring Fraction in Deep Potentials)")
        elif results['p_value'] < 0.05 and results['slope'] < 0:
             print("RESULT: CONTRADICTED (Lower Ring Fraction in Deep Potentials)")
        else:
             print("RESULT: NULL (No significant dependence)")

if __name__ == "__main__":
    main()
