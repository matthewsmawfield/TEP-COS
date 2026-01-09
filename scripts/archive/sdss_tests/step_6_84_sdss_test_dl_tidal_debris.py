#!/usr/bin/env python3
"""
Step 6.84: SDSS Test DL - Tidal Debris Lifetime (Phase Mixing)

Hypothesis:
Tidal tails and debris fade due to phase mixing (wrapping). This is a dynamical rate process. 
In deep potentials, phase mixing is time-dilated (slower). Tidal features should remain 
visible for longer after the interaction event. We expect a higher fraction of galaxies 
showing tidal features in high-sigma environments (at fixed merger rate).

Prediction:
Fraction of galaxies with Tidal Features increases with Sigma.

Data:
- zoo2MainSpecz: t08_odd_feature_a24_merger_debiased, t08_odd_feature_a22_disturbed_debiased
- emissionLinesPort: sigma_stars
- stellarMassFSPSGranWideDust: logMass

Method:
1. Download Zoo 2, emissionLines, and stellarMass tables separately (to avoid join timeouts).
2. Perform Client-Side Join on dr7objid / specObjID.
3. Define "Tidal Feature" galaxy: p_merger > threshold or p_disturbed > threshold.
4. Calculate Fraction(Tidal) in bins of Sigma.
5. Control for Mass.
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

def download_data():
    print("Downloading data for Test DL (Client-Side Join Strategy)...")
    
    # 1. Zoo 2
    # Limit to reliable votes
    print("  Querying zoo2MainSpecz...")
    sql_zoo = """
    SELECT TOP 2000
        dr7objid as specObjID,
        t08_odd_feature_a24_merger_debiased as p_merger,
        t08_odd_feature_a22_disturbed_debiased as p_disturbed,
        total_votes
    FROM zoo2MainSpecz
    WHERE total_votes > 20
    """
    df_zoo = query_sdss(sql_zoo)
    
    if df_zoo is None:
        print("  Failed to download Zoo data.")
        return None
        
    print(f"  Zoo data: {len(df_zoo)} rows.")
    
    # Get list of IDs for WHERE clause to optimize join? 
    # Or just download matching EmissionLines and StellarMass for a random set?
    # Better: Download a chunk of EmissionLines/Mass and join in memory.
    # Since we can't upload IDs, we have to rely on overlapping selection or download large chunks.
    # Strategy: Download 2000 random from emissionLinesPort with Sigma > 0 and Mass > 9
    # But probability of overlap with top 2000 Zoo is low if ordered differently.
    # Ordering by specObjID helps.
    
    # Re-query Zoo sorted by specObjID to sync
    sql_zoo = """
    SELECT TOP 5000
        dr7objid as specObjID,
        t08_odd_feature_a24_merger_debiased as p_merger,
        t08_odd_feature_a22_disturbed_debiased as p_disturbed
    FROM zoo2MainSpecz
    WHERE total_votes > 20
    ORDER BY dr7objid
    """
    df_zoo = query_sdss(sql_zoo)
    if df_zoo is None: return None
    
    min_id = df_zoo['specObjID'].min()
    max_id = df_zoo['specObjID'].max()
    
    # 2. EmissionLines (Sigma)
    print("  Querying emissionLinesPort...")
    sql_em = f"""
    SELECT specObjID, sigma_stars
    FROM emissionLinesPort
    WHERE specObjID BETWEEN {min_id} AND {max_id}
      AND sigma_stars > 0
    """
    df_em = query_sdss(sql_em)
    
    # 3. StellarMass
    print("  Querying stellarMass...")
    sql_mass = f"""
    SELECT specObjID, logMass
    FROM stellarMassFSPSGranWideDust
    WHERE specObjID BETWEEN {min_id} AND {max_id}
      AND logMass > 0
    """
    df_mass = query_sdss(sql_mass)
    
    # Join
    print("  Joining tables...")
    if df_em is not None and df_mass is not None:
        df = pd.merge(df_zoo, df_em, on='specObjID', how='inner')
        df = pd.merge(df, df_mass, on='specObjID', how='inner')
        print(f"  Final joined set: {len(df)} rows.")
        return df
    else:
        print("  Failed to download auxiliary tables.")
        return None

def analyze_tidal_debris(df):
    print("Analyzing Tidal Debris Lifetime...")
    
    if df is None or len(df) < 100:
        print("  Insufficient data.")
        return None
        
    df = df.dropna().copy()
    
    # Define Tidal Feature
    # Either merger or disturbed
    df['p_tidal'] = df[['p_merger', 'p_disturbed']].max(axis=1)
    
    # Threshold
    threshold = 0.4
    df['is_tidal'] = (df['p_tidal'] > threshold).astype(int)
    
    print("  Tidal Galaxy Counts:")
    print(df['is_tidal'].value_counts())
    
    n_tidal = df['is_tidal'].sum()
    if n_tidal < 10:
        print("  Not enough tidal galaxies.")
        return None
        
    # Bin by Sigma
    df['sigma_bin'] = pd.qcut(df['sigma_stars'], q=10, labels=False, duplicates='drop')
    
    grouped = df.groupby('sigma_bin').agg({
        'sigma_stars': 'mean',
        'is_tidal': ['mean', 'sem', 'count']
    }).reset_index()
    grouped.columns = ['bin', 'sigma', 'frac', 'sem', 'count']
    
    # Regression
    slope, intercept, r_val, p_val, std_err = stats.linregress(grouped['sigma'], grouped['frac'])
    
    print(f"  Correlation (Sigma vs Tidal Fraction): r={r_val:.3f}, p={p_val:.2e}, slope={slope:.5f}")
    
    # Plot
    plt.figure(figsize=(8, 6))
    plt.errorbar(grouped['sigma'], grouped['frac'], yerr=grouped['sem'], fmt='o-', color='magenta')
    
    x_range = np.linspace(grouped['sigma'].min(), grouped['sigma'].max(), 100)
    plt.plot(x_range, intercept + slope*x_range, 'k--', label=f'r={r_val:.2f}')
    
    plt.xlabel('Velocity Dispersion (km/s)')
    plt.ylabel('Fraction of Disturbed/Merger Galaxies')
    plt.title('Test DL: Tidal Debris Lifetime')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    out_path = os.path.join(FIGURES_DIR, 'sdss_test_dl_tidal_debris.png')
    plt.savefig(out_path, dpi=150)
    print(f"Figure saved to {out_path}")
    
    return {
        'correlation_r': r_val,
        'p_value': p_val,
        'slope': slope,
        'n_gal': int(len(df))
    }

def main():
    cache_path = os.path.join(DATA_DIR, 'sdss_tidal_debris.csv')
    
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

    results = analyze_tidal_debris(df)
    
    out_path = os.path.join(RESULTS_DIR, 'sdss_test_dl_results.json')
    if results:
        with open(out_path, 'w') as f:
            json.dump(results, f, indent=2)
            
        print("\nSUMMARY TEST DL:")
        print(f"Slope (Sigma vs Tidal Frac): {results['slope']:.5f}")
        
        if results['p_value'] < 0.05 and results['slope'] > 0:
             print("RESULT: SIGNAL (Higher Tidal Fraction in Deep Potentials)")
        elif results['p_value'] < 0.05 and results['slope'] < 0:
             print("RESULT: CONTRADICTED (Lower Tidal Fraction in Deep Potentials)")
        else:
             print("RESULT: NULL (No significant dependence)")

if __name__ == "__main__":
    main()
