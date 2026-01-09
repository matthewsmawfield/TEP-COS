#!/usr/bin/env python3
"""
Step 6.28: SDSS Test AU - Dust Survival (Sputtering Clock)

Hypothesis:
Dust grains are destroyed by sputtering in hot gas (rate process).
In deep potentials (massive elliptical halos), time dilation should slow down the sputtering rate.
Passive galaxies in high-sigma environments should retain a higher dust fraction than their low-sigma counterparts.
Dust fraction is measurable via Balmer decrement (Ha/Hb ratio) even in passive galaxies (weak lines).

Prediction:
For Passive Galaxies (D4000 > 1.8): H-alpha/H-beta ratio increases with sigma.
Standard: Passive galaxies are dust-free (Ratio ~ 2.86 or 3.1). Excess implies dust.

Data:
- emissionLinesPort: Ha, Hb fluxes, sigma_stars.
- galSpecIndx: d4000_n (Selection of passive).

Method:
1. Select Passive Galaxies (D4000 > 1.8, low EW(Ha)?).
2. Ensure S/N > 3 for Ha and Hb (might be hard for passive).
3. Compute Balmer Decrement: BD = Flux_Ha / Flux_Hb.
4. Correlate BD with sigma.
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

def download_data(limit=2000):
    print(f"Querying SDSS for Test AU (Limit: {limit})...")
    
    # Select passive galaxies with detectable emission lines
    # D4000 > 1.8 ensures old population.
    # We need emission lines to measure dust (from residual gas).
    # Many passive galaxies have LINER emission.
    
    sql = f"""
    SELECT TOP {limit}
        e.specObjID,
        e.sigmaStars as sigma,
        e.Flux_Ha_6562, e.Flux_Ha_6562_Err,
        e.Flux_Hb_4861, e.Flux_Hb_4861_Err,
        i.d4000_n
        
    FROM emissionLinesPort e
    JOIN galSpecIndx i ON e.specObjID = i.specObjID
    
    WHERE 
        e.sigmaStars > 50 AND e.sigmaStars < 400
        AND e.Flux_Ha_6562 > 0 AND e.Flux_Hb_4861 > 0
    """
    return query_sdss(sql)

def analyze_dust_survival(df):
    print("Analyzing Dust Survival (Balmer Decrement)...")
    
    # 1. Clean and Filter for Passive
    df_clean = df.dropna().copy()
    
    # Filter for Passive (D4000 > 1.8)
    df_clean = df_clean[df_clean['d4000_n'] > 1.8].copy()
    
    # S/N Cut
    df_clean = df_clean[
        (df_clean['Flux_Ha_6562']/df_clean['Flux_Ha_6562_Err'] > 3) &
        (df_clean['Flux_Hb_4861']/df_clean['Flux_Hb_4861_Err'] > 3)
    ].copy()
    
    df_clean['log_sigma'] = np.log10(df_clean['sigma'])
    
    # 2. Compute Balmer Decrement
    df_clean['balmer_dec'] = df_clean['Flux_Ha_6562'] / df_clean['Flux_Hb_4861']
    
    # Remove unphysical values
    # Theoretical min is ~2.86 (Case B). Allow some noise.
    df_clean = df_clean[(df_clean['balmer_dec'] > 2.0) & (df_clean['balmer_dec'] < 10.0)].copy()
    
    # 3. Correlation
    # TEP Prediction: BD increases with sigma (more dust survival)
    
    r_bd, p_bd = stats.pearsonr(df_clean['log_sigma'], df_clean['balmer_dec'])
    
    print(f"N = {len(df_clean)}")
    print(f"Mean Balmer Dec: {df_clean['balmer_dec'].mean():.2f}")
    print(f"Correlation r(BD, sigma): {r_bd:.4f} (p={p_bd:.2e})")
    
    # 4. Binning
    df_clean['sigma_bin'] = pd.qcut(df_clean['log_sigma'], 8)
    binned = df_clean.groupby('sigma_bin')['balmer_dec'].mean()
    print("\nMean Balmer Dec by Sigma Bin:")
    print(binned)
    
    return {
        'r_bd': float(r_bd),
        'p_bd': float(p_bd),
        'mean_bd': float(df_clean['balmer_dec'].mean()),
        'binned_means': binned.tolist(),
        'bin_centers': [mid.mid for mid in binned.index],
        'n_sample': int(len(df_clean))
    }, df_clean

def create_figure(df, results):
    print("Generating figure...")
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # Scatter
    ax.scatter(df['log_sigma'], df['balmer_dec'], alpha=0.1, s=2, c='k', label='Passive Galaxies')
    
    # Binned
    ax.plot(results['bin_centers'], results['binned_means'], 'r-o', lw=2, label='Mean Trend')
    
    ax.set_xlabel(r'$\log(\sigma)$')
    ax.set_ylabel(r'Balmer Decrement ($H\alpha / H\beta$)')
    ax.set_title(f"Test AU: Dust Survival (r={results['r_bd']:.3f})")
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.axhline(2.86, color='b', linestyle='--', label='Case B (Dust Free)')
    
    plt.tight_layout()
    out_path = os.path.join(FIGURES_DIR, 'sdss_test_au_dust_survival.png')
    plt.savefig(out_path, dpi=150)
    print(f"Figure saved to {out_path}")

def main():
    cache_path = os.path.join(DATA_DIR, 'sdss_dust_survival.csv')
    
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

    results, df_clean = analyze_dust_survival(df)
    
    out_path = os.path.join(RESULTS_DIR, 'sdss_test_au_results.json')
    with open(out_path, 'w') as f:
        results_json = results.copy()
        del results_json['bin_centers']
        json.dump(results_json, f, indent=2)
        
    create_figure(df_clean, results)
    
    print("\nSUMMARY TEST AU:")
    print("TEP Prediction: BD increases with sigma. r > 0.")
    print(f"Observed r: {results['r_bd']:.4f}")
    
    if results['r_bd'] > 0.05:
        print("RESULT: CONSISTENT (Higher dust content in deep potentials)")
    elif results['r_bd'] < -0.05:
        print("RESULT: CONTRADICTED (Less dust/cleaner in deep potentials)")
    else:
        print("RESULT: NULL/FLAT")

if __name__ == "__main__":
    main()
