#!/usr/bin/env python3
"""
Step 6.8: Test C - BPT Position Anomaly

TEP HYPOTHESIS:
BPT classification uses emission line ratios that depend on ionization equilibrium 
timescales. Under TEP, time dilation in deeper potentials could shift the equilibrium 
position, causing high-σ galaxies to appear offset in BPT space relative to their 
expected position from stellar metallicity alone.

TEP PREDICTION:
  Define: BPT_offset = observed log([NII]/Hα) - predicted log([NII]/Hα) from stellar [Z/H]
  At fixed stellar metallicity:
    r(BPT_offset, σ) ≠ 0 (Likely positive, shifting towards LINER/AGN region)

DATA:
  - emissionLinesPort
  - galSpecIndx (stellar metallicity)

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
    print(f"Querying SDSS for Test C (Limit: {limit})...")
    sql = f"""
    SELECT TOP {limit}
        p.specObjID,
        p.z AS redshift,
        p.sigmaStars AS sigma_stars,
        
        -- Emission Lines
        p.Flux_NII_6583,
        p.Flux_Ha_6562,
        p.Flux_OIII_5006,
        p.Flux_Hb_4861,
        
        -- Stellar Metallicity
        i.lick_mgb,
        i.lick_fe5270,
        i.lick_fe5335,
        
        -- Mass
        s.logMass
        
    FROM emissionLinesPort p
    JOIN galSpecIndx i ON p.specObjID = i.specObjID
    JOIN stellarMassFSPSGranWideDust s ON p.specObjID = s.specObjID
    
    WHERE 
        p.z BETWEEN 0.02 AND 0.20
        AND p.sigmaStars > 50 AND p.sigmaStars < 400
        AND p.Flux_NII_6583 > 0 AND p.Flux_Ha_6562 > 0
        AND p.Flux_OIII_5006 > 0 AND p.Flux_Hb_4861 > 0
        AND i.lick_mgb > 0
    """
    return query_sdss(sql)

def analyze_bpt_anomaly(df):
    print("Analyzing BPT anomalies...")
    
    # 1. Compute Ratios
    df['log_NII_Ha'] = np.log10(df['Flux_NII_6583'] / df['Flux_Ha_6562'])
    df['log_OIII_Hb'] = np.log10(df['Flux_OIII_5006'] / df['Flux_Hb_4861'])
    df['log_sigma'] = np.log10(df['sigma_stars'])
    
    # Stellar metallicity proxy: [MgFe]' or similar
    # [MgFe]' = sqrt(Mgb * (0.72*Fe5270 + 0.28*Fe5335)) is common, 
    # or just Mgb/<Fe> for alpha. Let's use Mgb as simple Z proxy for now.
    df['stellar_Z_proxy'] = df['lick_mgb'] 
    
    # 2. Predict log([NII]/Hα) from Stellar Metallicity and Mass
    # We expect nebular metallicity to track stellar metallicity.
    # We fit a relation to remove the "standard" dependence.
    
    # Only fit on Star Forming galaxies to establish the baseline?
    # SF criteria (Kauffmann et al. 2003): log([OIII]/Hbeta) < 0.61 / (log([NII]/Halpha) - 0.05) + 1.3
    # Approximate SF cut for simplicity in fitting: log([NII]/Hα) < -0.2
    
    sf_subset = df[df['log_NII_Ha'] < -0.2]
    if len(sf_subset) > 100:
        slope_Z, intercept_Z = np.polyfit(sf_subset['stellar_Z_proxy'], sf_subset['log_NII_Ha'], 1)
        print(f"Baseline Relation: log(NII/Ha) = {slope_Z:.3f} * Mgb + {intercept_Z:.3f}")
    else:
        slope_Z, intercept_Z = 0, -0.5 # Fallback
        
    df['predicted_NII_Ha'] = slope_Z * df['stellar_Z_proxy'] + intercept_Z
    df['bpt_offset'] = df['log_NII_Ha'] - df['predicted_NII_Ha']
    
    # 3. Correlate Offset with Sigma
    r_offset, p_offset = stats.pearsonr(df['log_sigma'], df['bpt_offset'])
    
    print(f"Correlation r(BPT Offset, σ): {r_offset:.4f} (p={p_offset:.2e})")
    
    return {
        'r_offset': float(r_offset),
        'p_offset': float(p_offset),
        'slope_Z': float(slope_Z),
        'intercept_Z': float(intercept_Z)
    }

def create_figure(df, results):
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    # BPT Diagram colored by Sigma
    ax = axes[0]
    sc = ax.scatter(df['log_NII_Ha'], df['log_OIII_Hb'], c=df['log_sigma'], 
                    cmap='viridis', s=2, alpha=0.5)
    plt.colorbar(sc, ax=ax, label=r'$\log(\sigma)$')
    ax.set_xlabel(r'$\log([NII]/H\alpha)$')
    ax.set_ylabel(r'$\log([OIII]/H\beta)$')
    ax.set_title('BPT Diagram (Sigma Colored)')
    
    # Offset vs Sigma
    ax = axes[1]
    ax.hexbin(df['log_sigma'], df['bpt_offset'], gridsize=40, cmap='Blues', mincnt=1)
    ax.set_xlabel(r'$\log(\sigma)$')
    ax.set_ylabel(r'BPT Offset (Obs - Pred from Z_star)')
    ax.set_title(f'BPT Offset vs $\sigma$ (r={results["r_offset"]:.3f})')
    
    plt.tight_layout()
    out_file = os.path.join(FIGURES_DIR, 'sdss_test_c_bpt_anomaly.png')
    plt.savefig(out_file, dpi=150)
    print(f"Figure saved to {out_file}")

def main():
    print("="*60)
    print("STEP 6.8: TEST C - BPT POSITION ANOMALY")
    print("="*60)
    
    cache_path = os.path.join(DATA_DIR, 'sdss_bpt_data.csv')
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
    
    results = analyze_bpt_anomaly(df)
    create_figure(df, results)
    
    out_path = os.path.join(RESULTS_DIR, 'sdss_test_c_results.json')
    with open(out_path, 'w') as f:
        json.dump(results, f, indent=2)
        
    print("\nSUMMARY:")
    print(f"Significant Correlation: {results['p_offset'] < 0.05}")
    print(f"Direction: {'Positive' if results['r_offset'] > 0 else 'Negative'}")

if __name__ == "__main__":
    main()
