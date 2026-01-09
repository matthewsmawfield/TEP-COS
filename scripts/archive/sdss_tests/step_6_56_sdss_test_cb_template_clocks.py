#!/usr/bin/env python3
"""
Step 6.56: SDSS Test CB - Template Clock Systematics (MaStar vs MILES)

Hypothesis:
Stellar Population models (Firefly) rely on empirical stellar libraries (templates).
MaStar uses stars from the SDSS footprint (including deep potential zones).
MILES uses local solar neighborhood stars (shallower potential).
If TEP is real, the "clocks" in these stars run at different rates.
Fitting the SAME galaxies with these two different libraries should yield systematic
age offsets that correlate with the target galaxy's sigma.

Prediction:
Delta_Age = Age(MaStar) - Age(MILES) correlates with Galaxy sigma.

Data:
- mangaFirefly_miles: LW_AGE_1RE
- mangaFirefly_mastar: LW_AGE_1RE
- mangaDAPall: stellar_sigma_1re

Method:
1. Join tables on mangaid.
2. Compute Delta_Age.
3. Correlate with sigma.
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

def download_data(limit=5000):
    print(f"Querying SDSS for Test CB (Limit: {limit})...")
    
    # We join on mangaid.
    # Note: Firefly tables might have multiple entries per mangaID if there are multiple plateifus?
    # Usually mangaid is unique per galaxy, but observations might differ.
    # We join on mangaid and IFUDSGN or just mangaid if usually 1-to-1.
    # Let's try simple join.
    
    sql = f"""
    SELECT TOP {limit}
        f1.mangaid,
        f1.LW_AGE_1RE as age_miles,
        f2.LW_AGE_1RE as age_mastar,
        d.stellar_sigma_1re as sigma,
        d.nsa_z as z
        
    FROM mangaFirefly_miles f1
    JOIN mangaFirefly_mastar f2 ON f1.mangaid = f2.mangaid AND f1.IFUDSGN = f2.IFUDSGN
    JOIN mangaDAPall d ON f1.mangaid = d.mangaid
    
    WHERE 
        d.drp3qual = 0
        AND f1.LW_AGE_1RE > 0
        AND f2.LW_AGE_1RE > 0
        AND d.stellar_sigma_1re > 0
    """
    return query_sdss(sql)

def analyze_templates(df):
    print("Analyzing Template Systematics...")
    
    df = df.dropna().copy()
    
    # Compute Delta Age
    df['delta_age'] = df['age_mastar'] - df['age_miles']
    df['log_sigma'] = np.log10(df['sigma'])
    
    print(f"  Sample size: {len(df)}")
    print(f"  Mean Delta Age (MaStar - MILES): {df['delta_age'].mean():.4f} Gyr")
    
    # Correlation
    r_val, p_val = stats.pearsonr(df['log_sigma'], df['delta_age'])
    print(f"  Correlation r(DeltaAge, logSigma): {r_val:.4f} (p={p_val:.2e})")
    
    slope, intercept, _, _, _ = stats.linregress(df['log_sigma'], df['delta_age'])
    print(f"  Slope: {slope:.4f} Gyr/dex")
    
    # Plot
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # Scatter
    sc = ax.scatter(df['log_sigma'], df['delta_age'], c=df['z'], cmap='viridis', alpha=0.5, s=10, label='Data')
    plt.colorbar(sc, label='Redshift')
    
    # Trend
    x_range = np.linspace(df['log_sigma'].min(), df['log_sigma'].max(), 100)
    ax.plot(x_range, slope * x_range + intercept, 'r-', label=f'Slope={slope:.2f}')
    
    ax.set_xlabel('log Sigma [km/s]')
    ax.set_ylabel('Age(MaStar) - Age(MILES) [Gyr]')
    ax.set_title(f'Test CB: Template Clock Systematics (r={r_val:.2f})')
    ax.axhline(0, color='k', linestyle='--', alpha=0.3)
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    out_path = os.path.join(FIGURES_DIR, 'sdss_test_cb_templates.png')
    plt.savefig(out_path, dpi=150)
    print(f"Figure saved to {out_path}")
    
    return {
        'r_val': r_val,
        'slope': slope,
        'mean_offset': df['delta_age'].mean(),
        'n_sample': int(len(df))
    }

def main():
    cache_path = os.path.join(DATA_DIR, 'sdss_templates.csv')
    
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

    results = analyze_templates(df)
    
    out_path = os.path.join(RESULTS_DIR, 'sdss_test_cb_results.json')
    with open(out_path, 'w') as f:
        json.dump(results, f, indent=2)
        
    print("\nSUMMARY TEST CB:")
    print("Prediction: Systematic age offset correlates with Sigma.")
    print(f"Observed r: {results['r_val']:.4f}")
    
    if abs(results['r_val']) > 0.1:
         print("RESULT: SIGNAL (Correlation observed)")
    else:
         print("RESULT: NULL (No significant correlation)")

if __name__ == "__main__":
    main()
