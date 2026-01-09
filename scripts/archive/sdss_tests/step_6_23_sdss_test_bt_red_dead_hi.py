#!/usr/bin/env python3
"""
Step 6.23: SDSS Test BT - Cold Gas in Quiescent Galaxies (Red & Dead HI)

Hypothesis:
Standard model: Quiescent galaxies are "red and dead" because they lack cold gas (strangulation).
TEP (RBH-1 Model): Quiescent galaxies may contain gas that is DYNAMICALLY STABILIZED by the soliton wake metric (reduced Jeans mass, prevents collapse) rather than thermally supported.
Prediction: A higher-than-expected detection rate of HI gas in high-sigma Quiescent galaxies.
Observable: HI Detection Rate in Red Sequence galaxies increases with sigma.

Data:
- mangaHIall: HI mass (logM_HI), Detection flag (conf_prob?).
- mangaDAPall: Stellar Sigma (stellar_sigma_1re).
- galSpecIndx: D4000 (to select Red Sequence).

Note: We need to link mangaHIall to DAP/Spec. HIall has 'mangaid'.
"""

import pandas as pd
import numpy as np
from scipy import stats
import os
import json
import matplotlib.pyplot as plt
import requests
import time
from sklearn.linear_model import LogisticRegression

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

def download_data(limit=1000):
    print(f"Querying SDSS for Test BT (Limit: {limit})...")
    
    # We select more galaxies and filter in Python to be safe with timeouts
    
    sql = f"""
    SELECT TOP {limit}
        h.mangaid,
        h.logMHI, 
        h.confprob,
        h.snr,
        d.stellar_sigma_1re as sigma,
        d.specindex_1re_d4000 as D4000
        
    FROM mangaHIall h
    JOIN mangaDAPall d ON h.mangaid = d.mangaid
    
    WHERE 
        d.stellar_sigma_1re > 50
    """
    return query_sdss(sql)

def analyze_hi_detection(df):
    print("Analyzing HI Detection in Red Sequence...")
    
    # 1. Clean and Filter for Red Sequence (D4000 > 1.6)
    df_clean = df.dropna(subset=['sigma', 'D4000', 'logMHI']).copy()
    df_clean = df_clean[df_clean['D4000'] > 1.6].copy()
    
    # 2. Define Detection
    # Usually confprob > 0.9 or SNR > 3
    # Let's use SNR > 3 as detection
    df_clean['detected'] = (df_clean['snr'] > 3.0).astype(int)
    
    print(f"N Total (Red Seq): {len(df_clean)}")
    print(f"N Detected HI: {df_clean['detected'].sum()}")
    print(f"Detection Rate: {df_clean['detected'].mean():.3f}")
    
    # 3. Bin by Sigma
    df_clean['log_sigma'] = np.log10(df_clean['sigma'])
    
    bins = np.linspace(df_clean['log_sigma'].min(), df_clean['log_sigma'].max(), 8)
    bin_centers = 0.5 * (bins[1:] + bins[:-1])
    rates = []
    
    for i in range(len(bins)-1):
        mask = (df_clean['log_sigma'] >= bins[i]) & (df_clean['log_sigma'] < bins[i+1])
        if mask.sum() > 20:
            rates.append(df_clean.loc[mask, 'detected'].mean())
        else:
            rates.append(np.nan)
            
    # 4. Logistic Regression
    # P(detected) ~ a*log_sigma
    # TEP predicts a > 0 (Higher detection in deep potentials)
    
    X = df_clean[['log_sigma']].values
    y = df_clean['detected'].values
    
    if len(np.unique(y)) > 1:
        clf = LogisticRegression(class_weight='balanced', solver='liblinear')
        clf.fit(X, y)
        coefs = clf.coef_[0]
        sigma_coef = coefs[0]
        mass_coef = 0.0 # Removed
    else:
        print("Warning: Only one class present (all detected or all non-detected).")
        sigma_coef = 0.0
        mass_coef = 0.0
        
    print(f"Logistic Coef (log_sigma): {sigma_coef:.4f}")
    
    return {
        'sigma_coef': float(sigma_coef),
        'mass_coef': float(mass_coef),
        'rates': rates,
        'bin_centers': bin_centers.tolist(),
        'n_sample': int(len(df_clean))
    }, df_clean

def create_figure(df, results):
    print("Generating figure...")
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # Plot Detection Rate vs Sigma
    ax.plot(results['bin_centers'], results['rates'], 'bo-', label='HI Detection Rate')
    
    ax.set_xlabel(r'$\log(\sigma)$')
    ax.set_ylabel('HI Detection Fraction (Red Sequence)')
    ax.set_title("Cold Gas Survival in Deep Potentials")
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    out_path = os.path.join(FIGURES_DIR, 'sdss_test_bt_red_hi.png')
    plt.savefig(out_path, dpi=150)
    print(f"Figure saved to {out_path}")

def main():
    cache_path = os.path.join(DATA_DIR, 'sdss_manga_hi_red.csv')
    
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

    results, df_clean = analyze_hi_detection(df)
    
    out_path = os.path.join(RESULTS_DIR, 'sdss_test_bt_results.json')
    with open(out_path, 'w') as f:
        results_json = results.copy()
        del results_json['bin_centers'] 
        json.dump(results_json, f, indent=2)
        
    create_figure(df_clean, results)
    
    print("\nSUMMARY TEST BT:")
    print("TEP Prediction: HI Detection increases with sigma (Coef > 0)")
    print(f"Logistic Coef (log_sigma): {results['sigma_coef']:.4f}")
    
    if results['sigma_coef'] > 0.2:
        print("RESULT: CONSISTENT (Higher gas fraction in deep potentials)")
    elif results['sigma_coef'] < -0.2:
        print("RESULT: CONTRADICTED (Gas stripped/depleted in deep potentials)")
    else:
        print("RESULT: NULL/FLAT")

if __name__ == "__main__":
    main()
