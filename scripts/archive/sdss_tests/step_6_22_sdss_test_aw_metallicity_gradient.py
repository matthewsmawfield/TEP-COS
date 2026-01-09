#!/usr/bin/env python3
"""
Step 6.22: SDSS Test AW - Metallicity Gradient Flattening

Hypothesis:
Standard inside-out formation creates steep negative metallicity gradients (centers enriched first).
TEP predicts centers experienced less proper time, potentially delaying enrichment relative to the outskirts (in cosmic time terms).
Prediction: This would FLATTEN the observed metallicity gradient in massive (deep potential) galaxies.
Observable: Metallicity Gradient (dex/Re) becomes flatter (less negative, closer to 0) as sigma increases.

Data:
- mangaFirefly: LW_Z_GRADIENT (Light-Weighted Metallicity Gradient).
- mangaDAPall: Stellar Sigma (stellar_sigma_1re).

Query:
SELECT f.mangaid, f.LW_Z_GRADIENT, d.stellar_sigma_1re
FROM mangaFirefly f JOIN mangaDAPall d ON f.PLATEIFU = d.plateifu
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
    print(f"Querying SDSS for Test AW (Limit: {limit})...")
    
    # Joining mangaTarget for mass
    sql = f"""
    SELECT TOP {limit}
        d.mangaid,
        f.LW_Z_GRADIENT as z_grad,
        d.stellar_sigma_1re as sigma
        
    FROM mangaFirefly_miles f
    JOIN mangaDAPall d ON f.PLATEIFU = d.plateifu
    
    WHERE 
        f.LW_Z_GRADIENT > -10 AND f.LW_Z_GRADIENT < 10 -- Reasonable bounds
        AND d.stellar_sigma_1re > 50
        AND d.drp3qual = 0
    """
    return query_sdss(sql)

def analyze_gradients(df):
    print("Analyzing Metallicity Gradients...")
    print(f"Columns: {df.columns.tolist()}")
    if not df.empty:
        print(df.head())
    
    # 1. Clean
    df_clean = df.dropna(subset=['z_grad', 'sigma']).copy()
    df_clean['log_sigma'] = np.log10(df_clean['sigma'])
    
    print(f"N = {len(df_clean)}")
    print(f"Mean Gradient: {df_clean['z_grad'].mean():.3f} dex/Re")
    
    # 2. Correlation
    # TEP Prediction: Gradient becomes flatter (closer to 0) as sigma increases.
    # Since gradients are typically negative, this means becoming LESS negative (increasing algebraically).
    # So we expect r(z_grad, sigma) > 0.
    
    r_simple, p_simple = stats.pearsonr(df_clean['log_sigma'], df_clean['z_grad'])
    
    print(f"r(Z_grad, log_sigma): {r_simple:.4f} (p={p_simple:.2e})")
    
    # 3. Bin Analysis
    df_clean['sigma_bin'] = pd.qcut(df_clean['log_sigma'], 5)
    binned = df_clean.groupby('sigma_bin')['z_grad'].mean()
    print("\nMean Z-Gradient by Sigma Bin:")
    print(binned)
    
    return {
        'r_sigma': float(r_simple),
        'p_sigma': float(p_simple),
        'mean_grad': float(df_clean['z_grad'].mean()),
        'n_sample': int(len(df_clean)),
        'binned_means': binned.tolist(),
        'bin_centers': [mid.mid for mid in binned.index]
    }, df_clean

def create_figure(df, results):
    print("Generating figure...")
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # Plot Grad vs Sigma
    ax.scatter(df['log_sigma'], df['z_grad'], alpha=0.1, s=2, c='purple')
    
    # Trend line
    m, b = np.polyfit(df['log_sigma'], df['z_grad'], 1)
    x = np.linspace(df['log_sigma'].min(), df['log_sigma'].max(), 100)
    ax.plot(x, m*x + b, 'r-', lw=2, label=f'r={results["r_sigma"]:.3f}')
    
    ax.set_xlabel(r'$\log(\sigma)$')
    ax.set_ylabel(r'Metallicity Gradient ($\nabla_Z$ [dex/$R_e$])')
    ax.set_title("Metallicity Gradient vs Potential Depth")
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.axhline(0, color='k', linestyle='--', alpha=0.5)
    
    plt.tight_layout()
    out_path = os.path.join(FIGURES_DIR, 'sdss_test_aw_metallicity_gradient.png')
    plt.savefig(out_path, dpi=150)
    print(f"Figure saved to {out_path}")

def main():
    cache_path = os.path.join(DATA_DIR, 'sdss_metallicity_gradients.csv')
    
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

    results, df_clean = analyze_gradients(df)
    
    out_path = os.path.join(RESULTS_DIR, 'sdss_test_aw_results.json')
    with open(out_path, 'w') as f:
        # Convert Intervals to string for JSON serialization
        results_json = results.copy()
        del results_json['bin_centers'] # Simplify or fix serialization if needed
        json.dump(results_json, f, indent=2)
        
    create_figure(df_clean, results)
    
    print("\nSUMMARY TEST AW:")
    print("TEP Prediction: Gradients flatten (become less negative) at high sigma. r > 0.")
    print(f"Observed r: {results['r_sigma']:.4f}")
    
    if results['r_sigma'] > 0.05:
        print("RESULT: CONSISTENT (Gradients flatten in deep potentials)")
    elif results['r_sigma'] < -0.05:
        print("RESULT: CONTRADICTED (Gradients steepen)")
    else:
        print("RESULT: NULL/FLAT")

if __name__ == "__main__":
    main()
