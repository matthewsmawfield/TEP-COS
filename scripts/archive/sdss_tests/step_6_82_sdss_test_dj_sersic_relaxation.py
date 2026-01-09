#!/usr/bin/env python3
"""
Step 6.82: SDSS Test DJ - Sersic Relaxation (Dynamical Clock)

Hypothesis:
Galaxies relax into high-Sersic index (n > 2) profiles via violent relaxation and mergers. 
These are dynamical processes. In deep potentials, relaxation is time-dilated. 
High-sigma galaxies might preserve lower-Sersic (disk-like) profiles for longer, 
or the n-sigma relation should show an offset compared to standard predictions.

Prediction:
Sersic Index n is lower than expected for high-sigma galaxies (at fixed Mass).

Data:
- mangaTarget: nsa_sersic_n, nsa_elpetro_mass
- mangaDAPall: stellar_sigma_1re

Method:
1. Join mangaTarget and mangaDAPall.
2. Control for Stellar Mass (Sersic index n strongly correlates with Mass).
3. Analyze residuals of n vs Sigma.
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

def download_data(limit=1000):
    print(f"Querying SDSS for Test DJ (Limit: {limit})...")
    
    sql = f"""
    SELECT TOP {limit}
        s.mangaid, 
        s.nsa_sersic_n as sersic_n,
        d.stellar_sigma_1re as sigma,
        s.nsa_elpetro_mass as mass
    FROM mangaTarget s
    JOIN mangaDAPall d ON s.mangaid = d.mangaid
    WHERE d.drp3qual = 0
      AND s.nsa_sersic_n > 0
      AND d.stellar_sigma_1re > 0
      AND s.nsa_elpetro_mass > 0
    """
    return query_sdss(sql)

def analyze_sersic_relaxation(df):
    print("Analyzing Sersic Relaxation...")
    
    if df is None or len(df) < 50:
        print("  Insufficient data.")
        return None
        
    df = df.dropna().copy()
    
    # Log quantities
    df['log_n'] = np.log10(df['sersic_n'])
    df['log_sigma'] = np.log10(df['sigma'])
    df['log_mass'] = np.log10(df['mass'])
    
    # 1. Control for Mass
    # Sersic n increases with Mass (bulge growth).
    from sklearn.linear_model import LinearRegression
    X = df[['log_mass']]
    y = df['log_n']
    
    reg = LinearRegression().fit(X, y)
    print(f"  Control Fit (Mass -> Sersic n) R2: {reg.score(X, y):.3f}")
    
    df['n_resid'] = y - reg.predict(X)
    
    # 2. Correlate residuals with Sigma
    # Prediction: High Sigma -> Lower n (Negative Slope)
    slope, intercept, r_val, p_val, std_err = stats.linregress(df['log_sigma'], df['n_resid'])
    
    print(f"  Correlation (log Sigma vs n Residual): r={r_val:.3f}, p={p_val:.2e}, slope={slope:.4f}")
    
    # Plot
    fig, ax = plt.subplots(1, 2, figsize=(12, 5))
    
    # Raw n vs Mass (color=sigma)
    sc = ax[0].scatter(df['log_mass'], df['log_n'], c=df['log_sigma'], cmap='viridis', s=10, alpha=0.6)
    plt.colorbar(sc, ax=ax[0], label='log(Sigma)')
    ax[0].set_xlabel('log(Stellar Mass)')
    ax[0].set_ylabel('log(Sersic n)')
    ax[0].set_title('Sersic Index vs Mass')
    
    # Residuals vs Sigma
    ax[1].scatter(df['log_sigma'], df['n_resid'], alpha=0.5, s=10, c='teal')
    
    x_range = np.linspace(df['log_sigma'].min(), df['log_sigma'].max(), 100)
    ax[1].plot(x_range, intercept + slope*x_range, 'k--', label=f'Slope={slope:.3f}')
    
    ax[1].set_xlabel('log(Velocity Dispersion)')
    ax[1].set_ylabel('Sersic n Residual (log)')
    ax[1].set_title('Test DJ: Sersic Relaxation vs Potential')
    ax[1].legend()
    ax[1].grid(True, alpha=0.3)
    
    out_path = os.path.join(FIGURES_DIR, 'sdss_test_dj_sersic.png')
    plt.savefig(out_path, dpi=150)
    print(f"Figure saved to {out_path}")
    
    return {
        'correlation_r': r_val,
        'p_value': p_val,
        'slope': slope,
        'n_gal': int(len(df))
    }

def main():
    cache_path = os.path.join(DATA_DIR, 'sdss_sersic_relax.csv')
    
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

    results = analyze_sersic_relaxation(df)
    
    out_path = os.path.join(RESULTS_DIR, 'sdss_test_dj_results.json')
    if results:
        with open(out_path, 'w') as f:
            json.dump(results, f, indent=2)
            
        print("\nSUMMARY TEST DJ:")
        print(f"Slope (Sigma vs n Resid): {results['slope']:.4f}")
        
        # Prediction: Lower n in high sigma -> Negative Slope
        if results['p_value'] < 0.05 and results['slope'] < -0.05:
             print("RESULT: SIGNAL (Lower Sersic n in Deep Potentials)")
        elif results['p_value'] < 0.05 and results['slope'] > 0.05:
             print("RESULT: CONTRADICTED (Higher Sersic n in Deep Potentials)")
        else:
             print("RESULT: NULL (No significant dependence)")

if __name__ == "__main__":
    main()
