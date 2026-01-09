#!/usr/bin/env python3
"""
Step 6.40: SDSS Test BI - Fundamental Metallicity Relation (FMR) Warp

Hypothesis:
The FMR relates Mass, SFR, and Metallicity. Both SFR (rate) and Metallicity (integrated rate) are affected by TEP.
The "surface" defined by the FMR should show a warp dependent on sigma, as the rate-based axes scale with sqrt(A(phi)).
Standard FMR: Z = f(M, SFR).
TEP Prediction: Residual Delta_Z = Z_obs - Z_pred(M, SFR) correlates with sigma.

Data:
- galSpecExtra: sfr_tot_p50 (log SFR), oh_p50 (Metallicity, 12+log(O/H)).
- emissionLinesPort: sigmaStars.
- stellarMassFSPSGranWideDust: logMass.

Method:
1. Select Star-Forming galaxies (SFR > -99, sigma > 0).
2. Fit FMR plane: Z = a*logM + b*logSFR + c*(logSFR)^2 + d?
   Or simpler: Mannucci et al. (2010) projection mu_alpha = logM - alpha * logSFR.
   Find alpha that minimizes scatter in Z vs mu_alpha.
   Then Fit Z vs mu_alpha.
3. Compute Residuals: Delta_Z = Z_obs - Z_fit.
4. Correlate Delta_Z with sigma (or log sigma).
"""

import pandas as pd
import numpy as np
from scipy import stats, optimize
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
    print(f"Querying SDSS for Test BI (Limit: {limit})...")
    
    sql = f"""
    SELECT TOP {limit}
        e.specObjID,
        e.sigmaStars as sigma,
        gx.oh_p50 as metallicity,
        gx.sfr_tot_p50 as log_sfr,
        s.logMass
        
    FROM galSpecExtra gx
    JOIN emissionLinesPort e ON gx.specObjID = e.specObjID
    JOIN stellarMassFSPSGranWideDust s ON gx.specObjID = s.specObjID
    
    WHERE 
        e.sigmaStars > 30
        AND gx.sfr_tot_p50 > -5
        AND gx.oh_p50 > 0
        AND s.logMass > 8.0
    """
    return query_sdss(sql)

def analyze_fmr_warp(df):
    print("Analyzing FMR Warp...")
    
    # 1. Clean
    df_clean = df.dropna().copy()
    
    # 2. Fit FMR
    # Mannucci projection: mu = logM - alpha * logSFR
    # Minimize scatter in Z vs mu
    
    def scatter_metric(alpha, data):
        mu = data['logMass'] - alpha * data['log_sfr']
        # Fit polynomial Z vs mu
        coeffs = np.polyfit(mu, data['metallicity'], 2)
        poly = np.poly1d(coeffs)
        resid = data['metallicity'] - poly(mu)
        return np.std(resid)
    
    res = optimize.minimize_scalar(scatter_metric, args=(df_clean,), bounds=(0.0, 1.0), method='bounded')
    best_alpha = res.x
    print(f"Best alpha for FMR: {best_alpha:.3f}")
    
    # Calculate residuals with best alpha
    df_clean['mu'] = df_clean['logMass'] - best_alpha * df_clean['log_sfr']
    coeffs = np.polyfit(df_clean['mu'], df_clean['metallicity'], 2)
    poly = np.poly1d(coeffs)
    df_clean['fmr_pred'] = poly(df_clean['mu'])
    df_clean['fmr_resid'] = df_clean['metallicity'] - df_clean['fmr_pred']
    
    # 3. Correlate with Sigma
    df_clean['log_sigma'] = np.log10(df_clean['sigma'])
    
    r_warp, p_warp = stats.pearsonr(df_clean['log_sigma'], df_clean['fmr_resid'])
    print(f"N = {len(df_clean)}")
    print(f"Correlation r(Residual, sigma): {r_warp:.4f} (p={p_warp:.2e})")
    
    # 4. Binning
    df_clean['sigma_bin'] = pd.qcut(df_clean['log_sigma'], 8)
    binned = df_clean.groupby('sigma_bin')['fmr_resid'].mean()
    print("\nMean FMR Residual by Sigma Bin:")
    print(binned)
    
    return {
        'best_alpha': float(best_alpha),
        'r_warp': float(r_warp),
        'p_warp': float(p_warp),
        'binned_means': binned.tolist(),
        'bin_centers': [mid.mid for mid in binned.index],
        'n_sample': int(len(df_clean))
    }, df_clean

def create_figure(df, results):
    print("Generating figure...")
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # Scatter
    ax.scatter(df['log_sigma'], df['fmr_resid'], alpha=0.2, s=5, c='k', label='Galaxies')
    
    # Binned
    ax.plot(results['bin_centers'], results['binned_means'], 'r-o', lw=2, label='Mean Residual')
    
    ax.set_xlabel(r'$\log(\sigma)$')
    ax.set_ylabel(r'FMR Residual $\Delta Z$')
    ax.set_title(f"Test BI: FMR Warp vs Potential (r={results['r_warp']:.3f})")
    ax.axhline(0, linestyle='--', color='b', alpha=0.5)
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    out_path = os.path.join(FIGURES_DIR, 'sdss_test_bi_fmr_warp.png')
    plt.savefig(out_path, dpi=150)
    print(f"Figure saved to {out_path}")

def main():
    cache_path = os.path.join(DATA_DIR, 'sdss_fmr_warp.csv')
    
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

    results, df_clean = analyze_fmr_warp(df)
    
    out_path = os.path.join(RESULTS_DIR, 'sdss_test_bi_results.json')
    with open(out_path, 'w') as f:
        results_json = results.copy()
        del results_json['bin_centers']
        json.dump(results_json, f, indent=2)
        
    create_figure(df_clean, results)
    
    print("\nSUMMARY TEST BI:")
    print("TEP Prediction: FMR Residual correlates with sigma. r != 0.")
    print(f"Observed r: {results['r_warp']:.4f}")
    
    if abs(results['r_warp']) > 0.1:
        print("RESULT: CONSISTENT (FMR Warped by potential)")
    else:
        print("RESULT: NULL (FMR is universal)")

if __name__ == "__main__":
    main()
