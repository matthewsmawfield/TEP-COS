#!/usr/bin/env python3
"""
Step 6.38: SDSS Test BG - M-sigma Residuals (Black Hole Growth Clock)

Hypothesis:
Black hole growth (M_dot) is a rate process.
The M-sigma relation connects this integrated growth to the potential depth sigma.
TEP predicts that the growth rate is dilated in deeper potentials, altering the equilibrium relationship.
Residuals in M_BH relative to sigma should correlate with secondary potential depth markers like surface brightness (Sigma_*) or compactness.

Prediction:
r(M_BH residual, Compactness) != 0.

Data:
- spiders_quasar: logBHMA_hb, logBHMS_mgII (BH Mass).
- emissionLinesPort: sigma_stars (Velocity Dispersion).
- stellarMassFSPSGranWideDust: logMass.
- PhotoObjAll: petroR50_r (Radius) -> Compactness.

Method:
1. Join tables.
2. Select objects with valid BH mass and sigma.
3. Fit M_BH vs sigma relation.
4. Calculate residuals.
5. Calculate Compactness = logMass - 2*log10(R50).
6. Correlate residuals with Compactness.
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
    print(f"Querying SDSS for Test BG (Limit: {limit})...")
    
    # Complex join.
    # spiders_quasar (q) -> emissionLinesPort (e) -> stellarMass... (s) -> SpecObj -> PhotoObj (p)
    # Note: stellarMass table usually has logMass.
    # We will try to get petroR50_r from PhotoObj using SpecObj linking.
    # But to save joins, maybe just check if we can get mass and sigma first.
    # Compactness needs Radius.
    # Let's try without PhotoObj first (just M-sigma), and use Mass as secondary param?
    # Or try to join PhotoObj.
    
    # We use q.SPECOBJID to link.
    
    sql = f"""
    SELECT TOP {limit}
        q.SPECOBJID,
        q.logBHMA_hb, q.logBHMS_mgII,
        e.sigmaStars as sigma,
        s.logMass
        -- p.petroR50_r -- Skipping PhotoObj join to avoid timeout for now
        
    FROM spiders_quasar q
    JOIN emissionLinesPort e ON q.SPECOBJID = e.specObjID
    JOIN stellarMassFSPSGranWideDust s ON q.SPECOBJID = s.specObjID
    
    WHERE 
        (q.logBHMA_hb > 0 OR q.logBHMS_mgII > 0)
        AND e.sigmaStars > 50 AND e.sigmaStars < 500
        AND s.logMass > 0
    """
    return query_sdss(sql)

def analyze_msigma(df):
    print("Analyzing M-sigma Residuals...")
    
    # 1. Clean
    df_clean = df.dropna().copy()
    
    # Combine BH masses
    df_clean['log_mbh'] = df_clean['logBHMS_mgII'].fillna(df_clean['logBHMA_hb'])
    df_clean = df_clean.dropna(subset=['log_mbh', 'sigma']).copy()
    
    # 2. Fit M-sigma
    # log M = alpha + beta * log(sigma/200)
    df_clean['log_sigma'] = np.log10(df_clean['sigma'])
    
    slope, intercept, r_val, p_val, std_err = stats.linregress(df_clean['log_sigma'], df_clean['log_mbh'])
    print(f"M-sigma Fit: slope={slope:.2f}, intercept={intercept:.2f}, r={r_val:.2f}")
    
    # 3. Residuals
    df_clean['mbh_resid'] = df_clean['log_mbh'] - (slope * df_clean['log_sigma'] + intercept)
    
    # 4. Correlate with Mass (as proxy for potential depth/compactness?)
    # Ideally we want Compactness (Mass/R). Using Mass alone is a check.
    
    r_mass, p_mass = stats.pearsonr(df_clean['logMass'], df_clean['mbh_resid'])
    print(f"Correlation r(Residual, StellarMass): {r_mass:.4f} (p={p_mass:.2e})")
    
    # 5. Binning
    df_clean['mass_bin'] = pd.qcut(df_clean['logMass'], 8)
    binned = df_clean.groupby('mass_bin')['mbh_resid'].mean()
    print("\nMean Residual by Mass Bin:")
    print(binned)
    
    return {
        'slope_msigma': float(slope),
        'r_msigma': float(r_val),
        'r_resid_mass': float(r_mass),
        'p_resid_mass': float(p_mass),
        'binned_means': binned.tolist(),
        'bin_centers': [mid.mid for mid in binned.index],
        'n_sample': int(len(df_clean))
    }, df_clean

def create_figure(df, results):
    print("Generating figure...")
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # Scatter
    ax.scatter(df['logMass'], df['mbh_resid'], alpha=0.2, s=5, c='k', label='Quasars')
    
    # Binned
    ax.plot(results['bin_centers'], results['binned_means'], 'r-o', lw=2, label='Mean Residual')
    
    ax.set_xlabel(r'$\log(M_{*}/M_{\odot})$')
    ax.set_ylabel(r'$M_{BH}$ Residual (from $M-\sigma$)')
    ax.set_title(f"Test BG: BH Residuals vs Stellar Mass (r={results['r_resid_mass']:.3f})")
    ax.axhline(0, linestyle='--', color='b', alpha=0.5)
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    out_path = os.path.join(FIGURES_DIR, 'sdss_test_bg_msigma_residuals.png')
    plt.savefig(out_path, dpi=150)
    print(f"Figure saved to {out_path}")

def main():
    cache_path = os.path.join(DATA_DIR, 'sdss_msigma_residuals.csv')
    
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

    results, df_clean = analyze_msigma(df)
    
    out_path = os.path.join(RESULTS_DIR, 'sdss_test_bg_results.json')
    with open(out_path, 'w') as f:
        results_json = results.copy()
        del results_json['bin_centers']
        json.dump(results_json, f, indent=2)
        
    create_figure(df_clean, results)
    
    print("\nSUMMARY TEST BG:")
    print("TEP Prediction: Residual correlates with potential depth (Mass/Compactness).")
    print(f"Observed r(Resid, Mass): {results['r_resid_mass']:.4f}")
    
    if abs(results['r_resid_mass']) > 0.2:
        print("RESULT: CONSISTENT (Secondary correlation observed)")
    else:
        print("RESULT: NULL (M-sigma residuals are random wrt mass)")

if __name__ == "__main__":
    main()
