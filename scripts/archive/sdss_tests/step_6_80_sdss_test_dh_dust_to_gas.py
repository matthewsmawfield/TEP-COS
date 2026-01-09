#!/usr/bin/env python3
"""
Step 6.80: SDSS Test DH - Dust-to-Gas Ratio (Survival Clock)

Hypothesis:
The Dust-to-Gas ratio (DGR) is a balance between dust formation (AGB/SN) and destruction 
(shocks/sputtering). Destruction timescales are often shorter than formation. 
If time dilation slows down destruction rates in high-sigma galaxies, the equilibrium 
DGR should be **higher**.

Prediction:
Dust-to-Gas Ratio increases with Velocity Dispersion.

Data:
- mangaHIall: logmhi (HI Mass)
- mangaDAPall: stellar_sigma_1re, emline_gflux_1re_ha_6564, emline_gflux_1re_hb_4862
- mangaTarget: nsa_elpetro_mass

Method:
1. Join MaNGA HI and DAP catalogs.
2. Calculate Balmer Decrement (Ha/Hb) as proxy for Dust Extinction (E(B-V)).
   - Dust Mass ~ E(B-V) * Area (roughly).
   - Or just use Extinction/Gas ratio.
3. Define DGR Proxy: log(Ha/Hb) / M_HI? 
   - Better: Analyze residuals of log(Ha/Hb) vs log(M_HI) against Sigma.
   - Or: log(Ha/Hb) vs Sigma (controlling for Mass and Gas Fraction).
4. Prediction: At fixed Gas Mass, High Sigma -> High Extinction.
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
    print(f"Querying SDSS for Test DH (Limit: {limit})...")
    
    # Note: mangaHIall logmhi is log10(M_sun).
    # Ha/Hb flux ratio.
    
    sql = f"""
    SELECT TOP {limit}
        h.mangaid, 
        h.logmhi,
        d.stellar_sigma_1re as sigma,
        s.nsa_elpetro_mass as mass,
        d.emline_gflux_1re_ha_6564 as flux_ha,
        d.emline_gflux_1re_hb_4862 as flux_hb
    FROM mangaHIall h
    JOIN mangaDAPall d ON h.mangaid = d.mangaid
    JOIN mangaTarget s ON h.mangaid = s.mangaid
    WHERE h.confprob > 0.9
      AND d.stellar_sigma_1re > 0
      AND d.emline_gflux_1re_hb_4862 > 0
      AND d.emline_gflux_1re_ha_6564 > 0
      AND h.logmhi > 0
    """
    return query_sdss(sql)

def analyze_dust_gas(df):
    print("Analyzing Dust-to-Gas Ratio...")
    
    if df is None or len(df) < 50:
        print("  Insufficient data.")
        return None
        
    df = df.dropna().copy()
    
    # Balmer Decrement
    df['balmer_dec'] = df['flux_ha'] / df['flux_hb']
    
    # Filter physical values
    # Theoretical Case B ~ 2.86. Observed can be lower due to errors, but dust implies > 2.86.
    df = df[(df['balmer_dec'] > 2.0) & (df['balmer_dec'] < 15.0)].copy()
    
    df['log_bd'] = np.log10(df['balmer_dec'])
    
    # DGR Proxy
    # We want to test if 'Dust per unit Gas' increases with Sigma.
    # Dust ~ log_bd
    # Gas ~ logmhi
    # Metric: log_bd - alpha * logmhi ?
    # Or simply: does log_bd correlate with Sigma at fixed Gas Mass?
    
    df['log_sigma'] = np.log10(df['sigma'])
    df['log_mass'] = np.log10(df['mass'])
    
    # Control for Gas Mass and Stellar Mass (Metallicity proxy)
    # Dust formation depends on Metallicity (Stellar Mass) and Gas content.
    
    from sklearn.linear_model import LinearRegression
    X = df[['logmhi', 'log_mass']]
    y = df['log_bd']
    
    reg = LinearRegression().fit(X, y)
    print(f"  Control Fit (Gas, Mass -> Dust) R2: {reg.score(X, y):.3f}")
    
    df['dgr_resid'] = y - reg.predict(X)
    
    # Correlate with Sigma
    slope, intercept, r_val, p_val, std_err = stats.linregress(df['log_sigma'], df['dgr_resid'])
    
    print(f"  Correlation (log Sigma vs DGR Residual): r={r_val:.3f}, p={p_val:.2e}, slope={slope:.4f}")
    
    # Plot
    fig, ax = plt.subplots(1, 2, figsize=(12, 5))
    
    # Raw
    sc = ax[0].scatter(df['logmhi'], df['log_bd'], c=df['log_sigma'], cmap='viridis', s=15, alpha=0.7)
    plt.colorbar(sc, ax=ax[0], label='log(Sigma)')
    ax[0].set_xlabel('log(HI Mass)')
    ax[0].set_ylabel('log(Balmer Decrement)')
    ax[0].set_title('Dust vs Gas (Color=Sigma)')
    
    # Residuals vs Sigma
    ax[1].scatter(df['log_sigma'], df['dgr_resid'], alpha=0.6, s=15, c='brown')
    
    x_range = np.linspace(df['log_sigma'].min(), df['log_sigma'].max(), 100)
    ax[1].plot(x_range, intercept + slope*x_range, 'k--', label=f'Slope={slope:.3f}')
    
    ax[1].set_xlabel('log(Velocity Dispersion)')
    ax[1].set_ylabel('DGR Residual (Dust excess)')
    ax[1].set_title('Test DH: Dust Survival vs Potential')
    ax[1].legend()
    ax[1].grid(True, alpha=0.3)
    
    out_path = os.path.join(FIGURES_DIR, 'sdss_test_dh_dust_gas.png')
    plt.savefig(out_path, dpi=150)
    print(f"Figure saved to {out_path}")
    
    return {
        'correlation_r': r_val,
        'p_value': p_val,
        'slope': slope,
        'n_gal': int(len(df))
    }

def main():
    cache_path = os.path.join(DATA_DIR, 'sdss_dust_gas.csv')
    
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

    results = analyze_dust_gas(df)
    
    out_path = os.path.join(RESULTS_DIR, 'sdss_test_dh_results.json')
    if results:
        with open(out_path, 'w') as f:
            json.dump(results, f, indent=2)
            
        print("\nSUMMARY TEST DH:")
        print(f"Slope (Sigma vs DGR): {results['slope']:.4f}")
        
        if results['p_value'] < 0.05 and results['slope'] > 0.05:
             print("RESULT: SIGNAL (Higher Dust-to-Gas in Deep Potentials)")
        elif results['p_value'] < 0.05 and results['slope'] < -0.05:
             print("RESULT: CONTRADICTED (Lower Dust-to-Gas in Deep Potentials)")
        else:
             print("RESULT: NULL (No significant dependence)")

if __name__ == "__main__":
    main()
