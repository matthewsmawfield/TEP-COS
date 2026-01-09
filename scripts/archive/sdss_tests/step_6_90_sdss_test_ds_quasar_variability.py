#!/usr/bin/env python3
"""
Step 6.90: SDSS Test DS - Quasar Variability Amplitude

Hypothesis:
Quasar variability amplitude (Structure Function A) is related to the driving mechanism stability. 
If timescales are dilated in deep potentials, the observed amplitude on a fixed terrestrial 
timescale (e.g., 10 years) might be suppressed because we are sampling a smaller fraction 
of the intrinsic characteristic timescale.

Prediction:
Variability Amplitude decreases as M_BH increases (stronger dilation).

Data:
- qsoVarStripe: VAR_A (Amplitude), RA, DEC
- mos_sdss_dr16_qso: logBH, ra, dec

Method:
1. Download Stripe 82 Variability catalog (qsoVarStripe).
2. Download DR16Q catalog (BH masses).
3. Cross-match on RA/DEC (1 arcsec).
4. Correlate VAR_A with logBH.
"""

import pandas as pd
import numpy as np
from scipy import stats
import os
import json
import matplotlib.pyplot as plt
import requests
import time
from astropy.coordinates import SkyCoord
from astropy import units as u

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

def download_data(limit=100):
    print(f"Downloading data for Test DS (Limit: {limit})...")
    
    # 1. Variability
    print("  Querying qsoVarStripe...")
    sql_var = f"""
    SELECT TOP {limit}
        VAR_OBJID, RA, DEC, VAR_A, VAR_GAMMA, VAR_CHI2
    FROM qsoVarStripe
    WHERE VAR_A > 0 AND VAR_CHI2 > 10
    """
    df_var = query_sdss(sql_var)
    if df_var is None: return None
    
    # 2. BH Masses (DR16Q) in Stripe 82
    # Stripe 82: RA -50 to 60, DEC -1.25 to 1.25
    # Simplest way: just get a bunch of DR16Q and cross-match locally.
    # Or try to join on server if possible? No, DR16Q tables often fail on join.
    # Let's get DR16Q in S82 region.
    
    print("  Querying DR16Q (Stripe 82)...")
    sql_qso = f"""
    SELECT TOP {limit}
        specObjID, ra, dec, logBH, logLbol, z
    FROM mos_sdss_dr16_qso
    WHERE (ra > 310 OR ra < 60)
      AND abs(dec) < 1.3
      AND logBH > 6
    """
    df_qso = query_sdss(sql_qso)
    
    if df_qso is None: 
        print("  Failed to download DR16Q.")
        return None
        
    # Cross-match
    print("  Cross-matching...")
    c_var = SkyCoord(ra=df_var['RA'].values*u.degree, dec=df_var['DEC'].values*u.degree)
    c_qso = SkyCoord(ra=df_qso['ra'].values*u.degree, dec=df_qso['dec'].values*u.degree)
    
    idx, d2d, d3d = c_var.match_to_catalog_sky(c_qso)
    
    # Match within 1 arcsec
    mask = d2d < 1.0 * u.arcsec
    
    df_matched = df_var[mask].copy()
    df_matched['logBH'] = df_qso.iloc[idx[mask]]['logBH'].values
    df_matched['z'] = df_qso.iloc[idx[mask]]['z'].values
    
    print(f"  Matched: {len(df_matched)} sources.")
    return df_matched

def analyze_variability(df):
    print("Analyzing Quasar Variability...")
    
    if df is None or len(df) < 50:
        print("  Insufficient data.")
        return None
        
    df = df.dropna().copy()
    
    # Correlation: VAR_A vs logBH
    # Prediction: Negative slope (High Mass -> Lower Amplitude / Slower timescales)
    
    slope, intercept, r_val, p_val, std_err = stats.linregress(df['logBH'], df['VAR_A'])
    
    print(f"  Correlation (logBH vs Amplitude): r={r_val:.3f}, p={p_val:.2e}, slope={slope:.4f}")
    
    # Plot
    plt.figure(figsize=(8, 6))
    sc = plt.scatter(df['logBH'], df['VAR_A'], c=df['z'], cmap='viridis', s=10, alpha=0.6)
    plt.colorbar(sc, label='Redshift')
    
    x_range = np.linspace(df['logBH'].min(), df['logBH'].max(), 100)
    plt.plot(x_range, intercept + slope*x_range, 'k--', label=f'Slope={slope:.3f}')
    
    plt.xlabel('log(Black Hole Mass)')
    plt.ylabel('Variability Amplitude (A)')
    plt.title('Test DS: Variability vs Potential')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # Limit y-axis if outliers
    plt.ylim(0, np.percentile(df['VAR_A'], 99))
    
    out_path = os.path.join(FIGURES_DIR, 'sdss_test_ds_variability.png')
    plt.savefig(out_path, dpi=150)
    print(f"Figure saved to {out_path}")
    
    return {
        'correlation_r': r_val,
        'p_value': p_val,
        'slope': slope,
        'n_qso': int(len(df))
    }

def main():
    cache_path = os.path.join(DATA_DIR, 'sdss_qso_var.csv')
    
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

    results = analyze_variability(df)
    
    out_path = os.path.join(RESULTS_DIR, 'sdss_test_ds_results.json')
    if results:
        with open(out_path, 'w') as f:
            json.dump(results, f, indent=2)
            
        print("\nSUMMARY TEST DS:")
        print(f"Slope (logBH vs Amplitude): {results['slope']:.4f}")
        
        # Prediction: Negative Slope
        if results['p_value'] < 0.05 and results['slope'] < 0:
             print("RESULT: SIGNAL (Variability suppressed in Deep Potentials)")
        elif results['p_value'] < 0.05 and results['slope'] > 0:
             print("RESULT: CONTRADICTED (Variability enhanced in Deep Potentials)")
        else:
             print("RESULT: NULL (No significant dependence)")

if __name__ == "__main__":
    main()
