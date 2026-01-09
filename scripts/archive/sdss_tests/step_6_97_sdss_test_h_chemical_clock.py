#!/usr/bin/env python3
"""
Step 6.97: Test H - Chemical Clock Discrepancy

Hypothesis:
[Mg/Fe] tracks the star formation timescale (alpha enrichment vs Fe delay).
Spectroscopic ages (D4000, Hb) track proper time.
Under TEP, high-sigma galaxies experience time dilation, meaning less proper time 
elapses for a given coordinate time interval. Type Ia SNe (coordinate time delay) 
should appear "delayed" relative to spectroscopic age (proper time).
High-sigma galaxies should show elevated [Mg/Fe] at fixed spectroscopic age.

Prediction:
At fixed Spectroscopic Age: [Mg/Fe] correlates positively with Sigma.

Data:
- galSpecIndx: Mgb, Fe5270, Fe5335, D4000, Hbeta
- galSpecInfo: Sigma
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

def query_sdss(limit=50000):
    print(f"Querying SDSS for Test H (Limit: {limit})...")
    
    # Use galSpec tables
    sql = f"""
    SELECT TOP {limit}
        g.specObjID,
        g.z,
        g.v_disp,
        i.lick_mgb,
        i.lick_fe5270,
        i.lick_fe5335,
        i.d4000_n,
        i.lick_hb,
        s.logMass
    FROM galSpecInfo g
    JOIN galSpecIndx i ON g.specObjID = i.specObjID
    JOIN stellarMassFSPSGranWideDust s ON g.specObjID = s.specObjID
    WHERE g.z BETWEEN 0.02 AND 0.20
      AND g.v_disp > 80 AND g.v_disp < 400
      AND g.reliable = 1
      AND i.lick_mgb > 0 AND i.lick_fe5270 > 0 AND i.lick_fe5335 > 0
      AND i.d4000_n > 1.2
      AND i.lick_mgb_err < 0.5
    """
    
    for attempt in range(3):
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
        time.sleep(2)
    return None

def analyze_chemical_clock(df):
    print("Analyzing Chemical Clock Discrepancy...")
    
    if df is None or len(df) < 100:
        print("  Insufficient data.")
        return None
    
    # Calculate Proxies
    # <Fe> = (Fe5270 + Fe5335) / 2
    df['Fe_avg'] = (df['lick_fe5270'] + df['lick_fe5335']) / 2
    
    # [Mg/Fe] proxy ~ log(Mgb / Fe_avg)
    # Calibrated relations exist, but raw ratio is sufficient for differential test
    df['MgFe_proxy'] = np.log10(df['lick_mgb'] / df['Fe_avg'])
    
    # Spectroscopic Age Proxy
    # D4000 is good, Hbeta adds leverage. 
    # Let's use D4000 primarily as it's robust.
    df['Age_proxy'] = df['d4000_n']
    
    # Clean
    df = df.replace([np.inf, -np.inf], np.nan).dropna()
    
    # Remove outliers
    df = df[(np.abs(stats.zscore(df['MgFe_proxy'])) < 3) & 
            (np.abs(stats.zscore(df['Age_proxy'])) < 3)].copy()
    
    # 1. Regress [Mg/Fe] on Age
    # Expected: [Mg/Fe] decreases with Age? 
    # Actually, [Mg/Fe] tracks formation timescale. Short timescale -> High [Mg/Fe].
    # Older galaxies (formed early) usually formed quickly -> High [Mg/Fe].
    # So [Mg/Fe] should correlate positively with Age (D4000).
    
    slope_age, intercept_age, r_age, _, _ = stats.linregress(df['Age_proxy'], df['MgFe_proxy'])
    print(f"  Mg/Fe vs Age Correlation: r={r_age:.3f}")
    
    # Calculate Residuals: Delta [Mg/Fe] at fixed Age
    df['MgFe_resid'] = df['MgFe_proxy'] - (intercept_age + slope_age * df['Age_proxy'])
    
    # 2. Correlate Residuals with Sigma
    df['log_sigma'] = np.log10(df['v_disp'])
    
    slope, intercept, r_val, p_val, std_err = stats.linregress(df['log_sigma'], df['MgFe_resid'])
    
    print(f"  Correlation (Delta[Mg/Fe] vs log_sigma): r={r_val:.3f}, p={p_val:.3e}, slope={slope:.4f}")
    
    # Plot
    plt.figure(figsize=(10, 5))
    
    # Left: Age vs Mg/Fe
    plt.subplot(1, 2, 1)
    plt.scatter(df['Age_proxy'], df['MgFe_proxy'], s=1, alpha=0.3, c='gray')
    x_range = np.linspace(df['Age_proxy'].min(), df['Age_proxy'].max(), 100)
    plt.plot(x_range, intercept_age + slope_age * x_range, 'r--', label='Mean Relation')
    plt.xlabel('Spectroscopic Age Proxy (D4000)')
    plt.ylabel('[Mg/Fe] Proxy')
    plt.title('Chemical Clock Calibration')
    plt.legend()
    
    # Right: Residuals vs Sigma
    plt.subplot(1, 2, 2)
    plt.scatter(df['log_sigma'], df['MgFe_resid'], s=1, alpha=0.3, c=df['logMass'], cmap='plasma')
    plt.colorbar(label='log Stellar Mass')
    
    x_range = np.linspace(df['log_sigma'].min(), df['log_sigma'].max(), 100)
    plt.plot(x_range, intercept + slope*x_range, 'k--', label=f'Slope={slope:.3f}')
    
    plt.xlabel('log(Velocity Dispersion)')
    plt.ylabel('$\Delta$[Mg/Fe] (at fixed Age)')
    plt.title('Test H: Clock Discrepancy')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    plt.tight_layout()
    out_path = os.path.join(FIGURES_DIR, 'sdss_test_h_chemical_clock.png')
    plt.savefig(out_path, dpi=150)
    print(f"Figure saved to {out_path}")
    
    return {
        'slope': slope,
        'p_val': p_val,
        'r_val': r_val,
        'n_gal': int(len(df))
    }

def main():
    cache_path = os.path.join(DATA_DIR, 'sdss_chemical_clock.csv')
    
    if os.path.exists(cache_path):
        print("Loading cached data...")
        df = pd.read_csv(cache_path)
    else:
        df = query_sdss()
        if df is not None:
            df.to_csv(cache_path, index=False)
        else:
            return

    results = analyze_chemical_clock(df)
    
    out_path = os.path.join(RESULTS_DIR, 'sdss_test_h_results.json')
    if results:
        with open(out_path, 'w') as f:
            json.dump(results, f, indent=2)
            
        print("\nSUMMARY TEST H:")
        print(f"Slope: {results['slope']:.4f}")
        
        if results['p_val'] < 0.05:
            if results['slope'] > 0:
                 print("RESULT: SIGNAL (High-sigma galaxies have enhanced Alpha at fixed Age)")
            else:
                 print("RESULT: CONTRADICTED (High-sigma galaxies have depleted Alpha)")
        else:
             print("RESULT: NULL (No discrepancy)")

if __name__ == "__main__":
    main()
