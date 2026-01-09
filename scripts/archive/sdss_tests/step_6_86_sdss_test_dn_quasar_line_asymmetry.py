#!/usr/bin/env python3
"""
Step 6.86: SDSS Test DN - Quasar Line Asymmetry (Inflow/Outflow Clock)

Hypothesis:
Broad emission lines often show asymmetry due to outflows or gravitational redshift. 
The kinematic structure of the BLR is governed by the potential. If TEP modifies 
the potential or the rates of inflow/outflow, the skewness or asymmetry of lines 
like H-beta or CIV should show a dependence on Black Hole Mass (potential depth) 
distinct from standard models.

Prediction:
Line Asymmetry (Skewness) scales with Black Hole Mass.

Data:
- mos_sdss_dr16_qso: logBH, z, FWHM_Hb, FWHM_MgII, FWHM_CIV?
  Note: DR16Q often has FWHM columns.
- spAll: detailed line fits (may be harder to join).

Method:
1. Use DR16Q catalog.
2. Define Asymmetry Proxy: FWHM difference between lines?
   - Or redshift difference: z_MgII - z_CIV (Blueshifts).
   - High blueshifts -> Outflows.
   - Does Outflow velocity depend on M_BH?
3. Calculate Delta_v = c * (z_MgII - z_CIV) / (1+z).
4. Correlate Delta_v with logBH.
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

def download_data(limit=100):
    print(f"Querying SDSS for Test DN (Limit: {limit})...")
    
    # Use redshift differences as asymmetry/outflow proxy
    # z_CIV is often blueshifted relative to z_MgII (systemic)
    
    sql = f"""
    SELECT TOP {limit}
        specObjID,
        logBH,
        logLbol,
        Z_CIV,
        Z_MGII,
        FWHM_CIV,
        FWHM_MGII
    FROM mos_sdss_dr16_qso
    WHERE logBH > 6
      AND Z_CIV > 0 AND Z_MGII > 0
      AND FWHM_CIV > 0 AND FWHM_MGII > 0
    """
    return query_sdss(sql)

def analyze_line_asymmetry(df):
    print("Analyzing Quasar Line Asymmetry...")
    
    if df is None or len(df) < 100:
        print("  Insufficient data.")
        return None
        
    df = df.dropna().copy()
    
    # Calculate Velocity Offset (Blueshift)
    # v = c * (z_sys - z_line) / (1+z_sys)
    # Use MgII as systemic proxy
    c = 300000.0 # km/s
    df['v_off_civ'] = c * (df['Z_MGII'] - df['Z_CIV']) / (1 + df['Z_MGII'])
    
    # Positive v_off -> CIV is blueshifted (outflow)
    
    # Calculate FWHM Ratio
    df['fwhm_ratio'] = df['FWHM_CIV'] / df['FWHM_MGII']
    
    # Correlate with Black Hole Mass
    # Prediction: In deep potentials (high M_BH), escape velocity is higher.
    # Outflows might be faster? Or suppressed?
    # TEP: Potential is deeper -> Slower outflows? Or time dilation effect?
    # If time runs slower, observed velocity v_obs = dl/dt_obs.
    # dt_obs > dt_local. So v_obs < v_local?
    # Or gravitational redshift term.
    
    # Let's check correlation
    slope, intercept, r_val, p_val, std_err = stats.linregress(df['logBH'], df['v_off_civ'])
    
    print(f"  Correlation (logBH vs CIV Offset): r={r_val:.3f}, p={p_val:.2e}, slope={slope:.4f}")
    
    # Plot
    fig, ax = plt.subplots(1, 2, figsize=(12, 5))
    
    # Offset vs BH Mass
    ax[0].scatter(df['logBH'], df['v_off_civ'], s=2, alpha=0.3, c='blue')
    
    # Fit line
    x_range = np.linspace(df['logBH'].min(), df['logBH'].max(), 100)
    ax[0].plot(x_range, intercept + slope*x_range, 'k--', label=f'Slope={slope:.2f}')
    
    ax[0].set_xlabel('log(Black Hole Mass)')
    ax[0].set_ylabel('CIV Velocity Offset (km/s)')
    ax[0].set_title('Test DN: Outflow Velocity vs Potential')
    ax[0].legend()
    
    # Width Ratio vs BH Mass
    slope2, int2, r2, p2, err2 = stats.linregress(df['logBH'], df['fwhm_ratio'])
    
    ax[1].scatter(df['logBH'], df['fwhm_ratio'], s=2, alpha=0.3, c='green')
    ax[1].plot(x_range, int2 + slope2*x_range, 'k--', label=f'Slope={slope2:.3f}')
    ax[1].set_xlabel('log(Black Hole Mass)')
    ax[1].set_ylabel('FWHM Ratio (CIV / MgII)')
    ax[1].set_title('Line Width Ratio')
    ax[1].legend()
    
    out_path = os.path.join(FIGURES_DIR, 'sdss_test_dn_line_asymmetry.png')
    plt.savefig(out_path, dpi=150)
    print(f"Figure saved to {out_path}")
    
    return {
        'slope_offset': slope,
        'p_value': p_val,
        'slope_width': slope2,
        'n_qso': int(len(df))
    }

def main():
    cache_path = os.path.join(DATA_DIR, 'sdss_qso_lines.csv')
    
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

    results = analyze_line_asymmetry(df)
    
    out_path = os.path.join(RESULTS_DIR, 'sdss_test_dn_results.json')
    if results:
        with open(out_path, 'w') as f:
            json.dump(results, f, indent=2)
            
        print("\nSUMMARY TEST DN:")
        print(f"Slope (logBH vs CIV Offset): {results['slope_offset']:.4f}")
        
        if results['p_value'] < 0.05 and abs(results['slope_offset']) > 50:
             print("RESULT: SIGNAL (Significant correlation)")
        else:
             print("RESULT: NULL (Weak or no correlation)")

if __name__ == "__main__":
    main()
