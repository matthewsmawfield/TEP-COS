#!/usr/bin/env python3
"""
Step 6.64: SDSS Test CK - Mass-Metallicity Relation (MZR) Residuals

Hypothesis:
The Mass-Metallicity Relation arises from the balance between gas inflow, star formation, and outflows. 
These are all rate-dependent processes. In TEP, these rates scale with the potential depth. 
If the scaling of inflows (gravitational) differs from outflows (feedback) due to time dilation, 
the equilibrium metallicity Z_eq should show residuals correlated with sigma at fixed Mass.

Prediction:
MZR Residual (Delta Z) correlates with velocity dispersion sigma.
Standard model also predicts correlation (Fundamental Metallicity Relation), 
typically Z depends on Mass and SFR. Sigma is a proxy for Mass/Density.
We need to see if Sigma adds information beyond Mass.

Data:
- galSpecExtra: oh_p50 (Oxygen Abundance / Metallicity)
- stellarMassFSPSGranWideDust: logMass
- galSpecInfo: v_disp (Sigma)

Method:
1. Fetch Mass, Metallicity, Sigma.
2. Fit MZR (Z vs Mass) using a polynomial.
3. Calculate Residuals: Delta Z = Z_obs - Z_model(Mass).
4. Correlate Delta Z with Sigma.
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
    print(f"Querying SDSS for Test CK (Limit: {limit})...")
    
    # Use galSpecInfo for v_disp (more robust than emissionLinesPort usually)
    
    sql = f"""
    SELECT TOP {limit}
        g.specObjID,
        s.logMass,
        g.oh_p50 as metallicity,
        i.v_disp as sigma,
        i.z
        
    FROM galSpecExtra g
    JOIN stellarMassFSPSGranWideDust s ON g.specObjID = s.specObjID
    JOIN galSpecInfo i ON g.specObjID = i.specObjID
    
    WHERE 
        g.oh_p50 > -9
        AND s.logMass > 8
        AND i.v_disp > 50 AND i.v_disp < 400
        AND i.z > 0.02 AND i.z < 0.2
        AND i.reliable = 1
    """
    return query_sdss(sql)

def analyze_mzr(df):
    print("Analyzing MZR Residuals...")
    
    # Clean
    df = df.dropna().copy()
    df = df[df['metallicity'] > -5] # Valid range
    
    print(f"  Sample size: {len(df)}")
    
    # 1. Fit MZR (Z vs Mass)
    # Usually a 2nd or 3rd order polynomial fits well in log-log
    # Z ~ 12 + log(O/H) usually. oh_p50 is usually 12+log(O/H) or relative?
    # Check range. 
    # Usually ranges 8.5 to 9.2.
    
    # Polynomial fit
    z_fit = np.polyfit(df['logMass'], df['metallicity'], 2)
    p = np.poly1d(z_fit)
    
    df['mzr_model'] = p(df['logMass'])
    df['mzr_resid'] = df['metallicity'] - df['mzr_model']
    
    print(f"  MZR Fit: {z_fit}")
    
    # 2. Correlate Residuals with Sigma
    df['log_sigma'] = np.log10(df['sigma'])
    
    r_val, p_val = stats.pearsonr(df['log_sigma'], df['mzr_resid'])
    print(f"  Correlation r(Delta Z, logSigma): {r_val:.4f} (p={p_val:.2e})")
    
    slope, intercept, _, _, _ = stats.linregress(df['log_sigma'], df['mzr_resid'])
    print(f"  Slope (Residual vs Sigma): {slope:.4f}")
    
    # Partial correlation controlling for Mass (redundant since we took residuals, but safe)
    # Technically we already removed mass dependence.
    
    # Plot
    fig, ax = plt.subplots(1, 2, figsize=(12, 5))
    
    # MZR
    ax[0].scatter(df['logMass'], df['metallicity'], alpha=0.1, s=2, c='k')
    x_range = np.linspace(df['logMass'].min(), df['logMass'].max(), 100)
    ax[0].plot(x_range, p(x_range), 'r-', label='MZR Fit')
    ax[0].set_xlabel('log Stellar Mass')
    ax[0].set_ylabel('Metallicity (12+log O/H)')
    ax[0].set_title('Mass-Metallicity Relation')
    ax[0].legend()
    ax[0].grid(True, alpha=0.3)
    
    # Residuals vs Sigma
    # Binning for clarity
    bins = np.linspace(df['log_sigma'].min(), df['log_sigma'].max(), 15)
    df['sig_bin'] = pd.cut(df['log_sigma'], bins=bins)
    binned = df.groupby('sig_bin')['mzr_resid'].agg(['mean', 'sem', 'count'])
    binned['sig_center'] = [i.mid for i in binned.index]
    
    ax[1].scatter(df['log_sigma'], df['mzr_resid'], alpha=0.05, s=2, c='gray')
    ax[1].errorbar(binned['sig_center'], binned['mean'], yerr=binned['sem'], fmt='o-', capsize=5, color='blue', label='Mean Residual')
    
    # Fit line
    ax[1].plot(x_range, np.zeros_like(x_range), 'k--', alpha=0.5) # Zero line? No x_range is mass.
    # Sigma range fit
    sig_range = np.linspace(df['log_sigma'].min(), df['log_sigma'].max(), 100)
    ax[1].plot(sig_range, slope * sig_range + intercept, 'r-', label=f'Slope={slope:.2f}')
    
    ax[1].set_xlabel('log Velocity Dispersion')
    ax[1].set_ylabel('MZR Residual (Delta Z)')
    ax[1].set_title(f'Test CK: MZR Residuals (r={r_val:.2f})')
    ax[1].legend()
    ax[1].grid(True, alpha=0.3)
    
    out_path = os.path.join(FIGURES_DIR, 'sdss_test_ck_mzr.png')
    plt.savefig(out_path, dpi=150)
    print(f"Figure saved to {out_path}")
    
    return {
        'r_val': r_val,
        'slope': slope,
        'n_sample': int(len(df))
    }

def main():
    cache_path = os.path.join(DATA_DIR, 'sdss_mzr.csv')
    
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

    results = analyze_mzr(df)
    
    out_path = os.path.join(RESULTS_DIR, 'sdss_test_ck_results.json')
    if results:
        def default(o):
            if isinstance(o, pd.Interval): return str(o)
            raise TypeError
            
        with open(out_path, 'w') as f:
            json.dump(results, f, indent=2, default=default)
            
        print("\nSUMMARY TEST CK:")
        print("Prediction: MZR Residuals correlate with Sigma (Fundamental Metallicity Relation).")
        print(f"Observed r: {results['r_val']:.4f}")
        
        if abs(results['r_val']) > 0.1:
             print("RESULT: SIGNAL (Correlation observed - FMR/TEP Consistent)")
        else:
             print("RESULT: NULL (No secondary dependence on sigma)")

if __name__ == "__main__":
    main()
