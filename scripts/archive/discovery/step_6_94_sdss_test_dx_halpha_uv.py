#!/usr/bin/env python3
"""
Step 6.94: SDSS Test DX - H-alpha vs UV Star Formation (Timescale Ratio)

Hypothesis:
H-alpha traces ionizing flux from O-stars (<10 Myr timescale). UV continuum traces 
B-stars (~100 Myr timescale). The ratio SFR(Ha)/SFR(UV) is sensitive to the 
"burstiness" of star formation. If time dilation stretches the duration of bursts 
or the lifetimes of massive stars in high-sigma galaxies, this ratio should show 
a systematic trend with sigma.

Prediction:
Flux Ratio F(Ha) / F(UV) correlates with Velocity Dispersion.

Data:
- galSpecLine: h_alpha_flux, h_beta_flux
- galSpecInfo: v_disp, z
- galSpecExtra: lgm_tot_p50 (Stellar Mass)
- SpecPhotoAll: modelMag_u (UV proxy)

Method:
1. Join galSpec tables and SpecPhotoAll.
2. Select Star Forming galaxies (Ha > 0).
3. Calculate F(UV) from u-band mag.
4. Correct for dust using Balmer decrement (Ha/Hb).
5. Define Ratio R = L(Ha) / L(u).
6. Correlate R with Sigma.
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
    print(f"Querying SDSS for Test DX (Using galSpec tables, Limit: {limit})...")
    
    # Using galSpec tables which might be more stable than emissionLinesPort
    
    sql = f"""
    SELECT TOP {limit}
        gl.specObjID, 
        gi.v_disp as sigma_stars,
        gl.h_alpha_flux, 
        gl.h_beta_flux,
        sp.modelMag_u, 
        gi.z,
        ge.lgm_tot_p50 as logMass
    FROM galSpecLine gl
    JOIN galSpecInfo gi ON gl.specObjID = gi.specObjID
    JOIN galSpecExtra ge ON gl.specObjID = ge.specObjID
    JOIN SpecPhotoAll sp ON gl.specObjID = sp.specObjID
    WHERE gi.v_disp > 0
      AND gl.h_alpha_flux > 10
      AND gl.h_beta_flux > 5
      AND sp.class = 'GALAXY'
      AND gi.z BETWEEN 0.02 AND 0.1
    """
    return query_sdss(sql)

def analyze_timescale_ratio(df):
    print("Analyzing H-alpha vs UV Star Formation...")
    
    if df is None or len(df) < 50:
        print("  Insufficient data.")
        return None
        
    df = df.dropna().copy()
    
    # 1. Dust Correction (Balmer Decrement)
    # Intrinsic Ha/Hb = 2.86 (Case B recombination, T=10000K)
    
    df['balmer_ratio'] = df['h_alpha_flux'] / df['h_beta_flux']
    df = df[df['balmer_ratio'] > 2.86].copy() # Physical decrement
    
    # E(B-V) approx
    # k_Ha = 2.53, k_Hb = 3.61. Delta k = 1.08.
    # E(B-V) = log10(R_obs/2.86) / 0.432
    
    df['ebv'] = np.log10(df['balmer_ratio'] / 2.86) / 0.432
    df['ebv'] = df['ebv'].clip(lower=0)
    
    # Correct Ha Flux
    # A_Ha = 2.53 * E(B-V)
    df['A_Ha'] = 2.53 * df['ebv']
    df['Ha_corr'] = df['h_alpha_flux'] * (10**(0.4 * df['A_Ha']))
    
    # Correct u-band Flux (UV Proxy)
    # A_u = 4.9 * E(B-V)
    df['A_u'] = 4.9 * df['ebv']
    df['u_corr'] = df['modelMag_u'] - df['A_u']
    
    # Convert u_corr to Flux density proxy
    # F_u ~ 10^(-0.4 * u_corr)
    df['F_u_corr'] = 10**(-0.4 * df['u_corr'])
    
    # Ratio: Ha / UV
    df['ratio'] = df['Ha_corr'] / df['F_u_corr']
    
    # Handle NaNs and Infs
    df = df.replace([np.inf, -np.inf], np.nan).dropna()
    print(f"  Records after cleaning NaNs/Infs: {len(df)}")
    
    if len(df) < 50:
        print("  Insufficient data after cleaning.")
        return None

    # Filter outliers
    df = df[df['ratio'] < np.percentile(df['ratio'], 99)].copy()
    
    # Correlate with Sigma
    df['log_sigma'] = np.log10(df['sigma_stars'])
    df['log_ratio'] = np.log10(df['ratio'])
    
    # Final check
    df = df.replace([np.inf, -np.inf], np.nan).dropna(subset=['log_sigma', 'log_ratio'])
    
    slope, intercept, r_val, p_val, std_err = stats.linregress(df['log_sigma'], df['log_ratio'])
    
    print(f"  Correlation (log Sigma vs log Ratio): r={r_val:.3f}, p={p_val:.2e}, slope={slope:.4f}")
    
    # Plot
    plt.figure(figsize=(8, 6))
    plt.scatter(df['log_sigma'], np.log10(df['ratio']), s=5, alpha=0.5, c=df['logMass'], cmap='inferno')
    plt.colorbar(label='log Stellar Mass')
    
    x_range = np.linspace(df['log_sigma'].min(), df['log_sigma'].max(), 100)
    plt.plot(x_range, intercept + slope*x_range, 'k--', label=f'Slope={slope:.3f}')
    
    plt.xlabel('log(Velocity Dispersion)')
    plt.ylabel('log(Ha / UV Flux Ratio)')
    plt.title('Test DX: SF Timescale Ratio vs Potential')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    out_path = os.path.join(FIGURES_DIR, 'sdss_test_dx_halpha_uv.png')
    plt.savefig(out_path, dpi=150)
    print(f"Figure saved to {out_path}")
    
    return {
        'correlation_r': r_val,
        'p_value': p_val,
        'slope': slope,
        'n_gal': int(len(df))
    }

def main():
    cache_path = os.path.join(DATA_DIR, 'sdss_halpha_uv.csv')
    
    # Force fresh download since we changed tables/columns
    if os.path.exists(cache_path):
        os.remove(cache_path)
        
    df = download_data()
    if df is not None:
        df.to_csv(cache_path, index=False)
    else:
        print("Download failed.")
        return

    results = analyze_timescale_ratio(df)
    
    out_path = os.path.join(RESULTS_DIR, 'sdss_test_dx_results.json')
    if results:
        with open(out_path, 'w') as f:
            json.dump(results, f, indent=2)
            
        print("\nSUMMARY TEST DX:")
        print(f"Slope (logSigma vs logRatio): {results['slope']:.4f}")
        
        if results['p_value'] < 0.05 and abs(results['slope']) > 0.1:
             print("RESULT: SIGNAL (Ratio depends on Sigma)")
        else:
             print("RESULT: NULL (No significant dependence)")

if __name__ == "__main__":
    main()
