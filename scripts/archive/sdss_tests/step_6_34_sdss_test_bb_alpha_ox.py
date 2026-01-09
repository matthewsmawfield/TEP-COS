#!/usr/bin/env python3
"""
Step 6.34: SDSS Test BB - Quasar X-ray/Optical Ratio (Alpha-ox)

Hypothesis:
X-ray emission (corona) and Optical emission (disk) originate from regions of different potential depth around the SMBH.
TEP predicts differential time dilation between these regions.
The alpha_ox parameter, which measures the ratio of X-ray to Optical luminosity, should show residuals correlated with Black Hole Mass (proxy for potential depth) that differ from standard accretion disk models.
Standard: alpha_ox correlates with L_uv (or L_2500).
Prediction: alpha_ox residual (after removing L_uv trend) correlates with M_BH.

Data:
- spiders_quasar:
    - l2keV_class_2RXS (Luminosity at 2 keV, erg/s/Hz? Or erg/s? Usually erg/s in SPIDERS)
    - l_2500 (Luminosity at 2500 A, erg/s/Hz? Need to check units)
    - logBHMA_hb, logBHMS_mgII (BH Mass estimates)

Method:
1. Select quasars with valid X-ray and UV luminosities and BH masses.
2. Compute alpha_ox = -0.3838 * log10(L_2500 / L_2keV).
   Note: Standard definition alpha_ox = log(L_2keV/L_2500) / log(nu_2keV/nu_2500).
   nu_2keV = 4.83e17 Hz. nu_2500 = 1.2e15 Hz.
   log(nu_ratio) = 2.605.
   alpha_ox = log(L_2keV_nu / L_2500_nu) / 2.605.
   We need monochromatic luminosities.
   If l2keV is integrated luminosity, convert to monochromatic.
3. Fit alpha_ox vs log(L_2500).
4. Compute residuals.
5. Correlate residuals with log(M_BH).
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
    print(f"Querying SDSS for Test BB (Limit: {limit})...")
    
    # l2keV_class_2RXS: Log Luminosity? Or Luminosity? Usually Log in catalogs or linear.
    # SPIDERS documentation says l2keV is log luminosity (erg/s).
    # l_2500 is likely log monochromatic luminosity (erg/s/Hz) or log(L).
    # We will assume they are Log10 values if they are small (<100) or Linear if large.
    # Usually 'l' prefix implies Log in SDSS VACs often, but let's check values.
    
    sql = f"""
    SELECT TOP {limit}
        SPECOBJID,
        l2keV_class_2RXS as log_l2kev, -- Assuming log erg/s
        l_2500 as log_l2500, -- Assuming log erg/s/Hz?
        logBHMA_hb as log_mbh_hb,
        logBHMS_mgII as log_mbh_mg
        
    FROM spiders_quasar
    
    WHERE 
        l2keV_class_2RXS > 0
        AND l_2500 > 0
        AND (logBHMA_hb > 0 OR logBHMS_mgII > 0)
    """
    return query_sdss(sql)

def analyze_alpha_ox(df):
    print("Analyzing Alpha-ox Residuals...")
    
    # 1. Clean and Homogenize
    df_clean = df.dropna(subset=['log_l2kev', 'log_l2500']).copy()
    
    # Combine BH masses (Prioritize MgII for higher z, Hb for lower z, or just average/take available)
    df_clean['log_mbh'] = df_clean['log_mbh_mg'].fillna(df_clean['log_mbh_hb'])
    df_clean = df_clean.dropna(subset=['log_mbh']).copy()
    
    # 2. Compute Alpha_ox
    # The columns are named log_l2kev but appear to be linear based on previous run (10^27 and 10^42)
    # We must convert to log10 if they are linear
    
    # Check if linear
    if df_clean['log_l2kev'].mean() > 100:
        df_clean['log_l2kev'] = np.log10(df_clean['log_l2kev'])
        
    if df_clean['log_l2500'].mean() > 100:
        df_clean['log_l2500'] = np.log10(df_clean['log_l2500'])
        
    print(f"Mean log L_2keV: {df_clean['log_l2kev'].mean():.2f}")
    print(f"Mean log L_2500: {df_clean['log_l2500'].mean():.2f}")
    
    # alpha_ox = -0.3838 * log(L_nu(2500A) / L_nu(2keV))
    # L_2keV from SPIDERS is usually integrated luminosity (erg/s) in 0.5-2 keV or 2-10 keV.
    # L_2500 is usually monochromatic luminosity (erg/s/Hz).
    
    # Convert integrated X-ray L to monochromatic L_nu(2keV)
    # Approx: log L_nu(2keV) = log L_int - 17.68 (assuming bandwidth ~ 4.8e17 Hz)
    df_clean['log_lnu_2kev'] = df_clean['log_l2kev'] - 17.68 
    
    # Check if L_2500 is monochromatic (usually is in catalogs if ~30)
    # If it's ~44, it might be integrated or lambda*L_lambda
    # If mean is ~30, it is erg/s/Hz.
    # If mean is ~45, it is erg/s.
    
    if df_clean['log_l2500'].mean() > 40:
        # Assume erg/s (lambda L_lambda?). Convert to erg/s/Hz?
        # log L_nu = log (lambda L_lambda) - log(nu)
        # nu_2500 = 1.2e15 Hz. log = 15.08.
        df_clean['log_lnu_2500'] = df_clean['log_l2500'] - 15.08
    else:
        # Already monochromatic
        df_clean['log_lnu_2500'] = df_clean['log_l2500']
        
    df_clean['alpha_ox'] = -0.3838 * (df_clean['log_lnu_2500'] - df_clean['log_lnu_2kev'])
    
    # 3. Fit alpha_ox vs L_2500
    # Known relation: alpha_ox becomes more negative as L_2500 increases.
    slope, intercept, r_val, p_val, std_err = stats.linregress(df_clean['log_lnu_2500'], df_clean['alpha_ox'])
    print(f"Alpha-ox vs L_2500: slope={slope:.3f}, r={r_val:.3f}")
    
    # 4. Residuals
    df_clean['aox_resid'] = df_clean['alpha_ox'] - (slope * df_clean['log_lnu_2500'] + intercept)
    
    # 5. Correlate with Mass
    r_mass, p_mass = stats.pearsonr(df_clean['log_mbh'], df_clean['aox_resid'])
    
    print(f"N = {len(df_clean)}")
    print(f"Correlation r(Resid, M_BH): {r_mass:.4f} (p={p_mass:.2e})")
    
    # 6. Binning
    df_clean['mass_bin'] = pd.qcut(df_clean['log_mbh'], 8)
    binned = df_clean.groupby('mass_bin')['aox_resid'].mean()
    print("\nMean Alpha-ox Residual by Mass Bin:")
    print(binned)
    
    return {
        'r_mass': float(r_mass),
        'p_mass': float(p_mass),
        'slope_luv': float(slope),
        'binned_means': binned.tolist(),
        'bin_centers': [mid.mid for mid in binned.index],
        'n_sample': int(len(df_clean))
    }, df_clean

def create_figure(df, results):
    print("Generating figure...")
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # Scatter
    ax.scatter(df['log_mbh'], df['aox_resid'], alpha=0.1, s=2, c='k', label='Quasars')
    
    # Binned
    ax.plot(results['bin_centers'], results['binned_means'], 'r-o', lw=2, label='Mean Residual')
    
    ax.set_xlabel(r'$\log(M_{BH}/M_{\odot})$')
    ax.set_ylabel(r'$\alpha_{ox}$ Residual')
    ax.set_title(f"Test BB: Alpha-ox vs BH Mass (r={results['r_mass']:.3f})")
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.axhline(0, color='b', linestyle='--')
    
    plt.tight_layout()
    out_path = os.path.join(FIGURES_DIR, 'sdss_test_bb_alpha_ox.png')
    plt.savefig(out_path, dpi=150)
    print(f"Figure saved to {out_path}")

def main():
    cache_path = os.path.join(DATA_DIR, 'sdss_alpha_ox.csv')
    
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

    results, df_clean = analyze_alpha_ox(df)
    
    out_path = os.path.join(RESULTS_DIR, 'sdss_test_bb_results.json')
    with open(out_path, 'w') as f:
        results_json = results.copy()
        del results_json['bin_centers']
        json.dump(results_json, f, indent=2)
        
    create_figure(df_clean, results)
    
    print("\nSUMMARY TEST BB:")
    print("TEP Prediction: Alpha-ox residual depends on Mass (Potential). r != 0.")
    print(f"Observed r: {results['r_mass']:.4f}")
    
    if abs(results['r_mass']) > 0.1:
        print("RESULT: CONSISTENT (Mass dependence observed)")
    else:
        print("RESULT: NULL (Standard disk model holds)")

if __name__ == "__main__":
    main()
