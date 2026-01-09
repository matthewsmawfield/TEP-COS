#!/usr/bin/env python3
"""
Step 6.60: SDSS Test CF - QSO Gravitational Redshift (Potential Mapping)

Hypothesis:
Broad emission lines in Quasars come from the BLR (deep in the BH potential), 
while narrow lines come from the NLR (further out). 
There is a gravitational redshift difference Delta_z = z_broad - z_narrow.
TEP modifies the potential Phi. The observed shift should scale with M_BH 
but with a normalization or slope distinct from GR if the metric is modified by the scalar field.

Prediction:
Gravitational Redshift Delta_v scales with M_BH differently than GR predicts.
GR: v ~ sqrt(M/R) ~ M^0.5 (if R ~ const?) or if R ~ L^0.5.
If R_BLR ~ L^0.5, then v ~ M^0.5 * L^-0.25.
We are looking for Delta_z (broad-narrow). z_broad > z_narrow (redshifted).
Delta_v = c * (z_broad - z_narrow) / (1 + z_narrow) approx.

Data:
- spAll: z (best fit, often dominated by narrow lines or template), z_err
- mos_sdss_dr16_qso: z_pca (PCA redshift, arguably systemic?), z_civ, z_mgii (Broad lines), logBH
  Actually, DR16Q has redshifts for different lines.
  z_sys usually from [OIII] or narrow Mg/C lines.
  z_broad from Broad components.
  DR16Q columns: Z_HW (Hewett & Wild, improved systemic), Z_PCA, Z_CIV, Z_CIII, Z_MgII.
  We want Delta_z = Z_Broad (MgII or CIV) - Z_Systemic (HW or PCA).

Method:
1. Fetch DR16Q data: Z_HW (Systemic), Z_MgII (Broad), logBH_MgII (Mass).
   Focus on Mg II region (0.4 < z < 2.0) as it's cleaner than CIV (outflows).
2. Calculate Delta_v = c * (Z_MgII - Z_HW) / (1 + Z_HW).
3. Analyze Delta_v vs Mass.
   Expect Positive Delta_v (Grav Redshift).
   Correlation with Mass?
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
    print(f"Querying SDSS for Test CF (Limit: {limit})...")
    
    # Use mos_sdss_dr16_qso directly as it has all Z estimates
    # Z_HW is generally considered the best systemic redshift (narrow lines/improved template)
    # Z_MgII is the broad line redshift
    
    sql = f"""
    SELECT TOP {limit}
        SPECOBJID,
        Z as z_pipe,
        Z_HW as z_sys,
        Z_MgII as z_broad,
        logBH_MgII as log_bh,
        LOGL1350 as log_lum
        
    FROM mos_sdss_dr16_qso
    WHERE 
        Z_HW > 0 AND Z_MgII > 0
        AND logBH_MgII > 6
    """
    return query_sdss(sql)

def analyze_grav_redshift(df):
    print("Analyzing Gravitational Redshift...")
    
    # Clean
    df = df.dropna().copy()
    
    # Calculate Velocity Shift
    # Delta_v = c * (z_broad - z_sys) / (1 + z_sys)
    c_kms = 299792.458
    df['delta_v'] = c_kms * (df['z_broad'] - df['z_sys']) / (1 + df['z_sys'])
    
    # Filter outliers (outflows can cause huge blueshifts - thousands of km/s)
    # Gravitational redshift should be small positive (< 1000 km/s?)
    # Broad lines are often blueshifted due to outflows (e.g. CIV). MgII is more stable.
    # But if outflows dominate, we might see negative delta_v.
    # We want to see if the *component* of grav redshift scales with Mass.
    # Outflows might also scale with L/Mass.
    
    print(f"  Sample size: {len(df)}")
    print(f"  Mean Delta_v: {df['delta_v'].mean():.2f} km/s")
    print(f"  Median Delta_v: {df['delta_v'].median():.2f} km/s")
    
    # Limit range to avoid extreme outflows
    df_clean = df[(df['delta_v'] > -2000) & (df['delta_v'] < 2000)].copy()
    print(f"  Clean Sample (|dv| < 2000): {len(df_clean)}")
    
    # Correlations
    r_mass, p_mass = stats.pearsonr(df_clean['log_bh'], df_clean['delta_v'])
    print(f"  Correlation r(Delta_v, logM): {r_mass:.4f} (p={p_mass:.2e})")
    
    slope, intercept, _, _, _ = stats.linregress(df_clean['log_bh'], df_clean['delta_v'])
    print(f"  Slope: {slope:.2f} km/s/dex")
    
    # Control for Luminosity (Outflow driver?)
    r_lum, p_lum = stats.pearsonr(df_clean['log_lum'], df_clean['delta_v'])
    print(f"  Correlation r(Delta_v, logL): {r_lum:.4f}")
    
    # Plot
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # Scatter
    # ax.scatter(df_clean['log_bh'], df_clean['delta_v'], alpha=0.1, s=2, c='gray')
    
    # Binning
    df_clean['mass_bin'] = pd.qcut(df_clean['log_bh'], 10)
    binned = df_clean.groupby('mass_bin')['delta_v'].agg(['mean', 'sem', 'count'])
    binned['mass_center'] = [i.mid for i in binned.index]
    
    ax.errorbar(binned['mass_center'], binned['mean'], yerr=binned['sem'], fmt='o-', capsize=5, label='Mean Shift')
    ax.set_xlabel('log(BH Mass) [M_sun]')
    ax.set_ylabel('Velocity Shift (Broad - Sys) [km/s]')
    ax.set_title(f'Test CF: QSO Gravitational Redshift (r={r_mass:.2f})')
    ax.axhline(0, color='k', linestyle='--', alpha=0.3)
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    out_path = os.path.join(FIGURES_DIR, 'sdss_test_cf_redshift.png')
    plt.savefig(out_path, dpi=150)
    print(f"Figure saved to {out_path}")
    
    return {
        'slope': slope,
        'r_val': r_mass,
        'mean_shift': df_clean['delta_v'].mean(),
        'n_sample': int(len(df_clean))
    }

def main():
    cache_path = os.path.join(DATA_DIR, 'sdss_grav_redshift.csv')
    
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

    results = analyze_grav_redshift(df)
    
    out_path = os.path.join(RESULTS_DIR, 'sdss_test_cf_results.json')
    with open(out_path, 'w') as f:
        json.dump(results, f, indent=2)
        
    print("\nSUMMARY TEST CF:")
    print("Prediction: Gravitational redshift (positive shift) scales with Mass.")
    print(f"Observed Mean Shift: {results['mean_shift']:.2f} km/s")
    print(f"Observed Slope vs Mass: {results['slope']:.2f} km/s/dex")
    
    if results['slope'] > 50: # Expect ~100-500 km/s range for massive BHs
         print("RESULT: CONSISTENT (Positive scaling observed)")
    elif results['slope'] < 0:
         print("RESULT: CONTRADICTED (Negative scaling - Outflows dominate?)")
    else:
         print("RESULT: NULL (Weak/No scaling)")

if __name__ == "__main__":
    main()
