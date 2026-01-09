#!/usr/bin/env python3
"""
Step 6.96: Test G - SN Ia Light Curve Stretch Anomaly

Hypothesis:
Type Ia supernova light curves are standardized using a "stretch" parameter (x1). 
Under TEP, supernovae in deeper gravitational potentials (high-sigma) should 
exhibit systematically longer observed timescales (higher x1) due to time dilation.

Prediction:
Correlation r(x1, sigma) > 0 at fixed color/redshift.

Data:
- Pantheon+ (Local CSV): x1, c, z, RA, DEC
- SDSS Hosts (galSpec tables): v_disp (sigma), coordinates
"""

import pandas as pd
import numpy as np
from scipy import stats, spatial
import os
import json
import matplotlib.pyplot as plt
import requests
import time
from astropy.coordinates import SkyCoord
import astropy.units as u

# Paths
DATA_DIR = os.path.join(os.path.dirname(__file__), '..', '..', 'data')
SN_DIR = os.path.join(DATA_DIR, 'supernovae')
RESULTS_DIR = os.path.join(os.path.dirname(__file__), '..', '..', 'results', 'outputs')
FIGURES_DIR = os.path.join(os.path.dirname(__file__), '..', '..', 'results', 'figures')

SDSS_URL = "https://skyserver.sdss.org/dr18/SkyServerWS/SearchTools/SqlSearch"

def query_sdss_hosts(limit=50000):
    print(f"Querying SDSS Hosts for Test G (Limit: {limit})...")
    
    # Use galSpec tables for robustness against HTTP 500
    sql = f"""
    SELECT TOP {limit}
        gi.specObjID,
        gi.ra, gi.dec, gi.z,
        gi.v_disp, gi.v_disp_err,
        ge.lgm_tot_p50 as logMass
    FROM galSpecInfo gi
    JOIN galSpecExtra ge ON gi.specObjID = ge.specObjID
    WHERE gi.z BETWEEN 0.01 AND 0.15
      AND gi.v_disp > 50 AND gi.v_disp < 400
      AND gi.v_disp_err < 30
      AND gi.reliable = 1
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

def load_pantheon():
    path = os.path.join(SN_DIR, 'pantheon_plus.csv')
    if not os.path.exists(path):
        print("Pantheon+ catalog not found.")
        return None
    return pd.read_csv(path)

def match_catalogs(sn_df, host_df):
    print("Cross-matching catalogs...")
    
    # Coordinates
    c_sn = SkyCoord(ra=sn_df['RA'].values*u.deg, dec=sn_df['DEC'].values*u.deg)
    c_host = SkyCoord(ra=host_df['ra'].values*u.deg, dec=host_df['dec'].values*u.deg)
    
    # Match
    max_sep = 5.0 * u.arcsec
    idx, d2d, d3d = c_sn.match_to_catalog_sky(c_host)
    
    # Filter matches
    mask = d2d < max_sep
    
    matched_sn = sn_df[mask].copy()
    matched_host = host_df.iloc[idx[mask]].copy()
    
    # Check redshift agreement
    z_diff = np.abs(matched_sn['zHD'].values - matched_host['z'].values)
    z_mask = z_diff < 0.01
    
    final_sn = matched_sn[z_mask].reset_index(drop=True)
    final_host = matched_host[z_mask].reset_index(drop=True)
    
    # Merge
    merged = pd.concat([final_sn, final_host.add_suffix('_host')], axis=1)
    
    print(f"  Matches found: {len(merged)}")
    return merged

def analyze_stretch(df):
    print("Analyzing Stretch vs Sigma...")
    
    if len(df) < 50:
        print("  Insufficient matches for robust analysis.")
        return None
    
    # Variables
    x1 = df['x1']
    sigma = df['v_disp']
    log_sigma = np.log10(sigma)
    
    # Remove outliers
    mask = (np.abs(stats.zscore(x1)) < 3) & (np.abs(stats.zscore(log_sigma)) < 3)
    df_clean = df[mask].copy()
    
    # Partial Correlation (Control for Color 'c' and Redshift 'zHD')
    # Simple approach: Residuals
    
    # 1. Regress x1 on color and z
    slope_c, intercept_c, r_c, p_c, _ = stats.linregress(df_clean['c'], df_clean['x1'])
    x1_resid = df_clean['x1'] - (intercept_c + slope_c * df_clean['c'])
    
    # 2. Correlate residuals with log_sigma
    slope, intercept, r_val, p_val, std_err = stats.linregress(df_clean['v_disp'].apply(np.log10), x1_resid)
    
    print(f"  Correlation (x1_resid vs log_sigma): r={r_val:.3f}, p={p_val:.3f}, slope={slope:.4f}")
    
    # Plot
    plt.figure(figsize=(8, 6))
    
    # Handle missing zHD for coloring
    if 'zHD' in df_clean.columns:
        c_vals = df_clean['zHD']
        label = 'Redshift'
    elif 'z' in df_clean.columns:
        c_vals = df_clean['z']
        label = 'Redshift'
    else:
        c_vals = 'blue' # Fallback color
        label = None
        
    scatter = plt.scatter(np.log10(df_clean['v_disp']), x1_resid, s=10, alpha=0.6, c=c_vals, cmap='viridis' if label else None)
    if label:
        plt.colorbar(scatter, label=label)
    
    x_range = np.linspace(np.log10(df_clean['v_disp'].min()), np.log10(df_clean['v_disp'].max()), 100)
    plt.plot(x_range, intercept + slope*x_range, 'k--', label=f'Slope={slope:.3f}')
    
    plt.xlabel('log(Velocity Dispersion)')
    plt.ylabel('Stretch (x1) Residuals (corrected for color)')
    plt.title('Test G: SN Ia Stretch vs Host Potential')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    out_path = os.path.join(FIGURES_DIR, 'sdss_test_g_sn_stretch.png')
    plt.savefig(out_path, dpi=150)
    print(f"Figure saved to {out_path}")
    
    return {
        'r_val': r_val,
        'p_val': p_val,
        'slope': slope,
        'n_sn': int(len(df_clean))
    }

def main():
    # 1. Load Data
    sn_df = load_pantheon()
    if sn_df is None: return
    
    # 2. Query Hosts (Cache if possible)
    cache_path = os.path.join(DATA_DIR, 'sdss_hosts_galspec.csv')
    if os.path.exists(cache_path):
        print("Loading cached hosts...")
        host_df = pd.read_csv(cache_path)
    else:
        host_df = query_sdss_hosts(limit=500000) # Get many to maximize overlap
        if host_df is not None:
            host_df.to_csv(cache_path, index=False)
        else:
            return

    # 3. Match
    merged_df = match_catalogs(sn_df, host_df)
    
    # Fallback if no matches (likely due to limited SDSS query overlap)
    if merged_df is None or len(merged_df) == 0:
        print("Live match failed. checking for pre-existing match file...")
        prematch_path = os.path.join(SN_DIR, 'sn_sdss_sigma_matches.csv')
        if os.path.exists(prematch_path):
             print(f"Loading pre-matched file: {prematch_path}")
             merged_df = pd.read_csv(prematch_path)
             # Rename columns to match expected format if needed
             # Expected: x1, v_disp (from sigma_host)
             if 'sigma_host' in merged_df.columns:
                 merged_df['v_disp'] = merged_df['sigma_host']
             if 'zHD' not in merged_df.columns and 'z' in merged_df.columns: # Pantheon+ usually has zHD
                 merged_df['zHD'] = merged_df['z'] # Proxy if needed
        else:
             print("No pre-existing match file found.")
             return

    # 4. Analyze
    if merged_df is not None:
        results = analyze_stretch(merged_df)
        
        out_path = os.path.join(RESULTS_DIR, 'sdss_test_g_results.json')
        if results:
            with open(out_path, 'w') as f:
                json.dump(results, f, indent=2)
                
            print("\nSUMMARY TEST G:")
            print(f"Slope: {results['slope']:.4f}")
            if results['p_val'] < 0.05:
                if results['slope'] > 0:
                    print("RESULT: SIGNAL (SNe slower in deep potentials)")
                else:
                    print("RESULT: CONTRADICTED (SNe faster in deep potentials)")
            else:
                print("RESULT: NULL (No correlation)")

if __name__ == "__main__":
    main()
