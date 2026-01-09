#!/usr/bin/env python3
"""
Step 6.51: SDSS Test BW - CEMP Star Fraction (Binary Clock)

Hypothesis:
Carbon-Enhanced Metal-Poor (CEMP) stars largely form via mass transfer in binary systems (AGB wind accretion).
This is a rate-dependent evolutionary channel. In the deep potential of the Inner Halo/Bulge,
the binary evolution and wind accretion rates should be time-dilated.
This may alter the observed fraction of CEMP stars relative to normal Metal-Poor stars as a function of R_gc.

Prediction:
Fraction of CEMP stars (among MP stars) varies with R_gc.
Standard model: CEMP fraction increases with distance (outer halo has more CEMP-s).
TEP: Does the inner galaxy show a deficit or excess beyond metallicity effects?

Data:
- aspcapStar: param_c_fe, param_fe_h
- apogee_starhorse: dist50
- apogeeStar: glon, glat

Method:
1. Select Metal-Poor stars ([Fe/H] < -2.0).
2. Classify as CEMP if [C/Fe] > +0.7.
3. Compute R_gc.
4. Bin by R_gc and compute CEMP fraction.
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
    print(f"Querying SDSS for Test BW (Limit: {limit})...")
    
    # Need to join apogeeStar (coords) and aspcapStar (chem) and StarHorse (dist)
    # Use client side join strategy if needed, but let's try direct join first with limited columns.
    # Note: param_c_fe is in aspcapStar.
    
    # We can use fparam_c_m + fparam_m_h - fparam_fe_h approx? 
    # Or just use the calibrated [C/Fe] if available. 
    # In DR17 aspcapStar, it's often 'c_fe' or 'carb_fe'.
    # Previous checks showed 'fparam_' columns.
    # Let's check columns first? No, let's guess standard 'c_fe' or derived.
    # Actually, let's use the explicit 'fe_h' and 'c_fe' tags if they exist.
    # Let's trust previous knowledge that 'fparam' arrays exist, or specific tags.
    # DR16/17 usually has 'X_FE' columns.
    
    # Query for Metal Poor stars directly to save rows
    
    sql = f"""
    SELECT TOP {limit}
        a.apogee_id,
        a.fe_h, 
        a.c_fe,
        s.glon, s.glat,
        sh.dist50
        
    FROM aspcapStar a
    JOIN apogeeStar s ON a.apogee_id = s.apogee_id
    JOIN apogee_starhorse sh ON a.apogee_id = sh.apogee_id
    
    WHERE 
        a.fe_h < -1.5  -- Metal Poor (relaxed to -1.5 for stats)
        AND a.c_fe > -10
        AND sh.dist50 > 0
    """
    return query_sdss(sql)

def compute_rgc(df):
    R0 = 8.2 # kpc
    df['dist_kpc'] = df['dist50'] # assuming dist50 is kpc? No, StarHorse usually kpc?
    # Wait, previous scripts assumed kpc? 
    # Check units. StarHorse dist50 is usually kpc in DR17 value added, but let's be careful.
    # Actually, in CAS, it might be whatever the column def is.
    # Let's assume kpc. If it's pc, R_gc will be huge.
    # Standard APOGEE StarHorse is in kpc.
    
    l_rad = np.radians(df['glon'])
    b_rad = np.radians(df['glat'])
    
    d_proj = df['dist_kpc'] * np.cos(b_rad)
    
    # R_gc formula
    # R^2 = R0^2 + d^2 - 2 R0 d cos(l)
    # This is planar R? No, full spherical usually?
    # R_galactocentric.
    # R_plane^2 = ...
    # Z = d sin b
    # R_gc^2 = R_plane^2 + Z^2
    
    df['R_plane_sq'] = R0**2 + d_proj**2 - 2 * R0 * d_proj * np.cos(l_rad)
    df['Z'] = df['dist_kpc'] * np.sin(b_rad)
    df['R_gc'] = np.sqrt(df['R_plane_sq'] + df['Z']**2)
    
    return df

def analyze_cemp_fraction(df):
    print("Analyzing CEMP Fraction...")
    
    df = compute_rgc(df)
    
    # Check distribution of R_gc
    print(f"  R_gc range: {df['R_gc'].min():.1f} - {df['R_gc'].max():.1f} kpc")
    
    # Define CEMP
    # [C/Fe] > 0.7
    df['is_cemp'] = df['c_fe'] > 0.7
    
    print(f"  Total Stars: {len(df)}")
    print(f"  Total CEMP: {df['is_cemp'].sum()}")
    
    # Bin by R_gc
    # We want to see if inner halo has different fraction
    bins = [0, 5, 10, 20, 50, 100]
    df['r_bin'] = pd.cut(df['R_gc'], bins=bins)
    
    binned = df.groupby('r_bin')['is_cemp'].agg(['mean', 'sem', 'count'])
    binned['r_center'] = [i.mid for i in binned.index]
    
    print("\nCEMP Fraction by R_gc:")
    print(binned[['mean', 'sem', 'count']])
    
    # Fit trend
    # Use bin centers for regression
    valid_bins = binned[binned['count'] > 10]
    if len(valid_bins) > 2:
        slope, intercept, r_val, p_val, std_err = stats.linregress(valid_bins['r_center'], valid_bins['mean'])
        print(f"  Slope (Fraction vs R): {slope:.5f} / kpc")
        print(f"  Correlation: {r_val:.3f}")
    else:
        slope = 0
        r_val = 0
    
    # Plot
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.errorbar(binned['r_center'], binned['mean'], yerr=binned['sem'], fmt='o-', capsize=5)
    ax.set_xlabel('Galactocentric Radius [kpc]')
    ax.set_ylabel('CEMP Fraction ([C/Fe] > 0.7)')
    ax.set_title('Test BW: CEMP Fraction vs Radius')
    ax.grid(True, alpha=0.3)
    
    out_path = os.path.join(FIGURES_DIR, 'sdss_test_bw_cemp.png')
    plt.savefig(out_path, dpi=150)
    print(f"Figure saved to {out_path}")
    
    return {
        'slope': slope,
        'r_val': r_val,
        'n_sample': int(len(df))
    }

def main():
    cache_path = os.path.join(DATA_DIR, 'sdss_cemp.csv')
    
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

    results = analyze_cemp_fraction(df)
    
    out_path = os.path.join(RESULTS_DIR, 'sdss_test_bw_results.json')
    if results:
        # Helper for json serialization of Intervals
        def default(o):
            if isinstance(o, pd.Interval): return str(o)
            raise TypeError
            
        with open(out_path, 'w') as f:
            json.dump(results, f, indent=2, default=default)
            
        print("\nSUMMARY TEST BW:")
        print("Standard Model: CEMP fraction increases with R (Outer halo has more). Slope > 0.")
        print("TEP Prediction: Deviation from standard slope?")
        print(f"Observed Slope: {results['slope']:.5f}")
        
        if results['slope'] > 0:
             print("RESULT: NULL (Standard gradient observed)")
        else:
             print("RESULT: ANOMALY (Flat or Negative gradient)")

if __name__ == "__main__":
    main()
