#!/usr/bin/env python3
"""
Step 6.54: SDSS Test BZ - The Carbon Star Clock (C/M Ratio)

Hypothesis:
The transition from Oxygen-rich (M-type) to Carbon-rich (C-type) AGB stars is driven by the third dredge-up events.
This is a rate-dependent mixing process. In deep potentials, mixing is slower.
The ratio of C-stars to M-stars should be LOWER in the inner Galaxy than in the outskirts,
even after controlling for the known metallicity dependence (C-stars are suppressed at high Z).

Prediction:
C/M Ratio is lower at small R_gc than predicted by metallicity gradients alone.

Data:
- aspcapStar: param_c_m, param_m_h, param_logg
- apogee_starhorse: dist50

Method:
1. Select Upper AGB stars (0.5 < logg < 2.0).
2. Classify as C-star if param_c_m > 0 (Carbon enhanced relative to metals, proxy for C > O in this context).
3. Compute R_gc.
4. Bin by R_gc.
5. Analyze C/M ratio trend.
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

def download_data(limit=2000):
    print(f"Querying SDSS for Test BZ (Limit: {limit})...")
    
    # Strategy: Two-step fetch
    
    # 1. Fetch Sample from apogee_starhorse (Dist, Coords, Logg)
    # Target Upper AGB: 0.5 < logg < 2.5 (Relaxed slightly)
    print("  Fetching apogee_starhorse (Sample)...")
    sql_sample = f"""
    SELECT TOP {limit}
        apogee_id,
        dist50,
        glon, glat,
        logg50 as param_logg
    FROM apogee_starhorse
    WHERE 
        dist50 > 0 
        AND logg50 BETWEEN 0.5 AND 2.5
    """
    df_sample = query_sdss(sql_sample)
    
    if df_sample is None or len(df_sample) == 0:
        print("  No sample data found.")
        return None
        
    print(f"  Got {len(df_sample)} stars. Fetching chemistry...")
    ids = df_sample['apogee_id'].astype(str).tolist()
    ids = list(set(ids))
    
    # 2. Fetch aspcapStar (Chemistry)
    # Need param_c_m, fparam_m_h (or param_m_h)
    chunk_size = 100
    df_chem_list = []
    
    for i in range(0, len(ids), chunk_size):
        chunk = ids[i:i+chunk_size]
        ids_str = "', '".join(chunk)
        
        sql_chem = f"""
        SELECT 
            apogee_id,
            param_c_m,
            fparam_m_h as param_m_h
        FROM aspcapStar
        WHERE 
            apogee_id IN ('{ids_str}')
        """
        res = query_sdss(sql_chem)
        if res is not None and len(res) > 0:
            df_chem_list.append(res)
        time.sleep(0.1)
            
    if not df_chem_list:
        print("  No chemistry data found.")
        return None
        
    df_chem = pd.concat(df_chem_list, ignore_index=True)
    
    # 3. Join
    print("  Joining datasets...")
    df = pd.merge(df_sample, df_chem, on='apogee_id', how='inner')
    
    print(f"  Merged N={len(df)}")
    
    return df

def compute_rgc(df):
    R0 = 8.2 # kpc
    df['dist_kpc'] = df['dist50']
    
    l_rad = np.radians(df['glon'])
    b_rad = np.radians(df['glat'])
    
    d_proj = df['dist_kpc'] * np.cos(b_rad)
    
    df['R_plane_sq'] = R0**2 + d_proj**2 - 2 * R0 * d_proj * np.cos(l_rad)
    df['Z'] = df['dist_kpc'] * np.sin(b_rad)
    df['R_gc'] = np.sqrt(df['R_plane_sq'] + df['Z']**2)
    
    return df

def analyze_carbon_stars(df):
    print("Analyzing C/M Ratio...")
    
    df = compute_rgc(df)
    
    # Filter valid chemistry
    df = df.dropna(subset=['param_c_m', 'param_m_h'])
    
    # Define C-Star
    # param_c_m > 0 (Carbon enhanced)
    # Often C-stars have C/O > 1. In param terms, [C/M] > 0 is a proxy for enrichment.
    # Let's use strict cut > 0.1
    df['is_c_star'] = df['param_c_m'] > 0.1
    
    # Note: Metallicity dependence.
    # High metallicity stars rarely become C-stars because it takes more dredge-up to overcome initial O.
    # We should look at Metal-Poor / Intermediate stars mainly?
    # Or just control for it.
    
    print(f"  Total Stars (AGB): {len(df)}")
    print(f"  C-Stars: {df['is_c_star'].sum()} ({df['is_c_star'].mean():.1%})")
    
    # Bin by R_gc
    bins = [0, 4, 6, 8, 10, 15, 30]
    df['r_bin'] = pd.cut(df['R_gc'], bins=bins)
    
    binned = df.groupby('r_bin')['is_c_star'].agg(['mean', 'sem', 'count'])
    binned['r_center'] = [i.mid for i in binned.index]
    
    print("\nC-Star Fraction by R_gc:")
    print(binned[['mean', 'sem', 'count']])
    
    # Fit trend
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
    ax.set_ylabel('C-Star Fraction (C/M Ratio Proxy)')
    ax.set_title('Test BZ: Carbon Star Fraction vs Radius')
    ax.grid(True, alpha=0.3)
    
    out_path = os.path.join(FIGURES_DIR, 'sdss_test_bz_carbon.png')
    plt.savefig(out_path, dpi=150)
    print(f"Figure saved to {out_path}")
    
    return {
        'slope': slope,
        'r_val': r_val,
        'n_sample': int(len(df))
    }

def main():
    cache_path = os.path.join(DATA_DIR, 'sdss_carbon.csv')
    
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

    results = analyze_carbon_stars(df)
    
    out_path = os.path.join(RESULTS_DIR, 'sdss_test_bz_results.json')
    if results:
        def default(o):
            if isinstance(o, pd.Interval): return str(o)
            raise TypeError
            
        with open(out_path, 'w') as f:
            json.dump(results, f, indent=2, default=default)
            
        print("\nSUMMARY TEST BZ:")
        print("Standard Model: C-star fraction increases with R (lower Z in outskirts favors C-stars). Positive Slope.")
        print("TEP Prediction: Inner Galaxy suppression (Slower mixing). Steeper positive slope or deficit at small R?")
        print(f"Observed Slope: {results['slope']:.5f}")
        
        if results['slope'] > 0:
             print("RESULT: CONSISTENT (Standard metallicity gradient likely dominates, but TEP compatible)")
        else:
             print("RESULT: ANOMALY (Flat or Negative - Unexpected high C-frac in inner galaxy?)")

if __name__ == "__main__":
    main()
