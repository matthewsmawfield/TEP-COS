#!/usr/bin/env python3
"""
Step 6.52: SDSS Test BX - Apogee Velocity Jitter (Close Binary Fraction)

Hypothesis:
Radial velocity scatter (`vscatter`) in multi-epoch APOGEE spectra is a proxy for close binary companions.
Binary orbital evolution (hardening) is a rate process. In deep potentials, binaries might harden slower or survive longer.
We expect the fraction of stars with high velocity jitter (close binaries) to show a gradient with potential depth (R_gc).

Prediction:
Mean Velocity Jitter (or Binary Fraction) correlates with Potential Depth.
Standard: Binary fraction might be constant or depend on formation density.
TEP: Potential-dependent rate of hardening/merging?
     If hardening is slower in deep potential -> More wide binaries, fewer close/hard binaries?
     Or survival of binaries is higher?
     Let's look for a gradient in High-Jitter fraction.

Data:
- apogeeStar: vscatter, nvisits, glon, glat, param_logg
- apogee_starhorse: dist50

Method:
1. Fetch vscatter, nvisits, coords.
2. Filter for Multi-epoch stars (nvisits > 2).
3. Compute R_gc.
4. Define "Binary" candidates as vscatter > 1 km/s (or some threshold).
5. Bin by R_gc.
6. Calculate Binary Fraction vs R_gc.
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
    print(f"Querying SDSS for Test BX (Limit: {limit})...")
    
    # Strategy: 
    # 1. Fetch Sample from apogee_starhorse (Dist, Coords, Logg). 
    #    StarHorse table seems robust.
    # 2. Fetch Jitter from apogeeStar for these IDs (Minimal columns).
    
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
        AND logg50 BETWEEN 0.5 AND 3.5 -- Giants
    """
    df_sample = query_sdss(sql_sample)
    
    if df_sample is None or len(df_sample) == 0:
        print("  No sample data found.")
        return None
        
    print(f"  Got {len(df_sample)} stars. Fetching jitter...")
    ids = df_sample['apogee_id'].astype(str).tolist()
    ids = list(set(ids))
    
    # 2. Fetch apogeeStar (Jitter)
    # Minimal columns: apogee_id, vscatter, nvisits
    chunk_size = 100
    df_jitter_list = []
    
    for i in range(0, len(ids), chunk_size):
        chunk = ids[i:i+chunk_size]
        ids_str = "', '".join(chunk)
        
        sql_jitter = f"""
        SELECT 
            apogee_id,
            vscatter,
            nvisits
        FROM apogeeStar
        WHERE 
            apogee_id IN ('{ids_str}')
        """
        res = query_sdss(sql_jitter)
        if res is not None and len(res) > 0:
            df_jitter_list.append(res)
        time.sleep(0.1)
            
    if not df_jitter_list:
        print("  No jitter data found.")
        return None
        
    df_jitter = pd.concat(df_jitter_list, ignore_index=True)
    
    # 3. Join
    print("  Joining datasets...")
    df = pd.merge(df_sample, df_jitter, on='apogee_id', how='inner')
    
    # Filter for valid jitter
    df = df[(df['nvisits'] > 2) & (df['vscatter'] > 0)]
    
    print(f"  Merged & Filtered N={len(df)}")
    
    return df

def compute_rgc(df):
    R0 = 8.2 # kpc
    df['dist_kpc'] = df['dist50'] # Assuming kpc (standard for StarHorse in these contexts)
    
    l_rad = np.radians(df['glon'])
    b_rad = np.radians(df['glat'])
    
    d_proj = df['dist_kpc'] * np.cos(b_rad)
    
    # R_gc formula
    df['R_plane_sq'] = R0**2 + d_proj**2 - 2 * R0 * d_proj * np.cos(l_rad)
    df['Z'] = df['dist_kpc'] * np.sin(b_rad)
    df['R_gc'] = np.sqrt(df['R_plane_sq'] + df['Z']**2)
    
    return df

def analyze_jitter(df):
    print("Analyzing Velocity Jitter...")
    
    df = compute_rgc(df)
    
    # Clean R_gc
    df = df[(df['R_gc'] > 0) & (df['R_gc'] < 30)]
    
    print(f"  R_gc range: {df['R_gc'].min():.1f} - {df['R_gc'].max():.1f} kpc")
    print(f"  Sample size: {len(df)}")
    
    # Define "Jittery" (Binary Candidate)
    # Threshold: vscatter > 3 km/s (typical for binaries vs intrinsic jitter < 0.5 km/s)
    threshold = 3.0
    df['is_binary'] = df['vscatter'] > threshold
    
    print(f"  Binary Candidates (vscatter > {threshold} km/s): {df['is_binary'].sum()} ({df['is_binary'].mean():.1%})")
    
    # Bin by R_gc
    bins = np.linspace(0, 20, 11)
    df['r_bin'] = pd.cut(df['R_gc'], bins=bins)
    
    binned = df.groupby('r_bin')['is_binary'].agg(['mean', 'sem', 'count'])
    binned['r_center'] = [i.mid for i in binned.index]
    
    print("\nBinary Fraction by R_gc:")
    print(binned[['mean', 'sem', 'count']])
    
    # Fit trend
    valid_bins = binned[binned['count'] > 50]
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
    ax.set_ylabel(f'High Jitter Fraction (> {threshold} km/s)')
    ax.set_title('Test BX: Binary Fraction (Velocity Jitter) vs Radius')
    ax.grid(True, alpha=0.3)
    
    out_path = os.path.join(FIGURES_DIR, 'sdss_test_bx_jitter.png')
    plt.savefig(out_path, dpi=150)
    print(f"Figure saved to {out_path}")
    
    return {
        'slope': slope,
        'r_val': r_val,
        'n_sample': int(len(df))
    }

def main():
    cache_path = os.path.join(DATA_DIR, 'sdss_jitter.csv')
    
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

    results = analyze_jitter(df)
    
    out_path = os.path.join(RESULTS_DIR, 'sdss_test_bx_results.json')
    if results:
        def default(o):
            if isinstance(o, pd.Interval): return str(o)
            raise TypeError
            
        with open(out_path, 'w') as f:
            json.dump(results, f, indent=2, default=default)
            
        print("\nSUMMARY TEST BX:")
        print("TEP Prediction: Gradient in Binary Fraction with Potential Depth.")
        print(f"Observed Slope: {results['slope']:.5f}")
        
        if abs(results['slope']) > 0.001: # Significant gradient
             print("RESULT: SIGNAL (Gradient observed)")
        else:
             print("RESULT: NULL (Flat profile)")

if __name__ == "__main__":
    main()
