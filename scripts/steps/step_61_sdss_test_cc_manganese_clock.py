#!/usr/bin/env python3
"""
Step 6.57: SDSS Test CC - The Manganese Clock (Type Ia Yield Physics)

Hypothesis:
Manganese (Mn) is produced in Type Ia supernovae, but unlike Iron, its yield is metallicity-dependent 
(production increases with progenitor metallicity). The [Mn/Fe] ratio is therefore a sensitive clock 
of chemical evolution "speed." If time dilation alters the rate of enrichment or the delay time 
distribution of SNe Ia in the inner Galaxy, the [Mn/Fe] vs [Fe/H] trends should diverge from the outer disk.

Prediction:
[Mn/Fe] at fixed [Fe/H] varies with R_gc.

Data:
- aspcapStar: mn_fe, fe_h (or param_fe_h)
- apogee_starhorse: dist50

Method:
1. Select stars with valid Mn, Fe.
2. Compute R_gc.
3. Select metallicity range (e.g., -1.5 < [Fe/H] < -0.5) where the trend is clear.
4. Analyze Mean [Mn/Fe] vs R_gc in this metallicity bin.
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
    print(f"Querying SDSS for Test CC (Limit: {limit})...")
    
    # Strategy: Two-step fetch to avoid timeouts and ensure good overlap
    
    # 1. Fetch Sample from apogee_starhorse
    print("  Fetching apogee_starhorse (Sample)...")
    sql_sample = f"""
    SELECT TOP {limit}
        apogee_id,
        dist50,
        glon, glat
    FROM apogee_starhorse
    WHERE dist50 > 0
    """
    df_sample = query_sdss(sql_sample)
    
    if df_sample is None or len(df_sample) == 0:
        print("  No sample data found.")
        return None
        
    print(f"  Got {len(df_sample)} stars. Fetching chemistry...")
    ids = df_sample['apogee_id'].astype(str).tolist()
    ids = list(set(ids))
    
    # 2. Fetch aspcapStar (Chemistry)
    # mn_fe, fe_h (using generic names, will map later if needed)
    # Using 'mn_fe' and 'fe_h' from previous check
    chunk_size = 100
    df_chem_list = []
    
    for i in range(0, len(ids), chunk_size):
        chunk = ids[i:i+chunk_size]
        ids_str = "', '".join(chunk)
        
        sql_chem = f"""
        SELECT 
            apogee_id,
            mn_fe,
            fe_h
        FROM aspcapStar
        WHERE 
            apogee_id IN ('{ids_str}')
            AND mn_fe > -10
            AND fe_h > -2.5
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

def analyze_manganese(df):
    print("Analyzing Manganese Clock...")
    
    df = compute_rgc(df)
    
    # Clean
    df = df.dropna(subset=['mn_fe', 'fe_h']).copy()
    
    # Select Metallicity Range
    # We want a range where [Mn/Fe] is evolving but not saturated
    # Typical range for disk/halo transition: -1.5 < [Fe/H] < -0.5
    # Let's use a slightly wider bin but control for Fe/H
    
    df_sub = df[(df['fe_h'] > -1.5) & (df['fe_h'] < -0.5)]
    print(f"  Subsample (-1.5 < [Fe/H] < -0.5): {len(df_sub)}")
    
    if len(df_sub) < 50:
        print("  Insufficient stars in metallicity bin.")
        return {'n_sample': 0, 'slope': 0}
        
    # Bin by R_gc
    bins = [0, 4, 6, 8, 10, 15, 30]
    df_sub['r_bin'] = pd.cut(df_sub['R_gc'], bins=bins)
    
    binned = df_sub.groupby('r_bin', observed=False)['mn_fe'].agg(['mean', 'sem', 'count'])
    binned['r_center'] = [i.mid for i in binned.index]
    
    print("\nMean [Mn/Fe] by R_gc (fixed [Fe/H]):")
    print(binned[['mean', 'sem', 'count']])
    
    # Fit trend
    valid_bins = binned[binned['count'] > 10]
    if len(valid_bins) > 2:
        slope, intercept, r_val, p_val, std_err = stats.linregress(valid_bins['r_center'], valid_bins['mean'])
        print(f"  Slope ([Mn/Fe] vs R): {slope:.5f} dex/kpc")
        print(f"  Correlation: {r_val:.3f}")
    else:
        slope = 0
        r_val = 0
        
    # Plot
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # Plot bins
    ax.errorbar(binned['r_center'], binned['mean'], yerr=binned['sem'], fmt='o-', capsize=5, label='Mean [Mn/Fe]')
    
    # Plot raw scatter (background)
    # ax.scatter(df_sub['R_gc'], df_sub['mn_fe'], alpha=0.1, s=2, c='gray')
    
    ax.set_xlabel('Galactocentric Radius [kpc]')
    ax.set_ylabel('[Mn/Fe] (dex)')
    ax.set_title(f'Test CC: Manganese Clock (-1.5 < [Fe/H] < -0.5)')
    ax.grid(True, alpha=0.3)
    ax.legend()
    
    out_path = os.path.join(FIGURES_DIR, 'step_61_sdss_test_cc_manganese.png')
    plt.savefig(out_path, dpi=150)
    print(f"Figure saved to {out_path}")
    
    return {
        'slope': slope,
        'r_val': r_val,
        'n_sample': int(len(df_sub))
    }

def main():
    cache_path = os.path.join(DATA_DIR, 'sdss_manganese.csv')
    
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

    results = analyze_manganese(df)
    
    out_path = os.path.join(RESULTS_DIR, 'step_61_sdss_test_cc_manganese_clock.json')
    if results:
        def default(o):
            if isinstance(o, pd.Interval): return str(o)
            raise TypeError
            
        with open(out_path, 'w') as f:
            json.dump(results, f, indent=2, default=default)
            
        print("\nSUMMARY TEST CC:")
        print("Prediction: [Mn/Fe] varies with R_gc at fixed [Fe/H].")
        print(f"Observed Slope: {results['slope']:.5f}")
        
        if abs(results['slope']) > 0.005:
             print("RESULT: SIGNAL (Gradient observed)")
        else:
             print("RESULT: NULL (Flat profile)")

if __name__ == "__main__":
    main()
