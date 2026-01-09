#!/usr/bin/env python3
"""
Step 6.69: SDSS Test CR - The Gravity Schism (Spec vs Isochrone log g)

Hypothesis:
Spectroscopic surface gravity (log g_spec) relies on line broadening (pressure/collisions). 
Isochrone surface gravity (log g_iso) relies on Luminosity/Temperature (Stefan-Boltzmann). 
TEP affects these physics differently (collisions vs radiation). 
We expect a systematic residual Delta log g = log g_spec - log g_iso that correlates with the background potential Phi.

Prediction:
Delta log g correlates with Galactocentric Radius R_gc.

Data:
- aspcapStar: logg (spec), glon, glat
- apogee_starhorse: logg50 (iso), dist50

Method:
1. Fetch spec logg and iso logg.
2. Compute Delta log g.
3. Compute R_gc.
4. Analyze Delta log g vs R_gc.
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

def download_data(limit=1000):
    print(f"Querying SDSS for Test CR (Limit: {limit})...")
    
    # 1. Fetch Sample from apogee_starhorse
    print("  Fetching apogee_starhorse (Sample)...")
    sql_sample = f"""
    SELECT TOP {limit}
        apogee_id,
        dist50,
        logg50 as logg_iso,
        glon, glat
    FROM apogee_starhorse
    WHERE dist50 > 0
    """
    df_sample = query_sdss(sql_sample)
    
    if df_sample is None or len(df_sample) == 0:
        print("  No sample data found.")
        return None
        
    print(f"  Got {len(df_sample)} stars. Fetching spec logg...")
    ids = df_sample['apogee_id'].astype(str).tolist()
    ids = list(set(ids))
    
    # 2. Fetch aspcapStar (Spec Gravity)
    # Using 'logg' from check script
    chunk_size = 50
    df_spec_list = []
    
    for i in range(0, len(ids), chunk_size):
        chunk = ids[i:i+chunk_size]
        ids_str = "', '".join(chunk)
        
        sql_spec = f"""
        SELECT 
            apogee_id,
            logg as logg_spec
        FROM aspcapStar
        WHERE 
            apogee_id IN ('{ids_str}')
            AND logg > -5
        """
        res = query_sdss(sql_spec)
        if res is not None and len(res) > 0:
            df_spec_list.append(res)
        time.sleep(0.2)
            
    if not df_spec_list:
        print("  No spectroscopic gravity data found.")
        return None
        
    df_spec = pd.concat(df_spec_list, ignore_index=True)
    
    # 3. Join
    print("  Joining datasets...")
    df = pd.merge(df_sample, df_spec, on='apogee_id', how='inner')
    
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

def analyze_gravity(df):
    print("Analyzing Gravity Schism...")
    
    df = compute_rgc(df)
    
    # Calculate Residual
    # Delta log g = Spec - Iso
    df['delta_logg'] = df['logg_spec'] - df['logg_iso']
    
    # Filter extreme outliers
    df_clean = df[(df['delta_logg'] > -1.0) & (df['delta_logg'] < 1.0)].copy()
    
    print(f"  Sample size: {len(df_clean)}")
    print(f"  Mean Delta logg: {df_clean['delta_logg'].mean():.4f} dex")
    
    # Bin by R_gc
    bins = [0, 4, 6, 8, 10, 15, 30]
    df_clean['r_bin'] = pd.cut(df_clean['R_gc'], bins=bins)
    
    binned = df_clean.groupby('r_bin')['delta_logg'].agg(['mean', 'sem', 'count'])
    binned['r_center'] = [i.mid for i in binned.index]
    
    print("\nDelta logg by R_gc:")
    print(binned[['mean', 'sem', 'count']])
    
    # Fit
    valid = binned[binned['count'] > 10]
    if len(valid) > 2:
        slope, intercept, r_val, p_val, std_err = stats.linregress(valid['r_center'], valid['mean'])
        print(f"  Slope (Delta logg vs R): {slope:.5f} dex/kpc")
        print(f"  Correlation: {r_val:.3f}")
    else:
        slope, r_val = 0, 0
        
    # Plot
    fig, ax = plt.subplots(figsize=(8, 6))
    
    ax.errorbar(binned['r_center'], binned['mean'], yerr=binned['sem'], fmt='o-', capsize=5)
    ax.set_xlabel('Galactocentric Radius [kpc]')
    ax.set_ylabel('Delta log g (Spec - Iso) [dex]')
    ax.set_title('Test CR: Gravity Schism')
    ax.axhline(0, color='k', linestyle='--', alpha=0.3)
    ax.grid(True, alpha=0.3)
    
    out_path = os.path.join(FIGURES_DIR, 'sdss_test_cr_gravity.png')
    plt.savefig(out_path, dpi=150)
    print(f"Figure saved to {out_path}")
    
    return {
        'slope': slope,
        'r_val': r_val,
        'mean_offset': df_clean['delta_logg'].mean(),
        'n_sample': int(len(df_clean))
    }

def main():
    cache_path = os.path.join(DATA_DIR, 'sdss_gravity.csv')
    
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

    results = analyze_gravity(df)
    
    out_path = os.path.join(RESULTS_DIR, 'sdss_test_cr_results.json')
    if results:
        def default(o):
            if isinstance(o, pd.Interval): return str(o)
            raise TypeError
            
        with open(out_path, 'w') as f:
            json.dump(results, f, indent=2, default=default)
            
        print("\nSUMMARY TEST CR:")
        print("Prediction: Delta log g (Spec - Iso) correlates with R_gc.")
        print(f"Observed Slope: {results['slope']:.5f}")
        
        if abs(results['slope']) > 0.002:
             print("RESULT: SIGNAL (Gradient observed)")
        else:
             print("RESULT: NULL (Flat profile)")

if __name__ == "__main__":
    main()
