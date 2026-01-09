#!/usr/bin/env python3
"""
Step 6.67: SDSS Test CP - The Aluminum Clock (Deep Mixing)

Hypothesis:
Aluminum is affected by the Mg-Al cycle in the deep interiors of Red Giants (deep mixing). 
The mixing rate depends on the thermal structure and rotation, which are potential-dependent (time dilated). 
We expect the surface abundance of [Al/Fe] or [Al/Mg] in giants to show residuals correlated with R_gc relative to standard chemical evolution.

Prediction:
[Al/Mg] (mixing indicator) shows an anomalous gradient with R_gc.
Standard expectation: [Al/Mg] might be constant or follow metallicity.
TEP: Deep mixing is slower in deep potential -> Less Al enhancement? Or less depletion? 
Mg-Al cycle converts Mg to Al. So high mixing -> High [Al/Mg].
TEP (slower mixing) -> Lower [Al/Mg] in Inner Galaxy?

Data:
- aspcapStar: al_fe, mg_fe, param_logg
- apogee_starhorse: dist50

Method:
1. Select Upper Red Giant Branch stars (logg < 2.0).
2. Compute R_gc.
3. Compute [Al/Mg] = [Al/Fe] - [Mg/Fe].
4. Analyze Mean [Al/Mg] vs R_gc.
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
    print(f"Querying SDSS for Test CP (Limit: {limit})...")
    
    # 1. Fetch Sample from apogee_starhorse (Giants)
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
        AND logg50 < 2.0 -- Upper RGB
    """
    df_sample = query_sdss(sql_sample)
    
    if df_sample is None or len(df_sample) == 0:
        print("  No sample data found.")
        return None
        
    print(f"  Got {len(df_sample)} stars. Fetching chemistry...")
    ids = df_sample['apogee_id'].astype(str).tolist()
    ids = list(set(ids))
    
    # 2. Fetch aspcapStar (Chemistry)
    # al_fe, mg_fe
    chunk_size = 50
    df_chem_list = []
    
    for i in range(0, len(ids), chunk_size):
        chunk = ids[i:i+chunk_size]
        ids_str = "', '".join(chunk)
        
        sql_chem = f"""
        SELECT 
            apogee_id,
            al_fe,
            mg_fe
        FROM aspcapStar
        WHERE 
            apogee_id IN ('{ids_str}')
            AND al_fe > -5 AND mg_fe > -5
        """
        res = query_sdss(sql_chem)
        if res is not None and len(res) > 0:
            df_chem_list.append(res)
        time.sleep(0.2)
            
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

def analyze_aluminum(df):
    print("Analyzing Aluminum Clock...")
    
    df = compute_rgc(df)
    
    # Compute [Al/Mg]
    df['al_mg'] = df['al_fe'] - df['mg_fe']
    
    print(f"  Sample size: {len(df)}")
    
    # Bin by R_gc
    bins = [0, 4, 6, 8, 10, 15, 30]
    df['r_bin'] = pd.cut(df['R_gc'], bins=bins)
    
    binned = df.groupby('r_bin')['al_mg'].agg(['mean', 'sem', 'count'])
    binned['r_center'] = [i.mid for i in binned.index]
    
    print("\nMean [Al/Mg] by R_gc:")
    print(binned[['mean', 'sem', 'count']])
    
    # Fit
    valid = binned[binned['count'] > 10]
    if len(valid) > 2:
        slope, intercept, r_val, p_val, std_err = stats.linregress(valid['r_center'], valid['mean'])
        print(f"  Slope ([Al/Mg] vs R): {slope:.5f} dex/kpc")
        print(f"  Correlation: {r_val:.3f}")
    else:
        slope, r_val = 0, 0
        
    # Plot
    fig, ax = plt.subplots(figsize=(8, 6))
    
    ax.errorbar(binned['r_center'], binned['mean'], yerr=binned['sem'], fmt='o-', capsize=5)
    ax.set_xlabel('Galactocentric Radius [kpc]')
    ax.set_ylabel('[Al/Mg] (Mixing Indicator)')
    ax.set_title('Test CP: Aluminum Clock (Deep Mixing)')
    ax.grid(True, alpha=0.3)
    
    out_path = os.path.join(FIGURES_DIR, 'sdss_test_cp_aluminum.png')
    plt.savefig(out_path, dpi=150)
    print(f"Figure saved to {out_path}")
    
    return {
        'slope': slope,
        'r_val': r_val,
        'n_sample': int(len(df))
    }

def main():
    cache_path = os.path.join(DATA_DIR, 'sdss_aluminum.csv')
    
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

    results = analyze_aluminum(df)
    
    out_path = os.path.join(RESULTS_DIR, 'sdss_test_cp_results.json')
    if results:
        def default(o):
            if isinstance(o, pd.Interval): return str(o)
            raise TypeError
            
        with open(out_path, 'w') as f:
            json.dump(results, f, indent=2, default=default)
            
        print("\nSUMMARY TEST CP:")
        print("Prediction: [Al/Mg] varies with R_gc (Mixing efficiency).")
        print(f"Observed Slope: {results['slope']:.5f}")
        
        if abs(results['slope']) > 0.005:
             print("RESULT: SIGNAL (Gradient observed)")
        else:
             print("RESULT: NULL (Flat profile)")

if __name__ == "__main__":
    main()
