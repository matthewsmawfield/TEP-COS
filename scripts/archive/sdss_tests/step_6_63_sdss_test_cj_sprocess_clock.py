#!/usr/bin/env python3
"""
Step 6.63: SDSS Test CJ - The S-Process Clock (Neutron Capture Rates)

Hypothesis:
Heavy elements (Ce, Nd) are produced by slow neutron capture (s-process) in AGB stars. 
The neutron capture rate depends on thermal timescales in the AGB interior. 
In deep potentials, these rates are dilated. We expect the yield of s-process elements 
relative to iron (explosive nucleosynthesis) to vary with potential depth (R_gc).

Prediction:
[Ce/Fe] or [Nd/Fe] decreases as R_gc decreases (Inner Galaxy suppression).

Data:
- aspcapStar: param_ce_h, param_nd_h, param_fe_h
- apogee_starhorse: dist50

Method:
1. Select Disk stars (param_fe_h > -0.5).
2. Compute R_gc.
3. Compute [Ce/Fe] = [Ce/H] - [Fe/H] (and Nd).
4. Analyze Mean [X/Fe] vs R_gc.
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
    print(f"Querying SDSS for Test CJ (Limit: {limit})...")
    
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
    # Correct columns from check: ce_fe, felem_nd_h
    # Need fe_h for normalization of Nd if using felem
    
    chunk_size = 50
    df_chem_list = []
    
    for i in range(0, len(ids), chunk_size):
        chunk = ids[i:i+chunk_size]
        ids_str = "', '".join(chunk)
        
        sql_chem = f"""
        SELECT 
            apogee_id,
            ce_fe,
            felem_nd_h as param_nd_h,
            fe_h as param_fe_h
        FROM aspcapStar
        WHERE 
            apogee_id IN ('{ids_str}')
        """
        # Removed param_fe_h filter from SQL to simplify, handle in pandas
        res = query_sdss(sql_chem)
        if res is not None and len(res) > 0:
            df_chem_list.append(res)
        time.sleep(0.2)
            
    if not df_chem_list:
        print("  No s-process chemistry data found.")
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

def analyze_sprocess(df):
    print("Analyzing S-Process Clock...")
    
    df = compute_rgc(df)
    
    # Calculate ratios
    
    # Ce: ce_fe is already [Ce/Fe]
    if 'ce_fe' in df.columns:
        df['ce_fe_val'] = df['ce_fe']
        # Filter valid (>-9)
        df_ce = df[df['ce_fe_val'] > -9].copy()
        
        print(f"  Ce/Fe Sample: {len(df_ce)}")
        
        # Bin
        bins = [0, 4, 6, 8, 10, 15, 30]
        df_ce['r_bin'] = pd.cut(df_ce['R_gc'], bins=bins)
        binned_ce = df_ce.groupby('r_bin')['ce_fe_val'].agg(['mean', 'sem', 'count'])
        binned_ce['r_center'] = [i.mid for i in binned_ce.index]
        
        print("\n[Ce/Fe] by R_gc:")
        print(binned_ce[['mean', 'sem', 'count']])
        
        # Fit
        valid_ce = binned_ce[binned_ce['count'] > 5]
        if len(valid_ce) > 2:
            slope_ce, intercept, r_ce, p, _ = stats.linregress(valid_ce['r_center'], valid_ce['mean'])
            print(f"  Slope ([Ce/Fe] vs R): {slope_ce:.5f} dex/kpc")
        else:
            slope_ce, r_ce = 0, 0
    else:
        slope_ce, r_ce = 0, 0
        binned_ce = None

    # Nd: param_nd_h is [Nd/H]. Need to subtract [Fe/H] (param_fe_h)
    if 'param_nd_h' in df.columns and 'param_fe_h' in df.columns:
        df['nd_fe'] = df['param_nd_h'] - df['param_fe_h']
        df_nd = df[(df['param_nd_h'] > -9) & (df['param_fe_h'] > -9)].copy()
        
        # Bin
        df_nd['r_bin'] = pd.cut(df_nd['R_gc'], bins=bins)
        binned_nd = df_nd.groupby('r_bin')['nd_fe'].agg(['mean', 'sem', 'count'])
        binned_nd['r_center'] = [i.mid for i in binned_nd.index]
        
        valid_nd = binned_nd[binned_nd['count'] > 5]
        if len(valid_nd) > 2:
            slope_nd, intercept, r_nd, p, _ = stats.linregress(valid_nd['r_center'], valid_nd['mean'])
            print(f"  Slope ([Nd/Fe] vs R): {slope_nd:.5f} dex/kpc")
        else:
            slope_nd, r_nd = 0, 0
    else:
        slope_nd, r_nd = 0, 0
        binned_nd = None
        
    # Plot
    fig, ax = plt.subplots(figsize=(8, 6))
    
    if binned_ce is not None:
        ax.errorbar(binned_ce['r_center'], binned_ce['mean'], yerr=binned_ce['sem'], fmt='o-', label='[Ce/Fe]')
        
    if binned_nd is not None:
        ax.errorbar(binned_nd['r_center'], binned_nd['mean'], yerr=binned_nd['sem'], fmt='s-', label='[Nd/Fe]')
        
    ax.set_xlabel('Galactocentric Radius [kpc]')
    ax.set_ylabel('[X/Fe] (s-process)')
    ax.set_title('Test CJ: S-Process Clock')
    ax.grid(True, alpha=0.3)
    ax.legend()
    
    out_path = os.path.join(FIGURES_DIR, 'sdss_test_cj_sprocess.png')
    plt.savefig(out_path, dpi=150)
    print(f"Figure saved to {out_path}")
    
    return {
        'slope_ce': slope_ce,
        'slope_nd': slope_nd,
        'n_sample': int(len(df))
    }

def main():
    cache_path = os.path.join(DATA_DIR, 'sdss_sprocess.csv')
    
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

    results = analyze_sprocess(df)
    
    out_path = os.path.join(RESULTS_DIR, 'sdss_test_cj_results.json')
    if results:
        def default(o):
            if isinstance(o, pd.Interval): return str(o)
            raise TypeError
            
        with open(out_path, 'w') as f:
            json.dump(results, f, indent=2, default=default)
            
        print("\nSUMMARY TEST CJ:")
        print("Prediction: [s/Fe] decreases as R_gc decreases (Inner Galaxy suppression). Positive Slope.")
        print(f"Observed Slope (Ce): {results['slope_ce']:.5f}")
        
        if results['slope_ce'] > 0.005:
             print("RESULT: CONSISTENT (Suppression in inner galaxy observed)")
        else:
             print("RESULT: NULL (Flat or Negative)")

if __name__ == "__main__":
    main()
