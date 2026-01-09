#!/usr/bin/env python3
"""
Step 6.36: SDSS Test BE - MDF Skewness (Enrichment Rates)

Hypothesis:
The Metallicity Distribution Function (MDF) skewness records the history of infall vs outflow rates.
If these rates are dilated in the inner galaxy relative to standard chemical evolution models, the MDF shape should show residuals in skewness/kurtosis vs R_gc that standard models (which assume universal rate constants) do not predict.
TEP: Rate dilation in deep potential might delay enrichment or alter the balance of infall/outflow.
Prediction: MDF Skewness varies with R_gc in a way not explained by simple inside-out formation (which predicts mean metallicity gradients but stable MDF shapes).

Data:
- apogeeStar: rv_feh (Metallicity), glon, glat, gaiaedr3_r_med_geo (Distance).

Method:
1. Select Disk stars (low alpha? or just all disk).
2. Compute R_gc.
3. Bin by R_gc.
4. Compute Skewness of [Fe/H] distribution in each bin.
   Standard Disk MDF is negatively skewed (tail towards low metallicity).
   Bulge MDF?
   We look for a gradient in skewness.
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
    print(f"Querying SDSS for Test BE (Limit: {limit})...")
    
    sql = f"""
    SELECT TOP {limit}
        apogee_id,
        rv_feh as fe_h,
        glon, glat,
        gaiaedr3_r_med_geo as dist_pc
        
    FROM apogeeStar
    
    WHERE 
        rv_feh > -2.5 AND rv_feh < 1.0
        AND gaiaedr3_r_med_geo > 0
        AND glat > -20 AND glat < 20 -- Disk focus
    """
    return query_sdss(sql)

def compute_rgc(d_kpc, l_deg, b_deg):
    R0 = 8.2 # kpc
    l = np.radians(l_deg)
    b = np.radians(b_deg)
    
    d_proj = d_kpc * np.cos(b)
    # R^2 = R0^2 + d^2 - 2 R0 d cos l
    R_plane_sq = R0**2 + d_proj**2 - 2 * R0 * d_proj * np.cos(l)
    Z = d_kpc * np.sin(b)
    R_gc = np.sqrt(R_plane_sq + Z**2)
    return R_gc

def analyze_mdf_skewness(df):
    print("Analyzing MDF Skewness...")
    
    # 1. Clean
    df_clean = df.dropna().copy()
    
    # 2. Compute R_gc
    df_clean['dist_kpc'] = df_clean['dist_pc'] / 1000.0
    df_clean['R_gc'] = compute_rgc(df_clean['dist_kpc'], df_clean['glon'], df_clean['glat'])
    
    # 3. Bin by R_gc
    bins = np.linspace(0, 15, 16) # 0 to 15 kpc
    
    results_list = []
    
    for i in range(len(bins)-1):
        r_min, r_max = bins[i], bins[i+1]
        sub = df_clean[(df_clean['R_gc'] >= r_min) & (df_clean['R_gc'] < r_max)]
        
        if len(sub) > 50:
            skew = stats.skew(sub['fe_h'])
            # Skewness error approx sqrt(6/N)
            err = np.sqrt(6/len(sub))
            mean_feh = sub['fe_h'].mean()
            
            results_list.append({
                'R_gc_bin': (r_min + r_max)/2,
                'skewness': skew,
                'error': err,
                'mean_feh': mean_feh,
                'n_stars': len(sub)
            })
            
    res_df = pd.DataFrame(results_list)
    print("\nMDF Skewness by R_gc:")
    print(res_df)
    
    # 4. Correlation
    r_skew, p_skew = stats.pearsonr(res_df['R_gc_bin'], res_df['skewness'])
    print(f"Correlation r(Skew, R_gc): {r_skew:.4f} (p={p_skew:.2e})")
    
    return {
        'r_skew': float(r_skew),
        'p_skew': float(p_skew),
        'binned_data': res_df.to_dict(orient='records'),
        'n_sample': int(len(df_clean))
    }, df_clean

def create_figure(df, results):
    print("Generating figure...")
    data = pd.DataFrame(results['binned_data'])
    
    fig, ax = plt.subplots(figsize=(8, 6))
    
    ax.errorbar(data['R_gc_bin'], data['skewness'], yerr=data['error'], fmt='o-', capsize=3, label='MDF Skewness')
    
    ax.set_xlabel(r'Galactocentric Radius $R_{GC}$ [kpc]')
    ax.set_ylabel(r'MDF Skewness')
    ax.set_title(f"Test BE: MDF Skewness vs Radius (r={results['r_skew']:.3f})")
    ax.axhline(0, linestyle='--', color='k', alpha=0.3)
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    out_path = os.path.join(FIGURES_DIR, 'sdss_test_be_mdf_skewness.png')
    plt.savefig(out_path, dpi=150)
    print(f"Figure saved to {out_path}")

def main():
    cache_path = os.path.join(DATA_DIR, 'sdss_mdf_skewness.csv')
    
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

    results, df_clean = analyze_mdf_skewness(df)
    
    out_path = os.path.join(RESULTS_DIR, 'sdss_test_be_results.json')
    with open(out_path, 'w') as f:
        json.dump(results, f, indent=2)
        
    create_figure(df_clean, results)
    
    print("\nSUMMARY TEST BE:")
    print("TEP Prediction: Skewness varies with R_gc (Potential depth).")
    print(f"Observed r: {results['r_skew']:.4f}")
    
    if abs(results['r_skew']) > 0.3:
        print("RESULT: CONSISTENT (MDF shape depends on potential)")
    else:
        print("RESULT: NULL/FLAT (MDF shape is universal/self-similar)")

if __name__ == "__main__":
    main()
