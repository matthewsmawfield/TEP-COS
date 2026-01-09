#!/usr/bin/env python3
"""
Step 6.32: SDSS Test AZ - HB/RGB Ratio (Helium Clock)

Hypothesis:
The ratio of Horizontal Branch (core He burning) to Red Giant Branch (H-shell burning) stars depends on the relative lifetimes of these phases.
If nuclear reaction rates differ due to time dilation in the deep potential of the inner Galaxy, this ratio should vary systematically with Galactocentric radius.
Standard theory: Ratio depends on Helium abundance (Y) and Age, but for disk populations these are relatively uniform or have known gradients.
TEP Prediction: Ratio N(HB) / N(RGB) varies with R_gc (at fixed Metallicity).

Data:
- aspcapStar: teff, logg, m_h (Metallicity).
- apogeeStar: gaiaedr3_r_med_geo (Distance), glon, glat.

Method:
1. Select Disk/Bulge stars (Metal rich > -0.5).
2. Classify:
   - HB (Red Clump): 2.3 < logg < 3.0 AND 4500 < teff < 5200.
   - RGB: 0.0 < logg < 3.0 AND teff < 4500.
3. Compute R_gc using Distance and Coordinates.
4. Bin by R_gc.
5. Compute Ratio = N_HB / N_RGB in each bin.
6. Correlate Ratio with R_gc.
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
    print(f"Querying SDSS for Test AZ (Limit: {limit})...")
    
    # Use apogeeStar only. It has rv_teff, rv_logg, rv_feh and coordinates.
    # This avoids the join with aspcapStar.
    
    sql = f"""
    SELECT TOP {limit}
        apogee_id,
        rv_teff as teff,
        rv_logg as logg,
        rv_feh as fe_h,
        glon, glat,
        gaiaedr3_r_med_geo as dist_pc
        
    FROM apogeeStar
    
    WHERE 
        rv_feh > -0.5 -- Disk/Bulge
        AND rv_logg < 3.5 -- Giants
        AND gaiaedr3_r_med_geo > 0
    """
    return query_sdss(sql)

def compute_rgc(d_kpc, l_deg, b_deg):
    R0 = 8.2 # kpc
    l = np.radians(l_deg)
    b = np.radians(b_deg)
    
    d_proj = d_kpc * np.cos(b)
    # X towards GC? Standard: Sun at -R0. GC at 0.
    # X = R0 - d cos l
    # R^2 = R0^2 + d^2 - 2 R0 d cos l
    R_plane_sq = R0**2 + d_proj**2 - 2 * R0 * d_proj * np.cos(l)
    Z = d_kpc * np.sin(b)
    R_gc = np.sqrt(R_plane_sq + Z**2)
    return R_gc

def analyze_hb_rgb_ratio(df):
    print("Analyzing HB/RGB Ratio...")
    
    # 1. Clean
    df_clean = df.dropna().copy()
    
    # 2. Compute R_gc
    df_clean['dist_kpc'] = df_clean['dist_pc'] / 1000.0
    df_clean['R_gc'] = compute_rgc(df_clean['dist_kpc'], df_clean['glon'], df_clean['glat'])
    
    # 3. Classify
    # HB (Red Clump): 2.3 < logg < 3.0 AND 4500 < teff < 5200 (Rough cut)
    # RGB: 0.0 < logg < 3.0 AND teff < 4500 (Cooler)
    
    def classify(row):
        if (2.3 < row['logg'] < 3.0) and (4500 < row['teff'] < 5200):
            return 'HB'
        elif (0.0 < row['logg'] < 3.0) and (row['teff'] < 4500):
            return 'RGB'
        else:
            return 'Other'
            
    df_clean['type'] = df_clean.apply(classify, axis=1)
    
    # Filter only HB and RGB
    df_clean = df_clean[df_clean['type'] != 'Other'].copy()
    
    # 4. Bin by R_gc
    # We want Ratio in bins.
    bins = np.linspace(0, 15, 16) # 0 to 15 kpc
    
    results_list = []
    
    for i in range(len(bins)-1):
        r_min, r_max = bins[i], bins[i+1]
        sub = df_clean[(df_clean['R_gc'] >= r_min) & (df_clean['R_gc'] < r_max)]
        
        n_hb = len(sub[sub['type'] == 'HB'])
        n_rgb = len(sub[sub['type'] == 'RGB'])
        
        if n_rgb > 10:
            ratio = n_hb / n_rgb
            # Error prop: sigma_R = R * sqrt(1/N_HB + 1/N_RGB)
            err = ratio * np.sqrt(1/n_hb + 1/n_rgb) if n_hb > 0 else 0
            
            results_list.append({
                'R_gc_bin': (r_min + r_max)/2,
                'ratio': ratio,
                'error': err,
                'n_hb': n_hb,
                'n_rgb': n_rgb
            })
            
    res_df = pd.DataFrame(results_list)
    print("\nHB/RGB Ratio by R_gc:")
    print(res_df)
    
    # 5. Correlation
    r_ratio, p_ratio = stats.pearsonr(res_df['R_gc_bin'], res_df['ratio'])
    print(f"Correlation r(Ratio, R_gc): {r_ratio:.4f} (p={p_ratio:.2e})")
    
    return {
        'r_ratio': float(r_ratio),
        'p_ratio': float(p_ratio),
        'binned_data': res_df.to_dict(orient='records'),
        'n_sample': int(len(df_clean))
    }, df_clean

def create_figure(res_df, results):
    print("Generating figure...")
    df = pd.DataFrame(results['binned_data'])
    
    fig, ax = plt.subplots(figsize=(8, 6))
    
    ax.errorbar(df['R_gc_bin'], df['ratio'], yerr=df['error'], fmt='o-', capsize=3, label='Observed Ratio')
    
    ax.set_xlabel(r'Galactocentric Radius $R_{GC}$ [kpc]')
    ax.set_ylabel(r'Ratio $N_{HB} / N_{RGB}$')
    ax.set_title(f"Test AZ: HB/RGB Ratio vs Radius (r={results['r_ratio']:.3f})")
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    out_path = os.path.join(FIGURES_DIR, 'sdss_test_az_hb_rgb_ratio.png')
    plt.savefig(out_path, dpi=150)
    print(f"Figure saved to {out_path}")

def main():
    cache_path = os.path.join(DATA_DIR, 'sdss_hb_rgb_ratio.csv')
    
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

    results, df_clean = analyze_hb_rgb_ratio(df)
    
    out_path = os.path.join(RESULTS_DIR, 'sdss_test_az_results.json')
    with open(out_path, 'w') as f:
        json.dump(results, f, indent=2)
        
    create_figure(None, results)
    
    print("\nSUMMARY TEST AZ:")
    print("TEP Prediction: Ratio varies with R_gc (Potential depth).")
    print(f"Observed r: {results['r_ratio']:.4f}")
    
    if abs(results['r_ratio']) > 0.3:
        print("RESULT: CONSISTENT (Significant variation with radius)")
    else:
        print("RESULT: NULL/FLAT (Ratio is uniform)")

if __name__ == "__main__":
    main()
