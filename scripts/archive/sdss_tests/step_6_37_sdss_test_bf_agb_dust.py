#!/usr/bin/env python3
"""
Step 6.37: SDSS Test BF - AGB Dust Production (Wind Clock)

Hypothesis:
AGB mass loss (M_dot) is driven by pulsations and radiation pressure on dust. These are rate processes.
In deep potentials (Inner Galaxy), M_dot should be suppressed (slower pulsations).
We expect lower IR excess (W3/W4 bands) for AGB stars in the bulge compared to the disk at fixed luminosity/temperature.

Prediction:
W3-W4 color (dust excess) is lower at small R_gc.

Data:
- mos_gaia_dr2_source: source_id, phot_g_mean_mag, bp_rp, l, b.
- mos_allwise: w3mpro, w4mpro.
- mos_geometric_distances_gaia_dr2: r_est (Distance).

Method:
1. Select AGB candidates:
   - Red: bp_rp > 1.5.
   - Luminous: Absolute G < 0 (approx). M_G = G - 5 log(d) + 5.
   - Quality: Good distance, valid W3/W4.
2. Compute R_gc.
3. Compute Dust Excess proxy: W3 - W4 (or W3 alone if K-W3 used, but we stick to WISE colors).
   Actually W3-W4 is a good proxy for cold dust/mass loss rate.
4. Bin by R_gc.
5. Correlate (W3-W4) with R_gc.
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

def download_data(limit=200):
    print(f"Querying SDSS for Test BF (Limit: {limit})...")
    
    # Use parallax from gaia source directly to estimate distance.
    # d ~ 1000 / parallax (mas) [pc].
    # Avoid triple join.
    
    sql = f"""
    SELECT TOP {limit}
        g.source_id,
        g.phot_g_mean_mag as g_mag,
        g.bp_rp,
        g.l, g.b,
        g.parallax,
        w.w3mpro, w.w4mpro
        
    FROM mos_gaia_dr2_source g
    JOIN mos_allwise w ON g.source_id = w.source_id
    
    WHERE 
        g.bp_rp > 1.5 -- Red
        AND g.parallax > 0.1 -- d < 10 kpc approx
        AND w.w3mpro > 0 AND w.w4mpro > 0
    """
    return query_sdss(sql)

def analyze_agb_dust(df):
    print("Analyzing AGB Dust Production...")
    
    # 1. Clean
    df_clean = df.dropna().copy()
    
    # 2. Compute Distance and Absolute Mag
    # Parallax is usually in mas? mos_gaia_dr2_source parallax unit?
    # Usually mas.
    df_clean['dist_pc'] = 1000.0 / df_clean['parallax']
    df_clean['dist_kpc'] = df_clean['dist_pc'] / 1000.0
    
    df_clean['abs_g'] = df_clean['g_mag'] - 5 * np.log10(df_clean['dist_pc']) + 5
    
    # Select Giants (Luminous)
    # AGBs typically M_G < 0 or < -1 depending on color.
    # Let's select M_G < 0.
    df_clean = df_clean[df_clean['abs_g'] < 0].copy()
    
    # 3. Compute R_gc
    df_clean['R_gc'] = compute_rgc(df_clean['dist_kpc'], df_clean['l'], df_clean['b'])
    
    # 4. Compute Dust Excess (W3 - W4)
    # Larger value = More excess (W4 brighter than W3 relative to Rayleigh-Jeans?)
    # AGBs with dust shells have red W3-W4.
    df_clean['w3_w4'] = df_clean['w3mpro'] - df_clean['w4mpro']
    
    # 5. Correlation
    # TEP: Lower excess (smaller W3-W4) at small R_gc.
    # Positive correlation: W3-W4 increases with R_gc.
    
    r_dust, p_dust = stats.pearsonr(df_clean['R_gc'], df_clean['w3_w4'])
    
    print(f"N = {len(df_clean)}")
    print(f"Correlation r(W3-W4, R_gc): {r_dust:.4f} (p={p_dust:.2e})")
    
    # 6. Binning
    df_clean['rgc_bin'] = pd.qcut(df_clean['R_gc'], 8)
    binned = df_clean.groupby('rgc_bin')['w3_w4'].mean()
    print("\nMean W3-W4 by R_gc:")
    print(binned)
    
    return {
        'r_dust': float(r_dust),
        'p_dust': float(p_dust),
        'binned_means': binned.tolist(),
        'bin_centers': [mid.mid for mid in binned.index],
        'n_sample': int(len(df_clean))
    }, df_clean

def create_figure(df, results):
    print("Generating figure...")
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # Scatter
    ax.scatter(df['R_gc'], df['w3_w4'], alpha=0.1, s=2, c='k', label='AGB Stars')
    
    # Binned
    ax.plot(results['bin_centers'], results['binned_means'], 'r-o', lw=2, label='Mean Color')
    
    ax.set_xlabel(r'Galactocentric Radius $R_{GC}$ [kpc]')
    ax.set_ylabel(r'Dust Excess $W3 - W4$ [mag]')
    ax.set_title(f"Test BF: AGB Dust vs Radius (r={results['r_dust']:.3f})")
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    out_path = os.path.join(FIGURES_DIR, 'sdss_test_bf_agb_dust.png')
    plt.savefig(out_path, dpi=150)
    print(f"Figure saved to {out_path}")

def main():
    cache_path = os.path.join(DATA_DIR, 'sdss_agb_dust.csv')
    
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

    results, df_clean = analyze_agb_dust(df)
    
    out_path = os.path.join(RESULTS_DIR, 'sdss_test_bf_results.json')
    with open(out_path, 'w') as f:
        results_json = results.copy()
        del results_json['bin_centers']
        json.dump(results_json, f, indent=2)
        
    create_figure(df_clean, results)
    
    print("\nSUMMARY TEST BF:")
    print("TEP Prediction: Dust excess (W3-W4) lower in inner galaxy (Small R_gc). r > 0.")
    print(f"Observed r: {results['r_dust']:.4f}")
    
    if results['r_dust'] > 0.05:
        print("RESULT: CONSISTENT (Dust suppressed in deep potential)")
    elif results['r_dust'] < -0.05:
        print("RESULT: CONTRADICTED (Dust enhanced in deep potential)")
    else:
        print("RESULT: NULL/FLAT")

if __name__ == "__main__":
    main()
