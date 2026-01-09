#!/usr/bin/env python3
"""
Step 6.25: SDSS Test AG - The Distance Schism (Spectro vs Geometric)

Hypothesis:
Spectrophotometric distance (d_L) relies on observed flux F = L/4pi d^2.
In TEP, time dilation reduces observed power (P_obs < P_emit) in deep potentials (Bulge).
Stars look dimmer -> inferred d_L is too large.
Geometric distance (d_G, parallax) is unaffected.

Prediction:
Residual Delta_d = (d_spectro - d_geo) is POSITIVE in deep potentials.
r(Delta_d/d, R_GC) < 0 (Excess distance decreases as we move out from Bulge).

Data:
- apogee_starhorse: dist50 (Spectro-photometric / Bayesian distance).
- apogeeStar: gaiaedr3_r_med_geo (Geometric distance from Gaia EDR3 parallax).
- Coordinates: glon, glat.

Analysis:
1. Join tables.
2. Filter for good parallax quality.
3. Compute R_GC using geometric distance.
4. Analyze Delta_ln_d vs R_GC.
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
    print(f"Querying SDSS for Test AG (Limit: {limit})...")
    
    # Use apogeeStar (a) and apogee_starhorse (s)
    # Filter for low extinction if possible? Or trust StarHorse.
    # We focus on Giants usually, but StarHorse handles all.
    
    sql = f"""
    SELECT TOP {limit}
        a.apogee_id,
        a.glon, a.glat,
        s.dist50 as d_spectro,
        a.gaiaedr3_r_med_geo as d_geo,
        a.gaiaedr3_parallax,
        a.gaiaedr3_parallax_error
        
    FROM apogeeStar a
    JOIN apogee_starhorse s ON a.apogee_id = s.apogee_id
    
    WHERE 
        a.gaiaedr3_r_med_geo > 0 
        AND s.dist50 > 0
        AND a.gaiaedr3_parallax > 0.1 -- d < 10 kpc approx
        AND (a.gaiaedr3_parallax_error / a.gaiaedr3_parallax) < 0.2 -- 20% precision
    """
    return query_sdss(sql)

def compute_rgc(d, l_deg, b_deg):
    # Galactic Center distance R0 = 8.2 kpc
    R0 = 8.2
    
    l = np.radians(l_deg)
    b = np.radians(b_deg)
    
    # Convert to Galactocentric Cylindrical/Spherical
    # Sun at (-R0, 0, 0) in standard frame? Or (R0, 0, 0)?
    # Standard convention: Sun is at X = -R0 (or +R0 depending on frame).
    # Distance from GC:
    # d_proj = d * cos(b)
    # X = R0 - d_proj * cos(l)  (Assuming Sun at R0, GC at 0)
    # Y = - d_proj * sin(l)
    # Z = d * sin(b)
    # R_GC = sqrt(X^2 + Y^2 + Z^2)
    # Using law of cosines on the plane:
    # R_plane^2 = R0^2 + (d cos b)^2 - 2 R0 d cos b cos l
    
    d_proj = d * np.cos(b)
    R_plane_sq = R0**2 + d_proj**2 - 2 * R0 * d_proj * np.cos(l)
    Z = d * np.sin(b)
    
    R_GC = np.sqrt(R_plane_sq + Z**2)
    return R_GC

def analyze_distance_schism(df):
    print("Analyzing Distance Schism...")
    
    # 1. Clean
    df_clean = df.dropna().copy()
    # Convert units: dist50 is kpc? gaiaedr3 is pc?
    # APOGEE StarHorse dist50 is usually in kpc.
    # Gaia geometric dist is usually in pc.
    # Let's verify units by ratio.
    
    # Quick check
    ratio_check = df_clean['d_spectro'].median() / df_clean['d_geo'].median()
    print(f"Median Ratio (Spectro/Geo) Raw: {ratio_check:.4f}")
    
    if ratio_check < 0.01:
        print("Detected unit mismatch: Spectro (kpc) vs Geo (pc). Scaling Geo to kpc.")
        df_clean['d_geo_kpc'] = df_clean['d_geo'] / 1000.0
        df_clean['d_spectro_kpc'] = df_clean['d_spectro']
    elif ratio_check > 100:
        print("Detected unit mismatch: Spectro (pc) vs Geo (kpc).")
        df_clean['d_geo_kpc'] = df_clean['d_geo']
        df_clean['d_spectro_kpc'] = df_clean['d_spectro'] / 1000.0
    else:
        print("Units appear consistent (likely both kpc or both pc). Assuming kpc based on standard APOGEE.")
        # Actually Gaia r_med_geo is in pc in the table usually?
        # StarHorse is kpc.
        # If ratio ~ 0.001, Spectro is kpc, Geo is pc.
        if df_clean['d_geo'].median() > 500: # likely pc
             df_clean['d_geo_kpc'] = df_clean['d_geo'] / 1000.0
             df_clean['d_spectro_kpc'] = df_clean['d_spectro'] # StarHorse is kpc
        else:
             df_clean['d_geo_kpc'] = df_clean['d_geo']
             df_clean['d_spectro_kpc'] = df_clean['d_spectro']
             
    # 2. Compute Log Ratio
    # Delta = ln(d_spectro) - ln(d_geo)
    # Positive means Spectro > Geo (Overestimated distance / Underestimated Flux)
    df_clean['delta_ln_d'] = np.log(df_clean['d_spectro_kpc']) - np.log(df_clean['d_geo_kpc'])
    
    # 3. Compute R_GC
    # Use Geometric distance for "True" position
    df_clean['R_GC'] = compute_rgc(df_clean['d_geo_kpc'], df_clean['glon'], df_clean['glat'])
    
    print(f"N = {len(df_clean)}")
    print(f"Mean Delta ln(d): {df_clean['delta_ln_d'].mean():.4f}")
    
    # 4. Correlation with R_GC
    # TEP: Delta decreases as R_GC increases (High in Bulge/Inner Galaxy)
    # r < 0
    
    r_rgc, p_rgc = stats.pearsonr(df_clean['R_GC'], df_clean['delta_ln_d'])
    print(f"Correlation r(Delta, R_GC): {r_rgc:.4f} (p={p_rgc:.2e})")
    
    # 5. Binning
    df_clean['rgc_bin'] = pd.qcut(df_clean['R_GC'], 8)
    binned = df_clean.groupby('rgc_bin')['delta_ln_d'].mean()
    print("\nMean Delta ln(d) by R_GC:")
    print(binned)
    
    return {
        'r_rgc': float(r_rgc),
        'p_rgc': float(p_rgc),
        'mean_delta': float(df_clean['delta_ln_d'].mean()),
        'binned_means': binned.tolist(),
        'bin_centers': [mid.mid for mid in binned.index]
    }, df_clean

def create_figure(df, results):
    print("Generating figure...")
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # Scatter
    ax.scatter(df['R_GC'], df['delta_ln_d'], alpha=0.1, s=1, c='k', label='Stars')
    
    # Binned
    ax.plot(results['bin_centers'], results['binned_means'], 'r-o', lw=2, label='Mean')
    
    ax.set_xlabel(r'Galactocentric Radius $R_{GC}$ [kpc]')
    ax.set_ylabel(r'Distance Schism $\ln(d_{spectro}/d_{geo})$')
    ax.set_title(f"Test AG: Distance Anomaly (r={results['r_rgc']:.3f})")
    ax.axhline(0, linestyle='--', color='b', alpha=0.5)
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    out_path = os.path.join(FIGURES_DIR, 'sdss_test_ag_distance_schism.png')
    plt.savefig(out_path, dpi=150)
    print(f"Figure saved to {out_path}")

def main():
    cache_path = os.path.join(DATA_DIR, 'sdss_apogee_distances.csv')
    
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

    results, df_clean = analyze_distance_schism(df)
    
    out_path = os.path.join(RESULTS_DIR, 'sdss_test_ag_results.json')
    with open(out_path, 'w') as f:
        results_json = results.copy()
        del results_json['bin_centers']
        json.dump(results_json, f, indent=2)
        
    create_figure(df_clean, results)
    
    print("\nSUMMARY TEST AG:")
    print("TEP Prediction: Delta > 0 in Inner Galaxy. r(Delta, R_GC) < 0.")
    print(f"Observed r: {results['r_rgc']:.4f}")
    
    if results['r_rgc'] < -0.05:
        print("RESULT: CONSISTENT (Spectro distances inflated in deep potential)")
    elif results['r_rgc'] > 0.05:
        print("RESULT: CONTRADICTED (Spectro distances suppressed in deep potential)")
    else:
        print("RESULT: NULL/FLAT")

if __name__ == "__main__":
    main()
