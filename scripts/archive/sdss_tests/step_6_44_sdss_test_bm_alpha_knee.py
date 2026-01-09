#!/usr/bin/env python3
"""
Step 6.44: SDSS Test BM - The Alpha Knee (SN Ia Delay Clock)

Hypothesis:
The "knee" in the [alpha/Fe] vs [Fe/H] diagram marks the onset of Type Ia supernovae (time delay ~1 Gyr).
In TEP, proper time flows slower in the Inner Galaxy. 1 Gyr of proper time (stellar evolution) corresponds to a longer cosmic interval.
Star formation continues during this extended interval, allowing [Fe/H] to reach higher levels before the Ia enrichment kicks in.
Prediction: The [Fe/H] position of the Alpha-Knee increases as R_gc decreases.
(Note: Standard chemical evolution also predicts this due to higher SFR efficiencies in the bulge, so this is a degenerate test, but TEP predicts an *enhancement* of this trend).

Data:
- apogeeStar: param_alpha_m ([alpha/M]), param_m_h ([M/H] ~ [Fe/H]), glon, glat.
- apogee_starhorse: dist50.

Method:
1. Select Red Clump or Giant stars (low surface gravity, high S/N).
2. Bin by R_gc.
3. In each bin, determine the "knee" position.
   - Simplified approach: Find the [M/H] where [alpha/M] drops below a certain threshold (e.g., 0.1) from the high-alpha plateau.
   - Or: Fit a broken line or hyperbolic tangent to the sequence.
   - Robust proxy: The mean [M/H] of stars in the transition region (0.05 < [alpha/M] < 0.15).
4. Correlate Knee Position with R_gc.
"""

import pandas as pd
import numpy as np
from scipy import stats, optimize
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
    print(f"Querying SDSS for Test BM (Client-side Join - Attempt 3)...")
    
    # Step 1: Fetch Chemistry from aspcapStar
    # Use fparam_alpha_m and fparam_m_h
    # Removed snr filter as column is missing in aspcapStar (it's in apogeeStar)
    print("  Fetching Chemistry from aspcapStar...")
    sql_chem = f"""
    SELECT TOP {limit}
        apogee_id,
        fparam_alpha_m as alpha_m,
        fparam_m_h as fe_h
    FROM aspcapStar
    WHERE 
        fparam_alpha_m > -1 
        AND fparam_m_h > -2.5
    """
    df_chem = query_sdss(sql_chem)
    
    if df_chem is None or len(df_chem) == 0:
        print("  No chemistry data found.")
        return None
        
    # Get ID list for IN clause
    # Clean IDs to ensure they are safe strings
    ids = df_chem['apogee_id'].astype(str).tolist()
    
    # Chunk the IDs if too many
    # Reduced chunk size to 50 to avoid URL length limits (HTTP 404 on GET)
    chunk_size = 50
    df_pos_list = []
    
    print(f"  Got {len(df_chem)} stars. Fetching positions in chunks of {chunk_size}...")
    
    for i in range(0, len(ids), chunk_size):
        chunk = ids[i:i+chunk_size]
        ids_str = "', '".join(chunk)
        
        # Filter by SNR here in apogeeStar
        sql_pos = f"""
        SELECT 
            apogee_id,
            glon, glat,
            gaiaedr3_r_med_photogeo as dist50
        FROM apogeeStar
        WHERE 
            apogee_id IN ('{ids_str}')
            AND gaiaedr3_r_med_photogeo > 0
            AND snr > 50
        """
        res = query_sdss(sql_pos)
        if res is not None and len(res) > 0:
            df_pos_list.append(res)
        time.sleep(0.5) # Be nice to the server
            
    if not df_pos_list:
        print("  No position data found.")
        return None
        
    df_pos = pd.concat(df_pos_list, ignore_index=True)
    print(f"  Got {len(df_pos)} position records.")
    
    # Step 3: Join
    print("  Joining datasets...")
    df = pd.merge(df_chem, df_pos, on='apogee_id', how='inner')
    print(f"  Merged N={len(df)}")
    
    return df

def compute_rgc(d_kpc, l_deg, b_deg):
    R0 = 8.2
    l = np.radians(l_deg)
    b = np.radians(b_deg)
    
    d_proj = d_kpc * np.cos(b)
    R_plane_sq = R0**2 + d_proj**2 - 2 * R0 * d_proj * np.cos(l)
    Z = d_kpc * np.sin(b)
    R_gc = np.sqrt(R_plane_sq + Z**2)
    return R_gc

def find_knee(df_bin):
    # Heuristic: The knee is the [Fe/H] where the high-alpha sequence turns over.
    # We select stars in the transition zone: 0.05 < alpha < 0.15
    # The mean [Fe/H] of these stars represents the "location" of the down-turn.
    
    transition = df_bin[(df_bin['alpha_m'] > 0.05) & (df_bin['alpha_m'] < 0.15)]
    if len(transition) < 10:
        return np.nan
    
    return transition['fe_h'].mean()

def analyze_alpha_knee(df):
    print("Analyzing Alpha Knee Position...")
    
    # 1. Clean
    df_clean = df.dropna().copy()
    
    # 2. Compute R_gc
    df_clean['R_gc'] = compute_rgc(df_clean['dist50'], df_clean['glon'], df_clean['glat'])
    
    # Filter valid R_gc
    df_clean = df_clean[(df_clean['R_gc'] > 0) & (df_clean['R_gc'] < 20)]
    
    # 3. Bin by R_gc
    bins = np.linspace(2, 16, 8)
    centers = []
    knees = []
    
    print("\nAlpha Knee [Fe/H] by R_gc:")
    for i in range(len(bins)-1):
        r_min, r_max = bins[i], bins[i+1]
        center = (r_min + r_max)/2
        sub = df_clean[(df_clean['R_gc'] >= r_min) & (df_clean['R_gc'] < r_max)]
        
        knee = find_knee(sub)
        print(f"  R={center:.1f} kpc: Knee [Fe/H] = {knee:.3f} (N={len(sub)})")
        
        if not np.isnan(knee):
            centers.append(center)
            knees.append(knee)
            
    # 4. Correlation
    if len(centers) > 2:
        slope, intercept, r_val, p_val, std_err = stats.linregress(centers, knees)
        print(f"\nCorrelation r(Knee, R_gc): {r_val:.4f} (p={p_val:.2e})")
        print(f"Slope: {slope:.4f} dex/kpc")
        
        return {
            'r_knee': float(r_val),
            'slope_knee': float(slope),
            'centers': centers,
            'knees': knees,
            'n_sample': int(len(df_clean))
        }, df_clean
    else:
        print("Not enough bins for correlation.")
        return {'r_knee': 0.0, 'slope_knee': 0.0, 'centers': [], 'knees': [], 'n_sample': int(len(df_clean))}, df_clean

def create_figure(df, results):
    print("Generating figure...")
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # Plot Knee trend
    ax.plot(results['centers'], results['knees'], 'o-', lw=2, color='red', label='Alpha Knee Position')
    
    # Fit line
    if len(results['centers']) > 1:
        x = np.array(results['centers'])
        y = results['slope_knee'] * x + (results['knees'][0] - results['slope_knee']*x[0]) # approx intercept
        ax.plot(x, y, 'k--', label=f"Slope: {results['slope_knee']:.3f} dex/kpc")

    ax.set_xlabel(r'Galactocentric Radius $R_{GC}$ [kpc]')
    ax.set_ylabel(r'Knee Metallicity [Fe/H]')
    ax.set_title(f"Test BM: Alpha Knee vs Radius (r={results['r_knee']:.3f})")
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    out_path = os.path.join(FIGURES_DIR, 'sdss_test_bm_alpha_knee.png')
    plt.savefig(out_path, dpi=150)
    print(f"Figure saved to {out_path}")

def main():
    cache_path = os.path.join(DATA_DIR, 'sdss_alpha_knee.csv')
    
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

    results, df_clean = analyze_alpha_knee(df)
    
    out_path = os.path.join(RESULTS_DIR, 'sdss_test_bm_results.json')
    with open(out_path, 'w') as f:
        json.dump(results, f, indent=2)
        
    create_figure(df_clean, results)
    
    print("\nSUMMARY TEST BM:")
    print("TEP Prediction: Knee [Fe/H] increases as R_gc decreases (Slope < 0).")
    print(f"Observed Slope: {results['slope_knee']:.4f} dex/kpc")
    
    if results['slope_knee'] < -0.01:
        print("RESULT: CONSISTENT (Knee shifts to higher metallicity in inner galaxy)")
    else:
        print("RESULT: NULL/CONTRADICTED")

if __name__ == "__main__":
    main()
