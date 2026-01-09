#!/usr/bin/env python3
"""
Step 6.43: SDSS Test BL - White Dwarf Mass Shift (AGB Wind Clock)

Hypothesis:
The Initial-Final Mass Relation (IFMR) depends on how much mass a star loses during the AGB phase.
Mass loss is driven by winds (rate process).
In deep potentials (Inner Galaxy), time dilation suppresses the wind rate.
Stars should retain more mass, resulting in more massive White Dwarfs for a given progenitor mass.
Prediction: Mean White Dwarf Mass (inferred from Mag/Color) decreases with R_gc.
Or: Inner WDs are more massive (fainter at fixed color/Teff).
Observable: Absolute Magnitude M_G vs R_gc at fixed color.
TEP: M_G should be fainter (larger) at small R_gc.

Data:
- mos_gaia_dr2_wd: source_id, g_mag_abs, bp_rp, pwd (probability WD).
- mos_geometric_distances_gaia_dr2: r_est (dist), l, b (coords via Gaia source match).
  Note: mos_gaia_dr2_wd might not have l, b directly. We need coordinates.
  Or calculate R_gc from r_est if l,b available.
  Let's assume we can get l,b from mos_gaia_dr2_source or if they are in the WD table.
  The query plan suggested joining geometric distances.

Method:
1. Select high-prob WDs (pwd > 0.9).
2. Compute R_gc.
3. Select a color range (e.g., -0.5 < bp_rp < 1.5) where WDs follow a cooling track.
4. Correct for Cooling Track: Fit M_G vs bp_rp. Compute Residuals Delta_M.
   Positive Delta_M = Fainter = More Massive (for WDs, Radius decreases with Mass).
5. Correlate Delta_M with R_gc.
   TEP: Delta_M > 0 at small R_gc. r < 0 (since Delta_M decreases as R_gc increases).
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
    print(f"Querying SDSS for Test BL (Limit: {limit})...")
    
    # Use mos_gaia_dr2_wd ONLY. 
    # Correct columns: g_gaia_mag, bpmag, rpmag, plx, glon, glat.
    
    sql = f"""
    SELECT TOP {limit}
        source_id,
        g_gaia_mag,
        bpmag, rpmag,
        plx as parallax,
        glon as l, 
        glat as b
        
    FROM mos_gaia_dr2_wd
    
    WHERE 
        pwd > 0.9 
        AND (bpmag - rpmag) BETWEEN -0.5 AND 1.5
        AND plx > 0.1 -- d < 10 kpc
    """
    return query_sdss(sql)

def compute_rgc(d_kpc, l_deg, b_deg):
    R0 = 8.2
    l = np.radians(l_deg)
    b = np.radians(b_deg)
    
    d_proj = d_kpc * np.cos(b)
    R_plane_sq = R0**2 + d_proj**2 - 2 * R0 * d_proj * np.cos(l)
    Z = d_kpc * np.sin(b)
    R_gc = np.sqrt(R_plane_sq + Z**2)
    return R_gc

def analyze_wd_mass_shift(df):
    print("Analyzing White Dwarf Mass Shift...")
    
    # 1. Clean
    df_clean = df.dropna().copy()
    
    # 2. Compute Derived Columns
    df_clean['bp_rp'] = df_clean['bpmag'] - df_clean['rpmag']
    
    # Dist = 1000 / parallax (mas)
    df_clean['dist_pc'] = 1000.0 / df_clean['parallax']
    df_clean['dist_kpc'] = df_clean['dist_pc'] / 1000.0
    
    # Absolute Mag M_G = m_g - 5 log(d/10)
    df_clean['g_mag_abs'] = df_clean['g_gaia_mag'] - 5 * np.log10(df_clean['dist_pc']) + 5
    
    df_clean['R_gc'] = compute_rgc(df_clean['dist_kpc'], df_clean['l'], df_clean['b'])
    
    # 3. Fit Cooling Sequence (M_G vs bp_rp)
    # WDs lie on a tight sequence.
    # Fit polynomial M_G = f(bp_rp).
    
    coeffs = np.polyfit(df_clean['bp_rp'], df_clean['g_mag_abs'], 3)
    poly = np.poly1d(coeffs)
    
    df_clean['mag_pred'] = poly(df_clean['bp_rp'])
    
    # Residual Delta_M = Obs - Pred.
    # Positive Residual -> Fainter than expected -> More Massive (smaller radius).
    # TEP: Expect Positive Residuals at small R_gc.
    df_clean['mag_resid'] = df_clean['g_mag_abs'] - df_clean['mag_pred']
    
    # 4. Correlation
    r_wd, p_wd = stats.pearsonr(df_clean['R_gc'], df_clean['mag_resid'])
    
    print(f"N = {len(df_clean)}")
    print(f"Correlation r(Delta_M, R_gc): {r_wd:.4f} (p={p_wd:.2e})")
    
    # 5. Binning
    df_clean['rgc_bin'] = pd.qcut(df_clean['R_gc'], 8)
    binned = df_clean.groupby('rgc_bin')['mag_resid'].mean()
    print("\nMean Mag Residual by R_gc:")
    print(binned)
    
    return {
        'r_wd': float(r_wd),
        'p_wd': float(p_wd),
        'binned_means': binned.tolist(),
        'bin_centers': [mid.mid for mid in binned.index],
        'n_sample': int(len(df_clean))
    }, df_clean

def create_figure(df, results):
    print("Generating figure...")
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # Scatter
    ax.scatter(df['R_gc'], df['mag_resid'], alpha=0.1, s=2, c='k', label='White Dwarfs')
    
    # Binned
    ax.plot(results['bin_centers'], results['binned_means'], 'r-o', lw=2, label='Mean Residual')
    
    ax.set_xlabel(r'Galactocentric Radius $R_{GC}$ [kpc]')
    ax.set_ylabel(r'WD Mag Residual $\Delta M_G$ [Fainter > 0]')
    ax.set_title(f"Test BL: WD Mass Shift (r={results['r_wd']:.3f})")
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.axhline(0, linestyle='--', color='b')
    
    plt.tight_layout()
    out_path = os.path.join(FIGURES_DIR, 'sdss_test_bl_wd_mass.png')
    plt.savefig(out_path, dpi=150)
    print(f"Figure saved to {out_path}")

def main():
    cache_path = os.path.join(DATA_DIR, 'sdss_wd_mass.csv')
    
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

    results, df_clean = analyze_wd_mass_shift(df)
    
    out_path = os.path.join(RESULTS_DIR, 'sdss_test_bl_results.json')
    with open(out_path, 'w') as f:
        results_json = results.copy()
        del results_json['bin_centers']
        json.dump(results_json, f, indent=2)
        
    create_figure(df_clean, results)
    
    print("\nSUMMARY TEST BL:")
    print("TEP Prediction: Inner WDs are massive/fainter (Delta_M > 0). r < 0.")
    print(f"Observed r: {results['r_wd']:.4f}")
    
    if results['r_wd'] < -0.05:
        print("RESULT: CONSISTENT (Massive WDs in inner galaxy)")
    elif results['r_wd'] > 0.05:
        print("RESULT: CONTRADICTED (Less massive WDs in inner galaxy)")
    else:
        print("RESULT: NULL/FLAT")

if __name__ == "__main__":
    main()
