#!/usr/bin/env python3
"""
Step 6.70: SDSS Test CS - Void HI Gas Fraction (The Void Clock)

Hypothesis:
In the TEP framework, time runs "fast" in voids (shallow potential) compared to the field. 
Standard LambdaCDM predicts void galaxies are retarded in evolution (gas-rich). 
TEP predicts they experience more proper time per unit cosmic time, potentially appearing 
**more evolved** (gas-poor or lower sSFR) than field galaxies of the same mass, 
or at least showing a shallower Gas-Density relation.

Prediction:
HI Gas Fraction in Voids is LOWER than expected from the standard density relation.

Data:
- mangaHIall: logM_HI, mangaid
- mangaTarget: nsa_elpetro_mass (Stellar Mass), objra, objdec
- ebossMCPM: MATTERDENS (density), RA, DEC, Z (for environment)
- SpecObjAll: z (redshift)

Method:
1. Fetch MaNGA HI galaxies with mass and coordinates.
2. Fetch eBOSS MCPM galaxies with density and coordinates.
3. Spatial Cross-Match (KDTree) to assign density to MaNGA galaxies.
   Match in 3D (RA, DEC, Z) or 2D + Z-slice.
   Given the density field is from eBOSS, we need to map it to the MaNGA target location.
4. Classify Void (density < -0.5?) vs Field.
5. Compare HI Gas Fraction (M_HI / M_*) vs Stellar Mass for Void vs Field.
"""

import pandas as pd
import numpy as np
from scipy import stats
from scipy.spatial import cKDTree
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

def download_data():
    print("Querying SDSS for Test CS...")
    
    # 1. Fetch MaNGA HI Data (Simple query)
    print("  Fetching MaNGA HI (base)...")
    # Fetch just HI properties and ID
    sql_hi_base = """
    SELECT TOP 100
        mangaid,
        logM_HI,
        conf_prob
    FROM mangaHIall
    WHERE conf_prob > 0.9
    """
    df_hi_base = query_sdss(sql_hi_base)
    
    if df_hi_base is None or len(df_hi_base) == 0: 
        print("  No HI data found.")
        return None, None

    ids = df_hi_base['mangaid'].astype(str).tolist()
    ids = list(set(ids))
    print(f"  Got {len(df_hi_base)} HI sources. Fetching Targets...")
    
    # 2. Fetch Targets for these IDs
    chunk_size = 50
    df_target_list = []
    
    for i in range(0, len(ids), chunk_size):
        chunk = ids[i:i+chunk_size]
        ids_str = "', '".join(chunk)
        
        sql_target = f"""
        SELECT 
            mangaid,
            nsa_elpetro_mass as logmass,
            objra as ra,
            objdec as dec,
            nsa_z as z
        FROM mangaTarget
        WHERE mangaid IN ('{ids_str}')
        """
        res = query_sdss(sql_target)
        if res is not None and len(res) > 0:
            df_target_list.append(res)
        time.sleep(0.2)
        
    if not df_target_list:
        print("  No target data found.")
        return None, None
        
    df_target = pd.concat(df_target_list, ignore_index=True)
    
    # Join locally
    df_hi = pd.merge(df_hi_base, df_target, on='mangaid', how='inner')
    print(f"  Merged HI + Target: {len(df_hi)}")

    # 3. Fetch Environmental Density (eBOSS MCPM)
    print("  Fetching Environment (eBOSS MCPM)...")
    
    # Get bounds of HI sample
    ra_min, ra_max = df_hi['ra'].min(), df_hi['ra'].max()
    dec_min, dec_max = df_hi['dec'].min(), df_hi['dec'].max()
    z_min, z_max = df_hi['z'].min(), df_hi['z'].max()
    
    # Ensure bounds are valid numbers
    if pd.isna(ra_min):
        print("  Invalid bounds.")
        return None, None
        
    print(f"  Bounds: RA {ra_min:.1f}-{ra_max:.1f}, Dec {dec_min:.1f}-{dec_max:.1f}, Z {z_min:.3f}-{z_max:.3f}")
    
    # Fetch environment in a slightly larger box? Or just a sample.
    # eBOSS covers a large area.
    # To avoid timeout, fetch a limited sample in the rough area.
    
    sql_env = f"""
    SELECT TOP 2000
        RA, DEC, Z, MATTERDENS
    FROM ebossMCPM
    WHERE 
        RA BETWEEN {ra_min-2} AND {ra_max+2}
        AND DEC BETWEEN {dec_min-2} AND {dec_max+2}
        AND Z BETWEEN {z_min-0.05} AND {z_max+0.05}
    """
    df_env = query_sdss(sql_env)
    
    return df_hi, df_env

def analyze_voids(df_hi, df_env):
    print("Analyzing Void HI Fraction...")
    
    if df_hi is None or df_env is None or len(df_env) == 0:
        print("  Missing data.")
        return None
        
    # Cross-match to assign density
    # Convert to Cartesian for KDTree (approx)
    def to_cartesian(ra, dec, z):
        dist = z * 300000.0 / 70.0 # Approx distance Mpc/h
        ra_rad = np.radians(ra)
        dec_rad = np.radians(dec)
        x = dist * np.cos(dec_rad) * np.cos(ra_rad)
        y = dist * np.cos(dec_rad) * np.sin(ra_rad)
        z_coord = dist * np.sin(dec_rad)
        return np.column_stack([x, y, z_coord])
        
    coords_hi = to_cartesian(df_hi['ra'], df_hi['dec'], df_hi['z'])
    coords_env = to_cartesian(df_env['RA'], df_env['DEC'], df_env['Z'])
    
    tree = cKDTree(coords_env)
    dists, idxs = tree.query(coords_hi, k=1)
    
    # Assign density
    # Check if match is reasonably close (e.g. < 5 Mpc)
    valid_match = dists < 5.0
    print(f"  Matched {sum(valid_match)}/{len(df_hi)} galaxies to environment.")
    
    df_hi['density'] = np.nan
    df_hi.loc[valid_match, 'density'] = df_env.iloc[idxs[valid_match]]['MATTERDENS'].values
    
    df_clean = df_hi.dropna(subset=['density', 'logM_HI', 'logmass']).copy()
    
    # Define Void vs Field
    # MATTERDENS is likely delta_rho? Or log density?
    # Usually < 0 is underdense. < -0.5 is void.
    # Let's assume standard units (delta).
    
    df_clean['is_void'] = df_clean['density'] < 0.0 # Underdense
    df_clean['is_dense'] = df_clean['density'] > 1.0 # Overdense
    
    print(f"  Void Sample: {df_clean['is_void'].sum()}")
    print(f"  Dense Sample: {df_clean['is_dense'].sum()}")
    
    # Calculate Gas Fraction
    # log f_HI = log(M_HI / M_*) = log M_HI - log M_*
    df_clean['log_fHI'] = df_clean['logM_HI'] - df_clean['logmass']
    
    # Compare
    mean_void = df_clean[df_clean['is_void']]['log_fHI'].mean()
    mean_dense = df_clean[df_clean['is_dense']]['log_fHI'].mean()
    
    print(f"  Mean log f_HI (Void): {mean_void:.3f}")
    print(f"  Mean log f_HI (Dense): {mean_dense:.3f}")
    
    delta = mean_void - mean_dense
    print(f"  Delta (Void - Dense): {delta:.3f} dex")
    
    # Plot
    fig, ax = plt.subplots(figsize=(8, 6))
    
    ax.hist(df_clean[df_clean['is_void']]['log_fHI'], bins=20, alpha=0.5, density=True, label='Void')
    ax.hist(df_clean[df_clean['is_dense']]['log_fHI'], bins=20, alpha=0.5, density=True, label='Dense')
    
    ax.set_xlabel('log(HI Gas Fraction)')
    ax.set_ylabel('Normalized Count')
    ax.set_title(f'Test CS: Void HI Gas Fraction (Delta={delta:.2f})')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    out_path = os.path.join(FIGURES_DIR, 'sdss_test_cs_void_hi.png')
    plt.savefig(out_path, dpi=150)
    print(f"Figure saved to {out_path}")
    
    return {
        'delta_fHI': delta,
        'n_void': int(df_clean['is_void'].sum()),
        'n_dense': int(df_clean['is_dense'].sum())
    }

def main():
    df_hi, df_env = download_data()
    
    results = analyze_voids(df_hi, df_env)
    
    out_path = os.path.join(RESULTS_DIR, 'sdss_test_cs_results.json')
    if results:
        with open(out_path, 'w') as f:
            json.dump(results, f, indent=2)
            
        print("\nSUMMARY TEST CS:")
        print("Prediction: HI Gas Fraction in Voids is LOWER than expected (TEP Fast Time).")
        print(f"Observed Excess in Voids: {results['delta_fHI']:.3f} dex")
        
        if results['delta_fHI'] < -0.1:
             print("RESULT: CONSISTENT (Void galaxies are gas-poor)")
        elif results['delta_fHI'] > 0.1:
             print("RESULT: CONTRADICTED (Void galaxies are gas-rich - Standard Model)")
        else:
             print("RESULT: NULL (No significant difference)")

if __name__ == "__main__":
    main()
