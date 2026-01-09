#!/usr/bin/env python3
"""
Step 6.11: Test F - Cluster Position Test

TEP HYPOTHESIS:
In TEP, the rate of proper time flow depends on the gravitational potential. 
Galaxies residing deeper in a cluster potential (closer to the center) experience 
greater time dilation relative to the cosmic mean than those at the outskirts.
Therefore, at fixed metallicity and formation epoch (controlled by [Mg/Fe]), 
galaxies in the center should appear "younger" (less evolved) spectroscopically 
than those in the outskirts, or show a systematic offset in age indicators 
that tracks the potential.

TEP PREDICTION:
  Within a cluster, at fixed [Mg/Fe] (formation timescale):
    r(Age_indicator, R_proj) > 0 
    (i.e. Age increases with distance from center; Center appears younger)

DATA:
  - Targeted query of massive clusters (Coma, Abell 85, etc.)
  - galSpecInfo, emissionLinesPort, galSpecIndx

Author: M. Smawfield
Date: January 2026
"""

import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
import json
import os
import requests
import time
from astropy.cosmology import Planck18
from astropy import units as u

# Configuration
SDSS_URL = "https://skyserver.sdss.org/dr18/SkyServerWS/SearchTools/SqlSearch"
DATA_DIR = os.path.join(os.path.dirname(__file__), '..', '..', 'data', 'sdss')
RESULTS_DIR = os.path.join(os.path.dirname(__file__), '..', '..', 'results', 'outputs')
FIGURES_DIR = os.path.join(os.path.dirname(__file__), '..', '..', 'results', 'figures')

os.makedirs(DATA_DIR, exist_ok=True)
os.makedirs(RESULTS_DIR, exist_ok=True)
os.makedirs(FIGURES_DIR, exist_ok=True)

# Target Clusters (Name, RA, Dec, z, R_search_deg)
# Selected for richness and SDSS coverage
CLUSTERS = [
    {"name": "Coma (A1656)", "ra": 194.9531, "dec": 27.9807, "z": 0.0231, "r_deg": 3.0},
    {"name": "Abell 85",     "ra": 10.4587,   "dec": -9.3019,  "z": 0.0551, "r_deg": 1.5},
    {"name": "Abell 119",    "ra": 14.0762,   "dec": -1.2536,  "z": 0.0442, "r_deg": 1.5},
    {"name": "Abell 168",    "ra": 18.7396,   "dec": 0.4306,   "z": 0.0450, "r_deg": 1.5},
    {"name": "Abell 2199",   "ra": 247.1594,  "dec": 39.5513,  "z": 0.0301, "r_deg": 2.0},
    {"name": "Abell 401",    "ra": 44.7408,   "dec": 13.5822,  "z": 0.0737, "r_deg": 1.0},
    {"name": "Abell 1795",   "ra": 207.2188,  "dec": 26.5916,  "z": 0.0625, "r_deg": 1.2},
    {"name": "Abell 2029",   "ra": 227.7342,  "dec": 5.7443,   "z": 0.0773, "r_deg": 1.0},
    {"name": "Abell 2142",   "ra": 239.5833,  "dec": 27.2334,  "z": 0.0909, "r_deg": 0.8},
    {"name": "Abell 2255",   "ra": 258.1292,  "dec": 64.0925,  "z": 0.0806, "r_deg": 1.0},
]

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
                    rows = data[0]["Rows"]
                    if len(rows) > 0:
                        return pd.DataFrame(rows)
                    else:
                        print("  No rows returned in JSON.")
                        return None
            else:
                print(f"  HTTP {response.status_code}")
                # Print snippet of error
                try:
                    print(f"  Response: {response.text[:200]}")
                except:
                    pass
        except Exception as e:
            print(f"  Error: {e}")
        if attempt < max_retries - 1:
            time.sleep(2)
    return None

def fetch_cluster_members(cluster):
    print(f"Querying {cluster['name']} (z={cluster['z']}, R={cluster['r_deg']} deg)...")
    
    # Define redshift range for membership (approx +/- 2000 km/s)
    dz = 0.015 * (1 + cluster['z']) # Broad cut, refine later
    z_min = cluster['z'] - dz
    z_max = cluster['z'] + dz
    
    # Radial search using fGetNearbyObjEq (or just box for speed if function fails, but WHERE box is better)
    # Using simple box cut first for SQL efficiency
    ra_min = cluster['ra'] - cluster['r_deg']
    ra_max = cluster['ra'] + cluster['r_deg']
    dec_min = cluster['dec'] - cluster['r_deg']
    dec_max = cluster['dec'] + cluster['r_deg']
    
    sql = f"""
    SELECT 
        g.specObjID,
        g.ra,
        g.dec,
        g.z,
        
        -- Age/Metallicity
        i.d4000_n,
        i.d4000_n_err,
        i.lick_hd_a,
        i.lick_mgb,
        i.lick_fe5270,
        i.lick_fe5335,
        
        -- Kinematics
        e.sigmaStars,
        e.sigmaStarsErr,
        
        -- Mass
        s.logMass
        
    FROM galSpecInfo g
    JOIN galSpecIndx i ON g.specObjID = i.specObjID
    JOIN emissionLinesPort e ON g.specObjID = e.specObjID
    JOIN stellarMassFSPSGranWideDust s ON g.specObjID = s.specObjID
    
    WHERE 
        g.z BETWEEN {z_min} AND {z_max}
        AND g.ra BETWEEN {ra_min} AND {ra_max}
        AND g.dec BETWEEN {dec_min} AND {dec_max}
        
        -- Quality Cuts
        AND g.reliable = 1
        AND i.d4000_n_err < 0.1
        AND e.sigmaStars > 50 AND e.sigmaStars < 400
        AND s.logMass > 9.0
    """
    
    df = query_sdss(sql)
    if df is not None:
        # Refine radial cut in Python (Great Circle distance)
        from astropy.coordinates import SkyCoord
        c_center = SkyCoord(ra=cluster['ra']*u.deg, dec=cluster['dec']*u.deg)
        c_gal = SkyCoord(ra=df['ra'].values*u.deg, dec=df['dec'].values*u.deg)
        sep = c_center.separation(c_gal).deg
        
        df['R_deg'] = sep
        df['cluster_name'] = cluster['name']
        df['cluster_z'] = cluster['z']
        
        # Keep only those within radius
        df = df[df['R_deg'] <= cluster['r_deg']].copy()
        
        # Convert R_deg to R_Mpc (approximate physical)
        kpc_per_deg = Planck18.kpc_proper_per_arcmin(cluster['z']).value * 60
        df['R_kpc'] = df['R_deg'] * kpc_per_deg
        df['R_Mpc'] = df['R_kpc'] / 1000.0
        
        print(f"  Found {len(df)} members.")
        return df
    else:
        print("  No data returned.")
        return None

def analyze_clusters(df):
    print("Analyzing cluster data...")
    
    # 1. Prepare Variables
    # [Mg/Fe] proxy
    # Avoid division by zero
    denom = 0.5 * (df['lick_fe5270'] + df['lick_fe5335'])
    df['MgFe'] = np.where(denom > 0, df['lick_mgb'] / denom, np.nan)
    
    # Avoid log of zero/negative
    df['log_MgFe'] = np.where(df['MgFe'] > 0, np.log10(df['MgFe']), np.nan)
    
    # Age proxy: D4000 (Higher = Older)
    # TEP Prediction: Higher D4000 at larger R (Positive Correlation)
    
    # Normalize R by some scale? Or just use Mpc.
    # Stacked analysis.
    
    # 2. Control for Metallicity and Mass
    # We want r(Age, R) | MgFe, Mass
    
    df_clean = df.dropna(subset=['d4000_n', 'R_Mpc', 'log_MgFe', 'logMass', 'sigmaStars']).copy()
    # Filter out any remaining infs
    df_clean = df_clean[np.isfinite(df_clean['log_MgFe']) & np.isfinite(df_clean['logMass']) & np.isfinite(df_clean['sigmaStars'])]
    
    from sklearn.linear_model import LinearRegression
    X = df_clean[['log_MgFe', 'logMass', 'sigmaStars']].values
    y = df_clean['d4000_n'].values
    
    # Regress out dependency on internal galaxy properties
    reg = LinearRegression().fit(X, y)
    df_clean['d4000_resid'] = y - reg.predict(X)
    
    # 3. Correlation with Radius
    r_simple, p_simple = stats.pearsonr(df_clean['R_Mpc'], df_clean['d4000_n'])
    r_resid, p_resid = stats.pearsonr(df_clean['R_Mpc'], df_clean['d4000_resid'])
    
    print(f"Simple r(D4000, R_Mpc): {r_simple:.4f} (p={p_simple:.2e})")
    print(f"Controlled r(D4000_resid, R_Mpc): {r_resid:.4f} (p={p_resid:.2e})")
    
    return {
        'r_simple': float(r_simple),
        'p_simple': float(p_simple),
        'r_controlled': float(r_resid),
        'p_controlled': float(p_resid),
        'n_sample': int(len(df_clean))
    }

def create_figure(df, results):
    # Need to reconstruct cleaned df or handle nans (same logic as above)
    denom = 0.5 * (df['lick_fe5270'] + df['lick_fe5335'])
    df['MgFe'] = np.where(denom > 0, df['lick_mgb'] / denom, np.nan)
    df['log_MgFe'] = np.where(df['MgFe'] > 0, np.log10(df['MgFe']), np.nan)

    df_clean = df.dropna(subset=['d4000_n', 'R_Mpc', 'log_MgFe', 'logMass', 'sigmaStars']).copy()
    df_clean = df_clean[np.isfinite(df_clean['log_MgFe']) & np.isfinite(df_clean['logMass']) & np.isfinite(df_clean['sigmaStars'])]

    from sklearn.linear_model import LinearRegression
    X = df_clean[['log_MgFe', 'logMass', 'sigmaStars']].values
    y = df_clean['d4000_n'].values
    reg = LinearRegression().fit(X, y)
    df_clean['d4000_resid'] = y - reg.predict(X)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    # D4000 vs R
    ax = axes[0]
    ax.hexbin(df_clean['R_Mpc'], df_clean['d4000_n'], gridsize=25, cmap='inferno', mincnt=1)
    ax.set_xlabel('Projected Radius (Mpc)')
    ax.set_ylabel('D4000 (Age Indicator)')
    ax.set_title(f'Age vs Radius (r={results["r_simple"]:.3f})')
    
    # Residual vs R
    ax = axes[1]
    ax.hexbin(df_clean['R_Mpc'], df_clean['d4000_resid'], gridsize=25, cmap='inferno', mincnt=1)
    ax.axhline(0, color='w', ls='--', alpha=0.5)
    ax.set_xlabel('Projected Radius (Mpc)')
    ax.set_ylabel('D4000 Residual (controlled)')
    ax.set_title(f'Controlled Age vs Radius (r={results["r_controlled"]:.3f})')
    
    plt.tight_layout()
    out_file = os.path.join(FIGURES_DIR, 'sdss_test_f_cluster_position.png')
    plt.savefig(out_file, dpi=150)
    print(f"Figure saved to {out_file}")

def main():
    print("="*60)
    print("STEP 6.11: TEST F - CLUSTER POSITION TEST")
    print("="*60)
    
    cache_path = os.path.join(DATA_DIR, 'sdss_cluster_members.csv')
    
    if os.path.exists(cache_path):
        print("Loading cached data...")
        df_all = pd.read_csv(cache_path)
    else:
        all_dfs = []
        for cluster in CLUSTERS:
            df = fetch_cluster_members(cluster)
            if df is not None:
                all_dfs.append(df)
        
        if all_dfs:
            df_all = pd.concat(all_dfs, ignore_index=True)
            df_all.to_csv(cache_path, index=False)
        else:
            print("No data collected.")
            return

    print(f"Total Members: {len(df_all)}")
    
    results = analyze_clusters(df_all)
    create_figure(df_all, results)
    
    out_path = os.path.join(RESULTS_DIR, 'sdss_test_f_results.json')
    with open(out_path, 'w') as f:
        json.dump(results, f, indent=2)
        
    print("\nSUMMARY:")
    print(f"TEP Consistent (r > 0): {results['r_controlled'] > 0}")
    print(f"Significance: p={results['p_controlled']:.2e}")

if __name__ == "__main__":
    main()
