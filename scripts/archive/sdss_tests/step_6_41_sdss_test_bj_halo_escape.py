#!/usr/bin/env python3
"""
Step 6.41: SDSS Test BJ - Halo Escape Velocity (Potential Mapper)

Hypothesis:
The maximum velocity of stars (escape velocity) directly traces the gravitational potential Phi(r).
TEP modifies the effective potential felt by matter (conformal factor).
The shape of the v_esc(r) curve derived from high-velocity halo stars should deviate from the NFW prediction in the inner galaxy if the potential depth is modified by time dilation.

Prediction:
v_esc(r) profile shape differs from NFW expectation.

Data:
- apogeeStar: vhelio_avg, ra, dec, pmra, pmdec.
- apogee_starhorse: dist50 (kpc).

Method:
1. Select high-velocity stars (|v_helio| > 150 km/s).
2. Compute full Galactocentric velocity (v_GC) using RV, PM, Distance.
   - Correct for Solar motion and LSR.
3. Compute Galactocentric radius (r_GC).
4. Bin by r_GC.
5. Estimate v_esc in each bin (e.g., 90th percentile or upper envelope fit).
6. Plot v_esc vs r_GC.
"""

import pandas as pd
import numpy as np
from scipy import stats
import os
import json
import matplotlib.pyplot as plt
import requests
import time
import astropy.units as u
from astropy.coordinates import SkyCoord, Galactocentric
import astropy.coordinates as coord

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
    print(f"Querying SDSS for Test BJ (Limit: {limit})...")
    
    # We need RV, PM, Dist
    
    sql = f"""
    SELECT TOP {limit}
        a.apogee_id,
        a.vhelio_avg,
        a.ra, a.dec,
        a.gaiaedr3_pmra as pmra,
        a.gaiaedr3_pmdec as pmdec,
        s.dist50 as dist_kpc
        
    FROM apogeeStar a
    JOIN apogee_starhorse s ON a.apogee_id = s.apogee_id
    
    WHERE 
        abs(a.vhelio_avg) > 150 -- High velocity candidate
        AND s.dist50 > 0
        AND a.gaiaedr3_pmra IS NOT NULL
        AND a.vhelio_avg > -9999
    """
    return query_sdss(sql)

def compute_galactocentric(df):
    print("Computing Galactocentric coordinates and velocities...")
    
    # Astropy coordinates
    # Assume standard Solar position/velocity
    # R0 = 8.122 kpc, z0 = 20.8 pc
    # v_sun = (12.9, 245.6, 7.78) km/s (Gravity 2018 or similar)
    
    c = SkyCoord(
        ra=df['ra'].values*u.deg, 
        dec=df['dec'].values*u.deg,
        distance=df['dist_kpc'].values*u.kpc,
        pm_ra_cosdec=df['pmra'].values*u.mas/u.yr,
        pm_dec=df['pmdec'].values*u.mas/u.yr,
        radial_velocity=df['vhelio_avg'].values*u.km/u.s
    )
    
    # Transform to Galactocentric
    gc = c.transform_to(coord.Galactocentric)
    
    # Extract R and V_tot
    # R_gc is sqrt(x^2 + y^2 + z^2)
    r_gc = np.sqrt(gc.x.value**2 + gc.y.value**2 + gc.z.value**2)
    
    # V_gc (total velocity in rest frame)
    v_gc = np.sqrt(gc.v_x.value**2 + gc.v_y.value**2 + gc.v_z.value**2)
    
    df['R_gc'] = r_gc
    df['V_gc'] = v_gc
    
    return df

def analyze_escape_velocity(df):
    print("Analyzing Escape Velocity Profile...")
    
    # 1. Clean
    df_clean = df.dropna(subset=['pmra', 'pmdec', 'dist_kpc', 'vhelio_avg']).copy()
    
    # 2. Compute Kinetics
    df_clean = compute_galactocentric(df_clean)
    
    # Filter reasonable bounds
    df_clean = df_clean[(df_clean['R_gc'] > 2) & (df_clean['R_gc'] < 20)].copy()
    df_clean = df_clean[df_clean['V_gc'] < 1000].copy() # Remove outliers
    
    # 3. Bin by Radius
    bins = np.linspace(2, 16, 8) # 2 kpc bins
    
    results_list = []
    
    for i in range(len(bins)-1):
        r_min, r_max = bins[i], bins[i+1]
        sub = df_clean[(df_clean['R_gc'] >= r_min) & (df_clean['R_gc'] < r_max)]
        
        if len(sub) > 10:
            # Estimate v_esc
            # Simple estimator: Max velocity? Or 95th percentile?
            # Leonard & Tremaine (1990): v_esc approx max(v). 
            # We use 90th and 95th percentile to be robust against outliers.
            v_90 = np.percentile(sub['V_gc'], 90)
            v_95 = np.percentile(sub['V_gc'], 95)
            v_max = np.max(sub['V_gc'])
            
            results_list.append({
                'R_gc_bin': (r_min + r_max)/2,
                'v_90': v_90,
                'v_95': v_95,
                'v_max': v_max,
                'n_stars': len(sub)
            })
            
    res_df = pd.DataFrame(results_list)
    print("\nEscape Velocity Profile:")
    print(res_df)
    
    # 4. Fit/Check Trend
    # NFW predicts v_esc decreases with R.
    # v_esc(r) = sqrt(2 * |Phi(r)|).
    # NFW Phi(r) ~ - ln(1+x)/x.
    # Check slope of v_95 vs R_gc.
    
    slope, intercept, r_val, p_val, std_err = stats.linregress(res_df['R_gc_bin'], res_df['v_95'])
    print(f"Slope of v_esc(95) vs R: {slope:.2f} km/s/kpc")
    
    return {
        'slope_vesc': float(slope),
        'r_vesc': float(r_val),
        'binned_data': res_df.to_dict(orient='records'),
        'n_sample': int(len(df_clean))
    }, df_clean

def create_figure(df, results):
    print("Generating figure...")
    data = pd.DataFrame(results['binned_data'])
    
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # Scatter all points
    ax.scatter(df['R_gc'], df['V_gc'], alpha=0.1, s=1, c='k', label='Stars')
    
    # Profile
    ax.plot(data['R_gc_bin'], data['v_95'], 'r-o', lw=2, label='v_esc (95%)')
    ax.plot(data['R_gc_bin'], data['v_max'], 'b--', lw=1, label='v_max')
    
    ax.set_xlabel(r'Galactocentric Radius $R_{GC}$ [kpc]')
    ax.set_ylabel(r'Total Velocity $v_{GC}$ [km/s]')
    ax.set_title(f"Test BJ: Escape Velocity Profile (Slope={results['slope_vesc']:.1f})")
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    out_path = os.path.join(FIGURES_DIR, 'sdss_test_bj_halo_escape.png')
    plt.savefig(out_path, dpi=150)
    print(f"Figure saved to {out_path}")

def main():
    cache_path = os.path.join(DATA_DIR, 'sdss_halo_escape.csv')
    
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

    results, df_clean = analyze_escape_velocity(df)
    
    out_path = os.path.join(RESULTS_DIR, 'sdss_test_bj_results.json')
    with open(out_path, 'w') as f:
        json.dump(results, f, indent=2)
        
    create_figure(df_clean, results)
    
    print("\nSUMMARY TEST BJ:")
    print("TEP Prediction: Deviation from NFW shape (Slope difference?).")
    print(f"Observed Slope: {results['slope_vesc']:.2f} km/s/kpc")
    
    # Standard NFW slope is approx -10 to -20 km/s/kpc in this range?
    # V_esc goes from ~550 at solar to ~400 at 20 kpc. Delta V ~ 150 / 12 ~ 12 km/s/kpc.
    # If slope is significantly different, or positive (?), that's interesting.
    
    if results['slope_vesc'] > -5: 
        print("RESULT: CONSISTENT? (Flat profile, deep potential?)")
    else:
        print("RESULT: NULL (Standard decline)")

if __name__ == "__main__":
    main()
