#!/usr/bin/env python3
"""
Step 6.95: SDSS Test DY - Vertical Phase Spirals (Dynamical Damping)

Hypothesis:
Gaia revealed "phase spirals" (snail shells) in the z-vz plane of the Galactic disk.
The winding of the spiral depends on the vertical frequency nu_z, determined by the 
local potential density rho. TEP modifies the effective density. The winding rate 
should differ from standard models.

Prediction:
Phase Spiral Winding (or existence) at different R_gc mismatches dynamical models.
(Here we characterize the vertical kinematics Z vs Vz as a function of R_gc).

Data:
- apogeeStar: vhelio_avg
- apogee_starhorse: dist50
- mos_gaia_dr2_source: pmra, pmdec
- Coordinates: ra, dec from apogeeStar

Method:
1. Join tables.
2. Calculate Galactocentric Cylindrical Coords (R, phi, Z) and Velocities (Vr, Vphi, Vz).
   - Use Astropy.
   - Solar position: R0 = 8.122 kpc, Z0 = 0.0208 kpc.
   - Solar velocity: (11.1, 242.0, 7.25) km/s (proper motion + LSR).
3. Select Disk stars (e.g. |Z| < 2 kpc).
4. Bin by R_gc (Inner: <7 kpc, Solar: 7-9 kpc, Outer: >9 kpc).
5. Analyze Phase Space (Z vs Vz).
   - Calculate vertical velocity dispersion sigma_z.
   - Visualize the phase space (save plot).
   - Check for "clumpiness" or "spirality" (difficult to quantify automatically, 
     will rely on visual inspection of the generated figure and sigma_z trends).
"""

import pandas as pd
import numpy as np
from scipy import stats
import os
import json
import matplotlib.pyplot as plt
import requests
import time
from astropy.coordinates import SkyCoord, Galactocentric
import astropy.units as u

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
    print(f"Querying SDSS for Test DY (Limit: {limit})...")
    
    # Need APOGEE RVs + Distances + Gaia PMs
    # utilizing gaiaedr3 columns directly in apogeeStar to avoid joins
    
    sql = f"""
    SELECT TOP {limit}
        a.apogee_id,
        a.ra, a.dec,
        a.vhelio_avg,
        s.dist50,
        a.gaiaedr3_pmra as pmra, 
        a.gaiaedr3_pmdec as pmdec
    FROM apogeeStar a
    JOIN apogee_starhorse s ON a.apogee_id = s.apogee_id
    WHERE a.vhelio_avg > -9999
      AND s.dist50 > 0
      AND abs(a.glat) > 10 -- Avoid midplane extinction/confusion, focus on vertical extent
      AND abs(a.glat) < 60 -- Stay within reasonable disk projection
      AND a.gaiaedr3_pmra != 0
    """
    return query_sdss(sql)

def process_coordinates(df):
    print("Converting to Galactocentric coordinates...")
    
    # Setup coordinates
    c = SkyCoord(
        ra=df['ra'].values*u.deg,
        dec=df['dec'].values*u.deg,
        distance=df['dist50'].values*u.kpc,
        pm_ra_cosdec=df['pmra'].values*u.mas/u.yr,
        pm_dec=df['pmdec'].values*u.mas/u.yr,
        radial_velocity=df['vhelio_avg'].values*u.km/u.s,
        frame='icrs'
    )
    
    # Galactocentric frame parameters
    # Using defaults for now: R_sun=8.122 kpc, z_sun=20.8 pc
    gc = c.transform_to(Galactocentric())
    
    # Cartesian components
    x = gc.x.to(u.kpc).value
    y = gc.y.to(u.kpc).value
    z = gc.z.to(u.kpc).value
    
    vx = gc.v_x.to(u.km/u.s).value
    vy = gc.v_y.to(u.km/u.s).value
    vz = gc.v_z.to(u.km/u.s).value
    
    # Calculate Cylindrical
    R = np.sqrt(x**2 + y**2)
    phi = np.degrees(np.arctan2(y, x))
    
    # v_R = (x*vx + y*vy) / R
    v_R = (x * vx + y * vy) / R
    
    # v_phi = (x*vy - y*vx) / R (Linear tangential velocity)
    v_phi = (x * vy - y * vx) / R
    
    df['R'] = R
    df['phi'] = phi
    df['z'] = z
    
    df['v_R'] = v_R
    df['v_phi'] = v_phi
    df['v_z'] = vz
    
    return df

def analyze_phase_space(df):
    print("Analyzing Vertical Phase Space...")
    
    if df is None or len(df) < 100:
        print("  Insufficient data.")
        return None
        
    df = process_coordinates(df)
    
    # Bin by R
    bins = [0, 7, 9, 15]
    labels = ['Inner (<7 kpc)', 'Solar (7-9 kpc)', 'Outer (>9 kpc)']
    
    df['R_bin'] = pd.cut(df['R'], bins=bins, labels=labels)
    
    results = {}
    
    plt.figure(figsize=(15, 5))
    
    for i, label in enumerate(labels):
        subset = df[df['R_bin'] == label]
        if len(subset) < 10:
            continue
            
        # Metrics
        sigma_z = subset['v_z'].std()
        sigma_dist = subset['z'].std()
        
        results[label] = {
            'sigma_z': float(sigma_z),
            'sigma_dist': float(sigma_dist),
            'count': int(len(subset))
        }
        
        # Plot Z vs Vz
        plt.subplot(1, 3, i+1)
        plt.scatter(subset['z'], subset['v_z'], s=1, alpha=0.3, color='black')
        # Density contours could be better but scatter is fine for "seeing" spirals
        
        plt.title(f"{label}\n$\sigma_z$={sigma_z:.1f} km/s")
        plt.xlabel('Z (kpc)')
        plt.ylabel('$V_Z$ (km/s)')
        plt.xlim(-2, 2)
        plt.ylim(-100, 100)
        plt.grid(True, alpha=0.3)
        
    plt.tight_layout()
    out_path = os.path.join(FIGURES_DIR, 'sdss_test_dy_phase_spirals.png')
    plt.savefig(out_path, dpi=150)
    print(f"Figure saved to {out_path}")
    
    # Check for trend in Sigma_Z vs R
    # In standard disk, Sigma_Z should decrease with R (exp scale length)
    # TEP might predict different scaling if potential is effectively deeper/shallower
    
    if 'Solar (7-9 kpc)' in results and 'Outer (>9 kpc)' in results:
        ratio = results['Outer (>9 kpc)']['sigma_z'] / results['Solar (7-9 kpc)']['sigma_z']
        print(f"  Outer/Solar Sigma_Z Ratio: {ratio:.3f}")
        results['outer_solar_ratio'] = ratio
    
    return results

def main():
    cache_path = os.path.join(DATA_DIR, 'sdss_phase_spirals.csv')
    
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

    results = analyze_phase_space(df)
    
    out_path = os.path.join(RESULTS_DIR, 'sdss_test_dy_results.json')
    if results:
        with open(out_path, 'w') as f:
            json.dump(results, f, indent=2)
            
        print("\nSUMMARY TEST DY:")
        if 'outer_solar_ratio' in results:
            print(f"Sigma_Z Ratio (Outer/Solar): {results['outer_solar_ratio']:.3f}")
            if results['outer_solar_ratio'] > 1.0:
                 print("RESULT: ANOMALY (Outer disk hotter than Solar?)")
            else:
                 print("RESULT: STANDARD (Outer disk cooler)")
        else:
            print("RESULT: INCONCLUSIVE (Insufficient coverage)")

if __name__ == "__main__":
    main()
