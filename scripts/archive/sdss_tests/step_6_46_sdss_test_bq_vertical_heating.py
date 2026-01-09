#!/usr/bin/env python3
"""
Step 6.46: SDSS Test BQ - Vertical Disk Heating (Kinematic Approach)

Hypothesis:
Stellar disks thicken over time due to dynamical scattering (heating). 
TEP predicts that in the deep potential of the inner Galaxy, heating rates are time-dilated.
Therefore, for a given stellar population (fixed age/abundance), the vertical velocity dispersion 
(sigma_z) should be LOWER at small R_gc than standard dynamic models predict (which scale with surface density).

Prediction:
Standard: sigma_z increases exponentially towards the center (following surface density).
TEP: sigma_z profile is flatter or depressed in the inner galaxy relative to the density scaling.

Data:
- aspcapStar: fparam_alpha_m, fparam_m_h (Chemistry)
- apogeeStar: ra, dec, gaiaedr3_pmra, gaiaedr3_pmdec, vhelio_avg, gaiaedr3_r_med_photogeo

Method:
1. Fetch Data (Chemistry + Kinematics).
2. Transform to Galactocentric Velocities (U, V, W) using astropy.
3. Select Mono-abundance populations (Thick/Thin disk).
4. Compute dispersion sigma_W (sigma_z) in bins of R_gc.
5. Analyze the slope of sigma_z vs R_gc.
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

def download_data(limit=1000):
    print(f"Querying SDSS for Test BQ (Kinematic) - Limit {limit}...")
    
    # Step 1: Chemistry
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
        
    ids = df_chem['apogee_id'].astype(str).tolist()
    print(f"  Got {len(df_chem)} stars. Fetching kinematics...")
    
    # Step 2: Kinematics
    # Need RA, Dec, PMs, RV, Dist
    # Reduce chunk size to avoid URL limit
    chunk_size = 50
    df_pos_list = []
    
    for i in range(0, len(ids), chunk_size):
        if i % 200 == 0:
             print(f"    Chunk {i} / {len(ids)}...")
        chunk = ids[i:i+chunk_size]
        ids_str = "', '".join(chunk)
        
        sql_pos = f"""
        SELECT 
            apogee_id,
            ra, dec,
            gaiaedr3_pmra as pmra,
            gaiaedr3_pmdec as pmdec,
            vhelio_avg as rv,
            gaiaedr3_r_med_photogeo as dist
        FROM apogeeStar
        WHERE 
            apogee_id IN ('{ids_str}')
            AND gaiaedr3_r_med_photogeo > 0
            AND vhelio_avg > -9000
        """
        res = query_sdss(sql_pos)
        if res is not None and len(res) > 0:
            df_pos_list.append(res)
        time.sleep(0.2)
            
    if not df_pos_list:
        print("  No kinematic data found.")
        return None
        
    df_pos = pd.concat(df_pos_list, ignore_index=True)
    
    # Step 3: Join
    print("  Joining datasets...")
    df = pd.merge(df_chem, df_pos, on='apogee_id', how='inner')
    print(f"  Merged N={len(df)}")
    
    return df

def compute_velocities(df):
    print("  Computing Galactocentric Velocities...")
    
    # Filter valid data
    df = df.dropna(subset=['ra', 'dec', 'pmra', 'pmdec', 'rv', 'dist']).copy()
    
    # Create SkyCoord object
    # Distance in pc (APOGEE uses pc usually? "r_med_photogeo" is pc in Gaia EDR3 catalog)
    # APOGEE query returns what? usually pc.
    
    c = SkyCoord(
        ra=df['ra'].values*u.degree,
        dec=df['dec'].values*u.degree,
        distance=df['dist'].values*u.pc,
        pm_ra_cosdec=df['pmra'].values*u.mas/u.yr,
        pm_dec=df['pmdec'].values*u.mas/u.yr,
        radial_velocity=df['rv'].values*u.km/u.s
    )
    
    # Transform to Galactocentric
    # Use default Astropy parameters (R0=8.122, z_sun=20.8, etc)
    gc = c.transform_to(Galactocentric())
    
    df['R_gc'] = gc.cylindrical.rho.to(u.kpc).value
    df['Z_gc'] = gc.cylindrical.z.to(u.kpc).value
    # Vertical velocity W is v_z (Cartesian Z velocity in Galactocentric frame)
    df['W_vel'] = gc.v_z.to(u.km/u.s).value
    
    # Also get azimuthal? Not strictly needed for BQ but good for context
    # V_phi approx
    
    return df

def analyze_vertical_heating(df):
    print("Analyzing Kinematics...")
    
    df = compute_velocities(df)
    
    # Clean spatial selection
    df = df[(df['R_gc'] > 3) & (df['R_gc'] < 15)].copy()
    
    # Define Populations
    pops = {
        'Thin': df[df['alpha_m'] < 0.10],
        'Thick': df[df['alpha_m'] > 0.15]
    }
    
    results = {}
    
    fig, ax = plt.subplots(figsize=(8, 6))
    
    for name, sub in pops.items():
        if len(sub) < 50:
            continue
            
        print(f"\nPopulation: {name} (N={len(sub)})")
        
        # Bin by R_gc
        sub['r_bin'] = pd.cut(sub['R_gc'], bins=np.linspace(3, 15, 7))
        
        # Calculate sigma_W (dispersion of vertical velocity)
        # Use robust estimator (std, or IQR)
        # std is fine if outliers removed.
        
        # Remove high velocity outliers (runaway stars)
        sub = sub[np.abs(sub['W_vel']) < 200]
        
        stats_df = sub.groupby('r_bin', observed=True)['W_vel'].apply(lambda x: np.nanstd(x)).reset_index()
        stats_df.columns = ['r_bin', 'sigma_z']
        
        # Get bin centers
        stats_df['r_center'] = stats_df['r_bin'].apply(lambda x: x.mid).astype(float)
        stats_df = stats_df.dropna()
        
        # Error bars: approx sigma / sqrt(2N)
        counts = sub.groupby('r_bin', observed=True)['W_vel'].count().values
        stats_df['err'] = stats_df['sigma_z'] / np.sqrt(2 * counts)
        
        print(stats_df)
        
        # Fit trend
        slope, intercept, r_val, p_val, std_err = stats.linregress(stats_df['r_center'], stats_df['sigma_z'])
        print(f"  Slope sigma_z vs R: {slope:.3f} km/s/kpc")
        
        results[name] = {
            'slope': slope,
            'r_val': r_val,
            'data': stats_df.to_dict(orient='records')
        }
        
        ax.errorbar(stats_df['r_center'], stats_df['sigma_z'], yerr=stats_df['err'], fmt='-o', label=f"{name} (Slope={slope:.3f})")
    
    ax.set_xlabel('Galactocentric Radius R [kpc]')
    ax.set_ylabel('Vertical Velocity Dispersion sigma_z [km/s]')
    ax.set_title('Test BQ: Vertical Disk Heating Profile')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    out_path = os.path.join(FIGURES_DIR, 'sdss_test_bq_kinematics.png')
    plt.savefig(out_path, dpi=150)
    print(f"Figure saved to {out_path}")
    
    return results

def main():
    cache_path = os.path.join(DATA_DIR, 'sdss_heating_kin.csv')
    
    # Always re-download for this new kinematic version or check if cache has velocity
    # Let's just force download to be safe/simple
    df = download_data()
    if df is not None:
        df.to_csv(cache_path, index=False)
    else:
        print("Download failed.")
        return

    results = analyze_vertical_heating(df)
    
    out_path = os.path.join(RESULTS_DIR, 'sdss_test_bq_results.json')
    
    def json_default(obj):
        if isinstance(obj, pd.Interval):
            return str(obj)
        raise TypeError
        
    with open(out_path, 'w') as f:
        json.dump(results, f, indent=2, default=json_default)
        
    print("\nSUMMARY TEST BQ (Kinematic):")
    # TEP Prediction: Inner disk heating is suppressed.
    # Standard: Sigma_z ~ exp(-R/2Rd). Slope should be NEGATIVE (high in center).
    # TEP: Slope should be LESS NEGATIVE (flatter) or even flat?
    # Actually, if TEP makes inner galaxy "younger" dynamically, sigma_z should be lower than standard.
    # Standard: Slope ~ -5 km/s/kpc (depends on scale length).
    # If we see a very flat profile, that might indicate suppression.
    
    if 'Thin' in results:
        slope = results['Thin']['slope']
        print(f"Thin Disk Slope: {slope:.3f}")
        if slope > -1.0: # Very flat or positive
             print("RESULT: CONSISTENT (Flat dispersion profile)")
        else:
             print("RESULT: NULL (Standard negative gradient)")

if __name__ == "__main__":
    main()
