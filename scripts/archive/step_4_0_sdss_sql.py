#!/usr/bin/env python3
"""
Step 4.0: Download SDSS Galaxy Data via SQL

Uses SDSS SkyServer SQL interface to download galaxy properties
for the temporal onion test. JSON format for reliable parsing.

Author: M. Smawfield
Date: January 2026
"""

import requests
import numpy as np
import pandas as pd
import json
import os
import time
from datetime import datetime
from astropy.cosmology import FlatLambdaCDM

DATA_DIR = os.path.join(os.path.dirname(__file__), '..', '..', 'data', 'sdss')
os.makedirs(DATA_DIR, exist_ok=True)

cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
SDSS_URL = "https://skyserver.sdss.org/dr18/SkyServerWS/SearchTools/SqlSearch"


def query_sdss_json(sql, max_retries=3):
    """Execute SQL query and return DataFrame."""
    for attempt in range(max_retries):
        try:
            params = {'cmd': sql, 'format': 'json'}
            response = requests.get(SDSS_URL, params=params, timeout=300)
            
            if response.status_code == 200:
                data = response.json()
                if data and len(data) > 0 and 'Rows' in data[0]:
                    rows = data[0]['Rows']
                    if len(rows) > 0:
                        return pd.DataFrame(rows)
                return None
            else:
                print(f"    HTTP {response.status_code}")
                
        except requests.exceptions.Timeout:
            print(f"    Timeout (attempt {attempt + 1})")
        except Exception as e:
            print(f"    Error: {e}")
        
        if attempt < max_retries - 1:
            time.sleep(3)
    
    return None


def download_galaxy_sample():
    """Download galaxies with spectroscopic properties."""
    print("=" * 70)
    print("DOWNLOADING SDSS DR18 GALAXIES VIA SQL")
    print("=" * 70)
    print(f"Started: {datetime.now().isoformat()}")
    
    # Query in redshift chunks to avoid timeout
    z_ranges = [
        (0.01, 0.05, 80000),
        (0.05, 0.10, 80000),
        (0.10, 0.20, 80000),
        (0.20, 0.35, 80000),
        (0.35, 0.55, 50000),
        (0.55, 0.75, 30000),
    ]
    
    all_data = []
    
    for z_min, z_max, limit in z_ranges:
        print(f"\nQuerying z = {z_min:.2f} - {z_max:.2f} (limit {limit})...")
        
        sql = f"""
        SELECT TOP {limit}
            s.specobjid,
            s.ra, s.dec,
            s.z as redshift,
            s.zerr as z_err,
            s.veldisp,
            s.veldisperr as veldisp_err,
            p.petroMag_g, p.petroMag_r, p.petroMag_i,
            p.petroR50_r, p.petroR90_r,
            p.fracDeV_r,
            p.expAB_r, p.deVAB_r,
            p.extinction_r
        FROM SpecObj s
        JOIN PhotoObj p ON s.bestobjid = p.objid
        WHERE s.class = 'GALAXY'
            AND s.z BETWEEN {z_min} AND {z_max}
            AND s.zWarning = 0
            AND s.veldisp > 20
            AND s.veldisp < 500
            AND s.veldisperr > 0
            AND s.veldisperr < 100
            AND p.petroMag_r BETWEEN 12 AND 22
        ORDER BY s.z
        """
        
        df = query_sdss_json(sql)
        
        if df is not None and len(df) > 0:
            print(f"  Retrieved {len(df)} galaxies")
            all_data.append(df)
        else:
            print(f"  No data retrieved")
    
    if len(all_data) == 0:
        print("\nERROR: No data retrieved")
        return None
    
    combined = pd.concat(all_data, ignore_index=True)
    print(f"\nTotal galaxies: {len(combined)}")
    
    # Remove duplicates
    combined = combined.drop_duplicates(subset=['specobjid'])
    print(f"After deduplication: {len(combined)}")
    
    return combined


def compute_properties(df):
    """Compute derived properties for temporal onion test."""
    print("\nComputing derived properties...")
    
    z = df['redshift'].values
    
    # Lookback time
    df['t_lookback'] = cosmo.lookback_time(z).value
    
    # Log velocity dispersion
    df['log_sigma'] = np.log10(df['veldisp'].clip(lower=20))
    
    # Stellar mass estimate from g-r color and velocity dispersion
    # Using Faber-Jackson + color correction
    g_r = df['petroMag_g'] - df['petroMag_r']
    df['log_mass'] = 10.3 + 4 * (df['log_sigma'] - 2.3) + 0.5 * (g_r - 0.7)
    df['log_mass'] = df['log_mass'].clip(8, 13)
    
    # SFR proxy from g-r color (bluer = more SF)
    df['log_sfr'] = 1.5 - 1.5 * g_r.clip(0, 1.5)
    
    # Concentration (morphology proxy)
    df['concentration'] = df['petroR90_r'] / df['petroR50_r'].replace(0, np.nan)
    
    # Sersic proxy from fracDeV
    df['sersic_proxy'] = 1 + 3 * df['fracDeV_r'].fillna(0.5).clip(0, 1)
    
    # Axis ratio
    df['axis_ratio'] = np.where(
        df['fracDeV_r'].fillna(0.5) > 0.5,
        df['deVAB_r'],
        df['expAB_r']
    ).clip(0.1, 1.0)
    
    # Filter valid
    valid = (
        np.isfinite(df['log_mass']) &
        np.isfinite(df['log_sigma']) &
        np.isfinite(df['concentration']) &
        (df['log_sigma'] > 1.3) & (df['log_sigma'] < 2.8) &
        (df['concentration'] > 1) & (df['concentration'] < 10)
    )
    
    df_valid = df[valid].copy()
    print(f"  Valid galaxies: {len(df_valid)}")
    
    return df_valid


def summarize(df):
    """Print summary."""
    print("\n" + "=" * 70)
    print("SAMPLE SUMMARY")
    print("=" * 70)
    
    print(f"\nTotal: {len(df):,} galaxies")
    
    print(f"\nRedshift distribution:")
    z_bins = [0.01, 0.05, 0.10, 0.20, 0.35, 0.55, 0.75]
    for i in range(len(z_bins) - 1):
        mask = (df['redshift'] >= z_bins[i]) & (df['redshift'] < z_bins[i+1])
        count = mask.sum()
        if count > 0:
            t = df.loc[mask, 't_lookback']
            print(f"  z = {z_bins[i]:.2f}-{z_bins[i+1]:.2f}: {count:,} "
                  f"(t = {t.min():.1f}-{t.max():.1f} Gyr)")
    
    print(f"\nLookback time: {df['t_lookback'].min():.1f} - {df['t_lookback'].max():.1f} Gyr")
    print(f"\nCompare to MaNGA: 10,000 galaxies, z<0.15, t<2 Gyr")
    print(f"SDSS advantage: {len(df)/10000:.0f}x galaxies, "
          f"{df['t_lookback'].max()/2:.1f}x deeper")


def main():
    """Main pipeline."""
    df_raw = download_galaxy_sample()
    
    if df_raw is None:
        return None
    
    # Save raw
    raw_path = os.path.join(DATA_DIR, 'sdss_galaxies_raw.csv')
    df_raw.to_csv(raw_path, index=False)
    print(f"\nRaw data saved: {raw_path}")
    
    # Process
    df = compute_properties(df_raw)
    
    # Summarize
    summarize(df)
    
    # Save processed
    output_path = os.path.join(DATA_DIR, 'sdss_galaxies.csv')
    df.to_csv(output_path, index=False)
    print(f"\nProcessed data saved: {output_path}")
    
    # Metadata
    meta = {
        'timestamp': datetime.now().isoformat(),
        'source': 'SDSS DR18 SQL',
        'n_galaxies': len(df),
        'z_range': [float(df['redshift'].min()), float(df['redshift'].max())],
        't_lookback_range': [float(df['t_lookback'].min()), float(df['t_lookback'].max())],
    }
    with open(os.path.join(DATA_DIR, 'sdss_metadata.json'), 'w') as f:
        json.dump(meta, f, indent=2)
    
    return df


if __name__ == '__main__':
    df = main()
