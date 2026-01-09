#!/usr/bin/env python3
"""
Step 4.0: Download SDSS DR18 Spectroscopic Galaxy Catalog

Downloads galaxy data from SDSS DR18 for the temporal onion test.
Uses the SDSS SkyServer SQL interface to query the catalog.

Target: ~100,000-500,000 galaxies with:
- Spectroscopic redshifts (z < 0.7)
- Stellar mass estimates
- Velocity dispersion measurements
- Star formation indicators
- Morphological parameters

Author: M. Smawfield
Date: January 2026
"""

import requests
import pandas as pd
import numpy as np
import json
import os
import time
from datetime import datetime

DATA_DIR = os.path.join(os.path.dirname(__file__), '..', '..', 'data', 'sdss')
os.makedirs(DATA_DIR, exist_ok=True)

# SDSS SkyServer SQL endpoint
SDSS_URL = "https://skyserver.sdss.org/dr18/SkyServerWS/SearchTools/SqlSearch"


def query_sdss(sql, max_retries=3):
    """
    Execute SQL query against SDSS SkyServer.
    
    Returns pandas DataFrame with results.
    """
    # Use the CasJobs-style endpoint with proper CSV format
    url = "https://skyserver.sdss.org/dr18/SkyServerWS/SearchTools/SqlSearch"
    
    for attempt in range(max_retries):
        try:
            print(f"  Executing query (attempt {attempt + 1})...")
            
            # POST request with form data
            data = {
                'cmd': sql,
                'format': 'csv'
            }
            response = requests.post(url, data=data, timeout=300)
            
            if response.status_code == 200:
                # Parse CSV response - skip comment lines
                from io import StringIO
                lines = response.text.strip().split('\n')
                # Filter out comment lines and empty lines
                data_lines = [l for l in lines if l and not l.startswith('#')]
                if len(data_lines) < 2:
                    print(f"  No data rows returned")
                    return None
                csv_text = '\n'.join(data_lines)
                df = pd.read_csv(StringIO(csv_text))
                return df
            else:
                print(f"  Error: HTTP {response.status_code}")
                print(f"  Response: {response.text[:500]}")
                
        except requests.exceptions.Timeout:
            print(f"  Timeout on attempt {attempt + 1}")
        except Exception as e:
            print(f"  Error: {e}")
        
        if attempt < max_retries - 1:
            time.sleep(5)
    
    return None


def download_galaxy_sample():
    """
    Download a sample of galaxies with spectroscopic data.
    
    We query the SpecObj and Galaxy tables to get:
    - Spectroscopic redshift
    - Stellar velocity dispersion
    - Photometric properties for mass/SFR estimation
    """
    print("\n" + "=" * 70)
    print("DOWNLOADING SDSS DR18 GALAXY CATALOG")
    print("=" * 70)
    print(f"Started: {datetime.now().isoformat()}")
    
    # Query for galaxies with spectroscopic data
    # We'll download in redshift chunks to avoid timeout
    
    z_ranges = [
        (0.01, 0.05),
        (0.05, 0.10),
        (0.10, 0.15),
        (0.15, 0.25),
        (0.25, 0.40),
        (0.40, 0.70),
    ]
    
    all_data = []
    
    for z_min, z_max in z_ranges:
        print(f"\nQuerying z = {z_min:.2f} - {z_max:.2f}...")
        
        # SQL query for galaxy properties
        # Using SpecObj for spectroscopy and PhotoObj for photometry
        sql = f"""
        SELECT TOP 100000
            s.specObjID,
            s.ra, s.dec,
            s.z as redshift,
            s.zErr as redshift_err,
            s.velDisp,
            s.velDispErr,
            s.class,
            s.subClass,
            p.petroMag_u, p.petroMag_g, p.petroMag_r, p.petroMag_i, p.petroMag_z,
            p.petroMagErr_u, p.petroMagErr_g, p.petroMagErr_r, p.petroMagErr_i, p.petroMagErr_z,
            p.petroR50_r,
            p.petroR90_r,
            p.expAB_r,
            p.deVAB_r,
            p.fracDeV_r,
            p.extinction_r
        FROM SpecObj s
        JOIN PhotoObj p ON s.bestObjID = p.objID
        WHERE s.class = 'GALAXY'
            AND s.z BETWEEN {z_min} AND {z_max}
            AND s.zWarning = 0
            AND s.velDisp > 0
            AND s.velDispErr > 0
            AND s.velDispErr < 50
            AND p.petroMag_r BETWEEN 10 AND 22
            AND p.petroMag_r > 0
        ORDER BY s.z
        """
        
        df = query_sdss(sql)
        
        if df is not None and len(df) > 0:
            print(f"  Retrieved {len(df)} galaxies")
            all_data.append(df)
        else:
            print(f"  No data retrieved for this range")
    
    if len(all_data) == 0:
        print("\nERROR: No data retrieved from SDSS")
        return None
    
    # Combine all chunks
    combined = pd.concat(all_data, ignore_index=True)
    print(f"\nTotal galaxies retrieved: {len(combined)}")
    
    # Check column names (SDSS may return lowercase)
    print(f"Columns: {list(combined.columns)[:5]}...")
    
    # Remove duplicates - use first column as ID
    id_col = combined.columns[0]
    combined = combined.drop_duplicates(subset=[id_col])
    print(f"After deduplication: {len(combined)}")
    
    # Save raw data
    output_path = os.path.join(DATA_DIR, 'sdss_dr18_galaxies_raw.csv')
    combined.to_csv(output_path, index=False)
    print(f"\nRaw data saved: {output_path}")
    
    return combined


def compute_derived_properties(df):
    """
    Compute derived properties from SDSS photometry.
    
    - Stellar mass from g-r color and r-band luminosity
    - Star formation rate proxy from u-r color
    - Concentration index (R90/R50) as morphology proxy
    """
    print("\nComputing derived properties...")
    
    # Cosmology for luminosity distance
    from astropy.cosmology import FlatLambdaCDM
    import astropy.units as u
    cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
    
    # Absolute magnitude
    z = df['redshift'].values
    dist_mod = cosmo.distmod(z).value
    
    # Extinction-corrected apparent magnitude
    app_mag_r = df['petroMag_r'] - df['extinction_r']
    
    # Absolute r-band magnitude
    abs_mag_r = app_mag_r - dist_mod
    
    # Stellar mass estimate using Bell et al. (2003) relation
    # log(M/L_r) = -0.306 + 1.097*(g-r)
    g_r = (df['petroMag_g'] - df['petroMag_r'])
    log_ml = -0.306 + 1.097 * g_r
    
    # Solar absolute magnitude in r-band
    M_sun_r = 4.65
    
    # Luminosity in solar units
    log_L = (M_sun_r - abs_mag_r) / 2.5
    
    # Stellar mass
    log_mass = log_ml + log_L
    
    # SFR proxy: u-r color (bluer = more star formation)
    u_r = df['petroMag_u'] - df['petroMag_r']
    # Convert to approximate log(SFR) using empirical relation
    # This is a rough proxy, not a precise measurement
    log_sfr = 2.0 - 0.5 * u_r  # Rough calibration
    log_sfr = np.clip(log_sfr, -3, 3)
    
    # Concentration index (morphology proxy)
    # R90/R50 > 2.6 typically indicates early-type (elliptical)
    concentration = df['petroR90_r'] / df['petroR50_r'].replace(0, np.nan)
    
    # Sersic-like index from fracDeV (fraction of light in de Vaucouleurs profile)
    # fracDeV ~ 1 means elliptical-like, fracDeV ~ 0 means disk-like
    sersic_proxy = 1 + 3 * df['fracDeV_r']  # Maps 0-1 to 1-4
    
    # Axis ratio (use exponential or deV depending on fracDeV)
    axis_ratio = np.where(
        df['fracDeV_r'] > 0.5,
        df['deVAB_r'],
        df['expAB_r']
    )
    
    # Add to dataframe
    df['log_mass'] = log_mass
    df['log_sigma'] = np.log10(df['velDisp'])
    df['log_sfr'] = log_sfr
    df['concentration'] = concentration
    df['sersic_proxy'] = sersic_proxy
    df['axis_ratio'] = axis_ratio
    df['abs_mag_r'] = abs_mag_r
    
    # Lookback time
    df['t_lookback'] = cosmo.lookback_time(z).value
    
    # Filter valid entries
    valid_mask = (
        (df['log_mass'] > 8) & (df['log_mass'] < 13) &
        (df['log_sigma'] > 1) & (df['log_sigma'] < 3) &
        np.isfinite(df['concentration']) &
        (df['concentration'] > 1) & (df['concentration'] < 10) &
        np.isfinite(df['axis_ratio']) &
        (df['axis_ratio'] > 0.1) & (df['axis_ratio'] < 1)
    )
    
    df_valid = df[valid_mask].copy()
    print(f"  Valid galaxies after filtering: {len(df_valid)}")
    
    return df_valid


def summarize_sample(df):
    """Print summary statistics of the sample."""
    print("\n" + "=" * 70)
    print("SAMPLE SUMMARY")
    print("=" * 70)
    
    print(f"\nTotal galaxies: {len(df)}")
    
    print(f"\nRedshift distribution:")
    z_bins = [0.01, 0.05, 0.10, 0.15, 0.25, 0.40, 0.70]
    for i in range(len(z_bins) - 1):
        count = ((df['redshift'] >= z_bins[i]) & (df['redshift'] < z_bins[i+1])).sum()
        t_lo = df[df['redshift'] >= z_bins[i]]['t_lookback'].min()
        t_hi = df[df['redshift'] < z_bins[i+1]]['t_lookback'].max()
        print(f"  z = {z_bins[i]:.2f} - {z_bins[i+1]:.2f}: {count:,} galaxies")
    
    print(f"\nProperty ranges:")
    print(f"  log(M*/M☉): {df['log_mass'].min():.1f} - {df['log_mass'].max():.1f}")
    print(f"  log(σ/km/s): {df['log_sigma'].min():.2f} - {df['log_sigma'].max():.2f}")
    print(f"  log(SFR): {df['log_sfr'].min():.1f} - {df['log_sfr'].max():.1f}")
    print(f"  Concentration: {df['concentration'].min():.1f} - {df['concentration'].max():.1f}")
    
    print(f"\nLookback time range: {df['t_lookback'].min():.1f} - {df['t_lookback'].max():.1f} Gyr")


def main():
    """Main download pipeline."""
    # Download raw data
    df_raw = download_galaxy_sample()
    
    if df_raw is None:
        print("Download failed. Please check network connection.")
        return None
    
    # Compute derived properties
    df = compute_derived_properties(df_raw)
    
    # Summarize
    summarize_sample(df)
    
    # Save processed data
    output_path = os.path.join(DATA_DIR, 'sdss_dr18_galaxies.csv')
    df.to_csv(output_path, index=False)
    print(f"\nProcessed data saved: {output_path}")
    
    # Save metadata
    metadata = {
        'timestamp': datetime.now().isoformat(),
        'source': 'SDSS DR18 SkyServer',
        'n_galaxies': len(df),
        'z_range': [float(df['redshift'].min()), float(df['redshift'].max())],
        't_lookback_range': [float(df['t_lookback'].min()), float(df['t_lookback'].max())],
        'columns': list(df.columns),
    }
    
    meta_path = os.path.join(DATA_DIR, 'sdss_dr18_metadata.json')
    with open(meta_path, 'w') as f:
        json.dump(metadata, f, indent=2)
    print(f"Metadata saved: {meta_path}")
    
    return df


if __name__ == '__main__':
    df = main()
