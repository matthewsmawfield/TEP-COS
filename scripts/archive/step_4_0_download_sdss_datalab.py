#!/usr/bin/env python3
"""
Step 4.0: Download SDSS Galaxy Data via NOIRLab Data Lab

Uses the NOIRLab Astro Data Lab to query SDSS DR17 eBOSS-DAP catalog,
which contains ~2 million galaxies with spectroscopic measurements.

Author: M. Smawfield
Date: January 2026
"""

import numpy as np
import pandas as pd
import json
import os
from datetime import datetime

DATA_DIR = os.path.join(os.path.dirname(__file__), '..', '..', 'data', 'sdss')
os.makedirs(DATA_DIR, exist_ok=True)

# Try to use Data Lab query client
try:
    from dl import queryClient as qc
    DATALAB_AVAILABLE = True
except ImportError:
    DATALAB_AVAILABLE = False
    print("NOIRLab Data Lab not installed. Installing...")


def install_datalab():
    """Install the Data Lab client."""
    import subprocess
    subprocess.check_call(['pip', 'install', 'astro-datalab'])


def query_datalab(sql):
    """Execute SQL query against NOIRLab Data Lab."""
    from dl import queryClient as qc
    
    print(f"  Executing query...")
    result = qc.query(sql, fmt='csv')
    
    from io import StringIO
    df = pd.read_csv(StringIO(result))
    return df


def download_eboss_sample():
    """
    Download eBOSS galaxy sample with derived properties.
    
    The eBOSS-DAP catalog contains:
    - Emission line measurements
    - Stellar kinematics (velocity dispersion)
    - Spectral indices
    - Redshifts up to z ~ 1.1
    """
    print("\n" + "=" * 70)
    print("DOWNLOADING SDSS eBOSS GALAXY CATALOG")
    print("=" * 70)
    print(f"Started: {datetime.now().isoformat()}")
    
    # Query for galaxies with spectroscopic data
    # Join eBOSS-DAP with main SpecObj for coordinates and redshift
    
    z_ranges = [
        (0.01, 0.10),
        (0.10, 0.25),
        (0.25, 0.50),
        (0.50, 0.80),
    ]
    
    all_data = []
    
    for z_min, z_max in z_ranges:
        print(f"\nQuerying z = {z_min:.2f} - {z_max:.2f}...")
        
        # Query eBOSS-DAP spectral indices for velocity dispersion
        # and join with main catalog for positions
        sql = f"""
        SELECT 
            s.specobjid,
            s.ra, s.dec,
            s.z as redshift,
            s.zerr as redshift_err,
            s.vdisp as veldisp,
            s.vdisp_err as veldisp_err,
            p.petroMag_u, p.petroMag_g, p.petroMag_r, p.petroMag_i, p.petroMag_z,
            p.petroR50_r,
            p.petroR90_r,
            p.fracdev_r,
            p.expab_r,
            p.devab_r,
            p.extinction_r
        FROM sdss_dr17.specobj s
        JOIN sdss_dr17.photoobj p ON s.bestobjid = p.objid
        WHERE s.class = 'GALAXY'
            AND s.z BETWEEN {z_min} AND {z_max}
            AND s.zwarning = 0
            AND s.vdisp > 30
            AND s.vdisp < 500
            AND s.vdisp_err > 0
            AND s.vdisp_err < 50
            AND p.petromag_r BETWEEN 12 AND 21
        LIMIT 150000
        """
        
        try:
            df = query_datalab(sql)
            if df is not None and len(df) > 0:
                print(f"  Retrieved {len(df)} galaxies")
                all_data.append(df)
            else:
                print(f"  No data retrieved")
        except Exception as e:
            print(f"  Error: {e}")
    
    if len(all_data) == 0:
        print("\nNo data retrieved. Trying alternative query...")
        return None
    
    # Combine all chunks
    combined = pd.concat(all_data, ignore_index=True)
    print(f"\nTotal galaxies retrieved: {len(combined)}")
    
    # Remove duplicates
    combined = combined.drop_duplicates(subset=[combined.columns[0]])
    print(f"After deduplication: {len(combined)}")
    
    return combined


def download_mpa_jhu_sample():
    """
    Alternative: Download MPA-JHU galaxy properties from DR8.
    
    This catalog has stellar masses and SFR estimates.
    """
    print("\n" + "=" * 70)
    print("DOWNLOADING SDSS MPA-JHU GALAXY CATALOG")
    print("=" * 70)
    print(f"Started: {datetime.now().isoformat()}")
    
    z_ranges = [
        (0.01, 0.08),
        (0.08, 0.15),
        (0.15, 0.25),
        (0.25, 0.40),
    ]
    
    all_data = []
    
    for z_min, z_max in z_ranges:
        print(f"\nQuerying z = {z_min:.2f} - {z_max:.2f}...")
        
        sql = f"""
        SELECT 
            g.specobjid,
            g.ra, g.dec,
            g.z as redshift,
            g.z_err as redshift_err,
            g.v_disp as veldisp,
            g.v_disp_err as veldisp_err,
            e.bptclass,
            e.lgm_tot_p50 as log_mass,
            e.sfr_tot_p50 as log_sfr,
            e.specsfr_tot_p50 as log_ssfr
        FROM sdss_dr8.galspecinfo g
        JOIN sdss_dr8.galspecextra e ON g.specobjid = e.specobjid
        WHERE g.z BETWEEN {z_min} AND {z_max}
            AND g.reliable = 1
            AND g.v_disp > 30
            AND g.v_disp < 500
            AND e.lgm_tot_p50 > 8
            AND e.lgm_tot_p50 < 13
        LIMIT 150000
        """
        
        try:
            df = query_datalab(sql)
            if df is not None and len(df) > 0:
                print(f"  Retrieved {len(df)} galaxies")
                all_data.append(df)
            else:
                print(f"  No data retrieved")
        except Exception as e:
            print(f"  Error: {e}")
    
    if len(all_data) == 0:
        return None
    
    combined = pd.concat(all_data, ignore_index=True)
    print(f"\nTotal galaxies retrieved: {len(combined)}")
    
    combined = combined.drop_duplicates(subset=[combined.columns[0]])
    print(f"After deduplication: {len(combined)}")
    
    return combined


def compute_derived_properties(df):
    """Compute derived properties."""
    print("\nComputing derived properties...")
    
    from astropy.cosmology import FlatLambdaCDM
    cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
    
    z = df['redshift'].values
    
    # Lookback time
    df['t_lookback'] = cosmo.lookback_time(z).value
    
    # Log velocity dispersion
    if 'veldisp' in df.columns:
        df['log_sigma'] = np.log10(df['veldisp'].clip(lower=10))
    
    # If we don't have stellar mass, estimate from velocity dispersion
    if 'log_mass' not in df.columns and 'log_sigma' in df.columns:
        # Faber-Jackson relation: M ~ sigma^4
        df['log_mass'] = 10.5 + 4 * (df['log_sigma'] - 2.3)
    
    # If we don't have SFR, use a placeholder
    if 'log_sfr' not in df.columns:
        # Assume main sequence relation as default
        df['log_sfr'] = 0.0  # Will be updated if we have color info
    
    # Morphology proxy from concentration if available
    if 'petroR90_r' in df.columns and 'petroR50_r' in df.columns:
        df['concentration'] = df['petroR90_r'] / df['petroR50_r'].replace(0, np.nan)
        df['sersic_proxy'] = np.where(df['concentration'] > 2.6, 4.0, 1.5)
    elif 'fracdev_r' in df.columns:
        df['sersic_proxy'] = 1 + 3 * df['fracdev_r'].fillna(0.5)
    else:
        df['sersic_proxy'] = 2.5  # Default
    
    # Axis ratio
    if 'expab_r' in df.columns and 'devab_r' in df.columns:
        df['axis_ratio'] = np.where(
            df.get('fracdev_r', 0.5) > 0.5,
            df['devab_r'],
            df['expab_r']
        ).clip(0.1, 1.0)
    else:
        df['axis_ratio'] = 0.7  # Default
    
    # Filter valid
    valid_mask = (
        np.isfinite(df['log_mass']) &
        np.isfinite(df['log_sigma']) &
        (df['log_mass'] > 8) & (df['log_mass'] < 13) &
        (df['log_sigma'] > 1) & (df['log_sigma'] < 3)
    )
    
    df_valid = df[valid_mask].copy()
    print(f"  Valid galaxies: {len(df_valid)}")
    
    return df_valid


def summarize_sample(df):
    """Print summary statistics."""
    print("\n" + "=" * 70)
    print("SAMPLE SUMMARY")
    print("=" * 70)
    
    print(f"\nTotal galaxies: {len(df)}")
    
    print(f"\nRedshift distribution:")
    z_bins = [0.01, 0.10, 0.25, 0.50, 0.80]
    for i in range(len(z_bins) - 1):
        count = ((df['redshift'] >= z_bins[i]) & (df['redshift'] < z_bins[i+1])).sum()
        print(f"  z = {z_bins[i]:.2f} - {z_bins[i+1]:.2f}: {count:,} galaxies")
    
    print(f"\nLookback time range: {df['t_lookback'].min():.1f} - {df['t_lookback'].max():.1f} Gyr")
    
    print(f"\nProperty ranges:")
    print(f"  log(M*/M☉): {df['log_mass'].min():.1f} - {df['log_mass'].max():.1f}")
    print(f"  log(σ/km/s): {df['log_sigma'].min():.2f} - {df['log_sigma'].max():.2f}")


def main():
    """Main download pipeline."""
    global DATALAB_AVAILABLE
    
    if not DATALAB_AVAILABLE:
        install_datalab()
        from dl import queryClient as qc
    
    # Try MPA-JHU catalog first (has stellar masses)
    df_raw = download_mpa_jhu_sample()
    
    if df_raw is None or len(df_raw) < 1000:
        print("\nMPA-JHU query failed, trying eBOSS...")
        df_raw = download_eboss_sample()
    
    if df_raw is None:
        print("All downloads failed.")
        return None
    
    # Save raw data
    raw_path = os.path.join(DATA_DIR, 'sdss_galaxies_raw.csv')
    df_raw.to_csv(raw_path, index=False)
    print(f"\nRaw data saved: {raw_path}")
    
    # Compute derived properties
    df = compute_derived_properties(df_raw)
    
    # Summarize
    summarize_sample(df)
    
    # Save processed data
    output_path = os.path.join(DATA_DIR, 'sdss_galaxies.csv')
    df.to_csv(output_path, index=False)
    print(f"\nProcessed data saved: {output_path}")
    
    # Save metadata
    metadata = {
        'timestamp': datetime.now().isoformat(),
        'source': 'SDSS via NOIRLab Data Lab',
        'n_galaxies': len(df),
        'z_range': [float(df['redshift'].min()), float(df['redshift'].max())],
        't_lookback_range': [float(df['t_lookback'].min()), float(df['t_lookback'].max())],
        'columns': list(df.columns),
    }
    
    meta_path = os.path.join(DATA_DIR, 'sdss_metadata.json')
    with open(meta_path, 'w') as f:
        json.dump(metadata, f, indent=2)
    
    return df


if __name__ == '__main__':
    df = main()
