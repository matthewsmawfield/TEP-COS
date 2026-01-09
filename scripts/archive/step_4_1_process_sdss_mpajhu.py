#!/usr/bin/env python3
"""
Step 4.1: Process SDSS MPA-JHU Galaxy Catalog

Processes the downloaded MPA-JHU DR7 catalog files to create a unified
galaxy sample for the temporal onion test.

Files used:
- gal_info_dr7_v5_2.fit.gz: Basic info (RA, Dec, z, velocity dispersion)
- totlgm_dr7_v5_2.fit.gz: Stellar masses
- gal_totsfr_dr7_v5_2.fits.gz: Star formation rates

Author: M. Smawfield
Date: January 2026
"""

import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.cosmology import FlatLambdaCDM
import json
import os
from datetime import datetime

DATA_DIR = os.path.join(os.path.dirname(__file__), '..', '..', 'data', 'sdss')
RESULTS_DIR = os.path.join(os.path.dirname(__file__), '..', '..', 'results', 'outputs', 'topology')

cosmo = FlatLambdaCDM(H0=70, Om0=0.3)


def load_mpajhu_catalog():
    """Load and merge MPA-JHU catalog files."""
    print("Loading MPA-JHU DR7 catalog files...")
    
    # Load galaxy info (positions, redshifts, velocity dispersions)
    info_path = os.path.join(DATA_DIR, 'gal_info_dr7_v5_2.fit.gz')
    with fits.open(info_path) as hdul:
        info = hdul[1].data
        print(f"  gal_info: {len(info)} galaxies")
        print(f"  Columns: {info.columns.names[:10]}...")
    
    # Load stellar masses
    mass_path = os.path.join(DATA_DIR, 'totlgm_dr7_v5_2.fit.gz')
    with fits.open(mass_path) as hdul:
        mass = hdul[1].data
        print(f"  totlgm: {len(mass)} entries")
        print(f"  Columns: {mass.columns.names}")
    
    # Load star formation rates
    sfr_path = os.path.join(DATA_DIR, 'gal_totsfr_dr7_v5_2.fits.gz')
    with fits.open(sfr_path) as hdul:
        sfr = hdul[1].data
        print(f"  gal_totsfr: {len(sfr)} entries")
        print(f"  Columns: {sfr.columns.names}")
    
    # Extract relevant columns
    n = len(info)
    
    galaxies = {
        'ra': info['RA'],
        'dec': info['DEC'],
        'z': info['Z'],
        'z_err': info['Z_ERR'],
        'veldisp': info['V_DISP'],
        'veldisp_err': info['V_DISP_ERR'],
        'sn_median': info['SN_MEDIAN'],
    }
    
    # Add stellar mass (median of posterior)
    # Column is typically 'AVG' or 'P50' for median
    if 'P50' in mass.columns.names:
        galaxies['log_mass'] = mass['P50']
    elif 'AVG' in mass.columns.names:
        galaxies['log_mass'] = mass['AVG']
    else:
        # Use first column after checking
        print(f"  Mass columns: {mass.columns.names}")
        galaxies['log_mass'] = mass[mass.columns.names[0]]
    
    # Add SFR (median of posterior)
    if 'P50' in sfr.columns.names:
        galaxies['log_sfr'] = sfr['P50']
    elif 'AVG' in sfr.columns.names:
        galaxies['log_sfr'] = sfr['AVG']
    else:
        print(f"  SFR columns: {sfr.columns.names}")
        galaxies['log_sfr'] = sfr[sfr.columns.names[0]]
    
    df = pd.DataFrame(galaxies)
    print(f"\nCombined catalog: {len(df)} galaxies")
    
    return df


def filter_and_process(df):
    """Filter to valid galaxies and compute derived properties."""
    print("\nFiltering and processing...")
    
    # Basic quality cuts
    valid_mask = (
        (df['z'] > 0.01) & (df['z'] < 0.7) &
        (df['z_err'] > 0) & (df['z_err'] < 0.01) &
        (df['veldisp'] > 30) & (df['veldisp'] < 500) &
        (df['veldisp_err'] > 0) & (df['veldisp_err'] < 50) &
        (df['log_mass'] > 8) & (df['log_mass'] < 13) &
        (df['log_sfr'] > -5) & (df['log_sfr'] < 5) &
        (df['sn_median'] > 5)
    )
    
    df_valid = df[valid_mask].copy()
    print(f"  After quality cuts: {len(df_valid)} galaxies")
    
    # Compute derived properties
    df_valid['log_sigma'] = np.log10(df_valid['veldisp'])
    df_valid['t_lookback'] = cosmo.lookback_time(df_valid['z'].values).value
    
    # Morphology proxy from mass-sigma relation residual
    # Ellipticals lie above the relation, disks below
    expected_sigma = 2.0 + 0.25 * (df_valid['log_mass'] - 10.5)
    sigma_residual = df_valid['log_sigma'] - expected_sigma
    df_valid['sersic_proxy'] = 2.5 + 2 * sigma_residual  # Maps to ~1-4
    df_valid['sersic_proxy'] = df_valid['sersic_proxy'].clip(0.5, 6)
    
    # Axis ratio proxy (not available, use placeholder)
    df_valid['axis_ratio'] = 0.7
    
    return df_valid


def summarize_sample(df):
    """Print summary statistics."""
    print("\n" + "=" * 70)
    print("SDSS MPA-JHU SAMPLE SUMMARY")
    print("=" * 70)
    
    print(f"\nTotal galaxies: {len(df):,}")
    
    print(f"\nRedshift distribution:")
    z_bins = [0.01, 0.05, 0.10, 0.15, 0.25, 0.40, 0.70]
    for i in range(len(z_bins) - 1):
        mask = (df['z'] >= z_bins[i]) & (df['z'] < z_bins[i+1])
        count = mask.sum()
        if count > 0:
            t_range = df.loc[mask, 't_lookback']
            print(f"  z = {z_bins[i]:.2f} - {z_bins[i+1]:.2f}: {count:,} galaxies "
                  f"(t = {t_range.min():.1f} - {t_range.max():.1f} Gyr)")
    
    print(f"\nLookback time range: {df['t_lookback'].min():.1f} - {df['t_lookback'].max():.1f} Gyr")
    
    print(f"\nProperty ranges:")
    print(f"  log(M*/M☉): {df['log_mass'].min():.1f} - {df['log_mass'].max():.1f}")
    print(f"  log(σ/km/s): {df['log_sigma'].min():.2f} - {df['log_sigma'].max():.2f}")
    print(f"  log(SFR): {df['log_sfr'].min():.1f} - {df['log_sfr'].max():.1f}")
    
    # Compare to MaNGA
    print(f"\n" + "=" * 70)
    print("COMPARISON TO MaNGA")
    print("=" * 70)
    print(f"  MaNGA: ~10,000 galaxies, z < 0.15, t < 2 Gyr")
    print(f"  SDSS:  {len(df):,} galaxies, z < 0.7, t < {df['t_lookback'].max():.1f} Gyr")
    print(f"  Improvement: {len(df)/10000:.0f}x more galaxies, "
          f"{df['t_lookback'].max()/2:.1f}x deeper in time")


def main():
    """Main processing pipeline."""
    print("=" * 70)
    print("PROCESSING SDSS MPA-JHU CATALOG FOR TEMPORAL ONION TEST")
    print("=" * 70)
    print(f"Started: {datetime.now().isoformat()}")
    
    # Load catalog
    df_raw = load_mpajhu_catalog()
    
    # Filter and process
    df = filter_and_process(df_raw)
    
    # Summarize
    summarize_sample(df)
    
    # Save processed data
    output_path = os.path.join(DATA_DIR, 'sdss_mpajhu_processed.csv')
    df.to_csv(output_path, index=False)
    print(f"\nProcessed data saved: {output_path}")
    
    # Save metadata
    metadata = {
        'timestamp': datetime.now().isoformat(),
        'source': 'SDSS MPA-JHU DR7',
        'n_galaxies': len(df),
        'z_range': [float(df['z'].min()), float(df['z'].max())],
        't_lookback_range': [float(df['t_lookback'].min()), float(df['t_lookback'].max())],
        'columns': list(df.columns),
    }
    
    meta_path = os.path.join(DATA_DIR, 'sdss_mpajhu_metadata.json')
    with open(meta_path, 'w') as f:
        json.dump(metadata, f, indent=2)
    print(f"Metadata saved: {meta_path}")
    
    return df


if __name__ == '__main__':
    df = main()
