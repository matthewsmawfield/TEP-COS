#!/usr/bin/env python3
"""
Step 6.73: SDSS Test CW - Stellar Twins (The Identical Twin Paradox)

Hypothesis:
Two stars with identical spectroscopic parameters (Teff, log g, [Fe/H]) are "twins" and should 
have the same intrinsic luminosity. In TEP, if one twin is in the Bulge (deep Phi) and the other is local, 
the locally measured photon rate (luminosity) may differ due to time dilation of the emitted flux 
or metric effects on the photosphere.

Prediction:
Inner Galaxy twins appear fainter (K_mag_abs) than Local twins.
Delta Mag = M_K(Inner) - M_K(Local) > 0 (Fainter = larger magnitude).

Data:
- aspcapStar: param_teff, param_logg, param_m_h (Spectroscopic params)
- apogeeObject: k (K-band magnitude), ak_targ (Extinction)
- apogee_starhorse: dist50 (Distance)

Method:
1. Select Red Giant Branch stars (reliable params, luminous enough to see in Bulge).
2. Bin stars in Teff, logg, [M/H] space (e.g., 50K, 0.1 dex, 0.1 dex).
3. In each bin, calculate Mean Absolute K Magnitude for Local (< 2 kpc) and Inner Galaxy (> 4 kpc, |l|<30) populations.
   M_K = k - 5 * log10(dist * 1000) + 5 - A_K
4. Calculate weighted mean difference across all bins.
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

def download_data(limit=100):
    print(f"Querying SDSS for Test CW (Limit: {limit})...")
    
    # We need to join aspcapStar, apogeeObject, and apogee_starhorse.
    # apogee_starhorse is the limiting factor usually (VAC).
    # apogeeObject contains photometry.
    
    sql = f"""
    SELECT TOP {limit}
        a.apogee_id,
        a.param_teff as teff,
        a.param_logg as logg,
        a.param_m_h as fe_h,
        a.glon, a.glat,
        o.k,
        o.ak_targ as ak,
        s.dist50 as dist
        
    FROM aspcapStar a
    JOIN apogeeObject o ON a.apogee_id = o.apogee_id
    JOIN apogee_starhorse s ON a.apogee_id = s.apogee_id
    
    WHERE 
        a.param_teff BETWEEN 3500 AND 5000 -- Red Giants
        AND a.param_logg BETWEEN 0.5 AND 3.5
        AND s.dist50 > 0
        AND o.k > 0 AND o.k < 14
        AND abs(a.glat) > 5 -- Avoid extreme extinction in plane if possible
        AND o.ak_targ > -0.1
    """
    return query_sdss(sql)

def analyze_twins(df):
    print("Analyzing Stellar Twins...")
    
    if df is None or len(df) < 100:
        print("  Insufficient data.")
        return None
        
    # Clean
    df = df.dropna().copy()
    
    # Calculate Absolute Magnitude
    # M = m - 5log(d) + 5 - A
    # dist is in kpc. log10(d*1000) = log10(d) + 3
    # 5(log d + 3) - 5 = 5log d + 15 - 5 = 5log d + 10
    # M = m - 5log10(dist_kpc) - 10 - A
    
    df['M_K'] = df['k'] - 5 * np.log10(df['dist']) - 10 - df['ak']
    
    print(f"  Sample size: {len(df)}")
    
    # Define Regions
    # Local: dist < 2 kpc
    # Inner: dist > 4 kpc AND |glon| < 40 (Towards Bulge)
    # Outer: dist > 4 kpc AND |glon| > 90? (For control)
    
    df['region'] = 'Other'
    df.loc[df['dist'] < 2.0, 'region'] = 'Local'
    df.loc[(df['dist'] > 4.0) & (df['glon'].abs() < 40), 'region'] = 'Inner'
    df.loc[(df['dist'] > 4.0) & (df['glon'].abs() > 90) & (df['glon'].abs() < 270), 'region'] = 'Outer' # Anti-center
    
    print("  Region counts:")
    print(df['region'].value_counts())
    
    local = df[df['region'] == 'Local'].copy()
    inner = df[df['region'] == 'Inner'].copy()
    
    if len(local) < 10 or len(inner) < 10:
        print("  Not enough stars in regions to compare.")
        return None
        
    # Binning Strategy
    # We want to compare M_K for stars with same (teff, logg, fe_h)
    # Create multi-dimensional bins
    
    teff_bins = np.arange(3500, 5100, 100) # 100K bins
    logg_bins = np.arange(0.5, 3.6, 0.2)   # 0.2 dex bins
    feh_bins = np.arange(-2.0, 0.6, 0.2)   # 0.2 dex bins
    
    local['teff_bin'] = pd.cut(local['teff'], teff_bins)
    local['logg_bin'] = pd.cut(local['logg'], logg_bins)
    local['feh_bin'] = pd.cut(local['fe_h'], feh_bins)
    
    inner['teff_bin'] = pd.cut(inner['teff'], teff_bins)
    inner['logg_bin'] = pd.cut(inner['logg'], logg_bins)
    inner['feh_bin'] = pd.cut(inner['fe_h'], feh_bins)
    
    # Group and Aggregate
    cols = ['teff_bin', 'logg_bin', 'feh_bin']
    
    local_grouped = local.groupby(cols, observed=True)['M_K'].agg(['mean', 'count', 'std'])
    inner_grouped = inner.groupby(cols, observed=True)['M_K'].agg(['mean', 'count', 'std'])
    
    # Merge
    merged = pd.merge(local_grouped, inner_grouped, on=cols, suffixes=('_loc', '_inn'))
    
    # Filter for robust comparison (at least N stars in each bin)
    valid = merged[(merged['count_loc'] >= 3) & (merged['count_inn'] >= 3)].copy()
    
    print(f"  Valid Twin Bins: {len(valid)}")
    
    if len(valid) == 0:
        print("  No overlapping bins with sufficient stars.")
        return None
        
    # Calculate Difference (Inner - Local)
    # Positive -> Inner is Fainter
    valid['delta_M'] = valid['mean_inn'] - valid['mean_loc']
    
    # Weighted Mean Delta
    # Weight by variance? Or just number?
    # Variance of diff = var1/n1 + var2/n2
    valid['var_diff'] = (valid['std_loc']**2 / valid['count_loc']) + (valid['std_inn']**2 / valid['count_inn'])
    valid['weight'] = 1.0 / (valid['var_diff'] + 1e-6)
    
    mean_delta = np.average(valid['delta_M'], weights=valid['weight'])
    sem_delta = np.sqrt(1.0 / np.sum(valid['weight'])) # Approx
    
    print(f"  Weighted Mean Delta M_K (Inner - Local): {mean_delta:.4f} +/- {sem_delta:.4f} mag")
    
    # Plot
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # Histogram of deltas
    ax.hist(valid['delta_M'], bins=20, weights=valid['weight'], alpha=0.7, density=True)
    ax.axvline(0, color='k', linestyle='--')
    ax.axvline(mean_delta, color='r', linestyle='-', label=f'Mean={mean_delta:.3f}')
    
    ax.set_xlabel('Delta M_K (Inner - Local) [mag]')
    ax.set_ylabel('Weighted Frequency')
    ax.set_title('Test CW: Stellar Twins Luminosity Difference')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    out_path = os.path.join(FIGURES_DIR, 'sdss_test_cw_twins.png')
    plt.savefig(out_path, dpi=150)
    print(f"Figure saved to {out_path}")
    
    return {
        'mean_delta_mag': mean_delta,
        'sem_delta_mag': sem_delta,
        'n_bins': int(len(valid)),
        'n_stars_local': int(valid['count_loc'].sum()),
        'n_stars_inner': int(valid['count_inn'].sum())
    }

def main():
    cache_path = os.path.join(DATA_DIR, 'sdss_twins.csv')
    
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

    results = analyze_twins(df)
    
    out_path = os.path.join(RESULTS_DIR, 'sdss_test_cw_results.json')
    if results:
        with open(out_path, 'w') as f:
            json.dump(results, f, indent=2)
            
        print("\nSUMMARY TEST CW:")
        print("Prediction: Inner Galaxy twins appear fainter (Delta M > 0).")
        print(f"Observed Delta M_K: {results['mean_delta_mag']:.4f} mag")
        
        if results['mean_delta_mag'] > 0.05:
             print("RESULT: SIGNAL (Inner stars are fainter)")
        elif results['mean_delta_mag'] < -0.05:
             print("RESULT: CONTRADICTED (Inner stars are brighter)")
        else:
             print("RESULT: NULL (Consistent luminosities)")

if __name__ == "__main__":
    main()
