#!/usr/bin/env python3
"""
Step 6.47: SDSS Test BR - IMF Variation (Jeans Mass Shift)

Hypothesis:
TEP predicts that the effective Jeans Mass scales with the scalar field A(phi).
In deep potentials, altered thermodynamics might favor a Bottom-Heavy IMF (more low mass stars).
The NaD absorption index is sensitive to the dwarf-to-giant ratio (surface gravity).
We expect a Sodium Excess (NaD strength) in high-sigma galaxies that tracks the potential depth,
even after controlling for metallicity.

Prediction:
NaD Index strength correlates with sigma (beyond metallicity effects).
Partial correlation r(NaD, sigma | [Mg/Fe], [Fe/H]) > 0.

Data:
- galSpecIndx: lick_nad, lick_mgb, lick_fe5270 (Indices)
- emissionLinesPort: sigma_stars (Velocity Dispersion)

Method:
1. Fetch NaD, Mgb, Fe5270, and Sigma.
2. Clean data (S/N cuts).
3. Compute Metallicity Proxy [Z/H] ~ Mgb + Fe.
4. Compute Partial Correlation between NaD and Sigma, controlling for Z.
5. Bin by Sigma and plot NaD excess.
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

def download_data(limit=5000):
    print(f"Querying SDSS for Test BR (Limit: {limit})...")
    
    # Check column names first? We know them from schema usually.
    # lick_nad, lick_mgb, lick_fe5270 are standard in galSpecIndx.
    # sigmaStars is correct in emissionLinesPort (from previous checks).
    
    sql = f"""
    SELECT TOP {limit}
        i.specObjID,
        i.lick_nad, i.lick_nad_err,
        i.lick_mgb, i.lick_mgb_err,
        i.lick_fe5270, i.lick_fe5270_err,
        e.sigmaStars as sigma_stars
        
    FROM galSpecIndx i
    JOIN emissionLinesPort e ON i.specObjID = e.specObjID
    
    WHERE 
        e.sigmaStars > 50 AND e.sigmaStars < 400
        AND i.lick_nad > 0
        AND i.lick_mgb > 0
        AND i.lick_fe5270 > 0
        AND i.lick_nad_err < 0.5
    """
    return query_sdss(sql)

def partial_corr(x, y, covar):
    """
    Returns the partial correlation of x and y, controlling for variables in covar (list of arrays).
    """
    data = np.column_stack([x, y] + covar)
    # Pandas correlation matrix inverse method
    df = pd.DataFrame(data)
    corr = df.corr()
    try:
        inv_corr = np.linalg.inv(corr.values)
        # Partial correlation r_xy.z = - P_xy / sqrt(P_xx * P_yy)
        r = -inv_corr[0, 1] / np.sqrt(inv_corr[0, 0] * inv_corr[1, 1])
        return r
    except:
        return np.nan

def analyze_imf_variation(df):
    print("Analyzing IMF Variation (NaD Excess)...")
    
    # Clean
    df = df.dropna().copy()
    
    # Define Metallicity Proxy
    # [Z/H] correlates with Mgb and Fe
    # Use simple sum or mean of metal lines as control
    df['metal_index'] = df['lick_mgb'] + df['lick_fe5270']
    
    # Variables
    # x = log(sigma)
    # y = NaD
    # z = Metal Index
    
    df['log_sigma'] = np.log10(df['sigma_stars'])
    
    # 1. Raw Correlation
    r_raw, p_raw = stats.pearsonr(df['log_sigma'], df['lick_nad'])
    print(f"  Raw Correlation r(NaD, sigma): {r_raw:.4f} (p={p_raw:.2e})")
    
    # 2. Correlation with Metallicity
    r_z, p_z = stats.pearsonr(df['metal_index'], df['lick_nad'])
    print(f"  Metallicity Correlation r(NaD, Z): {r_z:.4f}")
    
    # 3. Partial Correlation (controlling for Metallicity)
    r_partial = partial_corr(df['log_sigma'], df['lick_nad'], [df['metal_index']])
    print(f"  Partial Correlation r(NaD, sigma | Z): {r_partial:.4f}")
    
    # 4. Binning and Plotting
    # Compute NaD residual after removing metallicity trend
    slope_z, intercept_z, _, _, _ = stats.linregress(df['metal_index'], df['lick_nad'])
    df['nad_resid'] = df['lick_nad'] - (slope_z * df['metal_index'] + intercept_z)
    
    # Bin by Sigma
    df['sigma_bin'] = pd.qcut(df['sigma_stars'], 8)
    binned = df.groupby('sigma_bin')['nad_resid'].agg(['mean', 'std', 'count', 'median'])
    binned['sem'] = binned['std'] / np.sqrt(binned['count'])
    
    # Get bin centers (sigma)
    binned['sigma_center'] = [i.mid for i in binned.index]
    
    print("\nNaD Residual by Sigma Bin:")
    print(binned[['mean', 'sem', 'count']])
    
    # Plot
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    # Panel 1: NaD vs Sigma (Raw)
    axes[0].hexbin(df['log_sigma'], df['lick_nad'], gridsize=30, cmap='Greys', mincnt=1)
    axes[0].set_xlabel('log(Velocity Dispersion)')
    axes[0].set_ylabel('NaD Index [Angstrom]')
    axes[0].set_title(f'Raw Relation (r={r_raw:.2f})')
    
    # Panel 2: NaD Residual vs Sigma
    axes[1].errorbar(binned['sigma_center'], binned['mean'], yerr=binned['sem'], fmt='o-', capsize=3)
    axes[1].axhline(0, color='k', linestyle='--')
    axes[1].set_xlabel('Velocity Dispersion [km/s]')
    axes[1].set_ylabel('NaD Residual (Metallicity Subtracted)')
    axes[1].set_title(f'Excess NaD (Partial r={r_partial:.2f})')
    
    plt.tight_layout()
    out_path = os.path.join(FIGURES_DIR, 'sdss_test_br_imf.png')
    plt.savefig(out_path, dpi=150)
    print(f"Figure saved to {out_path}")
    
    return {
        'r_raw': r_raw,
        'r_partial': r_partial,
        'n_sample': int(len(df)),
        'slope_resid': float(stats.linregress(binned['sigma_center'], binned['mean'])[0])
    }

def main():
    cache_path = os.path.join(DATA_DIR, 'sdss_imf_nad.csv')
    
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

    results = analyze_imf_variation(df)
    
    out_path = os.path.join(RESULTS_DIR, 'sdss_test_br_results.json')
    with open(out_path, 'w') as f:
        json.dump(results, f, indent=2)
        
    print("\nSUMMARY TEST BR:")
    print("TEP Prediction: NaD Excess (positive residual) at high sigma.")
    print(f"Observed Partial r: {results['r_partial']:.4f}")
    
    if results['r_partial'] > 0.1:
        print("RESULT: CONSISTENT (NaD excess detected)")
    else:
        print("RESULT: NULL (NaD tracks metallicity only)")

if __name__ == "__main__":
    main()
