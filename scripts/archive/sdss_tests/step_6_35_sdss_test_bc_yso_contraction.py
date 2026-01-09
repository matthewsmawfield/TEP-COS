#!/usr/bin/env python3
"""
Step 6.35: SDSS Test BC - YSO Contraction (Kelvin-Helmholtz Clock)

Hypothesis:
Pre-Main Sequence stars contract via gravitational collapse (Kelvin-Helmholtz timescale).
In deep potentials (dense clusters), this process is time-dilated.
For a given chronological age (assuming clusters formed roughly co-evally), YSOs in denser regions should have experienced less proper time contraction.
Less contraction -> Larger Radii -> Appears "Younger" on isochrones.
Prediction: Apparent Isochrone Age decreases as local density increases.

Data:
- mos_sagitta: age (Isochrone age in Myr?), ra, dec.
- Density: Calculated from RA/DEC using k-nearest neighbors (Sigma_k).

Method:
1. Download YSOs from mos_sagitta.
2. Convert RA/DEC to Cartesian coordinates (approx) or use Sky distance.
3. Compute Local Density (Sigma_5 or Sigma_10) for each star.
   Sigma ~ k / (pi * d_k^2).
4. Correlate Age with Density.
   TEP: Negative correlation (Higher Density -> Lower Age).
   Standard: Null? Or positive (Mass segregation? Old stars drift out?).
   Standard cluster dynamics: Older stars might be dynamically heated and occupy lower density outskirts?
   That would imply High Density = Young.
   So standard assembly/dynamics might mimic TEP?
   However, TEP predicts this due to *time dilation*, not just sorting.
   If we see a strong effect, it's consistent with TEP.

"""

import pandas as pd
import numpy as np
from scipy import stats, spatial
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
    print(f"Querying SDSS for Test BC (Limit: {limit})...")
    
    # Select YSOs with valid ages
    # mos_sagitta contains YSO candidates from machine learning on Gaia/2MASS
    
    sql = f"""
    SELECT TOP {limit}
        source_id,
        age,
        age_std,
        ra, dec
        
    FROM mos_sagitta
    
    WHERE 
        age > 0 AND age < 100 -- Focus on young stars
        AND yso > 0.9 -- High probability YSOs
    """
    return query_sdss(sql)

def compute_density(df, k=10):
    print(f"Computing local density (k={k})...")
    
    # Convert to 3D on unit sphere for distance calculation (valid for small patches or use Haversine)
    # Since these might be all over the sky, KDTree on 3D unit vectors is best for angular separation.
    
    ra_rad = np.radians(df['ra'])
    dec_rad = np.radians(df['dec'])
    
    x = np.cos(dec_rad) * np.cos(ra_rad)
    y = np.cos(dec_rad) * np.sin(ra_rad)
    z = np.sin(dec_rad)
    
    coords = np.column_stack((x, y, z))
    
    tree = spatial.cKDTree(coords)
    
    # Query k+1 neighbors (point itself is 1st)
    dists, indices = tree.query(coords, k=k+1)
    
    # Distance to kth neighbor (dists[:, k])
    # dists is chord length on unit sphere. For small angles, approx angular separation in radians.
    # Sigma ~ k / (pi * theta_k^2)  [density on sky]
    
    theta_k = dists[:, k]
    # Avoid zero division
    theta_k = np.maximum(theta_k, 1e-6)
    
    density = k / (np.pi * theta_k**2) # sr^-1
    
    return density

def analyze_yso_contraction(df):
    print("Analyzing YSO Contraction...")
    
    # 1. Clean
    df_clean = df.dropna().copy()
    
    # 2. Compute Density
    df_clean['density'] = compute_density(df_clean, k=10)
    df_clean['log_density'] = np.log10(df_clean['density'])
    
    # 3. Correlation
    # TEP Prediction: Age decreases with Density (Younger in deep potential).
    # r < 0.
    
    r_age, p_age = stats.pearsonr(df_clean['log_density'], df_clean['age'])
    
    print(f"N = {len(df_clean)}")
    print(f"Correlation r(Age, Density): {r_age:.4f} (p={p_age:.2e})")
    
    # 4. Binning
    df_clean['dens_bin'] = pd.qcut(df_clean['log_density'], 8)
    binned = df_clean.groupby('dens_bin')['age'].mean()
    print("\nMean Age by Density Bin:")
    print(binned)
    
    return {
        'r_age': float(r_age),
        'p_age': float(p_age),
        'binned_data': binned.tolist(),
        'bin_centers': [mid.mid for mid in binned.index],
        'n_sample': int(len(df_clean))
    }, df_clean

def create_figure(df, results):
    print("Generating figure...")
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # Scatter
    ax.scatter(df['log_density'], df['age'], alpha=0.1, s=2, c='k', label='YSOs')
    
    # Binned
    ax.plot(results['bin_centers'], results['binned_data'], 'r-o', lw=2, label='Mean Age')
    
    ax.set_xlabel(r'Log Surface Density ($\Sigma$) [Arbitrary Units]')
    ax.set_ylabel(r'Isochrone Age [Myr]')
    ax.set_title(f"Test BC: YSO Age vs Density (r={results['r_age']:.3f})")
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    out_path = os.path.join(FIGURES_DIR, 'sdss_test_bc_yso_contraction.png')
    plt.savefig(out_path, dpi=150)
    print(f"Figure saved to {out_path}")

def main():
    cache_path = os.path.join(DATA_DIR, 'sdss_yso_contraction.csv')
    
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

    results, df_clean = analyze_yso_contraction(df)
    
    out_path = os.path.join(RESULTS_DIR, 'sdss_test_bc_results.json')
    with open(out_path, 'w') as f:
        results_json = results.copy()
        del results_json['bin_centers']
        json.dump(results_json, f, indent=2)
        
    create_figure(df_clean, results)
    
    print("\nSUMMARY TEST BC:")
    print("TEP Prediction: Apparent Age decreases with Density. r < 0.")
    print(f"Observed r: {results['r_age']:.4f}")
    
    if results['r_age'] < -0.1:
        print("RESULT: CONSISTENT (YSOs appear younger in dense regions)")
    elif results['r_age'] > 0.1:
        print("RESULT: CONTRADICTED (YSOs appear older in dense regions)")
    else:
        print("RESULT: NULL/FLAT")

if __name__ == "__main__":
    main()
