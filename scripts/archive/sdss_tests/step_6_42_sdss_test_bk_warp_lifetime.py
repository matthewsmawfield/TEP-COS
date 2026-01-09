#!/usr/bin/env python3
"""
Step 6.42: SDSS Test BK - Warp Lifetime (Relaxation Clock)

Hypothesis:
Galactic warps dissipate due to dynamical friction and torque relaxation (rate processes).
In deep potentials, relaxation should be time-dilated (slower).
High-sigma disk galaxies should exhibit a higher fraction of persistent warps than low-sigma disks.

Prediction:
Fraction of Warped Disks increases with sigma (at fixed mass/isolation).

Data:
- MaNGA_GZ2: Morphology classifications.
  - Look for 'Edge-on' -> 'Boxy/Peanut/Warped'?
  - Or 'Odd' -> 'Irregular/Disturbed/Warped'?
  - Specific column: t04_edgeon_a08_yes_fraction? Then t05_bulge_shape...?
  - Or t02_edgeon_a04_yes... ?
  - Let's use 'odd' feature fraction as proxy if explicit warp not found easily.
  - Actually, let's use the 't02_edgeon_a04_yes_fraction' to select edge-on disks,
    and then check if there is a warp classification?
    GZ2 Decision Tree:
    1. Smooth/Features? -> Features.
    2. Edge-on? -> Yes.
    3. Bulge Shape? -> Rounded/Boxy/No Bulge.
    4. Anything Odd? -> Yes -> Ring/Lens/Disturbed/Irregular/Other/Merge.
    Warp is hard to find directly in top-level columns.
    Let's use 'Irregular' or 'Disturbed' fraction in 'Odd' branch as a proxy for long-lived perturbations?
    Or better: Use 't09_bulge_shape_a27_irregular_fraction'? No.
    Let's use 't06_odd_a14_yes_fraction' as a general 'disturbed' proxy for disks.
    TEP: Disturbed states last longer -> Higher fraction.

- mangaDAPall: stellar_sigma_1re.
- mangaTarget: nsa_elpetro_mass.

Method:
1. Select Disk Galaxies (Features/Disk > 0.5).
2. Bin by Sigma.
3. Compute Fraction of 'Odd/Disturbed' galaxies.
4. Correlate with Sigma.
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

def download_data(limit=2000):
    print(f"Querying SDSS for Test BK (Limit: {limit})...")
    
    # Select Disks (t01_smooth_or_features_a02_features_or_disk_fraction > 0.5)
    # Check for 'Odd' fraction (t06_odd_a14_yes_fraction) or similar.
    # Note: GZ2 columns are verbose.
    # We will try to get relevant columns.
    
    sql = f"""
    SELECT TOP {limit}
        g.mangaid,
        g.t01_smooth_or_features_a02_features_or_disk_fraction as p_disk,
        g.t06_odd_a14_yes_fraction as p_odd,
        d.stellar_sigma_1re as sigma,
        t.nsa_elpetro_mass as logmass
        
    FROM MaNGA_GZ2 g
    JOIN mangaDAPall d ON g.mangaid = d.mangaid
    JOIN mangaTarget t ON g.mangaid = t.mangaid
    
    WHERE 
        d.drp3qual = 0
        AND d.stellar_sigma_1re > 0
        AND g.t01_smooth_or_features_a02_features_or_disk_fraction > 0.5
    """
    return query_sdss(sql)

def analyze_warp_lifetime(df):
    print("Analyzing Warp/Disturbance Lifetime...")
    
    # 1. Clean
    df_clean = df.dropna().copy()
    
    # 2. Define "Warped/Disturbed"
    # p_odd > 0.2? (Conservative cut)
    df_clean['is_warped'] = df_clean['p_odd'] > 0.2
    
    # 3. Bin by Sigma
    df_clean['log_sigma'] = np.log10(df_clean['sigma'])
    df_clean['sigma_bin'] = pd.qcut(df_clean['log_sigma'], 8)
    
    binned = df_clean.groupby('sigma_bin')['is_warped'].agg(['mean', 'count'])
    binned['sem'] = np.sqrt(binned['mean'] * (1 - binned['mean']) / binned['count'])
    
    print("\nWarp Fraction by Sigma Bin:")
    print(binned)
    
    # 4. Correlation
    r_warp, p_warp = stats.pearsonr(df_clean['log_sigma'], df_clean['is_warped'])
    print(f"Correlation r(Warped, Sigma): {r_warp:.4f} (p={p_warp:.2e})")
    
    return {
        'r_warp': float(r_warp),
        'p_warp': float(p_warp),
        'binned_data': binned.reset_index().to_dict(orient='records'),
        'n_sample': int(len(df_clean))
    }, df_clean

def create_figure(results):
    print("Generating figure...")
    # Plotting code omitted for brevity/simplicity in this step, focusing on analysis
    # ... (Standard binned plot) ...
    pass

def main():
    cache_path = os.path.join(DATA_DIR, 'sdss_warp_lifetime.csv')
    
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

    results, df_clean = analyze_warp_lifetime(df)
    
    out_path = os.path.join(RESULTS_DIR, 'sdss_test_bk_results.json')
    
    def json_default(obj):
        if isinstance(obj, pd.Interval):
            return str(obj)
        raise TypeError
        
    with open(out_path, 'w') as f:
        json.dump(results, f, indent=2, default=json_default)
        
    print("\nSUMMARY TEST BK:")
    print("TEP Prediction: Warp/Odd fraction increases with sigma (Slower relaxation). r > 0.")
    print(f"Observed r: {results['r_warp']:.4f}")
    
    if results['r_warp'] > 0.05:
        print("RESULT: CONSISTENT (Disturbances last longer in deep potentials)")
    else:
        print("RESULT: NULL/CONTRADICTED (Standard relaxation)")

if __name__ == "__main__":
    main()
