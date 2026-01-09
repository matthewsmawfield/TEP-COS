#!/usr/bin/env python3
"""
Step 6.78: SDSS Test DC - Pair Decay Ratio (Orbital Clock)

Hypothesis:
The ratio of close pairs (<20 kpc) to wide pairs (20-50 kpc) depends on the speed 
of orbital decay (dynamical friction). If this rate is dilated in high-sigma 
environments (deep potential), pairs spend longer in the "close" phase (or decay 
slower). We expect the ratio N_close / N_wide to vary with the pair's total 
mass/potential (proxied by sum of velocity dispersions).

Prediction:
Close/Wide Pair Ratio correlates with Pair Velocity Dispersion.

Data:
- Neighbors: objID, NeighborObjID, distance (arcmin?)
- SpecObjAll: z
- emissionLinesPort: sigma_stars

Method:
1. Select physical pairs:
   - Angular separation < 2 arcmin (~100 kpc at z=0.1)
   - dz < 0.005 (velocity separation < 1500 km/s)
2. Calculate projected separation (kpc).
3. Classify:
   - Close: < 20 kpc
   - Wide: 20 - 50 kpc
4. Calculate Pair Sigma = sigma1 + sigma2 (or sqrt(s1^2 + s2^2))
5. Bin by Pair Sigma.
6. Calculate Ratio = N_close / N_wide in each bin.
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

def download_data(limit=200):
    print(f"Querying SDSS for Test DC (Limit: {limit})...")
    
    # Note: Neighbors table distance is in arcmin
    
    sql = f"""
    SELECT TOP {limit}
        n.objID,
        n.distance as dist_arcmin,
        s1.z as z1, s2.z as z2,
        e1.sigma_stars as sig1, e2.sigma_stars as sig2
    FROM Neighbors n
    JOIN SpecObjAll s1 ON n.objID = s1.bestObjID
    JOIN SpecObjAll s2 ON n.NeighborObjID = s2.bestObjID
    LEFT JOIN emissionLinesPort e1 ON s1.specObjID = e1.specObjID
    LEFT JOIN emissionLinesPort e2 ON s2.specObjID = e2.specObjID
    WHERE n.distance < 2.0 
      AND abs(s1.z - s2.z) < 0.005
      AND s1.class = 'GALAXY' AND s2.class = 'GALAXY'
      AND s1.z BETWEEN 0.05 AND 0.15
      AND e1.sigma_stars > 0 AND e2.sigma_stars > 0
    """
    return query_sdss(sql)

def analyze_pair_decay(df):
    print("Analyzing Pair Decay Ratio...")
    
    if df is None or len(df) < 100:
        print("  Insufficient data.")
        return None
        
    df = df.dropna().copy()
    
    # Calculate Separation in kpc
    H0 = 70.0
    c = 300000.0
    
    # Mean redshift
    df['z_mean'] = (df['z1'] + df['z2']) / 2
    
    # Angular scale (kpc/arcmin)
    # 1 rad = D_A kpc
    # 1 arcmin = (1/60) * (pi/180) rad
    # D_A approx c*z/H0 (low z)
    df['dist_mpc'] = (c * df['z_mean']) / H0
    df['kpc_per_arcmin'] = (df['dist_mpc'] * 1000) * (np.pi / (180 * 60))
    
    df['sep_kpc'] = df['dist_arcmin'] * df['kpc_per_arcmin']
    
    # Classify
    # Close: < 20 kpc
    # Wide: 20 - 50 kpc (upper bound set by query radius ~100kpc, but let's restrict)
    
    df['type'] = 'Other'
    df.loc[df['sep_kpc'] < 20, 'type'] = 'Close'
    df.loc[(df['sep_kpc'] >= 20) & (df['sep_kpc'] < 50), 'type'] = 'Wide'
    
    df = df[df['type'] != 'Other']
    
    print("  Pair Counts:")
    print(df['type'].value_counts())
    
    # Pair Sigma
    df['pair_sigma'] = np.sqrt(df['sig1']**2 + df['sig2']**2)
    
    # Bin by Pair Sigma
    df['sigma_bin'] = pd.qcut(df['pair_sigma'], q=5, labels=False)
    
    results = []
    for bin_id in sorted(df['sigma_bin'].unique()):
        subset = df[df['sigma_bin'] == bin_id]
        n_close = len(subset[subset['type'] == 'Close'])
        n_wide = len(subset[subset['type'] == 'Wide'])
        mean_sigma = subset['pair_sigma'].mean()
        
        ratio = n_close / n_wide if n_wide > 0 else np.nan
        
        # Error on ratio? Poisson: err ~ ratio * sqrt(1/Nc + 1/Nw)
        if n_close > 0 and n_wide > 0:
            err = ratio * np.sqrt(1/n_close + 1/n_wide)
        else:
            err = 0
            
        results.append({
            'bin': bin_id,
            'sigma': mean_sigma,
            'ratio': ratio,
            'err': err,
            'n_close': n_close,
            'n_wide': n_wide
        })
        
    res_df = pd.DataFrame(results).dropna()
    print(res_df)
    
    if len(res_df) < 3:
        print("  Not enough bins.")
        return None
    
    # Regression
    slope, intercept, r_val, p_val, std_err = stats.linregress(res_df['sigma'], res_df['ratio'])
    
    print(f"  Correlation (Sigma vs Pair Ratio): r={r_val:.3f}, p={p_val:.2e}")
    
    # Plot
    plt.figure(figsize=(8, 6))
    plt.errorbar(res_df['sigma'], res_df['ratio'], yerr=res_df['err'], fmt='o-', color='green')
    
    x_range = np.linspace(res_df['sigma'].min(), res_df['sigma'].max(), 100)
    plt.plot(x_range, intercept + slope*x_range, 'k--', label=f'r={r_val:.2f}')
    
    plt.xlabel('Pair Velocity Dispersion $\sqrt{\sigma_1^2 + \sigma_2^2}$ (km/s)')
    plt.ylabel('Close/Wide Pair Ratio (<20kpc / 20-50kpc)')
    plt.title('Test DC: Pair Decay Rate vs Potential Depth')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    out_path = os.path.join(FIGURES_DIR, 'sdss_test_dc_pair_decay.png')
    plt.savefig(out_path, dpi=150)
    print(f"Figure saved to {out_path}")
    
    return {
        'correlation_r': r_val,
        'p_value': p_val,
        'slope': slope
    }

def main():
    cache_path = os.path.join(DATA_DIR, 'sdss_pairs.csv')
    
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

    results = analyze_pair_decay(df)
    
    out_path = os.path.join(RESULTS_DIR, 'sdss_test_dc_results.json')
    if results:
        with open(out_path, 'w') as f:
            json.dump(results, f, indent=2)
            
        print("\nSUMMARY TEST DC:")
        print(f"Correlation (Sigma vs Pair Ratio): {results['correlation_r']:.3f}")
        
        if results['p_value'] < 0.05 and abs(results['correlation_r']) > 0.3:
             print("RESULT: SIGNAL (Pair ratio depends on potential)")
        else:
             print("RESULT: NULL (No significant dependence)")

if __name__ == "__main__":
    main()
