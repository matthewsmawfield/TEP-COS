#!/usr/bin/env python3
"""
Step 6.65: SDSS Test CL - Spiral Chirality (Cosmic Coriolis)

Hypothesis:
If the TEP scalar field has a global vorticity or "Cosmic Coriolis" component (breaking parity on large scales), 
it might bias the formation of spiral arms, leading to an excess of S-wise vs Z-wise spirals 
(or vice versa) in a dipole pattern across the sky.

Prediction:
Dipole asymmetry in the fraction of S-wise vs Z-wise spirals (CW vs ACW).

Data:
- zooSpec (GZ1): p_cw, p_acw, ra, dec, z
- OR zoo2MainSpecz (GZ2): but GZ1 is simpler for chirality

Method:
1. Fetch GZ1 data.
2. Select Clean Spirals (p_cw + p_acw > 0.5? or spiral > 0.8?)
   Better: Use likelihoods directly or a clean cut.
   Likelihood: L = (p_cw - p_acw) / (p_cw + p_acw).
3. Analyze spatial variation of Mean(L) across the sky.
4. Fit Dipole: Mean(L) ~ A * cos(theta).
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

def download_data(limit=10000):
    print(f"Querying SDSS for Test CL (Limit: {limit})...")
    
    # Use zooSpec (GZ1) for simplicity
    
    sql = f"""
    SELECT TOP {limit}
        objid,
        ra, dec,
        p_cw,
        p_acw,
        p_mg as p_merger
        
    FROM zooSpec
    WHERE 
        (p_cw > 0.5 OR p_acw > 0.5)
        AND p_mg < 0.2
    """
    return query_sdss(sql)

def analyze_chirality(df):
    print("Analyzing Spiral Chirality...")
    
    # Clean
    df = df.dropna().copy()
    
    # Define Chirality Index
    # C = +1 for CW (Z-wise?), -1 for ACW (S-wise?)
    # Usually: Z-wise (CW on sky) vs S-wise (ACW on sky)
    # Let's use p_cw - p_acw
    # But normalized?
    # Since we selected p_cw > 0.5 OR p_acw > 0.5, the denominator is likely > 0.5
    
    df['chirality'] = (df['p_cw'] - df['p_acw']) / (df['p_cw'] + df['p_acw'])
    
    print(f"  Sample size: {len(df)}")
    print(f"  Mean Chirality (CW - ACW): {df['chirality'].mean():.4f}")
    print(f"  CW Fraction: {(df['chirality'] > 0).mean():.4f}")
    
    # Bin by RA (0 to 360)
    bins = np.linspace(0, 360, 13)
    df['ra_bin'] = pd.cut(df['ra'], bins=bins)
    
    binned = df.groupby('ra_bin')['chirality'].agg(['mean', 'sem', 'count'])
    binned['ra_center'] = [i.mid for i in binned.index]
    
    print("\nChirality by RA:")
    print(binned[['mean', 'sem', 'count']])
    
    # Fit Dipole: y = A * cos(ra - phase) + B
    valid = binned.dropna()
    x_rad = np.radians(valid['ra_center'])
    y = valid['mean']
    w = 1.0 / (valid['sem']**2 + 1e-6)
    
    X = np.column_stack([np.cos(x_rad), np.sin(x_rad), np.ones(len(x_rad))])
    W = np.diag(w)
    
    try:
        theta = np.linalg.inv(X.T @ W @ X) @ (X.T @ W @ y)
        A1, A2, B = theta
        
        amplitude = np.sqrt(A1**2 + A2**2)
        phase = np.degrees(np.arctan2(A2, A1))
        
        print(f"  Dipole Amplitude: {amplitude:.4f}")
        print(f"  Dipole Phase: {phase:.1f} deg")
        print(f"  Global Bias B: {B:.4f}")
        
    except Exception as e:
        print(f"  Fit failed: {e}")
        amplitude = 0
        phase = 0
        B = y.mean()
        
    # Plot
    fig, ax = plt.subplots(figsize=(10, 6))
    
    ax.errorbar(binned['ra_center'], binned['mean'], yerr=binned['sem'], fmt='o', label='Data')
    
    x_model = np.linspace(0, 360, 100)
    x_model_rad = np.radians(x_model)
    y_model = A1 * np.cos(x_model_rad) + A2 * np.sin(x_model_rad) + B
    
    ax.plot(x_model, y_model, 'r-', label=f'Dipole (Amp={amplitude:.3f})')
    
    ax.set_xlabel('Right Ascension [deg]')
    ax.set_ylabel('Mean Chirality (CW - ACW)')
    ax.set_title(f'Test CL: Spiral Chirality Dipole (Amp={amplitude:.3f})')
    ax.axhline(0, color='k', linestyle='--', alpha=0.3)
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    out_path = os.path.join(FIGURES_DIR, 'sdss_test_cl_chirality.png')
    plt.savefig(out_path, dpi=150)
    print(f"Figure saved to {out_path}")
    
    return {
        'amplitude': amplitude,
        'phase': phase,
        'bias': B,
        'n_sample': int(len(df))
    }

def main():
    cache_path = os.path.join(DATA_DIR, 'sdss_chirality.csv')
    
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

    results = analyze_chirality(df)
    
    out_path = os.path.join(RESULTS_DIR, 'sdss_test_cl_results.json')
    if results:
        with open(out_path, 'w') as f:
            json.dump(results, f, indent=2)
            
        print("\nSUMMARY TEST CL:")
        print("Prediction: Dipole asymmetry in chirality.")
        print(f"Observed Amplitude: {results['amplitude']:.4f}")
        
        if results['amplitude'] > 0.05:
             print("RESULT: SIGNAL (Dipole detected)")
        else:
             print("RESULT: NULL (No significant dipole)")

if __name__ == "__main__":
    main()
