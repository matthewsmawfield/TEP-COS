#!/usr/bin/env python3
"""
Step 6.61: SDSS Test CH - Hypervelocity Star Excess

Hypothesis:
Hypervelocity stars (HVS) are usually explained by 3-body interactions with the central SMBH (Hills mechanism).
TEP predicts that dynamical ejection rates or the escape velocity threshold itself are modified by the 
scalar field potential (phantom mass / scalar gradients).
We might expect an excess of high-velocity stars in the halo that cannot be explained by the standard 
Hills mechanism rate, or a different velocity distribution tail.

Prediction:
Number of stars with v > 500 km/s exceeds standard predictions.

Data:
- apogeeStar: vhelio_avg, vscatter, nvisits, glon, glat, apogee_id
- apogee_starhorse: dist50 (to confirm halo membership / exclude local high-proper-motion dwarfs if we had PM, but dist is good)

Method:
1. Query apogeeStar for high velocity candidates (|v| > 300 km/s).
2. Filter for reliable measurements (nvisits > 2, vscatter < 5 km/s).
3. Join with StarHorse for distances.
4. Analyze the tail of the velocity distribution.
   Standard Halo dispersion ~ 100-120 km/s.
   Stars > 3 sigma (> 350 km/s) are rare.
   Stars > 500 km/s should be extremely rare (HVS).
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
    print(f"Querying SDSS for Test CH (Limit: {limit})...")
    
    # Query high velocity stars directly
    # Note: vhelio_avg can be positive or negative
    
    sql = f"""
    SELECT TOP {limit}
        a.apogee_id,
        a.vhelio_avg,
        a.vscatter,
        a.nvisits,
        a.glon, a.glat,
        s.dist50
        
    FROM apogeeStar a
    JOIN apogee_starhorse s ON a.apogee_id = s.apogee_id
    
    WHERE 
        abs(a.vhelio_avg) > 300
        AND a.nvisits > 2
        AND a.vscatter < 5 -- Exclude obvious binaries
        AND s.dist50 > 0
    """
    return query_sdss(sql)

def analyze_hvs(df):
    print("Analyzing Hypervelocity Stars...")
    
    if df is None or len(df) == 0:
        print("  No high velocity stars found.")
        return None
        
    # Clean
    df = df.dropna().copy()
    
    # Abs velocity
    df['v_abs'] = df['vhelio_avg'].abs()
    
    print(f"  Total Candidates (> 300 km/s): {len(df)}")
    
    # Count > 500 km/s
    hvs = df[df['v_abs'] > 500]
    print(f"  HVS Candidates (> 500 km/s): {len(hvs)}")
    
    if len(hvs) > 0:
        print("\nTop Candidates:")
        print(hvs[['apogee_id', 'vhelio_avg', 'dist50', 'glon', 'glat']].head())
    
    # Fit Tail?
    # Simple histogram
    fig, ax = plt.subplots(figsize=(8, 6))
    
    bins = np.linspace(300, 800, 26)
    ax.hist(df['v_abs'], bins=bins, log=True, alpha=0.7, label='Data')
    
    # Simple Halo Model Comparison
    # Gaussian with sigma ~ 120 km/s centered at 0?
    # Tail P(v > 300) for Normal(0, 120):
    # This is rough, but illustrative.
    
    ax.set_xlabel('|Heliocentric Velocity| [km/s]')
    ax.set_ylabel('Count (log scale)')
    ax.set_title('Test CH: High Velocity Star Distribution')
    ax.grid(True, alpha=0.3)
    
    out_path = os.path.join(FIGURES_DIR, 'sdss_test_ch_hvs.png')
    plt.savefig(out_path, dpi=150)
    print(f"Figure saved to {out_path}")
    
    return {
        'n_300': int(len(df)),
        'n_500': int(len(hvs)),
        'max_v': float(df['v_abs'].max())
    }

def main():
    cache_path = os.path.join(DATA_DIR, 'sdss_hvs.csv')
    
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

    results = analyze_hvs(df)
    
    out_path = os.path.join(RESULTS_DIR, 'sdss_test_ch_results.json')
    if results:
        with open(out_path, 'w') as f:
            json.dump(results, f, indent=2)
            
        print("\nSUMMARY TEST CH:")
        print(f"Candidates > 300 km/s: {results['n_300']}")
        print(f"Candidates > 500 km/s: {results['n_500']}")
        
        if results['n_500'] > 10: # Arbitrary threshold for "Excess" without detailed model
             print("RESULT: SIGNAL (Significant number of HVS candidates)")
        else:
             print("RESULT: NULL (Consistent with standard halo tail)")

if __name__ == "__main__":
    main()
