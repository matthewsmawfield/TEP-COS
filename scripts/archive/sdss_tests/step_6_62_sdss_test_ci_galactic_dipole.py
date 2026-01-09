#!/usr/bin/env python3
"""
Step 6.62: SDSS Test CI - The Galactic Dipole (Cosmic Coriolis)

Hypothesis:
If the Solar System is moving through a gradient in the cosmic scalar field (or if the field has a dipole structure 
due to "Cosmic Coriolis" effects), then "clocks" (stellar ages, chemical clocks) should show a dipole asymmetry 
across the sky. Stars in the direction of motion might appear systematically younger/older than those in the anti-direction.

Prediction:
Dipole asymmetry in Mean Age (or [C/N] proxy) aligned with cosmic motion or scalar gradient.

Data:
- aspcapStar: param_c_m, param_n_m (to compute [C/N]), param_m_h
- apogee_starhorse: dist50, glon, glat, logg50

Method:
1. Select Red Giant stars (logg between 1 and 3.5) in the local volume (dist < 2 kpc) to minimize selection effects and extinction gradients.
   Ideally, we want a shell or a sphere around the sun.
2. Compute [C/N] = [C/M] - [N/M]. High [C/N] -> Old, Low [C/N] -> Young (typically).
   Check metallicity dependence?
3. Analyze spatial variation of mean [C/N] on the sky (glon, glat).
4. Fit a dipole: Mean([C/N]) ~ A * cos(glon - phase).
   Focus on Galactic Longitude trends for simplicity (planar dipole).
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
    print(f"Querying SDSS for Test CI (Limit: {limit})...")
    
    # 1. Fetch Sample from apogee_starhorse (Local Giants)
    print("  Fetching apogee_starhorse (Sample)...")
    sql_sample = f"""
    SELECT TOP {limit}
        apogee_id,
        dist50,
        glon, glat,
        logg50
    FROM apogee_starhorse
    WHERE 
        dist50 > 0 AND dist50 < 2.0 -- 2 kpc volume
        AND logg50 BETWEEN 1.0 AND 3.5 -- Red Giants
    """
    df_sample = query_sdss(sql_sample)
    
    if df_sample is None or len(df_sample) == 0:
        print("  No sample data found.")
        return None
        
    print(f"  Got {len(df_sample)} stars. Fetching chemistry...")
    ids = df_sample['apogee_id'].astype(str).tolist()
    ids = list(set(ids))
    
    # 2. Fetch aspcapStar (Chemistry)
    # param_c_m, param_n_m
    chunk_size = 20
    df_chem_list = []
    
    for i in range(0, len(ids), chunk_size):
        chunk = ids[i:i+chunk_size]
        ids_str = "', '".join(chunk)
        
        sql_chem = f"""
        SELECT 
            apogee_id,
            param_c_m,
            param_n_m,
            param_m_h
        FROM aspcapStar
        WHERE 
            apogee_id IN ('{ids_str}')
        """
        res = query_sdss(sql_chem)
        if res is not None and len(res) > 0:
            df_chem_list.append(res)
        time.sleep(0.2)
            
    if not df_chem_list:
        print("  No chemistry data found.")
        return None
        
    df_chem = pd.concat(df_chem_list, ignore_index=True)
    
    # 3. Join
    print("  Joining datasets...")
    df = pd.merge(df_sample, df_chem, on='apogee_id', how='inner')
    
    print(f"  Merged N={len(df)}")
    return df

def analyze_dipole(df):
    print("Analyzing Galactic Dipole...")
    
    # Compute [C/N]
    # [C/N] = [C/M] - [N/M]
    df['cn_ratio'] = df['param_c_m'] - df['param_n_m']
    
    # Remove outliers
    df = df[(df['cn_ratio'] > -1.0) & (df['cn_ratio'] < 1.0)]
    
    print(f"  Sample size: {len(df)}")
    print(f"  Mean [C/N]: {df['cn_ratio'].mean():.3f}")
    
    # Bin by Galactic Longitude (glon)
    # glon 0 to 360
    bins = np.linspace(0, 360, 13) # 30 deg bins
    df['l_bin'] = pd.cut(df['glon'], bins=bins)
    
    binned = df.groupby('l_bin')['cn_ratio'].agg(['mean', 'sem', 'count'])
    binned['l_center'] = [i.mid for i in binned.index]
    
    print("\n[C/N] by Galactic Longitude:")
    print(binned[['mean', 'sem', 'count']])
    
    # Fit Dipole: y = A * cos( (x - phase) * pi/180 ) + B
    # Linearize: y = A1 * cos(x) + A2 * sin(x) + B
    # A = sqrt(A1^2 + A2^2), phase = atan2(A2, A1)
    
    valid = binned.dropna()
    x_rad = np.radians(valid['l_center'])
    y = valid['mean']
    w = 1.0 / (valid['sem']**2 + 1e-6) # Inverse variance weighting
    
    # Design matrix
    X = np.column_stack([np.cos(x_rad), np.sin(x_rad), np.ones(len(x_rad))])
    
    # Weighted Least Squares
    # W = diag(w)
    # theta = (X.T W X)^-1 X.T W y
    W = np.diag(w)
    try:
        theta = np.linalg.inv(X.T @ W @ X) @ (X.T @ W @ y)
        A1, A2, B = theta
        
        amplitude = np.sqrt(A1**2 + A2**2)
        phase = np.degrees(np.arctan2(A2, A1))
        
        print(f"  Dipole Amplitude: {amplitude:.4f} dex")
        print(f"  Dipole Phase: {phase:.1f} deg")
        print(f"  Offset B: {B:.4f}")
        
    except Exception as e:
        print(f"  Fit failed: {e}")
        amplitude = 0
        phase = 0
        B = y.mean()
        
    # Plot
    fig, ax = plt.subplots(figsize=(10, 6))
    
    ax.errorbar(binned['l_center'], binned['mean'], yerr=binned['sem'], fmt='o', label='Data')
    
    x_model = np.linspace(0, 360, 100)
    x_model_rad = np.radians(x_model)
    y_model = A1 * np.cos(x_model_rad) + A2 * np.sin(x_model_rad) + B
    
    ax.plot(x_model, y_model, 'r-', label=f'Dipole (Amp={amplitude:.3f})')
    
    ax.set_xlabel('Galactic Longitude [deg]')
    ax.set_ylabel('Mean [C/N] (Age Proxy)')
    ax.set_title('Test CI: Galactic Dipole in Stellar Ages')
    ax.grid(True, alpha=0.3)
    ax.legend()
    
    out_path = os.path.join(FIGURES_DIR, 'sdss_test_ci_dipole.png')
    plt.savefig(out_path, dpi=150)
    print(f"Figure saved to {out_path}")
    
    return {
        'amplitude': amplitude,
        'phase': phase,
        'mean_cn': B,
        'n_sample': int(len(df))
    }

def main():
    cache_path = os.path.join(DATA_DIR, 'sdss_galactic_dipole.csv')
    
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

    results = analyze_dipole(df)
    
    out_path = os.path.join(RESULTS_DIR, 'sdss_test_ci_results.json')
    if results:
        with open(out_path, 'w') as f:
            json.dump(results, f, indent=2)
            
        print("\nSUMMARY TEST CI:")
        print("Prediction: Dipole asymmetry in [C/N] (Age).")
        print(f"Observed Amplitude: {results['amplitude']:.4f} dex")
        
        if results['amplitude'] > 0.05: # Significant gradient?
             print("RESULT: SIGNAL (Dipole detected)")
        else:
             print("RESULT: NULL (Isotropic / Weak dipole)")

if __name__ == "__main__":
    main()
