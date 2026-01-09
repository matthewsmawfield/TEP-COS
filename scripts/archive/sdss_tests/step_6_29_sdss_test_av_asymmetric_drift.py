#!/usr/bin/env python3
"""
Step 6.29: SDSS Test AV - Asymmetric Drift Anomaly

Hypothesis:
Asymmetric drift (v_gas - v_star) arises because stars have velocity dispersion (pressure support) which lags their rotation compared to the cold gas (circular velocity).
v_c^2 - v_star^2 = sigma_r^2 * (dln(rho)/dlnR + ...)
In TEP, if inner galaxy regions are dynamically "younger" (less heating time due to dilation), the stellar velocity dispersion might be lower than equilibrium, or the lag might be smaller?
Actually, TEP predicts "slower heating". So sigma should be lower for a given age?
Or, if we observe a system with high sigma (deep potential), TEP says it experienced LESS proper time than a low sigma system of the same coordinate age.
Less proper time = Less dynamical heating (scattering).
Less heating = Stellar velocity closer to Gas velocity (Cold disk).
Prediction: Asymmetric Drift (v_gas - v_star) is SMALLER in high-sigma galaxies than expected from standard Jeans modeling (or simply decreases with sigma at fixed mass/morphology).
Standard expectation: High sigma = High pressure support = Large drift.
TEP expectation: High sigma = Slower heating = "Colder" orbits relative to potential depth = Smaller drift (normalized).

Observable:
Drift = |v_gas| - |v_star| (at 1 Re).
Normalize by v_gas? Delta_v / v_gas.
Correlation with sigma.

Data:
- mangaDAPall: 
    - ha_gvel_1re (Gas velocity at 1 Re)
    - stellar_vel_1re (Stellar velocity at 1 Re)
    - stellar_sigma_1re (Velocity dispersion)
    - nsa_sersic_n (Morphology control)
    - nsa_sersic_ba (Inclination proxy)

Method:
1. Select rotating galaxies (v_gas > 50 km/s).
2. Compute Drift = abs(v_gas) - abs(v_star).
3. Normalize Drift: f_drift = Drift / abs(v_gas).
4. Correlate f_drift with sigma.
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

def download_data(limit=1000):
    print(f"Querying SDSS for Test AV (Limit: {limit})...")
    
    # We need galaxies with valid rotation curves in both gas and stars.
    # Using 'hi_clip' velocities as proxies for rotation amplitude at large radii
    
    sql = f"""
    SELECT TOP {limit}
        mangaid,
        ha_gvel_hi_clip as v_gas,
        stellar_vel_hi_clip as v_star,
        stellar_sigma_1re as sigma,
        nsa_sersic_n as sersic_n,
        nsa_sersic_ba as axis_ratio
        
    FROM mangaDAPall
    
    WHERE 
        drp3qual = 0
        AND stellar_sigma_1re > 0
    """
    return query_sdss(sql)

def analyze_asymmetric_drift(df):
    print("Analyzing Asymmetric Drift...")
    
    # 1. Clean
    df_clean = df.dropna().copy()
    
    # Filter inclination (axis ratio < 0.8)
    df_clean = df_clean[df_clean['axis_ratio'] < 0.8].copy()
    
    # Ensure significant rotation to avoid noise domination
    # v_gas > 50 km/s
    df_clean = df_clean[np.abs(df_clean['v_gas']) > 50].copy()
    
    df_clean['log_sigma'] = np.log10(df_clean['sigma'])
    
    # 2. Compute Drift
    # Drift = |v_gas| - |v_star|
    # Note: v_star should lag v_gas, so |v_gas| > |v_star| usually.
    # Check signs? Usually they rotate same direction.
    
    # We use absolute values because direction depends on orientation on sky
    df_clean['v_gas_abs'] = np.abs(df_clean['v_gas'])
    df_clean['v_star_abs'] = np.abs(df_clean['v_star'])
    
    df_clean['drift'] = df_clean['v_gas_abs'] - df_clean['v_star_abs']
    
    # Normalize
    df_clean['f_drift'] = df_clean['drift'] / df_clean['v_gas_abs']
    
    # Filter outliers
    df_clean = df_clean[(df_clean['f_drift'] > -0.5) & (df_clean['f_drift'] < 1.0)]
    
    # 3. Correlation with Sigma
    # Standard: High sigma -> High Drift (f_drift increases).
    # TEP: High sigma -> Less Heating -> Lower Drift (f_drift decreases or is lower than standard).
    
    r_drift, p_drift = stats.pearsonr(df_clean['log_sigma'], df_clean['f_drift'])
    
    print(f"N = {len(df_clean)}")
    print(f"Mean Drift Fraction: {df_clean['f_drift'].mean():.3f}")
    print(f"Correlation r(f_drift, sigma): {r_drift:.4f} (p={p_drift:.2e})")
    
    # 4. Binning
    df_clean['sigma_bin'] = pd.qcut(df_clean['log_sigma'], 8)
    binned = df_clean.groupby('sigma_bin')['f_drift'].mean()
    print("\nMean Drift Fraction by Sigma Bin:")
    print(binned)
    
    return {
        'r_drift': float(r_drift),
        'p_drift': float(p_drift),
        'mean_f_drift': float(df_clean['f_drift'].mean()),
        'binned_means': binned.tolist(),
        'bin_centers': [mid.mid for mid in binned.index],
        'n_sample': int(len(df_clean))
    }, df_clean

def create_figure(df, results):
    print("Generating figure...")
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # Scatter
    ax.scatter(df['log_sigma'], df['f_drift'], alpha=0.1, s=2, c='k', label='Galaxies')
    
    # Binned
    ax.plot(results['bin_centers'], results['binned_means'], 'r-o', lw=2, label='Mean Drift')
    
    ax.set_xlabel(r'$\log(\sigma)$')
    ax.set_ylabel(r'Asymmetric Drift Fraction ($v_{lag}/v_{circ}$)')
    ax.set_title(f"Test AV: Asymmetric Drift (r={results['r_drift']:.3f})")
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    out_path = os.path.join(FIGURES_DIR, 'sdss_test_av_asymmetric_drift.png')
    plt.savefig(out_path, dpi=150)
    print(f"Figure saved to {out_path}")

def main():
    cache_path = os.path.join(DATA_DIR, 'sdss_asymmetric_drift.csv')
    
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

    results, df_clean = analyze_asymmetric_drift(df)
    
    out_path = os.path.join(RESULTS_DIR, 'sdss_test_av_results.json')
    with open(out_path, 'w') as f:
        results_json = results.copy()
        del results_json['bin_centers']
        json.dump(results_json, f, indent=2)
        
    create_figure(df_clean, results)
    
    print("\nSUMMARY TEST AV:")
    print("TEP Prediction: Drift decreases with sigma (Slower heating). r < 0 (or less positive than standard).")
    print("Standard Prediction: Drift increases with sigma (Pressure support). r > 0.")
    print(f"Observed r: {results['r_drift']:.4f}")
    
    if results['r_drift'] < -0.1:
        print("RESULT: CONSISTENT (Drift decreases in deep potentials)")
    elif results['r_drift'] > 0.1:
        print("RESULT: CONTRADICTED (Standard pressure support dominates)")
    else:
        print("RESULT: NULL/FLAT")

if __name__ == "__main__":
    main()
