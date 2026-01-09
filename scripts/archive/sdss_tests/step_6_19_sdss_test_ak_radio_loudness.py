#!/usr/bin/env python3
"""
Step 6.19: SDSS Test AK - Radio Loudness Anomaly

Hypothesis:
Radio jet power depends on accretion rate and particle acceleration rates.
TEP predicts time dilation in deep potentials (high sigma).
If clock rates affect the observed power output or jet physics, we might expect a correlation between Radio Loudness and Host Sigma.
Prediction: Radio Loudness (L_radio / L_opt) decreases with host sigma? 
(Or increases if 'tired light' mechanism affects optical more? TEP suggests generalized slowing, so maybe power P = E/t is suppressed).
Let's test for ANY correlation r(R, sigma).

Data:
- FIRST: Radio Flux (peak, mJy).
- PhotoObjAll: Optical Magnitude (modelMag_r).
- emissionLinesPort: Sigma Stars.
- mos_sdss_dr16_qso: Bolometric Luminosity (optional control).

Radio Loudness R = log10(F_radio / F_optical)
"""

import pandas as pd
import numpy as np
from scipy import stats
import os
import json
import matplotlib.pyplot as plt
import requests
import time
from sklearn.linear_model import LinearRegression

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
    print(f"Querying SDSS for Test AK (Limit: {limit})...")
    
    # Use SpecPhoto for pre-joined spectroscopy and photometry
    sql = f"""
    SELECT TOP {limit}
        s.specObjID,
        s.z,
        p.sigmaStars,
        s.modelMag_r,
        f.peak as flux_radio_mJy
        
    FROM SpecPhoto s
    JOIN FIRST f ON s.objID = f.objID
    JOIN emissionLinesPort p ON s.specObjID = p.specObjID
    
    WHERE 
        s.z BETWEEN 0.05 AND 0.5
        AND p.sigmaStars > 70 AND p.sigmaStars < 400
        AND f.peak > 1.0 -- > 1 mJy detection
        AND s.modelMag_r > 10 AND s.modelMag_r < 22
        AND s.class IN ('QSO', 'GALAXY')
    """
    return query_sdss(sql)

def analyze_loudness(df):
    print("Analyzing Radio Loudness...")
    
    # 1. Calculate Optical Flux in mJy
    # SDSS AB Mag: F_nu (Jy) = 3631 * 10**(-0.4 * mag)
    # F_nu (mJy) = 3.631e6 * 10**(-0.4 * mag)
    
    df['flux_opt_mJy'] = 3.631e6 * 10**(-0.4 * df['modelMag_r'])
    
    # 2. Calculate Radio Loudness R
    # R = log10(F_radio / F_opt)
    df['R'] = np.log10(df['flux_radio_mJy'] / df['flux_opt_mJy'])
    
    # 3. Clean
    df_clean = df.dropna(subset=['R', 'sigmaStars']).copy()
    df_clean['log_sigma'] = np.log10(df_clean['sigmaStars'])
    
    print(f"N = {len(df_clean)}")
    
    # 4. Correlation: R vs Sigma
    r_simple, p_simple = stats.pearsonr(df_clean['log_sigma'], df_clean['R'])
    
    # 5. Control for Redshift (Selection effects) and L_bol (if available)
    # Simple control for z
    X = df_clean[['z']].values
    y = df_clean['R'].values
    reg = LinearRegression().fit(X, y)
    df_clean['R_resid'] = y - reg.predict(X)
    
    r_controlled, p_controlled = stats.pearsonr(df_clean['log_sigma'], df_clean['R_resid'])
    
    print(f"Simple r(R, log_sigma): {r_simple:.4f} (p={p_simple:.2e})")
    print(f"Controlled r(R_resid, log_sigma): {r_controlled:.4f} (p={p_controlled:.2e})")
    
    return {
        'r_simple': float(r_simple),
        'p_simple': float(p_simple),
        'r_controlled': float(r_controlled),
        'p_controlled': float(p_controlled),
        'n_sample': int(len(df_clean))
    }, df_clean

def create_figure(df, results):
    print("Generating figure...")
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # Plot 1: R vs Sigma
    ax = axes[0]
    ax.scatter(df['log_sigma'], df['R'], alpha=0.3, s=10, c='blue')
    m, b = np.polyfit(df['log_sigma'], df['R'], 1)
    x = np.linspace(df['log_sigma'].min(), df['log_sigma'].max(), 100)
    ax.plot(x, m*x + b, 'r-', label=f'r={results["r_simple"]:.3f}')
    ax.set_xlabel(r'$\log(\sigma)$')
    ax.set_ylabel(r'Radio Loudness $R = \log(F_{radio}/F_{opt})$')
    ax.set_title("Radio Loudness vs Potential Depth")
    ax.legend()
    
    # Plot 2: Residual R vs Sigma
    ax = axes[1]
    ax.scatter(df['log_sigma'], df['R_resid'], alpha=0.3, s=10, c='green')
    m, b = np.polyfit(df['log_sigma'], df['R_resid'], 1)
    ax.plot(x, m*x + b, 'r-', label=f'r={results["r_controlled"]:.3f}')
    ax.set_xlabel(r'$\log(\sigma)$')
    ax.set_ylabel(r'Residual Radio Loudness (z-corrected)')
    ax.set_title("Controlled Radio Loudness vs Potential")
    ax.legend()
    
    plt.tight_layout()
    out_path = os.path.join(FIGURES_DIR, 'sdss_test_ak_radio_loudness.png')
    plt.savefig(out_path, dpi=150)
    print(f"Figure saved to {out_path}")

def main():
    cache_path = os.path.join(DATA_DIR, 'sdss_radio_loudness.csv')
    
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

    results, df_clean = analyze_loudness(df)
    
    out_path = os.path.join(RESULTS_DIR, 'sdss_test_ak_results.json')
    with open(out_path, 'w') as f:
        json.dump(results, f, indent=2)
        
    create_figure(df_clean, results)
    
    print("\nSUMMARY TEST AK:")
    print("TEP Prediction: Radio Loudness suppressed in high sigma? (r < 0)")
    print(f"Observed Controlled r: {results['r_controlled']:.4f}")
    
    if results['r_controlled'] < -0.05:
        print("RESULT: CONSISTENT (Negative Correlation)")
    elif results['r_controlled'] > 0.05:
        print("RESULT: CONTRADICTED (Positive Correlation)")
    else:
        print("RESULT: NULL")

if __name__ == "__main__":
    main()
