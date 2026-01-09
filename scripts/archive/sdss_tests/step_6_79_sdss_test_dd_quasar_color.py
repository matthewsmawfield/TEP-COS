#!/usr/bin/env python3
"""
Step 6.79: SDSS Test DD - Quasar Color-Potential Relation

Hypothesis:
Quasar accretion disks emit thermal radiation. TEP modifies the effective potential 
and photon rates. For a given Luminosity and BH Mass, the spectral energy distribution (SED) 
should shift. Quasars in deeper potentials (high sigma_host) should appear **redder** 
(cooler effective temp) due to stronger time dilation of the emitting surface.

Prediction:
Quasar Color (u-g) correlates with Host Sigma (at fixed L, M_BH).

Data:
- mos_sdss_dr16_qso: logLbol, logBH, z
- SpecPhotoAll: modelMag_u, modelMag_g, modelMag_r
- emissionLinesPort: sigma_stars

Method:
1. Select low-z quasars (z < 0.8) where host sigma might be measurable.
2. Join QSO properties with photometry and host sigma.
3. Control for Redshift, Luminosity, and BH Mass (as color depends on accretion rate/mass).
4. Analyze residuals of Color vs Sigma.
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
    print(f"Querying SDSS for Test DD (Limit: {limit})...")
    
    sql = f"""
    SELECT TOP {limit}
        q.specObjID, 
        q.z,
        q.logLbol, 
        q.logBH,
        p.sigma_stars,
        s.modelMag_u - s.modelMag_g as u_g,
        s.modelMag_g - s.modelMag_r as g_r
    FROM mos_sdss_dr16_qso q
    JOIN SpecPhotoAll s ON q.specObjID = s.specObjID
    JOIN emissionLinesPort p ON q.specObjID = p.specObjID
    WHERE q.z < 0.8
      AND p.sigma_stars > 0
      AND q.logLbol > 0
    """
    return query_sdss(sql)

def analyze_quasar_color(df):
    print("Analyzing Quasar Color vs Potential...")
    
    if df is None or len(df) < 50:
        print("  Insufficient data.")
        return None
        
    df = df.dropna().copy()
    
    print(f"  Data points: {len(df)}")
    
    # We need to control for L, M_BH, and z, as they affect color.
    # Model: u-g = a*logL + b*logBH + c*z + d*log(sigma) + const
    
    df['log_sigma'] = np.log10(df['sigma_stars'])
    
    from sklearn.linear_model import LinearRegression
    
    # 1. Fit background model (L, M, z)
    X_control = df[['logLbol', 'logBH', 'z']]
    y = df['u_g']
    
    reg = LinearRegression().fit(X_control, y)
    print(f"  Control Fit R2: {reg.score(X_control, y):.3f}")
    
    # 2. Get residuals
    df['color_resid'] = y - reg.predict(X_control)
    
    # 3. Correlate residuals with Sigma
    slope, intercept, r_val, p_val, std_err = stats.linregress(df['log_sigma'], df['color_resid'])
    
    print(f"  Correlation (log Sigma vs Color Residual): r={r_val:.3f}, p={p_val:.2e}, slope={slope:.4f}")
    
    # Plot
    fig, ax = plt.subplots(1, 2, figsize=(12, 5))
    
    # Raw relation (just for context)
    ax[0].scatter(df['log_sigma'], df['u_g'], alpha=0.5, s=10)
    ax[0].set_xlabel('log(Sigma)')
    ax[0].set_ylabel('u-g Color')
    ax[0].set_title('Raw Data')
    
    # Residuals
    ax[1].scatter(df['log_sigma'], df['color_resid'], alpha=0.5, s=10, c='crimson')
    
    x_range = np.linspace(df['log_sigma'].min(), df['log_sigma'].max(), 100)
    ax[1].plot(x_range, intercept + slope*x_range, 'k--', label=f'Slope={slope:.3f}, p={p_val:.3f}')
    
    ax[1].set_xlabel('log(Host Velocity Dispersion)')
    ax[1].set_ylabel('Color Residual (u-g)')
    ax[1].set_title('Test DD: Quasar Reddening vs Potential')
    ax[1].legend()
    ax[1].grid(True, alpha=0.3)
    
    out_path = os.path.join(FIGURES_DIR, 'sdss_test_dd_quasar_color.png')
    plt.savefig(out_path, dpi=150)
    print(f"Figure saved to {out_path}")
    
    return {
        'correlation_r': r_val,
        'p_value': p_val,
        'slope': slope,
        'n_qso': int(len(df))
    }

def main():
    cache_path = os.path.join(DATA_DIR, 'sdss_qso_color.csv')
    
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

    results = analyze_quasar_color(df)
    
    out_path = os.path.join(RESULTS_DIR, 'sdss_test_dd_results.json')
    if results:
        with open(out_path, 'w') as f:
            json.dump(results, f, indent=2)
            
        print("\nSUMMARY TEST DD:")
        print(f"Slope (log Sigma vs Color Resid): {results['slope']:.4f}")
        
        # Prediction: Redder = Higher u-g (positive slope)
        if results['p_value'] < 0.05 and results['slope'] > 0.05:
             print("RESULT: SIGNAL (Quasars are redder in deep potentials)")
        elif results['p_value'] < 0.05 and results['slope'] < -0.05:
             print("RESULT: CONTRADICTED (Quasars are bluer in deep potentials)")
        else:
             print("RESULT: NULL (No significant dependence)")

if __name__ == "__main__":
    main()
