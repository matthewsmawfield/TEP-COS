#!/usr/bin/env python3
"""
Step 6.87: SDSS Test DP - TiO IMF Indicators (Bottom-Heavy IMF)

Hypothesis:
Titanium Oxide (TiO) absorption features are strong in low-mass stars. An enhancement 
of TiO indices relative to mean metallicity suggests a bottom-heavy Initial Mass 
Function (IMF). TEP thermodynamics (metric shock/Jeans mass) might favor low-mass 
star formation in deep potentials.

Prediction:
TiO indices are enhanced in high-sigma galaxies (IMF variation).

Data:
- galSpecIndx: lick_tio2, lick_mgb, lick_fe5270
- emissionLinesPort: sigma_stars

Method:
1. Join galSpecIndx with emissionLinesPort.
2. Control for Metallicity (Mg, Fe) as TiO also depends on Z.
3. Analyze residuals of TiO vs Sigma.
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
    print(f"Querying SDSS for Test DP (Limit: {limit})...")
    
    sql = f"""
    SELECT TOP {limit}
        i.specObjID, 
        i.lick_tio2 as tio2, 
        i.lick_tio2_err,
        e.sigma_stars as sigma,
        i.lick_mgb as mgb, 
        i.lick_fe5270 as fe5270
    FROM galSpecIndx i
    JOIN emissionLinesPort e ON i.specObjID = e.specObjID
    WHERE e.sigma_stars > 50
      AND i.lick_tio2 > 0
      AND i.lick_mgb > 0
    """
    return query_sdss(sql)

def analyze_tio_imf(df):
    print("Analyzing TiO IMF Indicators...")
    
    if df is None or len(df) < 50:
        print("  Insufficient data.")
        return None
        
    df = df.dropna().copy()
    
    # Metallicity Proxy: [MgFe]
    # [MgFe]' = sqrt(Mgb * (0.72*Fe5270 + 0.28*Fe5335))
    # We only have Fe5270 in query, let's use Mgb + Fe as proxy
    df['metal_proxy'] = df['mgb'] + df['fe5270']
    
    # 1. Control for Metallicity
    # TiO increases with Z. We want excess TiO at fixed Z.
    from sklearn.linear_model import LinearRegression
    X = df[['metal_proxy']]
    y = df['tio2']
    
    reg = LinearRegression().fit(X, y)
    print(f"  Control Fit (Z -> TiO) R2: {reg.score(X, y):.3f}")
    
    df['tio_resid'] = y - reg.predict(X)
    
    # 2. Correlate with Sigma
    df['log_sigma'] = np.log10(df['sigma'])
    
    slope, intercept, r_val, p_val, std_err = stats.linregress(df['log_sigma'], df['tio_resid'])
    
    print(f"  Correlation (log Sigma vs TiO Resid): r={r_val:.3f}, p={p_val:.2e}, slope={slope:.4f}")
    
    # Plot
    fig, ax = plt.subplots(1, 2, figsize=(12, 5))
    
    # Raw TiO vs Metal
    sc = ax[0].scatter(df['metal_proxy'], df['tio2'], c=df['log_sigma'], cmap='viridis', s=10, alpha=0.6)
    plt.colorbar(sc, ax=ax[0], label='log(Sigma)')
    ax[0].set_xlabel('Metallicity Proxy (Mg + Fe)')
    ax[0].set_ylabel('Lick TiO2 Index')
    ax[0].set_title('TiO vs Metallicity')
    
    # Residuals vs Sigma
    ax[1].scatter(df['log_sigma'], df['tio_resid'], alpha=0.5, s=10, c='maroon')
    
    x_range = np.linspace(df['log_sigma'].min(), df['log_sigma'].max(), 100)
    ax[1].plot(x_range, intercept + slope*x_range, 'k--', label=f'Slope={slope:.3f}')
    
    ax[1].set_xlabel('log(Velocity Dispersion)')
    ax[1].set_ylabel('TiO2 Residual')
    ax[1].set_title('Test DP: IMF (TiO) vs Potential')
    ax[1].legend()
    ax[1].grid(True, alpha=0.3)
    
    out_path = os.path.join(FIGURES_DIR, 'sdss_test_dp_tio_imf.png')
    plt.savefig(out_path, dpi=150)
    print(f"Figure saved to {out_path}")
    
    return {
        'correlation_r': r_val,
        'p_value': p_val,
        'slope': slope,
        'n_gal': int(len(df))
    }

def main():
    cache_path = os.path.join(DATA_DIR, 'sdss_tio_imf.csv')
    
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

    results = analyze_tio_imf(df)
    
    out_path = os.path.join(RESULTS_DIR, 'sdss_test_dp_results.json')
    if results:
        with open(out_path, 'w') as f:
            json.dump(results, f, indent=2)
            
        print("\nSUMMARY TEST DP:")
        print(f"Slope (Sigma vs TiO Resid): {results['slope']:.4f}")
        
        # Prediction: Positive Slope (Enhanced TiO in high sigma)
        if results['p_value'] < 0.05 and results['slope'] > 0.001:
             print("RESULT: SIGNAL (Bottom-heavy IMF in Deep Potentials)")
        elif results['p_value'] < 0.05 and results['slope'] < -0.001:
             print("RESULT: CONTRADICTED (Top-heavy IMF in Deep Potentials)")
        else:
             print("RESULT: NULL (No significant dependence)")

if __name__ == "__main__":
    main()
