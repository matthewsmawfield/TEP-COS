#!/usr/bin/env python3
"""
Step 6.81: SDSS Test DI - Cluster vs Field Stellar Spin

Hypothesis:
Stars spin down over time due to magnetic braking. If this process is time-dilated 
in deep cluster potentials, stars in clusters should appear "younger" (rapidly rotating) 
compared to field stars of the same age/type.

Prediction:
Mean v_sin(i) is higher in high-density environments (at fixed Type/Teff).

Data:
- sppParams: elodi_vsini, elodi_teff, elodi_fe_h
- ebossMCPM: mid_dens_1 (Density)

Method:
1. Select stars with reliable parameters (Teff > 5000, G/F dwarfs).
2. Match to environmental density.
3. Control for Teff and Metallicity (spin down depends on mass/structure).
4. Analyze residuals of v_sin(i) vs Density.
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
    print(f"Querying SDSS for Test DI (Limit: {limit})...")
    
    # Select stars (class check implicit via sppParams existence usually, but can check specObj)
    # Linking sppParams to ebossMCPM via specObjID
    
    sql = f"""
    SELECT TOP {limit}
        s.specObjID, 
        s.elodi_vsini as vsini, 
        s.elodi_teff as teff,
        s.elodi_fe_h as feh,
        e.mid_dens_1 as density
    FROM sppParams s
    JOIN ebossMCPM e ON s.specObjID = e.specObjID
    WHERE s.elodi_vsini > 0
      AND s.elodi_teff BETWEEN 5000 AND 7000 -- FG stars
      AND abs(e.mid_dens_1) < 10
    """
    return query_sdss(sql)

def analyze_stellar_spin(df):
    print("Analyzing Stellar Spin vs Density...")
    
    if df is None or len(df) < 100:
        print("  Insufficient data.")
        return None
        
    df = df.dropna().copy()
    
    print(f"  Data points: {len(df)}")
    
    # Control for Teff and FeH
    # v_rot depends strongly on Spectral Type (Mass)
    # Model: vsini = f(Teff, FeH)
    
    # Remove high outliers (very fast rotators) to avoid skewing by young stars?
    # Or keep them? 
    # Let's clip extreme values
    df = df[df['vsini'] < 100]
    
    from sklearn.linear_model import LinearRegression
    from sklearn.preprocessing import PolynomialFeatures
    
    X = df[['teff', 'feh']]
    y = df['vsini']
    
    # Use polynomial features for control (spin vs Teff is non-linear)
    poly = PolynomialFeatures(degree=2)
    X_poly = poly.fit_transform(X)
    
    reg = LinearRegression().fit(X_poly, y)
    print(f"  Control Fit (Teff, FeH -> vsini) R2: {reg.score(X_poly, y):.3f}")
    
    df['vsini_resid'] = y - reg.predict(X_poly)
    
    # Correlate with Density
    slope, intercept, r_val, p_val, std_err = stats.linregress(df['density'], df['vsini_resid'])
    
    print(f"  Correlation (Density vs vsini Residual): r={r_val:.3f}, p={p_val:.2e}, slope={slope:.4f}")
    
    # Plot
    fig, ax = plt.subplots(1, 2, figsize=(12, 5))
    
    # Raw Teff vs vsini
    ax[0].scatter(df['teff'], df['vsini'], c=df['density'], cmap='viridis', s=5, alpha=0.5)
    ax[0].set_xlabel('Teff (K)')
    ax[0].set_ylabel('v sin i (km/s)')
    ax[0].set_title('Stellar Rotation (Color=Density)')
    
    # Residuals vs Density
    # Bin density for clearer trend
    df['dens_bin'] = pd.qcut(df['density'], q=10)
    binned = df.groupby('dens_bin').agg({'density':'mean', 'vsini_resid':['mean', 'sem']}).reset_index()
    binned.columns = ['bin', 'density', 'resid_mean', 'resid_sem']
    
    ax[1].errorbar(binned['density'], binned['resid_mean'], yerr=binned['resid_sem'], fmt='o-', color='purple')
    ax[1].axhline(0, color='k', linestyle=':')
    ax[1].set_xlabel('Environment Density (delta)')
    ax[1].set_ylabel('vsini Residual (km/s)')
    ax[1].set_title('Test DI: Spin vs Environment')
    ax[1].grid(True, alpha=0.3)
    
    out_path = os.path.join(FIGURES_DIR, 'sdss_test_di_stellar_spin.png')
    plt.savefig(out_path, dpi=150)
    print(f"Figure saved to {out_path}")
    
    return {
        'correlation_r': r_val,
        'p_value': p_val,
        'slope': slope,
        'n_stars': int(len(df))
    }

def main():
    cache_path = os.path.join(DATA_DIR, 'sdss_stellar_spin.csv')
    
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

    results = analyze_stellar_spin(df)
    
    out_path = os.path.join(RESULTS_DIR, 'sdss_test_di_results.json')
    if results:
        with open(out_path, 'w') as f:
            json.dump(results, f, indent=2)
            
        print("\nSUMMARY TEST DI:")
        print(f"Slope (Density vs Spin): {results['slope']:.4f}")
        
        if results['p_value'] < 0.05 and results['slope'] > 0.1:
             print("RESULT: SIGNAL (Higher Spin in High Density)")
        elif results['p_value'] < 0.05 and results['slope'] < -0.1:
             print("RESULT: CONTRADICTED (Lower Spin in High Density)")
        else:
             print("RESULT: NULL (No significant dependence)")

if __name__ == "__main__":
    main()
