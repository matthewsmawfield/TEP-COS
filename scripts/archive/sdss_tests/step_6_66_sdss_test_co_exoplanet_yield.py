#!/usr/bin/env python3
"""
Step 6.66: SDSS Test CO - Exoplanet Yield (Migration Rates)

Hypothesis:
Planetary migration (Type I/II) is a rate-dependent process set by disk viscosity and torque timescales. 
In the deep potential of metal-rich stars (often Inner Galaxy, though MARVELS is local, metallicity correlates with kinematic origin), 
these rates are dilated. Planets might stall at different radii or survive migration differently. 
The distribution of minMass (M sin i) or Period for MARVELS companions should show dependence on host Fe/H (proxy for potential depth/origin).

Prediction:
Exoplanet Mass Function or Period Distribution varies with Galactic Potential (proxied by Metallicity).

Data:
- marvelsStar: starname, ra, dec, feh (Metallicity), logg
- marvelsVelocityCurveUF1D: minMass, period, eccentricity

Method:
1. Fetch MARVELS candidates with valid mass/period.
2. Join on starname.
3. Analyze Period vs Metallicity and Mass vs Metallicity.
   Are massive planets (Jupiters) more common around high-Z stars? (Known).
   Does the Period distribution shift? TEP predicts migration stalls earlier or later?
   Slower migration -> Planets stay further out? Or less efficient inward migration?
   Expect positive correlation between Period and Potential Depth (High Z)?
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
    print(f"Querying SDSS for Test CO (Limit: {limit})...")
    
    # MARVELS data
    
    sql = f"""
    SELECT TOP {limit}
        s.starname,
        s.ra, s.dec,
        s.FeH as fe_h,
        s.logg,
        v.minMass,
        v.period,
        v.eccentricity
        
    FROM marvelsStar s
    JOIN marvelsVelocityCurveUF1D v ON s.starname = v.starname
    WHERE 
        v.minMass > 0 
        AND v.period > 0
        AND s.FeH > -9
    """
    return query_sdss(sql)

def analyze_exoplanets(df):
    print("Analyzing Exoplanet Yield...")
    
    # Clean
    df = df.dropna().copy()
    
    # Log variables
    df['log_mass'] = np.log10(df['minMass']) # Jupiter masses? Usually M_Jup
    df['log_period'] = np.log10(df['period']) # Days?
    
    print(f"  Sample size: {len(df)}")
    print(f"  Mean Fe/H: {df['fe_h'].mean():.2f}")
    
    if len(df) < 10:
        print("  Insufficient sample size.")
        return {'n_sample': len(df), 'r_period': 0, 'r_mass': 0}
        
    # Correlations with Metallicity (Potential Proxy)
    # 1. Period vs Fe/H
    r_per, p_per = stats.pearsonr(df['fe_h'], df['log_period'])
    print(f"  Correlation r(Period, Fe/H): {r_per:.4f} (p={p_per:.2e})")
    
    # 2. Mass vs Fe/H
    r_mass, p_mass = stats.pearsonr(df['fe_h'], df['log_mass'])
    print(f"  Correlation r(Mass, Fe/H): {r_mass:.4f} (p={p_mass:.2e})")
    
    # Plot
    fig, ax = plt.subplots(1, 2, figsize=(12, 5))
    
    # Period
    ax[0].scatter(df['fe_h'], df['log_period'], alpha=0.5, s=20)
    # Fit line
    slope_p, int_p, _, _, _ = stats.linregress(df['fe_h'], df['log_period'])
    x_range = np.linspace(df['fe_h'].min(), df['fe_h'].max(), 100)
    ax[0].plot(x_range, slope_p * x_range + int_p, 'r-', label=f'Slope={slope_p:.2f}')
    ax[0].set_xlabel('[Fe/H]')
    ax[0].set_ylabel('log Period [days]')
    ax[0].set_title(f'Period vs Metallicity (r={r_per:.2f})')
    ax[0].legend()
    ax[0].grid(True, alpha=0.3)
    
    # Mass
    ax[1].scatter(df['fe_h'], df['log_mass'], alpha=0.5, s=20)
    slope_m, int_m, _, _, _ = stats.linregress(df['fe_h'], df['log_mass'])
    ax[1].plot(x_range, slope_m * x_range + int_m, 'r-', label=f'Slope={slope_m:.2f}')
    ax[1].set_xlabel('[Fe/H]')
    ax[1].set_ylabel('log Mass [M_Jup]')
    ax[1].set_title(f'Mass vs Metallicity (r={r_mass:.2f})')
    ax[1].legend()
    ax[1].grid(True, alpha=0.3)
    
    out_path = os.path.join(FIGURES_DIR, 'sdss_test_co_exoplanet.png')
    plt.savefig(out_path, dpi=150)
    print(f"Figure saved to {out_path}")
    
    return {
        'r_period': r_per,
        'slope_period': slope_p,
        'r_mass': r_mass,
        'n_sample': int(len(df))
    }

def main():
    cache_path = os.path.join(DATA_DIR, 'sdss_exoplanet.csv')
    
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

    results = analyze_exoplanets(df)
    
    out_path = os.path.join(RESULTS_DIR, 'sdss_test_co_results.json')
    if results:
        with open(out_path, 'w') as f:
            json.dump(results, f, indent=2)
            
        print("\nSUMMARY TEST CO:")
        print("Prediction: Exoplanet Period varies with Potential (Fe/H).")
        print(f"Observed r(Period, Fe/H): {results['r_period']:.4f}")
        
        if abs(results['r_period']) > 0.2:
             print("RESULT: SIGNAL (Correlation observed)")
        else:
             print("RESULT: NULL (Weak/No correlation)")

if __name__ == "__main__":
    main()
