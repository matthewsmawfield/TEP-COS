#!/usr/bin/env python3
"""
Step 6.18: SDSS Test AL - Cluster-Field Age Inversion

Hypothesis:
Standard hierarchical assembly predicts that galaxies in dense environments (clusters) formed earlier and quenched earlier than field galaxies of the same mass. 
Prediction (Standard): D4000 correlates POSITIVELY with Density at fixed Mass. (Cluster = Old).

TEP Hypothesis:
Galaxies in high-density regions reside in deeper gravitational potentials.
Time flows slower in these potentials.
Prediction (TEP): D4000 correlates NEGATIVELY with Density at fixed Mass (or the positive correlation is suppressed). 
(Cluster = Younger appearance due to time dilation).

Data:
- ebossMCPM: Environmental density (MATTERDENS).
- galSpecIndx: Age (D4000).
- stellarMassFSPSGranWideDust: Stellar Mass.
- emissionLinesPort: Sigma (for checking potential depth).

Join: ebossMCPM (Plate, MJD, FiberID) -> SpecObjAll -> others.
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

def download_data(limit=50000):
    print(f"Querying SDSS for Test AL (Limit: {limit})...")
    
    # MATTERDENS is log density? Or linear? 
    # Usually matter density in eBOSS MCPM is relative to mean? 
    # We will check distribution.
    
    sql = f"""
    SELECT TOP {limit}
        s.specObjID,
        m.Z as redshift,
        m.MATTERDENS as density,
        i.d4000_n as D4000,
        i.d4000_n_err,
        st.logMass,
        e.sigmaStars as sigma
        
    FROM ebossMCPM m
    JOIN SpecObjAll s ON m.PLATE = s.plate AND m.MJD = s.mjd AND m.FIBERID = s.fiberID
    JOIN galSpecIndx i ON s.specObjID = i.specObjID
    JOIN stellarMassFSPSGranWideDust st ON s.specObjID = st.specObjID
    JOIN emissionLinesPort e ON s.specObjID = e.specObjID
    
    WHERE 
        m.Z BETWEEN 0.02 AND 0.20
        AND i.d4000_n > 1.0 AND i.d4000_n < 3.0
        AND st.logMass > 9.0
        AND m.MATTERDENS > -90
    """
    return query_sdss(sql)

def analyze_inversion(df):
    print("Analyzing Cluster-Field Age Inversion...")
    
    # 1. Clean Data
    df_clean = df.dropna(subset=['D4000', 'density', 'logMass', 'sigma']).copy()
    print(f"N = {len(df_clean)}")
    
    # 2. Variables
    # Density: High MATTERDENS = Cluster, Low = Field/Void
    # Age: D4000
    # Mass: logMass
    
    # 3. Control for Mass
    # Standard Age-Mass relation (Downsizing): Massive galaxies are older.
    # We want to see the effect of ENVIRONMENT at fixed Mass.
    
    X = df_clean[['logMass']].values
    y = df_clean['D4000'].values
    
    reg = LinearRegression().fit(X, y)
    df_clean['D4000_resid'] = y - reg.predict(X)
    
    # 4. Correlation with Density
    # Standard: r(D4000_resid, density) > 0 (Assembly Bias)
    # TEP: r < 0 (Time Dilation)
    
    r_simple, p_simple = stats.pearsonr(df_clean['density'], df_clean['D4000'])
    r_controlled, p_controlled = stats.pearsonr(df_clean['density'], df_clean['D4000_resid'])
    
    print(f"Simple r(D4000, Density): {r_simple:.4f} (p={p_simple:.2e})")
    print(f"Controlled r(D4000_resid, Density): {r_controlled:.4f} (p={p_controlled:.2e})")
    
    # 5. Check against Sigma as well (Local Potential)
    r_sigma, p_sigma = stats.pearsonr(df_clean['sigma'], df_clean['D4000_resid'])
    print(f"Controlled r(D4000_resid, Sigma): {r_sigma:.4f} (p={p_sigma:.2e})")
    
    return {
        'r_simple': float(r_simple),
        'p_simple': float(p_simple),
        'r_controlled': float(r_controlled),
        'p_controlled': float(p_controlled),
        'r_sigma': float(r_sigma),
        'n_sample': int(len(df_clean))
    }, df_clean

def create_figure(df, results):
    print("Generating figure...")
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # Plot 1: D4000 vs Density (Raw)
    # Use log density for plotting if range is large
    # Check range
    min_d, max_d = df['density'].min(), df['density'].max()
    print(f"Density Range: {min_d} to {max_d}")
    
    x_col = df['density']
    if max_d > 100:
        # Likely linear overdensity, log it
        x_plot = np.log10(np.maximum(0.1, df['density']))
        xlabel = r'$\log(\rho)$'
    else:
        x_plot = df['density']
        xlabel = r'Density ($\rho$)'
        
    ax = axes[0]
    ax.scatter(x_plot, df['D4000'], alpha=0.1, s=2, c='blue')
    m, b = np.polyfit(x_plot, df['D4000'], 1)
    x_range = np.linspace(x_plot.min(), x_plot.max(), 100)
    ax.plot(x_range, m*x_range + b, 'r-', label=f'r={results["r_simple"]:.3f}')
    ax.set_xlabel(xlabel)
    ax.set_ylabel('D4000 (Age)')
    ax.set_title("Raw Age-Density Relation")
    ax.legend()
    
    # Plot 2: Residual D4000 vs Density
    ax = axes[1]
    ax.scatter(x_plot, df['D4000_resid'], alpha=0.1, s=2, c='green')
    m, b = np.polyfit(x_plot, df['D4000_resid'], 1)
    ax.plot(x_range, m*x_range + b, 'r-', label=f'r={results["r_controlled"]:.3f}')
    ax.set_xlabel(xlabel)
    ax.set_ylabel('Residual D4000 (Mass Controlled)')
    ax.set_title("Environment vs Age (Fixed Mass)")
    ax.legend()
    
    plt.tight_layout()
    out_path = os.path.join(FIGURES_DIR, 'sdss_test_al_cluster_inversion.png')
    plt.savefig(out_path, dpi=150)
    print(f"Figure saved to {out_path}")

def main():
    cache_path = os.path.join(DATA_DIR, 'sdss_eboss_env.csv')
    
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

    results, df_clean = analyze_inversion(df)
    
    out_path = os.path.join(RESULTS_DIR, 'sdss_test_al_results.json')
    with open(out_path, 'w') as f:
        json.dump(results, f, indent=2)
        
    create_figure(df_clean, results)
    
    print("\nSUMMARY TEST AL:")
    print("Standard Prediction: r > 0 (Dense = Old)")
    print("TEP Prediction: r < 0 (Dense = Young)")
    print(f"Observed Controlled r: {results['r_controlled']:.4f}")
    
    if results['r_controlled'] < -0.05:
        print("RESULT: CONSISTENT with TEP (Inversion Observed)")
    elif results['r_controlled'] > 0.05:
        print("RESULT: CONTRADICTED (Standard Assembly Bias Dominates)")
    else:
        print("RESULT: NULL/FLAT")

if __name__ == "__main__":
    main()
