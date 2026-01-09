#!/usr/bin/env python3
"""
Step 6.27: SDSS Test AT - Alpha Bimodality (The Valley of Time)

Hypothesis:
The "gap" between high-alpha (thick disk/early) and low-alpha (thin disk/late) populations depends on the timescale of the transition (gas accretion/infall).
If this timescale scales with sigma, the distinctness (depth) of the bimodal valley in the [alpha/Fe] vs [Fe/H] plane should vary with velocity dispersion.

Prediction:
Bimodality Gap Depth (or separation) varies with sigma.
TEP: Slower transition -> Gap might be wider or filled differently?
Usually, faster evolution = distinct populations. Slower = more mixing?
The prediction is simply that the morphology of the bimodality is potential-dependent beyond standard mass trends.

Data:
- galSpecIndx: Mgb, Fe5270, Fe5335 -> [Mg/Fe]
- emissionLinesPort: sigma_stars
- galSpecInfo: reliability
- stellarMass...: logMass

Method:
1. Calculate [Mg/Fe] and [Fe/H] (proxy from <Fe> index? Or use Mgb/<Fe> vs <Fe>).
2. Bin by sigma.
3. In each bin, fit a Double Gaussian to the [Mg/Fe] distribution.
4. Measure the "Valley Depth" (Peak-to-Valley ratio) or Separation.
5. Correlate with sigma.
"""

import pandas as pd
import numpy as np
from scipy import stats
from sklearn.mixture import GaussianMixture
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
    print(f"Querying SDSS for Test AT (Limit: {limit})...")
    
    sql = f"""
    SELECT TOP {limit}
        g.specObjID, 
        e.sigmaStars as sigma,
        i.lick_mgb, 
        i.lick_fe5270, 
        i.lick_fe5335,
        s.logMass
        
    FROM galSpecIndx i
    JOIN emissionLinesPort e ON i.specObjID = e.specObjID
    JOIN galSpecInfo g ON i.specObjID = g.specObjID
    JOIN stellarMassFSPSGranWideDust s ON i.specObjID = s.specObjID
    
    WHERE 
        e.sigmaStars > 50 AND e.sigmaStars < 400
        AND i.lick_mgb > 0 AND i.lick_fe5270 > 0 AND i.lick_fe5335 > 0
        AND i.lick_mgb_err < 0.5
    """
    return query_sdss(sql)

def analyze_bimodality(df):
    print("Analyzing Alpha Bimodality...")
    
    # 1. Clean
    df_clean = df.dropna().copy()
    df_clean['log_sigma'] = np.log10(df_clean['sigma'])
    
    # 2. Compute Proxies
    # <Fe>
    df_clean['avg_fe'] = (df_clean['lick_fe5270'] + df_clean['lick_fe5335']) / 2.0
    # [Mg/Fe] proxy
    df_clean['alpha_fe'] = np.log10(df_clean['lick_mgb'] / df_clean['avg_fe'])
    # [Fe/H] proxy (very rough, just <Fe> strength)
    df_clean['fe_h'] = np.log10(df_clean['avg_fe'])
    
    # Remove outliers
    df_clean = df_clean[(df_clean['alpha_fe'] > -0.5) & (df_clean['alpha_fe'] < 0.5)]
    
    # 3. Bin by Sigma
    # We want to see if the distribution of alpha_fe changes shape
    bins = np.linspace(df_clean['log_sigma'].min(), df_clean['log_sigma'].max(), 6)
    
    results_list = []
    
    plt.figure(figsize=(10, 6))
    
    for i in range(len(bins)-1):
        mask = (df_clean['log_sigma'] >= bins[i]) & (df_clean['log_sigma'] < bins[i+1])
        sub = df_clean[mask]
        
        if len(sub) < 100:
            continue
            
        # Fit GMM (2 components)
        X = sub[['alpha_fe']].values
        gmm = GaussianMixture(n_components=2, random_state=42)
        gmm.fit(X)
        
        means = gmm.means_.flatten()
        weights = gmm.weights_.flatten()
        
        # Sort by mean (Low Alpha, High Alpha)
        idx = np.argsort(means)
        means = means[idx]
        weights = weights[idx]
        
        separation = means[1] - means[0]
        
        # Valley Depth? 
        # Evaluate PDF at means and at midpoint
        x_grid = np.linspace(sub['alpha_fe'].min(), sub['alpha_fe'].max(), 100).reshape(-1, 1)
        log_prob = gmm.score_samples(x_grid)
        pdf = np.exp(log_prob)
        
        peak1 = np.max(pdf[x_grid.flatten() < (means[0]+means[1])/2])
        peak2 = np.max(pdf[x_grid.flatten() > (means[0]+means[1])/2])
        valley = np.min(pdf[(x_grid.flatten() > means[0]) & (x_grid.flatten() < means[1])])
        
        valley_ratio = valley / min(peak1, peak2) # Closer to 0 = Deep Valley, Closer to 1 = Flat
        
        results_list.append({
            'sigma_bin': (bins[i] + bins[i+1])/2,
            'separation': separation,
            'valley_ratio': valley_ratio,
            'weight_ratio': weights[1]/weights[0] # High/Low
        })
        
        # Plot for visual check (last bin)
        if i == len(bins)-2:
            plt.hist(sub['alpha_fe'], bins=30, density=True, alpha=0.3, label=f'Sigma Bin {i}')
            plt.plot(x_grid, pdf, label='GMM Fit')

    plt.xlabel('[Mg/Fe] Proxy')
    plt.legend()
    plt.savefig(os.path.join(FIGURES_DIR, 'sdss_test_at_dist_check.png'))
    plt.close()
    
    res_df = pd.DataFrame(results_list)
    print("\nBimodality Analysis:")
    print(res_df)
    
    # Correlation
    r_sep, p_sep = stats.pearsonr(res_df['sigma_bin'], res_df['separation'])
    r_val, p_val = stats.pearsonr(res_df['sigma_bin'], res_df['valley_ratio'])
    
    print(f"r(Separation, sigma): {r_sep:.4f}")
    print(f"r(ValleyRatio, sigma): {r_val:.4f}")
    
    return {
        'r_separation': float(r_sep),
        'r_valley': float(r_val),
        'binned_data': res_df.to_dict(orient='records')
    }, df_clean

def create_summary_figure(res_df, results):
    print("Generating summary figure...")
    fig, ax1 = plt.subplots(figsize=(8, 6))
    
    df = pd.DataFrame(results['binned_data'])
    
    color = 'tab:red'
    ax1.set_xlabel(r'$\log(\sigma)$')
    ax1.set_ylabel('Bimodality Separation', color=color)
    ax1.plot(df['sigma_bin'], df['separation'], color=color, marker='o')
    ax1.tick_params(axis='y', labelcolor=color)
    
    ax2 = ax1.twinx() 
    color = 'tab:blue'
    ax2.set_ylabel('Valley Ratio (Lower = Deeper)', color=color)  
    ax2.plot(df['sigma_bin'], df['valley_ratio'], color=color, marker='s')
    ax2.tick_params(axis='y', labelcolor=color)
    
    plt.title("Test AT: Alpha Bimodality vs Potential")
    fig.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, 'sdss_test_at_bimodality.png'), dpi=150)

def main():
    cache_path = os.path.join(DATA_DIR, 'sdss_alpha_bimodality.csv')
    
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

    results, df_clean = analyze_bimodality(df)
    
    out_path = os.path.join(RESULTS_DIR, 'sdss_test_at_results.json')
    with open(out_path, 'w') as f:
        json.dump(results, f, indent=2)
        
    create_summary_figure(None, results)
    
    print("\nSUMMARY TEST AT:")
    print(f"Separation Correlation: {results['r_separation']:.4f}")
    print(f"Valley Depth Correlation: {results['r_valley']:.4f}")
    
    if abs(results['r_separation']) > 0.5 or abs(results['r_valley']) > 0.5:
        print("RESULT: CONSISTENT (Bimodality morphology changes with potential)")
    else:
        print("RESULT: NULL/FLAT")

if __name__ == "__main__":
    main()
