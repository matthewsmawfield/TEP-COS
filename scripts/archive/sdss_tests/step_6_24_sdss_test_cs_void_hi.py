#!/usr/bin/env python3
"""
Step 6.24: SDSS Test CS - Void HI Gas Fraction (The Void Clock)

Hypothesis:
In TEP, time runs "fast" in voids (shallow potential) compared to the field.
Standard Lambda-CDM predicts void galaxies are retarded in evolution (gas-rich).
TEP predicts they experience more proper time per unit cosmic time, potentially appearing MORE EVOLVED (gas-poor) than field galaxies of the same mass.

Prediction:
HI Gas Fraction in Voids is LOWER than expected from the standard density relation.

Data:
- mangaHIall: logMHI (HI Mass).
- mangaDAPall: logmass (Stellar Mass).
- ebossMCPM: MATTERDENS (Environmental Density).

Method:
1. Join mangaHIall to mangaDAPall (by mangaid) and ebossMCPM (spatial).
2. Compute Gas Fraction: f_HI = M_HI / M_star.
3. Define Voids: Low percentile of MATTERDENS (e.g., < 20th percentile).
4. Compare mean f_HI in Voids vs Field at fixed Stellar Mass.
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

def download_data(limit=500):
    print(f"Querying SDSS for Test CS (Limit: {limit})...")
    
    # Optimizing spatial join:
    # Use BETWEEN for indices which is often faster than ABS
    
    sql = f"""
    SELECT TOP {limit}
        h.mangaid,
        h.logMHI,
        h.logmstars as logmass,
        m.MATTERDENS as density
        
    FROM mangaHIall h
    JOIN ebossMCPM m ON m.RA BETWEEN h.objra - 0.001 AND h.objra + 0.001
                    AND m.DEC BETWEEN h.objdec - 0.001 AND h.objdec + 0.001
    
    WHERE 
        h.logMHI > 0 
        AND h.logmstars > 0
    """
    return query_sdss(sql)

def analyze_void_fraction(df):
    print("Analyzing Void HI Fraction...")
    
    # 1. Clean
    df_clean = df.dropna().copy()
    print(f"N = {len(df_clean)}")
    
    # 2. Compute Gas Fraction (log)
    # log(f_HI) = log(M_HI) - log(M_star)
    df_clean['log_fHI'] = df_clean['logMHI'] - df_clean['logmass']
    
    # 3. Define Environment
    # Inspect density distribution
    print(f"Density Range: {df_clean['density'].min():.3f} - {df_clean['density'].max():.3f}")
    
    # Define Void as bottom 20%
    void_threshold = df_clean['density'].quantile(0.20)
    cluster_threshold = df_clean['density'].quantile(0.80)
    
    df_clean['env'] = 'Field'
    df_clean.loc[df_clean['density'] < void_threshold, 'env'] = 'Void'
    df_clean.loc[df_clean['density'] > cluster_threshold, 'env'] = 'Dense'
    
    print(f"Void Threshold (20%): {void_threshold:.3f}")
    
    # 4. Compare Means (controlling for mass?)
    # Simple comparison first
    means = df_clean.groupby('env')['log_fHI'].mean()
    print("\nMean log(f_HI) by Environment:")
    print(means)
    
    void_mean = means.get('Void', np.nan)
    field_mean = means.get('Field', np.nan)
    dense_mean = means.get('Dense', np.nan)
    
    diff_void_field = void_mean - field_mean
    
    # Ttest
    void_vals = df_clean[df_clean['env']=='Void']['log_fHI']
    field_vals = df_clean[df_clean['env']=='Field']['log_fHI']
    
    t_stat, p_val = stats.ttest_ind(void_vals, field_vals, equal_var=False)
    
    print(f"Void vs Field: Delta = {diff_void_field:.3f} dex")
    print(f"T-test: t={t_stat:.2f}, p={p_val:.2e}")
    
    # 5. Mass-Correction (Residuals)
    # Remove mass trend
    slope, intercept, _, _, _ = stats.linregress(df_clean['logmass'], df_clean['log_fHI'])
    df_clean['fHI_resid'] = df_clean['log_fHI'] - (slope * df_clean['logmass'] + intercept)
    
    means_resid = df_clean.groupby('env')['fHI_resid'].mean()
    print("\nMean Mass-Corrected Residuals:")
    print(means_resid)
    
    resid_diff = means_resid.get('Void', 0) - means_resid.get('Field', 0)
    
    return {
        'void_mean_logfHI': float(void_mean),
        'field_mean_logfHI': float(field_mean),
        'diff_raw': float(diff_void_field),
        'diff_resid': float(resid_diff),
        'p_value': float(p_val),
        'n_void': int(len(void_vals)),
        'n_field': int(len(field_vals)),
        'void_threshold': float(void_threshold)
    }, df_clean

def create_figure(df, results):
    print("Generating figure...")
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # Boxplot
    envs = ['Void', 'Field', 'Dense']
    data = [df[df['env']==e]['log_fHI'].values for e in envs]
    
    ax.boxplot(data, labels=envs, patch_artist=True)
    
    ax.set_ylabel(r'$\log(f_{HI}) = \log(M_{HI}/M_{*})$')
    ax.set_title("HI Gas Fraction by Environment")
    ax.grid(True, alpha=0.3)
    
    # Annotate TEP prediction
    # TEP: Void < Field (More evolved)
    # Standard: Void > Field (Less evolved)
    
    plt.tight_layout()
    out_path = os.path.join(FIGURES_DIR, 'sdss_test_cs_void_hi.png')
    plt.savefig(out_path, dpi=150)
    print(f"Figure saved to {out_path}")

def main():
    cache_path = os.path.join(DATA_DIR, 'sdss_void_hi.csv')
    
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

    results, df_clean = analyze_void_fraction(df)
    
    out_path = os.path.join(RESULTS_DIR, 'sdss_test_cs_results.json')
    with open(out_path, 'w') as f:
        json.dump(results, f, indent=2)
        
    create_figure(df_clean, results)
    
    print("\nSUMMARY TEST CS:")
    print("TEP Prediction: Voids have LOWER HI fraction (Fast time = Evolved).")
    print("Standard Prediction: Voids have HIGHER HI fraction (Retarded evolution).")
    print(f"Observed Diff (Void - Field): {results['diff_raw']:.3f} dex")
    
    if results['diff_raw'] < -0.1:
        print("RESULT: CONSISTENT (Void galaxies are gas-poor/evolved)")
    elif results['diff_raw'] > 0.1:
        print("RESULT: CONTRADICTED (Void galaxies are gas-rich)")
    else:
        print("RESULT: NULL/FLAT")

if __name__ == "__main__":
    main()
