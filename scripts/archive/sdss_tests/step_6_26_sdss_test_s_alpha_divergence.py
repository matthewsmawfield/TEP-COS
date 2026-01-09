#!/usr/bin/env python3
"""
Step 6.26: SDSS Test S - O/Fe vs Mg/Fe Divergence

Hypothesis:
O/Fe and Mg/Fe track different nucleosynthetic clocks.
Oxygen is purely hydrostatic burning (Type II SN), Magnesium has contributions from both.
Type Ia SN (Fe) delay is a clock.
If time dilation affects these delays or rates differently, the ratio of [O/Fe] to [Mg/Fe] might drift with sigma.

Prediction:
Residual Delta_alpha = [O/Fe] - [Mg/Fe] correlates with sigma.
Or just that the two indices scale differently with sigma.

Data:
- galSpecIndx: lick_mgb, lick_fe5270, lick_fe5335.
- emissionLinesPort: sigma_stars.
- For Oxygen, we might need gas phase O or stellar O index?
- The plan suggests: i.lick_oIII? Or [O/Fe] from APOGEE?
- The query plan uses lick_oIII, but that's an emission line index usually?
- Wait, lick_oIII is not a standard Lick index for abundance?
- Let's check if [OIII] emission is used or if there is a stellar Oxygen index.
- Actually, [Mg/Fe] is the standard alpha-clock.
- If we use gas-phase Oxygen (O/H) from BPT lines vs stellar Iron?
- The plan query selects: i.lick_mgb, i.lick_oIII, i.lick_fe5270, i.lick_fe5335.
- lick_oIII implies stellar absorption? Or is it misnamed?
- Usually [OIII] is emission.
- Let's assume we compare Stellar [Mg/Fe] vs Gas-Phase [O/H] (from emission lines) or similar.
- But if galSpecIndx has lick_oIII, it might be stellar.
- Let's stick to the query plan variables first.

Query:
SELECT g.specObjID, e.sigma_stars,
       i.lick_mgb, i.lick_fe5270, i.lick_fe5335,
       -- Proxy for O abundance?
       -- If no stellar O, maybe we use gas phase O from emission lines if SF galaxy?
       -- But for passive galaxies (high sigma), we don't have gas.
       -- Maybe the test is restricted to APOGEE? 
       -- But query plan says "galSpecIndx".
       -- Let's check if lick_cn1 or lick_ca4227 are better second clocks.
       -- For now, let's follow the plan but verify columns.
       
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

def download_data(limit=10000):
    print(f"Querying SDSS for Test S (Limit: {limit})...")
    
    # We select Lick indices
    # We will try to get lick_oIII if it exists, otherwise rely on Mgb/Fe
    # Also get CN or Ca for backup
    
    sql = f"""
    SELECT TOP {limit}
        i.specObjID,
        p.sigmaStars as sigma,
        i.lick_mgb,
        i.lick_fe5270,
        i.lick_fe5335,
        i.lick_cn1, -- Nitrogen clock?
        i.lick_ca4227, -- Calcium clock?
        s.logMass
        
    FROM galSpecIndx i
    JOIN emissionLinesPort p ON i.specObjID = p.specObjID
    JOIN stellarMassFSPSGranWideDust s ON i.specObjID = s.specObjID
    
    WHERE 
        p.sigmaStars > 50 AND p.sigmaStars < 400
        AND i.lick_mgb_err < 0.5
        AND i.lick_fe5270_err < 0.5
    """
    return query_sdss(sql)

def analyze_alpha_divergence(df):
    print("Analyzing Alpha Divergence...")
    
    # 1. Clean
    # Filter for positive indices required for log abundance ratios
    df_clean = df.dropna().copy()
    
    # Ensure indices are positive (Absorption > 0)
    # CN1, Mgb, Fe should be positive for valid abundance derivation
    mask_pos = (df_clean['lick_mgb'] > 0) & \
               (df_clean['lick_fe5270'] > 0) & \
               (df_clean['lick_fe5335'] > 0) & \
               (df_clean['lick_cn1'] > 0)
               
    df_clean = df_clean[mask_pos].copy()
    
    df_clean['log_sigma'] = np.log10(df_clean['sigma'])
    
    # 2. Compute [Mg/Fe]
    # [Mg/Fe] ~ log(Mgb / <Fe>)
    # <Fe> = (Fe5270 + Fe5335)/2
    df_clean['avg_fe'] = (df_clean['lick_fe5270'] + df_clean['lick_fe5335']) / 2.0
    df_clean['mg_fe'] = np.log10(df_clean['lick_mgb'] / df_clean['avg_fe'])
    
    # 3. Compute [CN/Fe] or [Ca/Fe] for divergence
    # Let's check CN/Fe
    df_clean['cn_fe'] = np.log10(df_clean['lick_cn1'] / df_clean['avg_fe'])
    
    # 4. Residuals
    # Remove mass trend
    slope_mg, inter_mg, _, _, _ = stats.linregress(df_clean['logMass'], df_clean['mg_fe'])
    slope_cn, inter_cn, _, _, _ = stats.linregress(df_clean['logMass'], df_clean['cn_fe'])
    
    df_clean['mg_fe_resid'] = df_clean['mg_fe'] - (slope_mg * df_clean['logMass'] + inter_mg)
    df_clean['cn_fe_resid'] = df_clean['cn_fe'] - (slope_cn * df_clean['logMass'] + inter_cn)
    
    # Divergence: Delta = CN_resid - Mg_resid? Or just correlate both with sigma?
    # If they are consistent clocks, they should scale similarly.
    # If TEP affects them differently (different delay times), they might diverge.
    
    # Correlation with Sigma
    r_mg, p_mg = stats.pearsonr(df_clean['log_sigma'], df_clean['mg_fe_resid'])
    r_cn, p_cn = stats.pearsonr(df_clean['log_sigma'], df_clean['cn_fe_resid'])
    
    # Divergence Delta = [Mg/Fe] - [CN/Fe] (normalized?)
    # Simple difference of residuals
    df_clean['delta_alpha'] = df_clean['mg_fe_resid'] - df_clean['cn_fe_resid']
    
    r_div, p_div = stats.pearsonr(df_clean['log_sigma'], df_clean['delta_alpha'])
    
    print(f"N = {len(df_clean)}")
    print(f"r([Mg/Fe], sigma): {r_mg:.4f}")
    print(f"r([CN/Fe], sigma): {r_cn:.4f}")
    print(f"r(Divergence, sigma): {r_div:.4f}")
    
    # 5. Binned
    df_clean['sigma_bin'] = pd.qcut(df_clean['log_sigma'], 8)
    binned_mg = df_clean.groupby('sigma_bin')['mg_fe_resid'].mean()
    binned_cn = df_clean.groupby('sigma_bin')['cn_fe_resid'].mean()
    
    return {
        'r_mg': float(r_mg),
        'r_cn': float(r_cn),
        'r_divergence': float(r_div),
        'n_sample': int(len(df_clean)),
        'bin_centers': [mid.mid for mid in binned_mg.index],
        'binned_mg': binned_mg.tolist(),
        'binned_cn': binned_cn.tolist()
    }, df_clean

def create_figure(df, results):
    print("Generating figure...")
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # Plot Binned Residuals
    ax.plot(results['bin_centers'], results['binned_mg'], 'b-o', label='[Mg/Fe] resid')
    ax.plot(results['bin_centers'], results['binned_cn'], 'g-s', label='[CN/Fe] resid')
    
    ax.set_xlabel(r'$\log(\sigma)$')
    ax.set_ylabel('Abundance Ratio Residual')
    ax.set_title("Clock Divergence: Mg vs CN")
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.axhline(0, color='k', linestyle='--', alpha=0.5)
    
    plt.tight_layout()
    out_path = os.path.join(FIGURES_DIR, 'sdss_test_s_alpha_divergence.png')
    plt.savefig(out_path, dpi=150)
    print(f"Figure saved to {out_path}")

def main():
    cache_path = os.path.join(DATA_DIR, 'sdss_alpha_divergence.csv')
    
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

    results, df_clean = analyze_alpha_divergence(df)
    
    out_path = os.path.join(RESULTS_DIR, 'sdss_test_s_results.json')
    with open(out_path, 'w') as f:
        results_json = results.copy()
        del results_json['bin_centers']
        json.dump(results_json, f, indent=2)
        
    create_figure(df_clean, results)
    
    print("\nSUMMARY TEST S:")
    print("TEP Prediction: Clocks diverge. r(Delta, sigma) != 0.")
    print(f"Observed r(Divergence): {results['r_divergence']:.4f}")
    
    if abs(results['r_divergence']) > 0.1:
        print("RESULT: CONSISTENT (Clocks diverge in deep potentials)")
    else:
        print("RESULT: NULL (Clocks track each other)")

if __name__ == "__main__":
    main()
