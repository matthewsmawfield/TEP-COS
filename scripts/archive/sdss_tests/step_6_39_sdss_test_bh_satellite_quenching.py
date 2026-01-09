#!/usr/bin/env python3
"""
Step 6.39: SDSS Test BH - Satellite Quenching Timescale

Hypothesis:
Satellite galaxies falling into clusters quench due to ram pressure and strangulation (rate processes).
In deep cluster potentials, these rates should be time-dilated.
Satellites should remain star-forming for *longer* (closer to the cluster core) than predicted by standard quenching models.
Prediction: Fraction of Star-Forming satellites at small R_proj is higher in high-mass clusters (deep potential) vs low-mass groups, or simply higher than standard models predict.
Since we lack a direct standard model comparison here, we look for the trend with potential depth.
Actually, standard theory predicts FASTER quenching in high-mass clusters (Ram pressure ~ rho v^2).
TEP predicts SLOWER quenching (Time dilation).
So if we find Quenched Fraction is LOWER (SF Fraction HIGHER) in deep potentials at fixed phase-space coordinates, that supports TEP.
However, density is higher in massive clusters.
Let's look at SF Fraction vs Projected Distance (R_proj).
TEP: SF fraction stays high until smaller R_proj.

Data:
- ebossMCPM: mid_dens_1 (Density), z.
- galSpecExtra: sfr_tot_p50.
- stellarMass...: logMass.
- SpecObjAll: z (for distance).

Method:
1. Select Cluster Members (High Density).
2. Compute Projected Distance (if cluster center known? We don't have cluster centers easily).
   Alternatively, use Local Density as proxy for R_proj.
   High Density ~ Small R_proj.
3. Compute SF Fraction vs Density.
   Standard: SF Fraction drops rapidly with Density.
   TEP: SF Fraction drops *more slowly* or stays higher at high density compared to standard?
   This is degenerate with efficiency of stripping.
   
   Alternative Plan: Compare SF Fraction at fixed Density for different Velocity Dispersions?
   If we can find cluster sigma.
   
   Let's stick to the Query Plan:
   "Fraction of Star-Forming satellites at small R_proj is higher in high-mass clusters."
   
   We need cluster mass/sigma. ebossMCPM has 'mid_dens_1'.
   Maybe we just look at SF fraction vs Density.
   If TEP is strong, maybe we see SF galaxies in very high density regions?
   
   Let's plot SF Fraction vs Density.
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

def download_data(limit=200):
    print(f"Querying SDSS for Test BH (Limit: {limit})...")
    
    # Select galaxies in dense environments
    # mid_dens_1 > 0 (Overdense)
    
    sql = f"""
    SELECT TOP {limit}
        m.specObjID,
        m.mid_dens_1 as density,
        e.sfr_tot_p50 as log_sfr,
        s.logMass,
        sp.z
        
    FROM ebossMCPM m
    JOIN galSpecExtra e ON m.specObjID = e.specObjID
    JOIN stellarMassFSPSGranWideDust s ON m.specObjID = s.specObjID
    JOIN SpecObjAll sp ON m.specObjID = sp.specObjID
    
    WHERE 
        m.mid_dens_1 > 0
        AND s.logMass > 9.0
        AND e.sfr_tot_p50 > -99
    """
    return query_sdss(sql)

def analyze_quenching(df):
    print("Analyzing Satellite Quenching...")
    
    # 1. Clean
    df_clean = df.dropna().copy()
    
    # 2. Define Star Forming
    # sSFR = SFR / Mass
    # log sSFR = log SFR - log Mass
    df_clean['log_ssfr'] = df_clean['log_sfr'] - df_clean['logMass']
    
    # Quenched definition: log sSFR < -11.0 (approx)
    df_clean['is_sf'] = df_clean['log_ssfr'] > -11.0
    
    # 3. Bin by Density
    # Density is log overdensity? mid_dens_1 description?
    # Usually ranges from -1 to 2?
    
    # Bin density
    df_clean['dens_bin'] = pd.qcut(df_clean['density'], 10)
    
    binned = df_clean.groupby('dens_bin')['is_sf'].agg(['mean', 'count'])
    binned['sem'] = np.sqrt(binned['mean'] * (1 - binned['mean']) / binned['count'])
    
    print("\nSF Fraction by Density:")
    print(binned)
    
    # 4. Correlation
    # Standard: SF Fraction decreases with Density (r < 0).
    # TEP: If suppression is weaker, r might be less negative? 
    # But we check the slope.
    
    r_sf, p_sf = stats.pearsonr(df_clean['density'], df_clean['is_sf'])
    print(f"Correlation r(SF, Density): {r_sf:.4f} (p={p_sf:.2e})")
    
    return {
        'r_sf': float(r_sf),
        'p_sf': float(p_sf),
        'binned_data': binned.reset_index().to_dict(orient='records'),
        'n_sample': int(len(df_clean))
    }, df_clean

def create_figure(df, results):
    print("Generating figure...")
    
    # Re-bin for plotting
    bins = np.linspace(df['density'].min(), df['density'].max(), 10)
    centers = (bins[:-1] + bins[1:]) / 2
    
    fracs = []
    errs = []
    
    for i in range(len(bins)-1):
        sub = df[(df['density'] >= bins[i]) & (df['density'] < bins[i+1])]
        if len(sub) > 10:
            frac = sub['is_sf'].mean()
            err = np.sqrt(frac * (1-frac) / len(sub))
            fracs.append(frac)
            errs.append(err)
        else:
            fracs.append(np.nan)
            errs.append(np.nan)
            
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.errorbar(centers, fracs, yerr=errs, fmt='o-', capsize=3, label='SF Fraction')
    
    ax.set_xlabel(r'Environment Density (log $\delta$)')
    ax.set_ylabel(r'Star Forming Fraction')
    ax.set_title(f"Test BH: Quenching vs Density (r={results['r_sf']:.3f})")
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    out_path = os.path.join(FIGURES_DIR, 'sdss_test_bh_satellite_quenching.png')
    plt.savefig(out_path, dpi=150)
    print(f"Figure saved to {out_path}")

def main():
    cache_path = os.path.join(DATA_DIR, 'sdss_satellite_quenching.csv')
    
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

    results, df_clean = analyze_quenching(df)
    
    out_path = os.path.join(RESULTS_DIR, 'sdss_test_bh_results.json')
    
    # Helper to serialize Intervals
    def json_default(obj):
        if isinstance(obj, pd.Interval):
            return str(obj)
        raise TypeError
        
    with open(out_path, 'w') as f:
        json.dump(results, f, indent=2, default=json_default)
        
    create_figure(df_clean, results)
    
    print("\nSUMMARY TEST BH:")
    print("TEP Prediction: Slower quenching (SF Fraction higher than expected at high density).")
    print("Standard: Strong negative correlation.")
    print(f"Observed r: {results['r_sf']:.4f}")
    
    if results['r_sf'] > -0.1:
        print("RESULT: CONSISTENT (Quenching suppressed in high density)")
    else:
        print("RESULT: CONTRADICTED (Standard environmental quenching dominates)")

if __name__ == "__main__":
    main()
