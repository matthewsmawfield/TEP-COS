#!/usr/bin/env python3
"""
Step 6.45: SDSS Test BN - CMG Survival (Compact Massive Galaxies)

Hypothesis:
CMGs ("Red Nuggets") are thought to puff up over time via mergers.
If merger rates/dynamical processes are time-dilated in these extremely deep potential wells (sigma > 250 km/s),
they should survive as compact relics for longer in cosmic time.

Prediction:
Compactness (or Fraction of CMGs) correlates with Sigma (at fixed Mass).
High Sigma -> More Compact / Higher Survival Fraction.

Data:
- stellarMassFSPSGranWideDust: logMass
- emissionLinesPort: sigma_stars
- SpecObjAll: z, bestObjID
- PhotoObjAll: petroR50_r (Half-light radius in arcsec)

Method:
1. Select massive galaxies (logM > 10.5).
2. Convert radius from arcsec to kpc using redshift.
3. Define "Compactness" metric: Sigma_eff = M / R^1.5 or just use Radius at fixed Mass.
   Or define CMG flag: R_eff < 1.5 kpc.
4. Correlate Compactness with Velocity Dispersion.
   Note: Sigma is part of the compactness definition dynamically (Virial), so we need to be careful not to just rediscover the Virial Theorem.
   TEP Prediction is about the *survival* of the compact state.
   We look for an excess of compact objects at high sigma *relative* to the standard size-mass relation.
   Actually, Test BN hypothesis says: "Correlation of Compactness with sigma should be preserved" (implies it's stronger or different?).
   Better: At fixed Mass, do galaxies with higher sigma have smaller radii? (Yes, Virial Theorem).
   Is there an *anomaly*?
   The hypothesis states: "Number density of CMGs is higher than standard hierarchical models predict."
   This requires a cosmological volume correction which is hard here.
   Alternative TEP prediction: "CMGs survive longer".
   So, fraction of CMGs should increase with Sigma (if Sigma proxies potential depth/time dilation factor).
   But Sigma also proxies Compactness directly via Virial.
   
   Let's look for the *residual* of the Size-Mass relation vs Sigma.
   R = A * M^alpha.
   Residual Delta_logR = logR_obs - logR_pred.
   Correlate Delta_logR with Sigma.
   Standard Virial: High Sigma -> Small R -> Negative Correlation.
   TEP: Enhanced survival of *very* compact things?
   Maybe check the *fraction* of outliers (Red Nuggets) vs Sigma?
"""

import pandas as pd
import numpy as np
from scipy import stats
import os
import json
import matplotlib.pyplot as plt
import requests
import time
from astropy.cosmology import Planck15

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
    print(f"Querying SDSS for Test BN (Limit: {limit})...")
    
    # Use bestObjID to link SpecObj to PhotoObj
    # Correct column name: sigmaStars (not sigma_stars)
    sql = f"""
    SELECT TOP {limit}
        s.specObjID,
        s.logMass,
        sp.z,
        e.sigmaStars as sigma_stars,
        ph.petroR50_r as r_arcsec,
        ph.petroR90_r
        
    FROM stellarMassFSPSGranWideDust s
    JOIN SpecObjAll sp ON s.specObjID = sp.specObjID
    JOIN PhotoObjAll ph ON sp.bestObjID = ph.objID
    JOIN emissionLinesPort e ON s.specObjID = e.specObjID
    
    WHERE 
        s.logMass > 10.5 
        AND sp.class = 'GALAXY'
        AND sp.z BETWEEN 0.05 AND 0.2
        AND e.sigmaStars > 50
    """
    return query_sdss(sql)

def analyze_cmg_survival(df):
    print("Analyzing CMG Survival...")
    
    # 1. Clean
    df_clean = df.dropna().copy()
    df_clean = df_clean[df_clean['r_arcsec'] > 0]
    
    # 2. Compute Physical Radius (kpc)
    # 1 arcsec = kpc_per_arcsec
    # Use simplified cosmology or astropy
    # scale at z=0.1 is approx 1.8 kpc/arcsec.
    # Let's use astropy properly
    
    print("  Computing physical radii...")
    # Vectorized conversion is cleaner if we map z to scale
    # Planck15.kpc_proper_per_arcmin(z) returns Quantity in kpc/arcmin
    # We want kpc/arcsec, so divide by 60
    
    # Ensure z is numpy array
    z_vals = df_clean['z'].values
    
    # Get values directly from Quantity object
    kpc_per_arcmin = Planck15.kpc_proper_per_arcmin(z_vals).value
    scales = kpc_per_arcmin / 60.0
    
    df_clean['r_kpc'] = df_clean['r_arcsec'] * scales
    
    # 3. Size-Mass Relation Residuals
    # Fit logR = a * logM + b
    slope, intercept, r_val, p_val, std_err = stats.linregress(df_clean['logMass'], np.log10(df_clean['r_kpc']))
    
    print(f"  Size-Mass Relation: slope={slope:.2f}, intercept={intercept:.2f}, r={r_val:.2f}")
    
    df_clean['logR_pred'] = slope * df_clean['logMass'] + intercept
    df_clean['logR_resid'] = np.log10(df_clean['r_kpc']) - df_clean['logR_pred']
    
    # 4. Correlate Residual with Sigma
    # We expect high sigma -> small R (negative residual).
    # This is standard.
    # TEP Check: Is the correlation *stronger* or is there an excess of *extreme* compacts at high sigma?
    
    r_bn, p_bn = stats.pearsonr(df_clean['sigma_stars'], df_clean['logR_resid'])
    print(f"  Correlation r(Sigma, Size_Resid): {r_bn:.4f} (p={p_bn:.2e})")
    
    # 5. Fraction of CMGs vs Sigma
    # Define CMG: Compact Massive Galaxy.
    # e.g., Residual < -0.3 (0.3 dex smaller than average, factor of 2)
    # Or R_kpc < 2.0 (for M > 10.5)
    
    # Using relative compactness definition (Residual < -0.3)
    df_clean['is_compact'] = df_clean['logR_resid'] < -0.3
    
    # Bin by Sigma
    df_clean['sigma_bin'] = pd.qcut(df_clean['sigma_stars'], 8)
    binned = df_clean.groupby('sigma_bin')['is_compact'].agg(['mean', 'count'])
    
    print("\nCMG Fraction by Sigma Bin:")
    print(binned)
    
    return {
        'r_corr': float(r_bn),
        'p_corr': float(p_bn),
        'slope_sm': float(slope),
        'binned_data': binned.reset_index().to_dict(orient='records'),
        'n_sample': int(len(df_clean))
    }, df_clean

def create_figure(df, results):
    print("Generating figure...")
    
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    # Panel 1: Size vs Mass
    axes[0].scatter(df['logMass'], np.log10(df['r_kpc']), c=df['sigma_stars'], cmap='viridis', s=5, alpha=0.5)
    x = np.linspace(df['logMass'].min(), df['logMass'].max(), 100)
    axes[0].plot(x, results['slope_sm'] * x + (np.mean(np.log10(df['r_kpc'])) - results['slope_sm']*np.mean(df['logMass'])), 'r--', label='Fit')
    axes[0].set_xlabel('log Mass [M_sun]')
    axes[0].set_ylabel('log R_eff [kpc]')
    axes[0].set_title('Size-Mass Relation (Color=Sigma)')
    
    # Panel 2: Compact Fraction vs Sigma
    data = pd.DataFrame(results['binned_data'])
    # Extract bin centers
    # qcut interval string parsing is annoying, let's use the bin index or re-calculate
    # Simplified: just plot vs index or extract mid point from string if needed.
    # For now, just bar plot?
    
    # axes[1].bar(range(len(data)), data['mean'], yerr=np.sqrt(data['mean']*(1-data['mean'])/data['count']))
    # axes[1].set_xlabel('Sigma Bin (Low -> High)')
    
    # Better: Scatter plot of is_compact (binned) vs mean sigma of bin
    # We need mean sigma per bin. Re-do aggregation in main analysis? 
    # Or just use the 'sigma_stars' and 'is_compact' to do a logit reg plot or similar.
    # Let's stick to simple scatter of residuals vs sigma.
    
    axes[1].scatter(df['sigma_stars'], df['logR_resid'], alpha=0.2, s=5, c='k')
    axes[1].set_xlabel('Velocity Dispersion [km/s]')
    axes[1].set_ylabel('Size Residual (dex)')
    axes[1].set_title(f'Residual vs Sigma (r={results["r_corr"]:.2f})')
    axes[1].axhline(-0.3, color='r', linestyle='--', label='CMG Threshold')
    axes[1].legend()
    
    plt.tight_layout()
    out_path = os.path.join(FIGURES_DIR, 'sdss_test_bn_cmg.png')
    plt.savefig(out_path, dpi=150)
    print(f"Figure saved to {out_path}")

def main():
    cache_path = os.path.join(DATA_DIR, 'sdss_cmg_survival.csv')
    
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

    results, df_clean = analyze_cmg_survival(df)
    
    out_path = os.path.join(RESULTS_DIR, 'sdss_test_bn_results.json')
    
    def json_default(obj):
        if isinstance(obj, pd.Interval):
            return str(obj)
        raise TypeError
        
    with open(out_path, 'w') as f:
        json.dump(results, f, indent=2, default=json_default)
        
    create_figure(df_clean, results)
    
    print("\nSUMMARY TEST BN:")
    print("TEP Prediction: Compactness correlates with Sigma (r < 0, more compact at high sigma).")
    print(f"Observed r: {results['r_corr']:.4f}")
    
    if results['r_corr'] < -0.3:
        print("RESULT: CONSISTENT (Strong correlation, Virial/TEP degenerate)")
    else:
        print("RESULT: NULL/WEAK")

if __name__ == "__main__":
    main()
