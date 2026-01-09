#!/usr/bin/env python3
"""
Step 6.16: SDSS Test J - Ionization Equilibrium Timescale

Hypothesis:
Emission line ratios depend on the ionization equilibrium of the gas.
High-ionization species (e.g., OIII, NeIII) might equilibrate on different timescales compared to low-ionization species (NII, SII) or recombination lines.
Under TEP, if proper time flows slower in deeper potentials (high sigma), the gas might appear "out of equilibrium" or show systematic shifts in ionization parameter U relative to standard photoionization models at fixed metallicity/SED.

TEP Prediction:
r(Ionization Parameter, sigma) != 0 at fixed stellar metallicity.
Specifically, if time is slower, maybe recombination is delayed? Or ionization?
The prediction in the plan is: High-sigma galaxies appear "out of equilibrium" -> deviations in [OIII]/[NII] or [OIII]/[OII].

Data:
- emissionLinesPort: Fluxes for OII, OIII, NII, SII, NeIII, Ha, Hb.
- galSpecIndx: Stellar metallicity (Mgb).
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

def download_data(limit=10000):
    print(f"Querying SDSS for Test J (Limit: {limit})...")
    # Note: Using Flux_OII_3728 based on introspection, likely sum or doublet component
    sql = f"""
    SELECT TOP {limit}
        p.specObjID,
        p.z,
        p.sigmaStars,
        
        -- High Ionization
        p.Flux_OIII_5006,
        p.Flux_OIII_5006_Err,
        p.Flux_OII_3726,
        p.Flux_OII_3728, -- Using 3728 as seen in check
        p.Flux_NeIII_3868,
        
        -- Low Ionization
        p.Flux_NII_6583,
        p.Flux_SII_6716,
        p.Flux_SII_6730,
        p.Flux_OI_6300,
        
        -- Balmer
        p.Flux_Ha_6562,
        p.Flux_Hb_4861,
        
        -- Stellar Metallicity
        i.lick_mgb,
        i.lick_fe5270,
        i.lick_fe5335,
        
        -- BPT Class from galSpecExtra
        e.bptclass
        
    FROM emissionLinesPort p
    JOIN galSpecIndx i ON p.specObjID = i.specObjID
    JOIN stellarMassFSPSGranWideDust s ON p.specObjID = s.specObjID
    JOIN galSpecExtra e ON p.specObjID = e.specObjID
    
    WHERE 
        p.z BETWEEN 0.02 AND 0.15
        AND p.sigmaStars > 50 AND p.sigmaStars < 400
        -- Strong detections for ratios
        AND p.Flux_OIII_5006 > 0 
        AND p.Flux_NII_6583 > 0
        AND p.Flux_Ha_6562 > 0 
        AND p.Flux_Hb_4861 > 0
        -- Star forming only to avoid AGN physics dominating (bptclass 1 = Star Forming)
        AND e.bptclass = 1
    """
    return query_sdss(sql)

def analyze_ionization(df):
    print("Analyzing Ionization Equilibrium...")
    
    # 1. Prepare Data
    # Sum OII doublet if possible, or use what we have. 
    # If 3726 and 3728 are both present, sum them.
    df['Flux_OII'] = df['Flux_OII_3726'] + df['Flux_OII_3728']
    
    # Calculate Ratios (log10)
    # [OIII]/[OII] - Ionization Parameter proxy (U)
    # Avoid div zero
    df = df[df['Flux_OII'] > 0].copy()
    df['log_O32'] = np.log10(df['Flux_OIII_5006'] / df['Flux_OII'])
    
    # [OIII]/[NII] - Hardness? Metallicity?
    df['log_O3N2'] = np.log10(df['Flux_OIII_5006'] / df['Flux_NII_6583'])
    
    # [NII]/[SII] - Metallicity sensitive, ionization insensitive?
    df['Flux_SII'] = df['Flux_SII_6716'] + df['Flux_SII_6730']
    df = df[df['Flux_SII'] > 0].copy()
    df['log_N2S2'] = np.log10(df['Flux_NII_6583'] / df['Flux_SII'])
    
    # Stellar Metallicity
    df['Fe_avg'] = (df['lick_fe5270'] + df['lick_fe5335']) / 2
    df['MgFe'] = np.sqrt(np.maximum(0, df['lick_mgb'] * df['Fe_avg']))
    
    df['log_sigma'] = np.log10(df['sigmaStars'])
    
    # 2. Analysis: Ionization Parameter (O32) vs Sigma
    # Control for Metallicity (Gas phase or Stellar)
    # O32 depends on U and Z. 
    # Stellar MgFe is a proxy for Z.
    
    df_clean = df.dropna(subset=['log_O32', 'log_sigma', 'MgFe']).copy()
    
    X = df_clean[['MgFe']].values
    y = df_clean['log_O32'].values
    
    reg = LinearRegression().fit(X, y)
    df_clean['log_O32_resid'] = y - reg.predict(X)
    
    r_simple, p_simple = stats.pearsonr(df_clean['log_sigma'], df_clean['log_O32'])
    r_resid, p_resid = stats.pearsonr(df_clean['log_sigma'], df_clean['log_O32_resid'])
    
    print(f"N = {len(df_clean)}")
    print(f"Simple r(log O32, sigma): {r_simple:.4f} (p={p_simple:.2e})")
    print(f"Controlled r(log O32 resid, sigma): {r_resid:.4f} (p={p_resid:.2e})")
    
    # 3. Analysis: O3N2 vs Sigma (Check consistency)
    df_clean2 = df.dropna(subset=['log_O3N2', 'log_sigma', 'MgFe']).copy()
    X2 = df_clean2[['MgFe']].values
    y2 = df_clean2['log_O3N2'].values
    reg2 = LinearRegression().fit(X2, y2)
    df_clean2['log_O3N2_resid'] = y2 - reg2.predict(X2)
    
    r_o3n2, p_o3n2 = stats.pearsonr(df_clean2['log_sigma'], df_clean2['log_O3N2_resid'])
    print(f"Controlled r(log O3N2 resid, sigma): {r_o3n2:.4f} (p={p_o3n2:.2e})")

    return {
        'r_simple': float(r_simple),
        'p_simple': float(p_simple),
        'r_controlled': float(r_resid),
        'p_controlled': float(p_resid),
        'r_o3n2_resid': float(r_o3n2),
        'p_o3n2_resid': float(p_o3n2),
        'n_sample': int(len(df_clean))
    }, df_clean

def create_figure(df, results):
    print("Generating figure...")
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # Plot 1: O32 vs Sigma
    ax = axes[0]
    ax.scatter(df['log_sigma'], df['log_O32'], alpha=0.1, s=2, c='blue')
    m, b = np.polyfit(df['log_sigma'], df['log_O32'], 1)
    x = np.linspace(df['log_sigma'].min(), df['log_sigma'].max(), 100)
    ax.plot(x, m*x + b, 'r-', label=f'r={results["r_simple"]:.3f}')
    ax.set_xlabel(r'$\log(\sigma)$')
    ax.set_ylabel(r'$\log([OIII]/[OII])$ (Ionization Param)')
    ax.set_title("Raw Ionization vs Sigma")
    ax.legend()
    
    # Plot 2: Residual O32 vs Sigma
    ax = axes[1]
    ax.scatter(df['log_sigma'], df['log_O32_resid'], alpha=0.1, s=2, c='green')
    m, b = np.polyfit(df['log_sigma'], df['log_O32_resid'], 1)
    ax.plot(x, m*x + b, 'r-', label=f'r={results["r_controlled"]:.3f}')
    ax.set_xlabel(r'$\log(\sigma)$')
    ax.set_ylabel(r'Residual $\log([OIII]/[OII])$ (Fixed Stellar Z)')
    ax.set_title("Controlled Ionization vs Sigma")
    ax.legend()
    
    plt.tight_layout()
    out_path = os.path.join(FIGURES_DIR, 'sdss_test_j_ionization.png')
    plt.savefig(out_path, dpi=150)
    print(f"Figure saved to {out_path}")

def main():
    cache_path = os.path.join(DATA_DIR, 'sdss_ionization_data.csv')
    
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

    results, df_clean = analyze_ionization(df)
    
    out_path = os.path.join(RESULTS_DIR, 'sdss_test_j_results.json')
    with open(out_path, 'w') as f:
        json.dump(results, f, indent=2)
        
    create_figure(df_clean, results)
    
    print("\nSUMMARY TEST J:")
    print("TEP Prediction: Deviations in Ionization Equilibrium (r != 0)")
    print(f"Observed Controlled r: {results['r_controlled']:.4f}")
    
    if abs(results['r_controlled']) > 0.05:
        print("RESULT: DETECTED (Correlation exists)")
    else:
        print("RESULT: NULL")

if __name__ == "__main__":
    main()
