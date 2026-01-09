#!/usr/bin/env python3
"""
Step 6.17: SDSS Test N - QSO Variability Timescale vs BH Mass

Hypothesis:
QSO variability is driven by accretion disk instabilities. The characteristic timescale correlates with BH mass.
Under TEP, time dilation near the SMBH is extreme.
If time flows slower in deeper potentials (higher BH mass or host sigma), the OBSERVED variability should be:
1. Slower (longer timescales) than predicted by standard accretion physics?
2. Or: Variability amplitude at a fixed rest-frame timescale might be suppressed?

Actually, the TEP prediction in the plan is:
"Variability timescales should be LONGER than expected from pure BH mass scaling.
The 'excess' should correlate with host sigma (total potential depth)."

Since we don't have explicit timescales (tau) in qsoVarPTF (only VAR_A, VAR_GAMMA), we can test:
1. Variability Amplitude (VAR_A) vs BH Mass / Host Sigma
2. Structure Function slope (VAR_GAMMA) vs BH Mass / Host Sigma

Hypothesis Refinement for PTF Data:
VAR_A is the amplitude of variability.
VAR_GAMMA is the structure function power law index.
If time is dilated, processes look "frozen". 
Maybe amplitude is suppressed at fixed observed cadence? 
Or maybe characteristic timescale is stretched -> variability looks "redder" (steeper slope)?

Let's test:
r(VAR_A, sigma_host) | M_BH
r(VAR_GAMMA, sigma_host) | M_BH

Data:
- spiders_quasar: BH Mass (logBHMA_hb), OIII Width (proxy for host sigma)
- qsoVarPTF: Variability (VAR_A, VAR_GAMMA)
- Overlap: ~848 objects

Join path: spiders_quasar -> SpecObjAll -> qsoVarPTF
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

def download_data():
    print(f"Querying SDSS for Test N (Spiders+PTF)...")
    
    # We use OIII width as a proxy for host velocity dispersion (sigma_stars)
    # sigma ~ FWHM_OIII / 2.35
    
    sql = """
    SELECT TOP 5000
        s.SPECOBJID,
        sp.z AS redshift,
        
        -- BH Mass (H-beta based)
        s.logBHMA_hb,
        s.errlogBHMA_hb,
        
        -- Host Sigma Proxy (OIII width)
        s.width2_OIII5007 AS OIII_width, -- 2nd component often better
        s.errwidth2_OIII5007,
        
        -- Variability from PTF
        v.VAR_A,      -- Amplitude
        v.VAR_GAMMA,  -- Slope
        v.VAR_CHI2    -- Quality
        
    FROM spiders_quasar s
    JOIN SpecObjAll sp ON s.SPECOBJID = sp.specObjID
    JOIN qsoVarPTF v ON sp.bestObjID = v.VAR_OBJID
    
    WHERE 
        s.width2_OIII5007 > 0
        AND s.logBHMA_hb > 0
        AND v.VAR_A > 0
    """
    return query_sdss(sql)

def analyze_variability(df):
    print("Analyzing Variability...")
    
    # 1. Prepare Variables
    # Sigma proxy: sigma = FWHM / 2.35. Width in spiders is likely sigma or FWHM. 
    # Usually 'width' implies sigma in Gaussian fits, but let's check. 
    # Assuming it's sigma (velocity dispersion of the line).
    df['log_sigma_host'] = np.log10(df['OIII_width'])
    
    df['log_MBH'] = df['logBHMA_hb']
    
    # Filter quality
    # High variability signal? Or just valid fits.
    # VAR_CHI2?
    
    df_clean = df.dropna(subset=['log_MBH', 'log_sigma_host', 'VAR_A']).copy()
    
    # 2. Hypothesis: Time dilation affects variability amplitude/slope?
    # TEP: High sigma host -> Slower time -> ??
    # If time is slower, variability timescales stretch. 
    # For a fixed window (PTF), stretching time might reduce the observed amplitude (red noise).
    # So r(VAR_A, sigma_host) < 0?
    # Or correlation with BH mass is standard.
    
    # First, let's control for BH Mass (standard scaling).
    X = df_clean[['log_MBH', 'redshift']].values
    y = df_clean['VAR_A'].values
    
    reg = LinearRegression().fit(X, y)
    df_clean['VAR_A_resid'] = y - reg.predict(X)
    
    # Correlation with Host Sigma
    r_simple, p_simple = stats.pearsonr(df_clean['log_sigma_host'], df_clean['VAR_A'])
    r_resid, p_resid = stats.pearsonr(df_clean['log_sigma_host'], df_clean['VAR_A_resid'])
    
    print(f"N = {len(df_clean)}")
    print(f"Simple r(Amp, Sigma_host): {r_simple:.4f} (p={p_simple:.2e})")
    print(f"Controlled r(Amp_resid, Sigma_host): {r_resid:.4f} (p={p_resid:.2e})")
    
    return {
        'r_simple': float(r_simple),
        'p_simple': float(p_simple),
        'r_controlled': float(r_resid),
        'p_controlled': float(p_resid),
        'n_sample': int(len(df_clean))
    }, df_clean

def create_figure(df, results):
    print("Generating figure...")
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # Plot 1: Amp vs Host Sigma
    ax = axes[0]
    ax.scatter(df['log_sigma_host'], df['VAR_A'], alpha=0.3, s=10, c='blue')
    m, b = np.polyfit(df['log_sigma_host'], df['VAR_A'], 1)
    x = np.linspace(df['log_sigma_host'].min(), df['log_sigma_host'].max(), 100)
    ax.plot(x, m*x + b, 'r-', label=f'r={results["r_simple"]:.3f}')
    ax.set_xlabel(r'$\log(\sigma_{host})$ (OIII Width)')
    ax.set_ylabel('Variability Amplitude (A)')
    ax.set_title("Raw Variability vs Host Potential")
    ax.legend()
    
    # Plot 2: Residual Amp vs Host Sigma
    ax = axes[1]
    ax.scatter(df['log_sigma_host'], df['VAR_A_resid'], alpha=0.3, s=10, c='green')
    m, b = np.polyfit(df['log_sigma_host'], df['VAR_A_resid'], 1)
    ax.plot(x, m*x + b, 'r-', label=f'r={results["r_controlled"]:.3f}')
    ax.set_xlabel(r'$\log(\sigma_{host})$ (OIII Width)')
    ax.set_ylabel(r'Residual Amplitude (Fixed $M_{BH}, z$)')
    ax.set_title("Controlled Variability vs Host Potential")
    ax.legend()
    
    plt.tight_layout()
    out_path = os.path.join(FIGURES_DIR, 'sdss_test_n_variability.png')
    plt.savefig(out_path, dpi=150)
    print(f"Figure saved to {out_path}")

def main():
    cache_path = os.path.join(DATA_DIR, 'sdss_qso_ptf.csv')
    
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

    results, df_clean = analyze_variability(df)
    
    out_path = os.path.join(RESULTS_DIR, 'sdss_test_n_results.json')
    with open(out_path, 'w') as f:
        json.dump(results, f, indent=2)
        
    create_figure(df_clean, results)
    
    print("\nSUMMARY TEST N:")
    print("TEP Prediction: Variability correlates with Host Potential (r != 0)")
    print(f"Observed Controlled r: {results['r_controlled']:.4f}")
    
    if abs(results['r_controlled']) > 0.1 and results['p_controlled'] < 0.05:
        print("RESULT: DETECTED (Significant Correlation)")
    else:
        print("RESULT: NULL/INCONCLUSIVE")

if __name__ == "__main__":
    main()
