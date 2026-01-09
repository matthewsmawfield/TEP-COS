#!/usr/bin/env python3
"""
Step 6.20: SDSS Test AJ - Green Valley Dwell Time

Hypothesis:
The transition from Star Forming (Blue Cloud) to Quiescent (Red Sequence) is a rate-limited process (quenching).
TEP predicts that time flows slower in deep potentials (high sigma).
If quenching is a physical process with a characteristic timescale (e.g. gas depletion), then in deep potentials, this process should appear SLOWER to an outside observer.
Prediction: Galaxies in high-sigma environments spend LONGER traversing the Green Valley.
Observable: The fraction of galaxies in the Green Valley (f_GV) should be HIGHER at high sigma (controlling for Mass).

Data:
- galSpecExtra: sSFR (specsfr_tot_p50).
- emissionLinesPort: Sigma Stars.
- stellarMassFSPSGranWideDust: Stellar Mass.

Definition:
Green Valley defined by specific SFR (sSFR).
Blue: log sSFR > -10.5
Red: log sSFR < -11.5
GV: -11.5 < log sSFR < -10.5
(Values approximate, will adjust based on distribution).
"""

import pandas as pd
import numpy as np
from scipy import stats
import os
import json
import matplotlib.pyplot as plt
import requests
import time
from sklearn.linear_model import LogisticRegression

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
    print(f"Querying SDSS for Test AJ (Limit: {limit})...")
    
    sql = f"""
    SELECT TOP {limit}
        e.specObjID,
        e.specsfr_tot_p50 as log_sSFR, -- This is usually log10(sSFR)
        p.sigmaStars,
        s.logMass,
        sp.z
        
    FROM galSpecExtra e
    JOIN emissionLinesPort p ON e.specObjID = p.specObjID
    JOIN stellarMassFSPSGranWideDust s ON e.specObjID = s.specObjID
    JOIN SpecObjAll sp ON e.specObjID = sp.specObjID
    
    WHERE 
        sp.z BETWEEN 0.02 AND 0.15
        AND p.sigmaStars > 50 AND p.sigmaStars < 400
        AND s.logMass > 9.0
        AND e.specsfr_tot_p50 > -14 -- Valid measurements
    """
    return query_sdss(sql)

def analyze_green_valley(df):
    print("Analyzing Green Valley Dwell Time...")
    
    # 1. Clean
    df_clean = df.dropna(subset=['log_sSFR', 'sigmaStars', 'logMass']).copy()
    
    # 2. Define Green Valley
    # Let's inspect distribution briefly
    # Typically: Blue > -10.8, Red < -11.8 roughly?
    # Or -11 to -10?
    
    # We define GV broadly: -11.8 < log sSFR < -10.8
    # (Checking standard literature definitions for SDSS)
    GV_LOW = -11.8
    GV_HIGH = -10.8
    
    df_clean['is_GV'] = ((df_clean['log_sSFR'] > GV_LOW) & (df_clean['log_sSFR'] < GV_HIGH)).astype(int)
    
    print(f"N Total: {len(df_clean)}")
    print(f"N GV: {df_clean['is_GV'].sum()}")
    print(f"Global Fraction: {df_clean['is_GV'].mean():.3f}")
    
    # 3. Bin by Sigma and Calculate Fraction
    df_clean['log_sigma'] = np.log10(df_clean['sigmaStars'])
    
    bins = np.linspace(df_clean['log_sigma'].min(), df_clean['log_sigma'].max(), 10)
    bin_centers = 0.5 * (bins[1:] + bins[:-1])
    fractions = []
    
    for i in range(len(bins)-1):
        mask = (df_clean['log_sigma'] >= bins[i]) & (df_clean['log_sigma'] < bins[i+1])
        if mask.sum() > 50:
            fractions.append(df_clean.loc[mask, 'is_GV'].mean())
        else:
            fractions.append(np.nan)
            
    # 4. Logistic Regression to control for Mass
    # P(is_GV) ~ a*log_sigma + b*logMass + c*z
    # We want coefficient of log_sigma to be positive (Higher sigma -> Higher probability of being in GV)
    
    X = df_clean[['log_sigma', 'logMass', 'z']].values
    y = df_clean['is_GV'].values
    
    clf = LogisticRegression(class_weight='balanced', solver='liblinear')
    clf.fit(X, y)
    
    coefs = clf.coef_[0]
    names = ['log_sigma', 'logMass', 'z']
    
    print("Logistic Regression Coefs:")
    for n, c in zip(names, coefs):
        print(f"  {n}: {c:.4f}")
        
    sigma_coef = coefs[0]
    
    # Also simple correlation on binned data?
    # Let's use the coefficient significance roughly or just the sign/magnitude.
    
    return {
        'sigma_coef': float(sigma_coef),
        'fractions': fractions, # list of floats
        'bin_centers': bin_centers.tolist(),
        'n_sample': int(len(df_clean))
    }, df_clean

def create_figure(df, results):
    print("Generating figure...")
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # Plot GV Fraction vs Sigma
    ax.plot(results['bin_centers'], results['fractions'], 'bo-', label='Observed Fraction')
    
    ax.set_xlabel(r'$\log(\sigma)$')
    ax.set_ylabel('Green Valley Fraction ($f_{GV}$)')
    ax.set_title("Green Valley Dwell Time vs Potential Depth")
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    out_path = os.path.join(FIGURES_DIR, 'sdss_test_aj_green_valley.png')
    plt.savefig(out_path, dpi=150)
    print(f"Figure saved to {out_path}")

def main():
    cache_path = os.path.join(DATA_DIR, 'sdss_green_valley.csv')
    
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

    results, df_clean = analyze_green_valley(df)
    
    out_path = os.path.join(RESULTS_DIR, 'sdss_test_aj_results.json')
    with open(out_path, 'w') as f:
        json.dump(results, f, indent=2)
        
    create_figure(df_clean, results)
    
    print("\nSUMMARY TEST AJ:")
    print("TEP Prediction: f_GV increases with sigma (Coef > 0)")
    print(f"Logistic Coef (log_sigma): {results['sigma_coef']:.4f}")
    
    if results['sigma_coef'] > 0.1:
        print("RESULT: CONSISTENT (Higher GV fraction in deep potentials)")
    elif results['sigma_coef'] < -0.1:
        print("RESULT: CONTRADICTED (Lower GV fraction in deep potentials)")
    else:
        print("RESULT: NULL/FLAT")

if __name__ == "__main__":
    main()
