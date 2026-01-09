#!/usr/bin/env python3
"""
Step 6.49: SDSS Test BU - Spiral Arm Winding (The Winding Clock)

Hypothesis:
Spiral arms wind up over time due to differential rotation. In deep potentials, dynamical time is slower.
For a given age, spiral arms in high-sigma galaxies should appear "less wound" (looser, larger pitch angle)
than in low-sigma galaxies, contradicting the standard trend where massive disks usually have tighter arms.

Prediction:
Fraction of "Tight" arms decreases as sigma increases (at fixed mass/color).
Or: Mean winding score (0=Tight, 1=Medium, 2=Loose) increases with sigma.

Data:
- MaNGA_GZ2: Spiral winding probabilities (t10_arms_winding_a28_tight_fraction, etc.)
- mangaDAPall: stellar_sigma_1re
- mangaTarget: nsa_elpetro_mass (Stellar Mass)

Method:
1. Fetch Winding Fractions and Sigma.
2. Clean data (Confirmed spirals only).
3. Compute Winding Score: W = P(Loose) * 2 + P(Medium) * 1 + P(Tight) * 0.
   High W = Looser Arms.
4. Analyze Correlation r(W, sigma).
5. Bin by Sigma and plot.
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

def download_data(limit=5000):
    print(f"Querying SDSS for Test BU (Limit: {limit})...")
    
    # MaNGA_GZ2 uses MANGAID. mangaDAPall uses mangaid.
    # Join on MANGAID.
    
    sql = f"""
    SELECT TOP {limit}
        g.MANGAID,
        g.t10_arms_winding_a28_tight_fraction as f_tight,
        g.t10_arms_winding_a29_medium_fraction as f_medium,
        g.t10_arms_winding_a30_loose_fraction as f_loose,
        g.t04_spiral_a08_spiral_fraction as f_spiral,
        
        d.stellar_sigma_1re as sigma,
        s.nsa_elpetro_mass as logmass
        
    FROM MaNGA_GZ2 g
    JOIN mangaDAPall d ON g.MANGAID = d.mangaid
    JOIN mangaTarget s ON g.MANGAID = s.mangaid
    
    WHERE 
        d.stellar_sigma_1re > 50 AND d.stellar_sigma_1re < 400
        AND g.t04_spiral_a08_spiral_fraction > 0.5
        AND d.drp3qual = 0
    """
    return query_sdss(sql)

def analyze_winding(df):
    print("Analyzing Spiral Winding...")
    
    # Clean
    df = df.dropna().copy()
    
    # Compute Winding Score (0=Tight, 1=Medium, 2=Loose)
    # Using probability weighted average
    df['winding_score'] = (df['f_tight'] * 0.0 + df['f_medium'] * 1.0 + df['f_loose'] * 2.0) / (df['f_tight'] + df['f_medium'] + df['f_loose'])
    
    # Clean up any division by zero (though spiral frac > 0.5 should prevent this)
    df = df.dropna(subset=['winding_score'])
    
    # Variables
    df['log_sigma'] = np.log10(df['sigma'])
    
    # 1. Raw Correlation
    r_raw, p_raw = stats.pearsonr(df['log_sigma'], df['winding_score'])
    print(f"  Correlation r(Winding, sigma): {r_raw:.4f} (p={p_raw:.2e})")
    print("  (Positive r = Looser arms at high sigma)")
    
    # 2. Control for Mass?
    # Massive galaxies tend to have tighter arms (lower score).
    # Expected correlation is NEGATIVE in standard model.
    # TEP predicts POSITIVE (or less negative) - "Looser than expected".
    
    r_mass, p_mass = stats.pearsonr(df['logmass'], df['winding_score'])
    print(f"  Correlation r(Winding, Mass): {r_mass:.4f}")
    
    # Partial correlation
    # r_xy.z
    def partial_corr(x, y, z):
        c_xy = stats.pearsonr(x, y)[0]
        c_xz = stats.pearsonr(x, z)[0]
        c_yz = stats.pearsonr(y, z)[0]
        return (c_xy - c_xz * c_yz) / np.sqrt((1 - c_xz**2) * (1 - c_yz**2))
        
    r_partial = partial_corr(df['log_sigma'], df['winding_score'], df['logmass'])
    print(f"  Partial Correlation r(Winding, sigma | Mass): {r_partial:.4f}")
    
    # 3. Binning
    df['sigma_bin'] = pd.qcut(df['sigma'], 6)
    binned = df.groupby('sigma_bin')['winding_score'].agg(['mean', 'sem', 'count'])
    binned['sigma_center'] = [i.mid for i in binned.index]
    
    print("\nWinding Score by Sigma Bin:")
    print(binned[['mean', 'sem', 'count']])
    
    # Plot
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # Scatter
    # ax.scatter(df['sigma'], df['winding_score'], alpha=0.1, s=5, c='gray')
    
    # Binned
    ax.errorbar(binned['sigma_center'], binned['mean'], yerr=binned['sem'], fmt='o-', capsize=5, label='Mean Winding Score')
    
    ax.set_xlabel('Velocity Dispersion [km/s]')
    ax.set_ylabel('Spiral Winding Score (0=Tight, 2=Loose)')
    ax.set_title(f'Test BU: Spiral Winding vs Sigma (r={r_raw:.2f})')
    ax.grid(True, alpha=0.3)
    ax.legend()
    
    out_path = os.path.join(FIGURES_DIR, 'sdss_test_bu_winding.png')
    plt.savefig(out_path, dpi=150)
    print(f"Figure saved to {out_path}")
    
    return {
        'r_raw': r_raw,
        'r_partial': r_partial,
        'n_sample': int(len(df))
    }

def main():
    cache_path = os.path.join(DATA_DIR, 'sdss_spiral_winding.csv')
    
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

    results = analyze_winding(df)
    
    out_path = os.path.join(RESULTS_DIR, 'sdss_test_bu_results.json')
    with open(out_path, 'w') as f:
        json.dump(results, f, indent=2)
        
    print("\nSUMMARY TEST BU:")
    print("Standard Model: High Mass/Sigma -> Tighter Arms (Winding Score decreases). r < 0.")
    print("TEP Prediction: High Sigma -> Looser Arms (Winding Score increases relative to standard). r > r_standard.")
    print(f"Observed r: {results['r_raw']:.4f}")
    
    if results['r_raw'] > 0:
        print("RESULT: CONSISTENT (Looser arms in deep potentials)")
    elif results['r_raw'] > -0.1:
        print("RESULT: WEAK/NULL")
    else:
        print("RESULT: CONTRADICTED (Standard tightening observed)")

if __name__ == "__main__":
    main()
