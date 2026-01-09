#!/usr/bin/env python3
"""
Step 6.33: SDSS Test BA - Massive Fast Rotators (Kinematic Preservation)

Hypothesis:
In standard Lambda-CDM, massive galaxies are slow rotators (formed by dry mergers).
TEP allows for "phantom mass" which might support rotation without the dynamical heating associated with dark matter halos.
Alternatively, time dilation might slow the dynamical friction process that spins down galaxies.
We expect a higher fraction of massive fast rotators than standard predictions.

Prediction:
Fraction of Fast Rotators (lambda_R > 0.1 or V/sigma > X) at high Mass is higher than standard predictions.
Standard: Fraction drops significantly above logM ~ 11.0.

Data:
- mangaDAPall: stellar_vel, stellar_sigma.
- mangaTarget: nsa_elpetro_mass.

Method:
1. Join DAPall and Target.
2. Calculate kinematic proxy:
   V_rot = stellar_vel_hi_clip (Max rotation).
   Sigma = stellar_sigma_1re (Central dispersion).
   Proxy Lambda ~ V / sqrt(V^2 + Sigma^2).
   Or just V/Sigma.
3. Define Fast Rotator: Lambda_proxy > 0.2 (approx).
4. Bin by Stellar Mass.
5. Calculate Fraction of Fast Rotators in each bin.
6. Check trend at high mass.
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

def download_data(limit=2000):
    print(f"Querying SDSS for Test BA (Limit: {limit})...")
    
    sql = f"""
    SELECT TOP {limit}
        d.mangaid,
        d.stellar_vel_hi_clip as v_star,
        d.stellar_sigma_1re as sigma,
        t.nsa_elpetro_mass as logmass,
        d.nsa_sersic_ba as axis_ratio
        
    FROM mangaDAPall d
    JOIN mangaTarget t ON d.mangaid = t.mangaid
    
    WHERE 
        d.drp3qual = 0
        AND d.stellar_sigma_1re > 0
        AND t.nsa_elpetro_mass > 9.0 -- Focus on moderate to high mass
    """
    return query_sdss(sql)

def analyze_fast_rotators(df):
    print("Analyzing Fast Rotators...")
    
    # 1. Clean
    df_clean = df.dropna().copy()
    
    # 2. Compute Kinematic Proxy
    # Lambda_R proxy ~ V / sqrt(V^2 + sigma^2)
    # Use absolute V
    v_abs = np.abs(df_clean['v_star'])
    df_clean['lambda_proxy'] = v_abs / np.sqrt(v_abs**2 + df_clean['sigma']**2)
    
    # Define Fast Rotator
    # Emsellem et al (2007): Lambda_R > 0.1 (Slow/Fast boundary is roughly 0.1 to 0.2 depending on ellipticity)
    # Let's use 0.2 as a safe cut for "Fast".
    cut = 0.2
    df_clean['is_fast'] = df_clean['lambda_proxy'] > cut
    
    # 3. Bin by Mass
    bins = np.linspace(9.0, 12.0, 10)
    df_clean['mass_bin'] = pd.cut(df_clean['logmass'], bins)
    
    binned = df_clean.groupby('mass_bin')['is_fast'].agg(['mean', 'count'])
    binned['sem'] = np.sqrt(binned['mean'] * (1 - binned['mean']) / binned['count'])
    
    print("\nFast Rotator Fraction by Mass:")
    print(binned)
    
    # 4. Check Trend
    # Correlation of lambda_proxy with Mass
    r_lam, p_lam = stats.pearsonr(df_clean['logmass'], df_clean['lambda_proxy'])
    print(f"Correlation r(Lambda, Mass): {r_lam:.4f} (p={p_lam:.2e})")
    
    return {
        'r_lambda': float(r_lam),
        'p_lambda': float(p_lam),
        'binned_data': binned.reset_index().to_dict(orient='records'),
        'n_sample': int(len(df_clean))
    }, df_clean

def create_figure(results):
    print("Generating figure...")
    data = results['binned_data']
    
    m_centers = []
    fracs = []
    errs = []
    
    for bin_data in data:
        interval = bin_data['mass_bin']
        # Interval is a string or Interval object?
        # In dict it might be preserved or str
        # If it's from cut, it's an Interval.
        # Let's reconstruct centers roughly
        if isinstance(interval, pd.Interval):
            center = interval.mid
        else:
            # Fallback if serialization changed it (should handle json serialization carefully)
            # For now, let's just use the bin index in plotting or assume linear
            pass
            
    # Actually, let's use the dataframe we printed
    # Re-extract from results not easy if Interval is serialized to str in JSON
    # Better to plot inside analysis or use simple centers
    
    # Let's just create a quick plot using matplotlib directly
    pass

def create_figure_direct(df, results):
    print("Generating figure...")
    
    # Re-bin for plotting
    bins = np.linspace(9.0, 12.0, 10)
    centers = (bins[:-1] + bins[1:]) / 2
    
    fracs = []
    errs = []
    
    for i in range(len(bins)-1):
        sub = df[(df['logmass'] >= bins[i]) & (df['logmass'] < bins[i+1])]
        if len(sub) > 10:
            frac = sub['is_fast'].mean()
            err = np.sqrt(frac * (1-frac) / len(sub))
            fracs.append(frac)
            errs.append(err)
        else:
            fracs.append(np.nan)
            errs.append(np.nan)
            
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.errorbar(centers, fracs, yerr=errs, fmt='o-', capsize=3, label='Observed Fraction')
    
    ax.set_xlabel(r'$\log(M_{*}/M_{\odot})$')
    ax.set_ylabel(r'Fraction of Fast Rotators ($\lambda_{proxy} > 0.2$)')
    ax.set_title(f"Test BA: Fast Rotators vs Mass")
    ax.grid(True, alpha=0.3)
    ax.set_ylim(0, 1.0)
    
    # Standard expectation (Atlas3D): Fraction drops above 11.0
    ax.axvline(11.0, color='r', linestyle='--', label='Transition Mass')
    ax.legend()
    
    plt.tight_layout()
    out_path = os.path.join(FIGURES_DIR, 'sdss_test_ba_fast_rotators.png')
    plt.savefig(out_path, dpi=150)
    print(f"Figure saved to {out_path}")

def main():
    cache_path = os.path.join(DATA_DIR, 'sdss_fast_rotators.csv')
    
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

    results, df_clean = analyze_fast_rotators(df)
    
    out_path = os.path.join(RESULTS_DIR, 'sdss_test_ba_results.json')
    
    # Helper to serialize Intervals
    def json_default(obj):
        if isinstance(obj, pd.Interval):
            return str(obj)
        raise TypeError
        
    with open(out_path, 'w') as f:
        json.dump(results, f, indent=2, default=json_default)
        
    create_figure_direct(df_clean, results)
    
    print("\nSUMMARY TEST BA:")
    print("TEP Prediction: Fraction remains high at high mass.")
    print("Standard Prediction: Fraction drops at high mass (Slow rotators dominate).")
    print(f"Correlation r(Lambda, Mass): {results['r_lambda']:.4f}")
    
    if results['r_lambda'] > -0.1:
         print("RESULT: CONSISTENT? (Weak dependence on mass?)")
    else:
         print("RESULT: CONTRADICTED (Massive galaxies are slow rotators, r < 0)")

if __name__ == "__main__":
    main()
