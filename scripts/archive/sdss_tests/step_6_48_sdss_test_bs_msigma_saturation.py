#!/usr/bin/env python3
"""
Step 6.48: SDSS Test BS - M-sigma Saturation (The UCD Limit)

Hypothesis:
If there is a Universal Critical Density (UCD) or saturation scale in the TEP scalar field, 
the M_BH - sigma relation might show a break or curvature at the high-mass end (M_BH > 10^9 M_sun), 
where the horizon density of the BH drops below the critical density, or where the host potential saturates.

Prediction:
Slope of M_BH vs sigma changes (flattens or steepens) at high Mass/Sigma.
Standard GR predicts a single power law (roughly M ~ sigma^4-5).
TEP saturation implies a deviation from power law.

Data:
- spiders_quasar: BH mass estimates (logBHMS_mgII, logBHMA_hb)
- emissionLinesPort: sigmaStars

Method:
1. Fetch BH Mass and Sigma.
2. Clean data (ensure valid mass and sigma).
3. Analyze M_BH vs Sigma.
4. Fit single power law vs broken power law (or quadratic).
5. Check for residuals at high mass.
"""

import pandas as pd
import numpy as np
from scipy import stats
from scipy.optimize import curve_fit
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
    print(f"Querying SDSS for Test BS (Limit: {limit})...")
    
    # Use SPIDERS quasar catalog joined with galSpecInfo (MPA-JHU)
    # emissionLinesPort had poor overlap. galSpecInfo has v_disp.
    
    sql = f"""
    SELECT TOP {limit}
        q.SPECOBJID as specObjID,
        q.logBHMS_mgII, -- Magnesium based mass
        q.logBHMA_hb,   -- H-beta based mass
        g.v_disp as sigmaStars,
        g.v_disp_err as sigmaStarsErr
        
    FROM spiders_quasar q
    JOIN galSpecInfo g ON q.SPECOBJID = g.specObjID
    
    WHERE 
        g.v_disp > 50 AND g.v_disp < 400
        AND (
            (q.logBHMS_mgII > 0 AND q.logBHMS_mgII < 99) 
            OR 
            (q.logBHMA_hb > 0 AND q.logBHMA_hb < 99)
        )
    """
    return query_sdss(sql)

def power_law(x, a, b):
    return a * x + b

def analyze_msigma(df):
    print("Analyzing M-sigma Relation...")
    
    # 1. Consolidate BH Mass
    # Prefer MgII, fallback to Hb
    df['logBH'] = df['logBHMS_mgII']
    mask_nan = df['logBH'].isna() | (df['logBH'] == -99) | (df['logBH'] == 0)
    df.loc[mask_nan, 'logBH'] = df.loc[mask_nan, 'logBHMA_hb']
    
    # Clean
    df_clean = df.dropna(subset=['logBH', 'sigmaStars']).copy()
    df_clean = df_clean[(df_clean['logBH'] > 6) & (df_clean['sigmaStars'] > 0)]
    
    print(f"  Sample size: {len(df_clean)}")
    
    if len(df_clean) < 50:
        print("  Insufficient sample size.")
        return None
        
    df_clean['log_sigma'] = np.log10(df_clean['sigmaStars'])
    
    # 2. Fit Linear Relation (Power Law in log-log)
    slope, intercept, r_val, p_val, std_err = stats.linregress(df_clean['log_sigma'], df_clean['logBH'])
    print(f"  Linear Fit: log(M) = {slope:.2f} * log(sigma) + {intercept:.2f}")
    print(f"  Correlation: r={r_val:.3f}")
    
    # 3. Analyze Residuals vs Sigma
    df_clean['logBH_pred'] = slope * df_clean['log_sigma'] + intercept
    df_clean['resid'] = df_clean['logBH'] - df_clean['logBH_pred']
    
    # Check for curvature: Fit quadratic to residuals
    # resid = a * log_sigma^2 + b * log_sigma + c
    quad_coeffs = np.polyfit(df_clean['log_sigma'], df_clean['resid'], 2)
    curvature = quad_coeffs[0]
    print(f"  Residual Curvature (coeff of x^2): {curvature:.3f}")
    
    # 4. Binning
    df_clean['sigma_bin'] = pd.qcut(df_clean['sigmaStars'], 6)
    binned = df_clean.groupby('sigma_bin')['logBH'].agg(['mean', 'std', 'count'])
    binned['sem'] = binned['std'] / np.sqrt(binned['count'])
    binned['sigma_center'] = [i.mid for i in binned.index]
    binned['log_sigma_center'] = np.log10(binned['sigma_center'])
    
    print("\nM-sigma Binned:")
    print(binned[['mean', 'sem', 'count']])
    
    # Plot
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    # Panel 1: Relation
    axes[0].scatter(df_clean['log_sigma'], df_clean['logBH'], alpha=0.3, s=10)
    x_range = np.linspace(df_clean['log_sigma'].min(), df_clean['log_sigma'].max(), 100)
    axes[0].plot(x_range, slope*x_range + intercept, 'r-', label=f'Slope={slope:.2f}')
    
    # Plot binned
    axes[0].errorbar(binned['log_sigma_center'], binned['mean'], yerr=binned['sem'], fmt='ko', label='Binned')
    
    axes[0].set_xlabel('log(Velocity Dispersion)')
    axes[0].set_ylabel('log(BH Mass)')
    axes[0].set_title('M-sigma Relation')
    axes[0].legend()
    
    # Panel 2: Residuals
    axes[1].scatter(df_clean['log_sigma'], df_clean['resid'], alpha=0.3, s=10)
    axes[1].axhline(0, color='k', linestyle='--')
    
    # Plot quadratic fit to residuals
    axes[1].plot(x_range, np.polyval(quad_coeffs, x_range), 'g--', label=f'Curv={curvature:.2f}')
    
    axes[1].set_xlabel('log(Velocity Dispersion)')
    axes[1].set_ylabel('log(M) Residual')
    axes[1].set_title('Residuals from Linear Fit')
    axes[1].legend()
    
    plt.tight_layout()
    out_path = os.path.join(FIGURES_DIR, 'sdss_test_bs_msigma.png')
    plt.savefig(out_path, dpi=150)
    print(f"Figure saved to {out_path}")
    
    return {
        'slope': slope,
        'intercept': intercept,
        'r_val': r_val,
        'curvature': curvature,
        'n_sample': int(len(df_clean))
    }

def main():
    cache_path = os.path.join(DATA_DIR, 'sdss_msigma.csv')
    
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

    results = analyze_msigma(df)
    
    out_path = os.path.join(RESULTS_DIR, 'sdss_test_bs_results.json')
    if results:
        with open(out_path, 'w') as f:
            json.dump(results, f, indent=2)
            
        print("\nSUMMARY TEST BS:")
        print(f"Observed Slope: {results['slope']:.2f}")
        print(f"Residual Curvature: {results['curvature']:.3f}")
        
        # TEP Prediction: Curvature/Break at high sigma?
        # If curvature is significantly negative (flattening) or positive (steepening)
        if abs(results['curvature']) > 1.0:
             print("RESULT: DEVIATION (Significant curvature)")
        else:
             print("RESULT: NULL (Consistent with power law)")

if __name__ == "__main__":
    main()
