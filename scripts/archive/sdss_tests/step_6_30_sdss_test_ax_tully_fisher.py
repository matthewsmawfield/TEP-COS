#!/usr/bin/env python3
"""
Step 6.30: SDSS Test AX - Tully-Fisher Residuals (Metric Potential)

Hypothesis:
The Tully-Fisher relation connects Luminosity (Mass) to Rotation Velocity (Potential Gradient).
TEP modifies the effective potential. 
Residuals from the standard TF relation should correlate with the depth of the potential (approximated by central velocity dispersion or surface brightness).
Ideally, we compare rotation velocity (potential gradient) with luminosity (mass).
Standard: L ~ v^alpha.
TEP: Effective potential is deeper? Or time dilation dims surface brightness?
If time dilation makes L appear lower for a given Mass (dimming), then for a given v_rot, L is lower.
So galaxy falls below the TF line.
Or if v_rot is boosted?
Prediction: TF Residual (Delta logM or Delta logL) correlates with sigma_central.
Define Residual = logM_obs - logM_pred(v_rot).
If L is dimmed (M_obs is lower), Residual < 0 at high sigma.

Data:
- mangaDAPall:
    - ha_gvel_hi_clip (v_rot proxy)
    - nsa_elpetro_mass (Stellar Mass)
    - stellar_sigma_1re (Central sigma)
    - nsa_sersic_ba (Inclination)

Method:
1. Select rotating disks (v_rot > 50, b/a < 0.8).
2. Fit TF relation: logM = a * log(v_rot) + b.
3. Compute residuals: Delta = logM_obs - (a * log(v_rot) + b).
4. Correlate Delta with sigma.
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

def download_data(limit=1000):
    print(f"Querying SDSS for Test AX (Limit: {limit})...")
    
    sql = f"""
    SELECT TOP {limit}
        d.mangaid,
        d.ha_gvel_hi_clip as v_rot_raw,
        t.nsa_elpetro_mass as logmass,
        d.stellar_sigma_1re as sigma,
        d.nsa_sersic_ba as axis_ratio
        
    FROM mangaDAPall d
    JOIN mangaTarget t ON d.mangaid = t.mangaid
    
    WHERE 
        d.drp3qual = 0
        AND t.nsa_elpetro_mass > 0
    """
    return query_sdss(sql)

def analyze_tully_fisher(df):
    print("Analyzing Tully-Fisher Residuals...")
    
    # 1. Clean
    df_clean = df.dropna().copy()
    
    # 2. Inclination Correction
    # sin(i) ~ sqrt(1 - (b/a)^2) for thin disk?
    # Or just standard approximation if q0=0.2
    q = df_clean['axis_ratio']
    q0 = 0.2
    # If q < q0, set to q0
    q = np.maximum(q, q0)
    
    cos2i = (q**2 - q0**2) / (1 - q0**2)
    cos2i = np.clip(cos2i, 0, 1)
    sin_i = np.sqrt(1 - cos2i)
    
    # Avoid face-on
    df_clean['sin_i'] = sin_i
    df_clean = df_clean[df_clean['sin_i'] > 0.5].copy() # i > 30 deg
    
    # Correct velocity: v_rot = v_obs / sin(i)
    df_clean['v_rot'] = np.abs(df_clean['v_rot_raw']) / df_clean['sin_i']
    
    # Filter valid rotation
    df_clean = df_clean[df_clean['v_rot'] > 50].copy()
    
    # 3. Fit TF Relation
    # logM = a * logV + b
    df_clean['log_v'] = np.log10(df_clean['v_rot'])
    
    slope, intercept, r_val, p_val, std_err = stats.linregress(df_clean['log_v'], df_clean['logmass'])
    print(f"TF Fit: slope={slope:.2f}, intercept={intercept:.2f}, r={r_val:.2f}")
    
    # 4. Residuals
    df_clean['tf_resid'] = df_clean['logmass'] - (slope * df_clean['log_v'] + intercept)
    
    # 5. Correlation with Sigma
    # sigma is partial proxy for mass, so we must control for mass?
    # Or is sigma the "potential depth" variable we test against?
    # Residual is "Mass at fixed Velocity".
    # If we correlate with sigma, we might just recover the fact that sigma and v_rot are correlated?
    # But we subtracted v_rot dependence.
    # So we are asking: At fixed v_rot, does having higher sigma (bulge?) change Mass?
    
    # Filter valid sigma
    df_clean = df_clean[df_clean['sigma'] > 0].copy()
    
    df_clean['log_sigma'] = np.log10(df_clean['sigma'])
    
    # Drop any remaining NaNs
    df_clean = df_clean.dropna(subset=['log_sigma', 'tf_resid'])
    
    r_resid, p_resid = stats.pearsonr(df_clean['log_sigma'], df_clean['tf_resid'])
    
    print(f"N = {len(df_clean)}")
    print(f"Correlation r(Residual, sigma): {r_resid:.4f} (p={p_resid:.2e})")
    
    # 6. Binning
    df_clean['sigma_bin'] = pd.qcut(df_clean['log_sigma'], 8)
    binned = df_clean.groupby('sigma_bin')['tf_resid'].mean()
    print("\nMean TF Residual by Sigma Bin:")
    print(binned)
    
    return {
        'tf_slope': float(slope),
        'r_resid': float(r_resid),
        'p_resid': float(p_resid),
        'mean_resid': float(df_clean['tf_resid'].mean()),
        'binned_means': binned.tolist(),
        'bin_centers': [mid.mid for mid in binned.index],
        'n_sample': int(len(df_clean))
    }, df_clean

def create_figure(df, results):
    print("Generating figure...")
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # Scatter
    ax.scatter(df['log_sigma'], df['tf_resid'], alpha=0.1, s=2, c='k', label='Galaxies')
    
    # Binned
    ax.plot(results['bin_centers'], results['binned_means'], 'r-o', lw=2, label='Mean Residual')
    
    ax.set_xlabel(r'$\log(\sigma)$')
    ax.set_ylabel(r'TF Residual ($\Delta \log M_{*}$)')
    ax.set_title(f"Test AX: TF Residual vs Potential (r={results['r_resid']:.3f})")
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.axhline(0, color='b', linestyle='--')
    
    plt.tight_layout()
    out_path = os.path.join(FIGURES_DIR, 'sdss_test_ax_tully_fisher.png')
    plt.savefig(out_path, dpi=150)
    print(f"Figure saved to {out_path}")

def main():
    cache_path = os.path.join(DATA_DIR, 'sdss_tully_fisher.csv')
    
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

    results, df_clean = analyze_tully_fisher(df)
    
    out_path = os.path.join(RESULTS_DIR, 'sdss_test_ax_results.json')
    with open(out_path, 'w') as f:
        results_json = results.copy()
        del results_json['bin_centers']
        json.dump(results_json, f, indent=2)
        
    create_figure(df_clean, results)
    
    print("\nSUMMARY TEST AX:")
    print("TEP Prediction: Residual < 0 at high sigma (Dimming/Lower Mass inferred). r < 0.")
    print(f"Observed r: {results['r_resid']:.4f}")
    
    if results['r_resid'] < -0.1:
        print("RESULT: CONSISTENT (Mass/Luminosity suppressed in deep potentials)")
    elif results['r_resid'] > 0.1:
        print("RESULT: CONTRADICTED (Mass excess in deep potentials)")
    else:
        print("RESULT: NULL/FLAT")

if __name__ == "__main__":
    main()
