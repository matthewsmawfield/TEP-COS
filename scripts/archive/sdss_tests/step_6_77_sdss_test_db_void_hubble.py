#!/usr/bin/env python3
"""
Step 6.77: SDSS Test DB - Void Hubble Drift (Expansion Rate Variance)

Hypothesis:
In TEP, proper time flows faster in voids (shallow potential) than in filaments/clusters. 
The observed expansion rate H_local ~ 1/dt should appear higher in voids. 
Galaxies in underdense regions should have larger redshifts for a given true distance 
(Fundamental Plane distance) than those in overdensities.

Prediction:
Hubble Residual (v_peculiar) is positive (outward/fast) in Voids.
Delta v = c*z_obs - c*z_flow(D_FP) > 0 in Voids.

Data:
- ebossMCPM: mid_dens_1 (Density)
- galSpecInfo: z
- emissionLinesPort: sigma_stars
- PhotoObjAll: deVRad_r, deVAB_r, modelMag_r

Method:
1. Select Early-Type Galaxies (ETGs) suitable for Fundamental Plane (FP).
   - High sigma (> 100 km/s), red color (u-r > 2.2), concentration.
   - Use simplified selection: Sigma > 70.
2. Calibrate FP using the full sample (or Field sample).
   - log(R_e) = a * log(sigma) + b * SB_e + c
   - Residuals in Radius <-> Residuals in Distance <-> Residuals in Velocity.
3. Compare FP residuals in Voids (density < -0.5) vs Field (density ~ 0).
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

def download_data(limit=500):
    print(f"Querying SDSS for Test DB (Limit: {limit})...")
    
    sql = f"""
    SELECT TOP {limit}
        m.mid_dens_1 as density,
        g.z,
        e.sigma_stars,
        ph.deVRad_r, -- Effective Radius (arcsec)
        ph.modelMag_r,
        ph.deVAB_r -- Axis ratio
    FROM ebossMCPM m
    JOIN galSpecInfo g ON m.specObjID = g.specObjID
    JOIN emissionLinesPort e ON m.specObjID = e.specObjID
    JOIN PhotoObjAll ph ON m.specObjID = ph.objID
    WHERE e.sigma_stars > 70 
      AND g.z BETWEEN 0.05 AND 0.15
      AND abs(m.mid_dens_1) < 10
      AND ph.deVRad_r > 0
    """
    return query_sdss(sql)

def analyze_hubble_drift(df):
    print("Analyzing Void Hubble Drift...")
    
    if df is None or len(df) < 100:
        print("  Insufficient data.")
        return None
        
    df = df.dropna().copy()
    
    # Fundamental Plane Parameters
    # X = log(sigma)
    # Y = Surface Brightness
    # Z = log(R_physical)
    
    # Convert arcsec to physical kpc? 
    # For FP residuals, we often look at offset in log(R) at fixed sigma and SB.
    # Prediction: In Voids, z is "too high" for the distance.
    # If z is high, we infer a larger distance D(z).
    # But the object is actually closer (or standard size).
    # If we use z to calculate Physical Radius, the Void object will appear LARGER than it is.
    # So, Residual = log(R_measured_via_z) - log(R_predicted_by_FP).
    # Void objects should have Positive Residuals (appear larger due to z-boost).
    
    # 1. Calculate Surface Brightness (mu_e)
    # mu_e = m_r + 2.5 log(2 pi R_e^2)  (approx)
    # Use simpler: log(I_e) ~ -0.4 * SB
    # Standard FP: log(R_e) = a log(sigma) + b mu_e + c
    
    # We need angular diameter distance to convert deVRad_r (arcsec) to R_kpc
    # Approximation in local universe: D ~ c*z/H0
    H0 = 70.0
    c = 300000.0
    df['dist_mpc'] = (c * df['z']) / H0
    df['r_kpc'] = df['deVRad_r'] * (df['dist_mpc'] * 1000) / 206265.0
    
    df['log_r'] = np.log10(df['r_kpc'])
    df['log_sigma'] = np.log10(df['sigma_stars'])
    
    # Surface Brightness (mag/arcsec^2)
    # mu = m + 2.5 log(2 pi r^2) - 2.5 log(b/a) ?
    # Simpler: mean surface brightness within effective radius <mu>_e
    # <mu>_e = m + 2.5 log(2 pi r_eff^2)
    df['mu_e'] = df['modelMag_r'] + 2.5 * np.log10(2 * np.pi * df['deVRad_r']**2)
    
    # Fit FP to Field population (Density > 0)
    field = df[(df['density'] > 0) & (df['density'] < 2)] # Non-void, non-extreme cluster
    
    if len(field) < 50:
        print("  Not enough field galaxies for calibration.")
        return None
    
    # Fit: log_r = a * log_sigma + b * mu_e + c
    from sklearn.linear_model import LinearRegression
    X = field[['log_sigma', 'mu_e']]
    y = field['log_r']
    
    reg = LinearRegression().fit(X, y)
    print(f"  FP Fit: log(R) = {reg.coef_[0]:.3f} log(sigma) + {reg.coef_[1]:.3f} mu + {reg.intercept_:.3f}")
    
    # Apply to all
    df['log_r_pred'] = reg.predict(df[['log_sigma', 'mu_e']])
    df['fp_resid'] = df['log_r'] - df['log_r_pred']
    
    # Check Voids
    voids = df[df['density'] < -0.5]
    print(f"  N_Field: {len(field)}, N_Void: {len(voids)}")
    
    if len(voids) < 10:
        print("  Insufficient void galaxies.")
        return None
        
    void_mean_resid = voids['fp_resid'].mean()
    void_sem = voids['fp_resid'].sem()
    
    print(f"  Void Mean FP Residual: {void_mean_resid:.4f} +/- {void_sem:.4f}")
    
    # T-test
    t_stat, p_val = stats.ttest_ind(voids['fp_resid'], field['fp_resid'], equal_var=False)
    print(f"  T-test (Void vs Field): t={t_stat:.2f}, p={p_val:.2e}")
    
    # Plot
    fig, ax = plt.subplots(1, 2, figsize=(12, 5))
    
    # Histogram of residuals
    ax[0].hist(field['fp_resid'], bins=30, density=True, alpha=0.5, color='gray', label='Field')
    ax[0].hist(voids['fp_resid'], bins=30, density=True, alpha=0.5, color='blue', label='Void')
    ax[0].axvline(void_mean_resid, color='blue', linestyle='--', label=f'Void Mean={void_mean_resid:.3f}')
    ax[0].set_xlabel('FP Residual (log R_obs - log R_pred)')
    ax[0].set_ylabel('Density')
    ax[0].set_title('Fundamental Plane Residuals')
    ax[0].legend()
    
    # Residual vs Density
    # Bin by density
    df['dens_bin'] = pd.qcut(df['density'], q=10)
    binned = df.groupby('dens_bin').agg({'density':'mean', 'fp_resid':['mean', 'sem']}).reset_index()
    binned.columns = ['bin', 'density', 'resid_mean', 'resid_sem']
    
    ax[1].errorbar(binned['density'], binned['resid_mean'], yerr=binned['resid_sem'], fmt='o-', color='purple')
    ax[1].axhline(0, color='k', linestyle=':')
    ax[1].set_xlabel('Environment Density (delta)')
    ax[1].set_ylabel('Mean FP Residual')
    ax[1].set_title('Test DB: Hubble Drift Check')
    ax[1].grid(True, alpha=0.3)
    
    out_path = os.path.join(FIGURES_DIR, 'sdss_test_db_void_hubble.png')
    plt.savefig(out_path, dpi=150)
    print(f"Figure saved to {out_path}")
    
    return {
        'void_offset': void_mean_resid,
        'p_value': p_val,
        'n_void': int(len(voids))
    }

def main():
    cache_path = os.path.join(DATA_DIR, 'sdss_void_hubble.csv')
    
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

    results = analyze_hubble_drift(df)
    
    out_path = os.path.join(RESULTS_DIR, 'sdss_test_db_results.json')
    if results:
        with open(out_path, 'w') as f:
            json.dump(results, f, indent=2)
            
        print("\nSUMMARY TEST DB:")
        print(f"Void Residual Offset: {results['void_offset']:.4f} (Expected > 0 for Faster Void Time)")
        
        if results['p_value'] < 0.05 and results['void_offset'] > 0.01:
             print("RESULT: SIGNAL (Void galaxies appear larger/faster)")
        elif results['p_value'] < 0.05 and results['void_offset'] < -0.01:
             print("RESULT: CONTRADICTED (Void galaxies appear smaller/slower)")
        else:
             print("RESULT: NULL (No significant difference)")

if __name__ == "__main__":
    main()
