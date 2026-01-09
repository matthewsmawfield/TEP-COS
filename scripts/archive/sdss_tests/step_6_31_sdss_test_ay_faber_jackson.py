#!/usr/bin/env python3
"""
Step 6.31: SDSS Test AY - Faber-Jackson Curvature (Surface Brightness Dimming)

Hypothesis:
The Faber-Jackson relation (L ~ sigma^4) connects luminosity to potential depth.
In TEP, time dilation in deep potentials reduces the observed photon arrival rate (surface brightness dimming) beyond the standard redshift effect.
This should cause the FJ relation to curve or steepen at the high-sigma end (galaxies appear dimmer than expected for their sigma).

Prediction:
Luminosity falls below the standard FJ power law at high sigma.
Residual = logL_obs - logL_pred(sigma).
Residual < 0 at high sigma.

Data:
- emissionLinesPort: sigmaStars.
- SpecPhotoAll: modelMag_r, z (Redshift).
- stellarMassFSPSGranWideDust: logMass (Optional check).

Method:
1. Select Elliptical Galaxies (Class=GALAXY, High Sigma, maybe fracDeV > 0.8?).
   - SpecPhotoAll has 'type' but no fracDeV. PhotoObjAll has fracDeV.
   - We'll stick to High Sigma (> 150 km/s) usually implies ETGs, or just checking all galaxies.
   - FJ is for ETGs.
2. Compute Absolute Magnitude M_r.
   - D_L approx cz/H0 (H0=70).
   - K-correction? Approx 2.5 log(1+z).
3. Fit FJ: M_r = a * log(sigma) + b.
4. Compute Residuals.
5. Check for Curvature (Correlation of Residual vs Sigma? Or Quadratic term?).
   - If linear fit removes the main trend, a quadratic deviation will show up as curvature in residuals.
   - Or simpler: Residual vs Sigma slope at high sigma end?
   - The prediction says "steepen", so L is lower (M_r is fainter/more positive) at high sigma.
   - Standard FJ: L ~ sigma^4 -> M ~ -10 log sigma.
   - TEP: L_obs = L_int * A^2? A < 1.
   - So L_obs is lower. M_obs is higher (fainter).
   - Residual (M_obs - M_pred) should be POSITIVE (Fainter) at high sigma?
   - Wait. If L falls below, M (magnitude) is above (fainter).
   - Let's check sign convention carefully.

   Fit: M_pred = -10 log(sigma) + C.
   Residual = M_obs - M_pred.
   If M_obs is fainter (larger) than M_pred at high sigma -> Residual > 0.
   So we expect Positive Residual correlation or upward curvature.

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
    print(f"Querying SDSS for Test AY (Limit: {limit})...")
    
    # Select galaxies.
    # Join SpecPhotoAll for mags and z.
    # Join emissionLinesPort for sigma.
    
    sql = f"""
    SELECT TOP {limit}
        sp.specObjID,
        sp.z,
        sp.modelMag_r,
        e.sigmaStars as sigma
        
    FROM SpecPhotoAll sp
    JOIN emissionLinesPort e ON sp.specObjID = e.specObjID
    
    WHERE 
        sp.class = 'GALAXY'
        AND sp.z > 0.01 AND sp.z < 0.2
        AND e.sigmaStars > 70 AND e.sigmaStars < 450
        AND sp.modelMag_r > 0 AND sp.modelMag_r < 25
    """
    return query_sdss(sql)

def analyze_faber_jackson(df):
    print("Analyzing Faber-Jackson Curvature...")
    
    # 1. Clean
    df_clean = df.dropna().copy()
    
    # 2. Compute Absolute Magnitude
    # Distance Modulus mu = 5 log10(D_L) + 25
    # D_L = (c z / H0) * (1 + z/2) approx for low z?
    # Simple Hubble Law: D = v/H0 = cz/H0. D_L = D(1+z).
    # D_L = cz(1+z)/H0.
    # c = 3e5 km/s. H0 = 70 km/s/Mpc.
    # D_L in Mpc.
    
    H0 = 70.0
    c = 300000.0
    df_clean['D_L_Mpc'] = (c * df_clean['z'] * (1 + df_clean['z'])) / H0
    df_clean['dist_mod'] = 5 * np.log10(df_clean['D_L_Mpc'] * 1e6 / 10.0)
    
    # K-correction approximation for r-band (Chilingarian et al 2010 approx or just 2.5 log(1+z))
    # Simple K-corr:
    df_clean['k_corr'] = 2.5 * np.log10(1 + df_clean['z']) # Very rough
    
    df_clean['abs_mag_r'] = df_clean['modelMag_r'] - df_clean['dist_mod'] - df_clean['k_corr']
    
    # 3. Fit Faber-Jackson
    # M_r = a * log10(sigma) + b
    # Expect a ~ -10 (L ~ sigma^4 => -2.5 log L ~ -10 log sigma)
    
    df_clean['log_sigma'] = np.log10(df_clean['sigma'])
    
    slope, intercept, r_val, p_val, std_err = stats.linregress(df_clean['log_sigma'], df_clean['abs_mag_r'])
    print(f"FJ Fit: slope={slope:.2f} (Expected ~ -10), intercept={intercept:.2f}, r={r_val:.2f}")
    
    # 4. Residuals
    # Residual = Obs - Pred
    df_clean['fj_resid'] = df_clean['abs_mag_r'] - (slope * df_clean['log_sigma'] + intercept)
    
    # 5. Check for Curvature / Correlation with Sigma
    # If standard FJ holds, residuals should be flat vs sigma.
    # If TEP holds (Dimming at high sigma), M_obs > M_pred (Fainter).
    # Residual > 0 at high sigma.
    # Positive correlation?
    
    r_resid, p_resid = stats.pearsonr(df_clean['log_sigma'], df_clean['fj_resid'])
    
    print(f"N = {len(df_clean)}")
    print(f"Correlation r(Residual, sigma): {r_resid:.4f} (p={p_resid:.2e})")
    
    # 6. Binning
    df_clean['sigma_bin'] = pd.qcut(df_clean['log_sigma'], 8)
    binned = df_clean.groupby('sigma_bin')['fj_resid'].mean()
    print("\nMean FJ Residual by Sigma Bin:")
    print(binned)
    
    return {
        'fj_slope': float(slope),
        'r_resid': float(r_resid),
        'p_resid': float(p_resid),
        'mean_resid': float(df_clean['fj_resid'].mean()),
        'binned_means': binned.tolist(),
        'bin_centers': [mid.mid for mid in binned.index],
        'n_sample': int(len(df_clean))
    }, df_clean

def create_figure(df, results):
    print("Generating figure...")
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # Scatter
    ax.scatter(df['log_sigma'], df['fj_resid'], alpha=0.1, s=2, c='k', label='Galaxies')
    
    # Binned
    ax.plot(results['bin_centers'], results['binned_means'], 'r-o', lw=2, label='Mean Residual')
    
    ax.set_xlabel(r'$\log(\sigma)$')
    ax.set_ylabel(r'FJ Residual ($\Delta M_r$) [Fainter > 0]')
    ax.set_title(f"Test AY: FJ Residuals (r={results['r_resid']:.3f})")
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.axhline(0, color='b', linestyle='--')
    
    plt.tight_layout()
    out_path = os.path.join(FIGURES_DIR, 'sdss_test_ay_faber_jackson.png')
    plt.savefig(out_path, dpi=150)
    print(f"Figure saved to {out_path}")

def main():
    cache_path = os.path.join(DATA_DIR, 'sdss_faber_jackson.csv')
    
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

    results, df_clean = analyze_faber_jackson(df)
    
    out_path = os.path.join(RESULTS_DIR, 'sdss_test_ay_results.json')
    with open(out_path, 'w') as f:
        results_json = results.copy()
        del results_json['bin_centers']
        json.dump(results_json, f, indent=2)
        
    create_figure(df_clean, results)
    
    print("\nSUMMARY TEST AY:")
    print("TEP Prediction: Residual > 0 at high sigma (Fainter). r > 0.")
    print(f"Observed r: {results['r_resid']:.4f}")
    
    if results['r_resid'] > 0.05:
        print("RESULT: CONSISTENT (Galaxies appear fainter in deep potentials)")
    elif results['r_resid'] < -0.05:
        print("RESULT: CONTRADICTED (Galaxies appear brighter in deep potentials)")
    else:
        print("RESULT: NULL/FLAT")

if __name__ == "__main__":
    main()
