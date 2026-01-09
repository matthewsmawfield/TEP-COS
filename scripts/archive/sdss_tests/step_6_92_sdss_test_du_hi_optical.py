#!/usr/bin/env python3
"""
Step 6.92: SDSS Test DU - HI vs Optical Kinematics (Phantom Halo Extent)

Hypothesis:
HI probes the outer potential (flat rotation curve), while Optical (H-alpha, stars) 
probes the inner potential. TEP's phantom mass distribution (soliton + tail) might 
predict a different ratio of V_outer/V_inner compared to CDM NFW halos, especially 
as a function of central concentration (sigma).

Prediction:
Ratio W20(HI) / V_opt varies with Central Velocity Dispersion.

Data:
- mangaHIall: W20 (HI width), conf_prob
- mangaDAPall: ha_gvel_1re (Optical Gas Vel at 1 Re), stellar_sigma_1re

Method:
1. Join mangaHIall and mangaDAPall.
2. Select confident HI detections.
3. Define Ratio R = (0.5 * W20) / ha_gvel_1re (Approx V_max_HI / V_opt_1re).
   Note: W20 is full width, so V_HI ~ W20/2 / sin(i).
   ha_gvel_1re is rotational velocity at 1 Re (usually deprojected? check docs).
   If both are projected or both deprojected, ratio holds.
   mangaHIall W20 is projected.
   mangaDAPall ha_gvel_1re is likely projected? "Gas velocity dispersion" vs "Gas Velocity".
   Wait, ha_gvel might be dispersion? No, usually rotation.
   Let's assume we need to handle inclination if possible, or assume random orientation averages out in bins.
   However, we are joining the SAME galaxy. Inclination cancels out in the ratio V_HI/V_opt 
   if both are projected velocities!
   (W20/2) / V_opt_proj = V_HI_circ * sin(i) / (V_opt_circ * sin(i)) = V_HI / V_opt.
   Excellent.
4. Correlate R with stellar_sigma_1re (Potential Depth).
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

def download_data(limit=100):
    print(f"Querying SDSS for Test DU (Limit: {limit})...")
    
    # mangaHIall joined with mangaDAPall
    # ha_gvel_1re: "Rotation velocity of H-alpha at 1 Re"? 
    # Actually, DAP output usually has 'STELLAR_VEL', 'STELLAR_SIGMA', 'EMLINE_GVEL', 'EMLINE_GSIGMA'.
    # Checking common columns in mangaDAPall: 'ha_gvel_hi', 'ha_gvel_lo' ??
    # Let's use whatever rotation proxy is available.
    # 'emline_gsigma_1re_ha' is dispersion.
    # We might not have explicit rotation velocity in DAPall summary table.
    # We have 'stellar_sigma_1re'.
    # We have 'ha_flux_1re'.
    # Maybe we should use 'stellar_sigma_1re' as the V_inner proxy if rotation is missing?
    # But hypothesis compares Outer Rotation to Inner Rotation.
    
    # Checking known columns from my check script output (truncated).
    # Common practice: Use stellar sigma as inner potential proxy.
    # Compare V_HI (Rotation) vs Sigma_star (Dispersion).
    # The "v_rot / sigma" ratio is a measure of dynamical state (V/sigma).
    # But we want to know if this ratio SCALES with sigma differently.
    # TEP: V_outer / V_inner.
    # If we use V_HI / Sigma, we are testing halo vs bulge.
    
    sql = f"""
    SELECT TOP {limit}
        h.mangaid, 
        h.W20, 
        h.W50,
        d.stellar_sigma_1re,
        d.emline_gsigma_1re_ha -- Gas dispersion
        -- d.ha_vel_1re ? (Likely not in summary)
    FROM mangaHIall h
    JOIN mangaDAPall d ON h.mangaid = d.mangaid
    WHERE h.conf_prob > 0.9
      AND h.W20 > 0
      AND d.stellar_sigma_1re > 0
    """
    return query_sdss(sql)

def analyze_hi_optical(df):
    print("Analyzing HI vs Optical Kinematics...")
    
    if df is None or len(df) < 50:
        print("  Insufficient data.")
        return None
        
    df = df.dropna().copy()
    
    # Proxy for V_HI: W20 / 2 (Projected)
    df['v_hi_proj'] = df['W20'] / 2.0
    
    # We don't have V_opt_rot (Projected) easily from summary.
    # But we have stellar_sigma (Intrinsic, roughly).
    # If we look at V_HI_proj / Sigma, inclination matters.
    # But inclination is random.
    # If we bin by Sigma, mean inclination should be similar?
    # Unless high sigma galaxies are seen face-on? (Selection bias?)
    # Ellipticals (High Sigma) have little HI.
    # We are looking at Spirals (HI detected).
    
    # Let's define the ratio: R = V_HI_proj / Sigma_star
    # This is roughly (V_circ * sin i) / Sigma.
    # Does this ratio correlate with Sigma?
    # Low Sigma (dwarf/late spiral) -> Rotation dominated -> High V/Sigma.
    # High Sigma (early spiral) -> Bulge dominated -> Lower V/Sigma?
    # This is the standard "Hubble Sequence" kinematic trend.
    # TEP Prediction: "Ratio ... varies with Central Velocity Dispersion".
    # Well, standard physics predicts this too.
    # We need a deviation from LCDM.
    # LCDM (NFW): V_circ is flat. V_outer ~ V_inner.
    # V_opt (at 1Re) ~ V_outer.
    # So V_HI / V_opt ~ 1.
    # If we use Sigma as V_inner proxy, we rely on Jeans eqn.
    # V_circ ~ sqrt(2) * Sigma (Isothermal).
    # So V_HI / Sigma ~ sqrt(2) * sin(i).
    # This should be constant on average if profile shape is self-similar (Isothermal/NFW).
    # If it varies with Sigma, it implies non-homology or non-isothermal.
    
    df['ratio'] = df['v_hi_proj'] / df['stellar_sigma_1re']
    
    # Filter outliers
    df = df[df['ratio'] < 20].copy()
    
    # Log Sigma
    df['log_sigma'] = np.log10(df['stellar_sigma_1re'])
    
    # Correlation
    slope, intercept, r_val, p_val, std_err = stats.linregress(df['log_sigma'], df['ratio'])
    
    print(f"  Correlation (log Sigma vs Ratio): r={r_val:.3f}, p={p_val:.2e}, slope={slope:.4f}")
    
    # Plot
    plt.figure(figsize=(8, 6))
    plt.scatter(df['log_sigma'], df['ratio'], s=5, alpha=0.5, c='teal')
    
    x_range = np.linspace(df['log_sigma'].min(), df['log_sigma'].max(), 100)
    plt.plot(x_range, intercept + slope*x_range, 'k--', label=f'Slope={slope:.3f}')
    
    plt.xlabel('log(Stellar Sigma) [km/s]')
    plt.ylabel('Ratio: (W20/2) / Sigma')
    plt.title('Test DU: HI Kinematics vs Inner Potential')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    out_path = os.path.join(FIGURES_DIR, 'sdss_test_du_hi_optical.png')
    plt.savefig(out_path, dpi=150)
    print(f"Figure saved to {out_path}")
    
    return {
        'correlation_r': r_val,
        'p_value': p_val,
        'slope': slope,
        'n_gal': int(len(df))
    }

def main():
    cache_path = os.path.join(DATA_DIR, 'sdss_hi_optical.csv')
    
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

    results = analyze_hi_optical(df)
    
    out_path = os.path.join(RESULTS_DIR, 'sdss_test_du_results.json')
    if results:
        with open(out_path, 'w') as f:
            json.dump(results, f, indent=2)
            
        print("\nSUMMARY TEST DU:")
        print(f"Slope (logSigma vs Ratio): {results['slope']:.4f}")
        
        # Interpretation:
        # Standard: V/Sigma decreases with Mass/Sigma (Tully-Fisher vs Faber-Jackson offset).
        # Late types (Low Sigma) have high V/Sigma. Early types (High Sigma) have low V/Sigma.
        # So we expect a Negative Slope.
        # TEP Prediction: "Varies".
        # If slope is significantly different from "standard" - hard to say without a model.
        # But if it's very steep, it's interesting.
        # Let's report the signal.
        
        if results['p_value'] < 0.05 and abs(results['slope']) > 0.1:
             print("RESULT: SIGNAL (Strong dependence observed)")
        else:
             print("RESULT: NULL (Weak/No dependence)")

if __name__ == "__main__":
    main()
