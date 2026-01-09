#!/usr/bin/env python3
"""
Step 6.55: SDSS Test CA - BLR Kinematics (Virial Factor Anomaly)

Hypothesis:
In AGN Broad Line Regions (BLR), cloud velocities v_FWHM balance gravity. v^2 ~ M_BH/R_BLR.
The radius R_BLR is set by photoionization physics (lag tau ~ L^0.5), a rate-based distance.
TEP modifies both gravity (metric) and rates (c_eff).
The relationship between v_FWHM and Luminosity L should show residuals correlated with the host potential depth (sigma_stars).

Prediction:
BLR FWHM (at fixed L) correlates with Host sigma.
Standard Virial Relation: v_FWHM \propto L^{-0.25} * M_BH^0.5
We are looking for deviations in v_FWHM vs L that correlate with sigma_host (independent of M_BH? No, M_BH correlates with sigma).
Better approach: 
The virial factor f in M_BH = f * R * v^2 / G might depend on inclination/geometry.
But if TEP is real, for a fixed M_BH (proxy sigma) and fixed L (proxy R), v should match standard gravity.
If TEP modifies gravity or R definition, v might deviate.
Specifically, if R_BLR (reverberation) is rate-based dist, and Gravity is metric based.
Let's check if FWHM residuals vs L correlate with sigma.

Data:
- spiders_quasar: fwhm1_hb, l_bol1, SPECOBJID
- galSpecInfo: v_disp (sigma), specObjID

Method:
1. Fetch Quasar properties (FWHM H-beta, L_bol).
2. Fetch Host properties (v_disp).
3. Join.
4. Fit FWHM vs L relation (Slope approx -0.25 expected).
5. Analyze residuals vs sigma.
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
    print(f"Querying SDSS for Test CA (Limit: {limit})...")
    
    # Join SPIDERS (Quasar) and galSpecInfo (Host)
    # Use Plate, MJD, FiberID for robust join
    # spiders: Plate, MJD, FiberID
    # galSpecInfo: plateid, mjd, fiberid
    
    sql = f"""
    SELECT TOP {limit}
        q.SPECOBJID,
        q.fwhm1_hb as fwhm_hb,
        q.l_bol1 as l_bol,
        g.v_disp as sigma,
        g.v_disp_err as sigma_err,
        g.z
        
    FROM spiders_quasar q
    JOIN galSpecInfo g ON q.Plate = g.plateid AND q.MJD = g.mjd AND q.FiberID = g.fiberid
    
    WHERE 
        q.fwhm1_hb > 0
        AND q.l_bol1 > 0
        AND g.v_disp > 50 
        AND g.v_disp < 400
        AND g.v_disp_err > 0
        AND g.v_disp_err < 50
    """
    return query_sdss(sql)

def analyze_blr(df):
    print("Analyzing BLR Kinematics...")
    
    if df is None or len(df) == 0:
        print("  No data to analyze.")
        return None
        
    print(f"  Columns: {df.columns.tolist()}")
    
    # Clean
    df = df.dropna().copy()
    
    # Handle potentially different column names (case sensitivity)
    # The SQL aliased as fwhm_hb, but sometimes it comes back as uppercase or original
    if 'fwhm_hb' not in df.columns:
        # Try to find match
        for col in df.columns:
            if col.lower() == 'fwhm_hb':
                df['fwhm_hb'] = df[col]
                break
            if col.lower() == 'fwhm1_hb': # Original name
                df['fwhm_hb'] = df[col]
                break
                
    if 'l_bol' not in df.columns:
        for col in df.columns:
            if col.lower() == 'l_bol':
                df['l_bol'] = df[col]
                break
            if col.lower() == 'l_bol1':
                df['l_bol'] = df[col]
                break
                
    # Log variables
    try:
        df['log_fwhm'] = np.log10(df['fwhm_hb'])
        df['log_lbol'] = np.log10(df['l_bol'])
        df['log_sigma'] = np.log10(df['sigma'])
    except KeyError as e:
        print(f"  KeyError: {e} - Columns missing after mapping attempt.")
        return None
    
    print(f"  Sample size: {len(df)}")
    
    if len(df) < 50:
        print("  Sample too small for meaningful analysis.")
        return {'n_sample': len(df), 'r_resid': 0}

    # 1. Fit FWHM vs Luminosity (The R-L relation inverse)
    # v \propto R^-0.5 \propto L^-0.25
    # log v = -0.25 log L + C
    
    slope, intercept, r_val, p_val, std_err = stats.linregress(df['log_lbol'], df['log_fwhm'])
    print(f"  FWHM-L Relation: Slope={slope:.4f} (Expected -0.25), r={r_val:.4f}")
    
    # Calculate Residuals (Excess Velocity at fixed L)
    df['fwhm_resid'] = df['log_fwhm'] - (slope * df['log_lbol'] + intercept)
    
    # 2. Correlate Residuals with Host Sigma
    # Does the BLR move faster/slower than expected in deep potentials?
    # Note: M_BH is correlated with sigma (M-sigma). 
    # And v^2 ~ M/R. So v ~ M^0.5. 
    # So v should correlate with sigma positively.
    # Residuals effectively control for L (R_BLR).
    # So we expect Residuals ~ 0.5 * log(M_BH).
    # Since M_BH ~ sigma^4, we expect Residuals ~ 0.5 * 4 * log(sigma) ~ 2 * log(sigma).
    # Standard Model prediction: Strong Positive Correlation.
    
    r_resid, p_resid = stats.pearsonr(df['log_sigma'], df['fwhm_resid'])
    print(f"  Residuals vs Sigma: r={r_resid:.4f} (p={p_resid:.2e})")
    
    slope_resid, intercept_resid, _, _, _ = stats.linregress(df['log_sigma'], df['fwhm_resid'])
    print(f"  Slope (Resid vs logSigma): {slope_resid:.4f}")
    
    # TEP Check: Is the slope consistent with standard M-sigma?
    # Standard: Slope ~ 2 (if M ~ sigma^4).
    # If TEP modifies gravity/clocks, we might see deviation?
    # Actually, this is a consistency check.
    
    # Plot
    fig, ax = plt.subplots(1, 2, figsize=(12, 5))
    
    # FWHM-L
    ax[0].scatter(df['log_lbol'], df['log_fwhm'], alpha=0.3, s=10)
    ax[0].plot(df['log_lbol'], slope * df['log_lbol'] + intercept, 'r-', label=f'Slope={slope:.2f}')
    ax[0].set_xlabel('log L_bol')
    ax[0].set_ylabel('log FWHM (H-beta)')
    ax[0].set_title('BLR Kinematics')
    ax[0].legend()
    ax[0].grid(True, alpha=0.3)
    
    # Resid-Sigma
    ax[1].scatter(df['log_sigma'], df['fwhm_resid'], alpha=0.3, s=10)
    ax[1].plot(df['log_sigma'], slope_resid * df['log_sigma'] + intercept_resid, 'r-', label=f'Slope={slope_resid:.2f}')
    ax[1].set_xlabel('log Sigma (Host)')
    ax[1].set_ylabel('FWHM Residuals (Fixed L)')
    ax[1].set_title(f'Residuals vs Potential (r={r_resid:.2f})')
    ax[1].legend()
    ax[1].grid(True, alpha=0.3)
    
    out_path = os.path.join(FIGURES_DIR, 'sdss_test_ca_blr.png')
    plt.savefig(out_path, dpi=150)
    print(f"Figure saved to {out_path}")
    
    return {
        'slope_vl': slope,
        'r_resid': r_resid,
        'slope_resid': slope_resid,
        'n_sample': int(len(df))
    }

def main():
    cache_path = os.path.join(DATA_DIR, 'sdss_blr_kinematics.csv')
    
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

    results = analyze_blr(df)
    
    out_path = os.path.join(RESULTS_DIR, 'sdss_test_ca_results.json')
    with open(out_path, 'w') as f:
        json.dump(results, f, indent=2)
        
    print("\nSUMMARY TEST CA:")
    print("Standard Model: v ~ M^0.5 ~ sigma^2. Residuals vs Sigma should have slope ~2.")
    print(f"Observed Slope (Resid vs logSigma): {results['slope_resid']:.4f}")
    
    if abs(results['slope_resid'] - 2.0) < 1.0:
         print("RESULT: CONSISTENT (Standard M-sigma scaling)")
    else:
         print("RESULT: ANOMALY (Slope deviates from M-sigma expectation)")

if __name__ == "__main__":
    main()
