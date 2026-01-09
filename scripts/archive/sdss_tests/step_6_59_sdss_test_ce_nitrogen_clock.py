#!/usr/bin/env python3
"""
Step 6.59: SDSS Test CE - Nitrogen/Oxygen Clock (Secondary Enrichment)

Hypothesis:
Nitrogen is largely a "secondary" element produced from Oxygen and Carbon in the CNO cycle.
Its abundance scales as (O/H)^2. The N/O ratio tracks the integrated star formation history.
In galaxies with deep potentials (slow time), the enrichment cycles proceed slower.
At fixed O/H, we might expect different N/O ratios compared to galaxies where time runs "fast" (shallow potential).

Prediction:
N/O ratio at fixed O/H varies with sigma.

Data:
- emissionLinesPort: Flux_NII_6583, Flux_OII_3726, Flux_OIII_5006, Flux_Ha_6562, Flux_Hb_4861, sigma_stars

Method:
1. Select galaxies with reliable emission lines (S/N > 3).
2. Calculate R3, R23, N2 metallicity indicators.
3. Estimate O/H and N/O.
   Simplified Proxy: N2O2 = log(NII/OII) correlates with N/O.
   Or use N2 = log(NII/Ha) and O3N2 = log((OIII/Hb)/(NII/Ha)).
   Let's use N/O proxy directly if possible, or standard calibrations.
   A robust empirical proxy for N/O is N2S2 (NII/SII) or just NII/OII (if reddening corrected).
   Let's use NII/OII (flux ratio) but we need to worry about reddening.
   Ideally, use NII/Ha vs OIII/Hb (BPT) to ensure Star Forming.
   Then look at N/O at fixed O/H.
   
   Proxy for N/O: log([NII]6583 / [OII]3726) - but this is sensitive to reddening.
   Better: log([NII]6583 / [SII]6717,6731) - close in wavelength.
   But we checked columns and SII wasn't explicitly in the short list (Flux_SII_6716 exists usually).
   Let's use the available columns: NII, OII, OIII, Ha, Hb.
   We can correct for reddening using Ha/Hb (Balmer decrement).
   
   Steps:
   a. Compute E(B-V) from Ha/Hb (intrinsic 2.86).
   b. De-redden fluxes.
   c. Calculate 12+log(O/H) using O3N2 calibration.
   d. Calculate log(N/O) using N2O2 calibration? Or just look at NII/OII trend.
   e. Analyze Residuals of N/O vs O/H against Sigma.

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
    print(f"Querying SDSS for Test CE (Limit: {limit})...")
    
    # We need emission lines and sigma
    # Try minimal columns first to debug
    
    sql = f"""
    SELECT TOP {limit}
        specObjID,
        sigma_stars,
        Flux_NII_6583 as nii,
        Flux_OII_3726 as oii,
        Flux_OIII_5006 as oiii,
        Flux_Ha_6562 as ha,
        Flux_Hb_4861 as hb
        
    FROM emissionLinesPort
    WHERE 
        sigma_stars > 50
    """
    return query_sdss(sql)

def analyze_nitrogen(df):
    print("Analyzing Nitrogen Clock...")
    
    # Clean
    df = df.dropna().copy()
    
    # 1. BPT Selection (Star Forming only)
    # log([OIII]/Hb) vs log([NII]/Ha)
    # Kauffmann et al. 2003 line: 0.61 / (log(NII/Ha) - 0.05) + 1.3
    
    df['log_n2ha'] = np.log10(df['nii'] / df['ha'])
    df['log_o3hb'] = np.log10(df['oiii'] / df['hb'])
    
    # Keep SF galaxies (below Ka03 line)
    # Also NII/Ha < 0 to avoid LINERs/AGN
    mask_sf = (df['log_n2ha'] < 0) & (df['log_o3hb'] < (0.61 / (df['log_n2ha'] - 0.05) + 1.3))
    df_sf = df[mask_sf].copy()
    
    print(f"  Total: {len(df)}")
    print(f"  Star Forming: {len(df_sf)}")
    
    if len(df_sf) < 100:
        print("  Insufficient SF galaxies.")
        return None
        
    # 2. Reddening Correction (Cardelli et al 1989 / Calzetti)
    # Ha/Hb intrinsic = 2.86
    # E(B-V) = log10((Ha/Hb)/2.86) / 0.4 / (k_Ha - k_Hb)
    # k_Ha ~ 2.53, k_Hb ~ 3.61 (Calzetti) -> delta_k ~ -1.08? 
    # Simplified: A_V = E(B-V) * R_V.
    # Flux_corr = Flux_obs * 10^(0.4 * A_lambda)
    # Let's use a simple approximation for NII/OII relative reddening.
    # Lambda NII ~ 6583, OII ~ 3726. Large range. Sensitive.
    
    # Calculate E(B-V)
    # balmer_ratio = df_sf['ha'] / df_sf['hb']
    # ebv = 1.97 * np.log10(balmer_ratio / 2.86) # Approx coefficient
    # ebv = np.maximum(ebv, 0)
    
    # Instead of full correction, let's use O3N2 for Metallicity (immune to reddening largely)
    # O3N2 = log((OIII/Hb) / (NII/Ha))
    # 12 + log(O/H) = 8.73 - 0.32 * O3N2 (Pettini & Pagel 2004)
    
    df_sf['o3n2'] = df_sf['log_o3hb'] - df_sf['log_n2ha']
    df_sf['oh'] = 8.73 - 0.32 * df_sf['o3n2']
    
    # N/O Proxy
    # N2 = log(NII/Ha). This tracks N/O strongly but also O/H.
    # N2O2 = log(NII/OII). Tracks N/O better but needs reddening correction.
    # Let's use N2 (log NII/Ha) as the N-abundance proxy, and control for O/H (via O3N2).
    # Since O3N2 involves N2, they are correlated by definition.
    # But N2 is N/H * H/Ha? No, NII/Ha ~ N/H if T_e is fixed.
    # Let's use the Residuals of N2 vs O3N2?
    # Actually, let's look at NII/Ha at fixed O3N2.
    # If N/O is enhanced, NII/Ha should be higher for a given O3N2 (which is driven by excitation + Z).
    
    # Better: Use N/O calibration from N2O2 (with reddening correction).
    # k_NII ~ 2.5, k_OII ~ 4.7. diff ~ 2.2.
    # A_NII - A_OII = E(B-V) * (k_NII - k_OII) ~ -2.2 * E(B-V)
    # log(NII/OII)_int = log(NII/OII)_obs + 0.4 * (A_NII - A_OII)
    
    balmer = df_sf['ha'] / df_sf['hb']
    ebv = np.log10(balmer / 2.86) / (0.4 * (3.7 - 2.5)) # k approx diff ~ 1.2?
    # Let's just use raw flux ratios and hope scatter averages out, or simple correction.
    
    # Let's use the variable: log(NII/Ha) vs log(OIII/Hb).
    # TEP Prediction: At fixed excitation (OIII/Hb), NII/Ha (Nitrogen) varies with Sigma.
    
    # Fit NII/Ha vs OIII/Hb (The SF sequence)
    slope_bpt, intercept_bpt, _, _, _ = stats.linregress(df_sf['log_o3hb'], df_sf['log_n2ha'])
    print(f"  BPT Slope (NII/Ha vs OIII/Hb): {slope_bpt:.4f}")
    
    # Residuals of NII/Ha
    df_sf['n2_resid'] = df_sf['log_n2ha'] - (slope_bpt * df_sf['log_o3hb'] + intercept_bpt)
    
    # Correlate with Sigma
    df_sf['log_sigma'] = np.log10(df_sf['sigma_stars'])
    
    r_val, p_val = stats.pearsonr(df_sf['log_sigma'], df_sf['n2_resid'])
    print(f"  Correlation r(N2 Residual, logSigma): {r_val:.4f} (p={p_val:.2e})")
    
    slope_sig, intercept_sig, _, _, _ = stats.linregress(df_sf['log_sigma'], df_sf['n2_resid'])
    print(f"  Slope (Resid vs logSigma): {slope_sig:.4f}")
    
    # Plot
    fig, ax = plt.subplots(1, 2, figsize=(12, 5))
    
    # BPT
    ax[0].scatter(df_sf['log_n2ha'], df_sf['log_o3hb'], alpha=0.3, s=5, c='gray')
    ax[0].plot(df_sf['log_n2ha'], 0.61 / (df_sf['log_n2ha'] - 0.05) + 1.3, 'k--', label='Ka03')
    ax[0].set_xlabel('log([NII]/Ha)')
    ax[0].set_ylabel('log([OIII]/Hb)')
    ax[0].set_title('BPT Diagram (SF Selection)')
    ax[0].set_xlim(-2.5, 0.5)
    ax[0].set_ylim(-1.5, 1.0)
    
    # Resid vs Sigma
    # ax[1].scatter(df_sf['log_sigma'], df_sf['n2_resid'], alpha=0.3, s=5)
    # Binning
    bins = np.linspace(df_sf['log_sigma'].min(), df_sf['log_sigma'].max(), 10)
    df_sf['sig_bin'] = pd.cut(df_sf['log_sigma'], bins=bins)
    binned = df_sf.groupby('sig_bin')['n2_resid'].agg(['mean', 'sem', 'count'])
    binned['sig_center'] = [i.mid for i in binned.index]
    
    ax[1].errorbar(binned['sig_center'], binned['mean'], yerr=binned['sem'], fmt='o-', capsize=5)
    ax[1].set_xlabel('log Sigma [km/s]')
    ax[1].set_ylabel('Delta log(NII/Ha) (at fixed OIII/Hb)')
    ax[1].set_title(f'Test CE: Nitrogen Clock (r={r_val:.2f})')
    ax[1].grid(True, alpha=0.3)
    
    out_path = os.path.join(FIGURES_DIR, 'sdss_test_ce_nitrogen.png')
    plt.savefig(out_path, dpi=150)
    print(f"Figure saved to {out_path}")
    
    return {
        'r_val': r_val,
        'slope': slope_sig,
        'n_sample': int(len(df_sf))
    }

def main():
    cache_path = os.path.join(DATA_DIR, 'sdss_nitrogen.csv')
    
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

    results = analyze_nitrogen(df)
    
    out_path = os.path.join(RESULTS_DIR, 'sdss_test_ce_results.json')
    if results:
        # Fix interval serialization
        def default(o):
            if isinstance(o, pd.Interval): return str(o)
            raise TypeError
            
        with open(out_path, 'w') as f:
            json.dump(results, f, indent=2, default=default)
            
        print("\nSUMMARY TEST CE:")
        print("Prediction: N/O ratio (proxied by NII/Ha at fixed OIII/Hb) varies with Sigma.")
        print(f"Observed r: {results['r_val']:.4f}")
        
        if abs(results['r_val']) > 0.1:
             print("RESULT: SIGNAL (Correlation observed)")
        else:
             print("RESULT: NULL (No significant correlation)")

if __name__ == "__main__":
    main()
