#!/usr/bin/env python3
"""
Step 6.15: SDSS Test I - Post-Starburst Timing Anomaly

Hypothesis:
Post-Starburst (E+A) galaxies are "clocks" marking a specific recent event (starburst shutoff).
The Balmer absorption (HdeltaA) fades over ~1 Gyr.
Under TEP, time flows slower in deep potentials. 
Therefore, the fading of HdeltaA should proceed slower in high-sigma galaxies.
We should see a systematic "younger" appearance (stronger HdeltaA) or an excess of strong-HdeltaA systems 
in high-sigma environments compared to standard expectations?
Prediction: r(HdeltaA, sigma) > 0 for PSB galaxies.
(High sigma -> Slower fading -> Observed HdeltaA stays high longer -> Appears younger/stronger)

Data:
- sdss_spectral_indices.csv (HdeltaA, D4000, Sigma)
- sdss_bpt_data.csv (Halpha Flux for selection)

Selection:
- Strong HdeltaA (> 3 or 4)
- Weak Halpha emission (EW > -3 or low S/N)
"""

import pandas as pd
import numpy as np
from scipy import stats
import os
import json
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression

# Paths
DATA_DIR = os.path.join(os.path.dirname(__file__), '..', '..', 'data', 'sdss')
RESULTS_DIR = os.path.join(os.path.dirname(__file__), '..', '..', 'results', 'outputs')
FIGURES_DIR = os.path.join(os.path.dirname(__file__), '..', '..', 'results', 'figures')

def load_data():
    indices_path = os.path.join(DATA_DIR, 'sdss_spectral_indices.csv')
    bpt_path = os.path.join(DATA_DIR, 'sdss_bpt_data.csv')
    
    print(f"Loading {indices_path}...")
    df_idx = pd.read_csv(indices_path)
    
    print(f"Loading {bpt_path}...")
    df_bpt = pd.read_csv(bpt_path)
    
    # Merge
    print("Merging datasets...")
    # Ensure columns don't overlap except specObjID
    cols_to_use = df_bpt.columns.difference(df_idx.columns).tolist()
    cols_to_use.append('specObjID')
    
    # Check casing of specObjID
    if 'specobjid' in df_idx.columns:
        df_idx = df_idx.rename(columns={'specobjid': 'specObjID'})
    
    df = pd.merge(df_idx, df_bpt, on='specObjID', how='inner')
    print(f"Merged size: {len(df)}")
    return df

def analyze_psb(df):
    print("Analyzing Post-Starbursts...")
    
    # 1. Variables
    # HdeltaA
    # Usually in indices csv as 'hdelta_a'
    
    # Halpha
    # In bpt csv as 'Flux_Ha_6562'
    # We ideally want EW, but Flux is what we have in BPT file. 
    # Wait, the prompt plan mentioned EW. We only have Flux in bpt_data.csv.
    # However, sdss_spectral_indices.csv usually comes from galSpecIndx which has indices.
    # Is EW in there? 'hbeta' is Lick index.
    # Let's use Flux/Continuum proxy or just Flux cuts?
    # Or rely on 'bptclass'.
    # PSB usually means NO emission.
    # bptclass: -1=Inactive, 1=SF, 2=Comp, 3=AGN, 4=Seyfert, 5=Liner
    # So bptclass = -1 is good for "No Emission".
    # Or just low Halpha flux.
    
    # Selection Criteria for E+A (Goto et al. 2003):
    # HdeltaA > 4.0
    # Halpha EW > -3.0 (Emission is negative in some conventions, but usually Absorption positive)
    # SDSS emission lines: Positive flux = emission.
    # We want LOW emission.
    
    # Let's check column names in merged df
    # print(df.columns)
    
    # Define proxies
    df['HdeltaA'] = df['hdelta_a']
    df['log_sigma'] = np.log10(df['veldisp'])
    
    # Select PSBs
    # 1. Strong Absorption
    hdelta_cut = 3.0 # Slightly relaxed from 4.0 to get stats
    
    # 2. Weak Emission
    # Flux Halpha < threshold?
    # Or use bptclass == -1 (Passive) or maybe composite?
    # Strictly, E+A has no SF.
    # Let's look at Flux_Ha_6562.
    # Normalize by something?
    # Let's take the bottom 25% of Halpha emitters or just a raw flux cut if continuum is unavailable.
    # Or use the 'log_sfr' from indices file if available.
    # indices file has 'log_sfr'.
    # E+A should have low sSFR.
    
    # Let's use:
    # HdeltaA > 3.0
    # log_sfr < -1.0 (or some low value)
    # And valid sigma
    
    mask = (
        (df['HdeltaA'] > 3.0) & (df['HdeltaA'] < 15.0) & # Physical range
        (df['log_sfr'] < -0.5) & # Quiescent/Green valley
        (df['veldisp'] > 50) & (df['veldisp'] < 450) &
        (df['z_err'] < 0.001)
    )
    
    df_psb = df[mask].copy()
    print(f"  Selected {len(df_psb)} Post-Starburst candidates (Hdelta > 3, low SFR)")
    
    if len(df_psb) < 50:
        print("  WARNING: Low sample size. Relaxing cuts...")
        mask = (
            (df['HdeltaA'] > 2.0) & # Relaxed
            (df['log_sfr'] < 0.0) & # Relaxed
            (df['veldisp'] > 50) & (df['veldisp'] < 450)
        )
        df_psb = df[mask].copy()
        print(f"  Selected {len(df_psb)} relaxed candidates")
    
    # 3. Correlation
    # r(HdeltaA, sigma)
    
    r_simple, p_simple = stats.pearsonr(df_psb['log_sigma'], df_psb['HdeltaA'])
    
    # 4. Control for D4000 (Age within the PSB phase?)
    # If D4000 tracks the fading too, controlling for it removes the signal?
    # TEP prediction: HdeltaA is stronger at fixed "evolution"?
    # Actually, D4000 and HdeltaA are both evolutionary indicators.
    # If they are tight, controlling for one kills the other.
    # We check if the *sample* has higher HdeltaA at high sigma.
    
    print(f"  r(HdeltaA, sigma): {r_simple:.4f} (p={p_simple:.2e})")
    
    return {
        'r_simple': float(r_simple),
        'p_simple': float(p_simple),
        'n_sample': int(len(df_psb))
    }, df_psb

def create_figure(df, results):
    print("Generating figure...")
    fig, ax = plt.subplots(figsize=(8, 6))
    
    ax.scatter(df['log_sigma'], df['HdeltaA'], alpha=0.3, s=5, c='purple', label='PSB Candidates')
    
    # Fit
    m, b = np.polyfit(df['log_sigma'], df['HdeltaA'], 1)
    x = np.linspace(df['log_sigma'].min(), df['log_sigma'].max(), 100)
    ax.plot(x, m*x + b, 'k--', lw=2, label=f'Fit (r={results["r_simple"]:.3f})')
    
    ax.set_xlabel(r'$\log(\sigma)$')
    ax.set_ylabel(r'H$\delta_A$ Index (Absorption)')
    ax.set_title(f"Test I: Post-Starburst Timing\nStronger Absorption = 'Younger' Feature")
    ax.legend()
    
    plt.tight_layout()
    out_path = os.path.join(FIGURES_DIR, 'sdss_test_i_post_starburst.png')
    plt.savefig(out_path, dpi=150)
    print(f"Figure saved to {out_path}")

def main():
    df = load_data()
    results, df_psb = analyze_psb(df)
    
    out_path = os.path.join(RESULTS_DIR, 'sdss_test_i_results.json')
    with open(out_path, 'w') as f:
        json.dump(results, f, indent=2)
        
    create_figure(df_psb, results)
    
    print("\nSUMMARY TEST I:")
    print("TEP Prediction: r > 0 (High sigma -> Stronger HdeltaA / Younger Appearance)")
    print(f"Observed: r = {results['r_simple']:.4f}")
    
    if results['r_simple'] > 0.05 and results['p_simple'] < 0.05:
        print("RESULT: CONSISTENT with TEP.")
    elif results['r_simple'] < -0.05 and results['p_simple'] < 0.05:
        print("RESULT: CONTRADICTED (High sigma -> Weaker HdeltaA / Older).")
    else:
        print("RESULT: NULL/INCONCLUSIVE.")

if __name__ == "__main__":
    main()
