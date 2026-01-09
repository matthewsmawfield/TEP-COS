#!/usr/bin/env python3
"""
Step 6.13: SDSS Test K - Size-Age Sign Test

Hypothesis:
Standard hierarchical formation: Compact galaxies (at fixed mass) formed earlier at high redshift.
Prediction: Compactness correlates POSITIVELY with Age (D4000).

TEP Hypothesis:
Compact galaxies have deeper gravitational potentials.
Time flows slower in deep potentials.
Prediction: Compactness correlates NEGATIVELY with apparent Age (D4000) at fixed formation timescale ([Mg/Fe]).
(Or at least, the positive correlation is suppressed).

Variables:
- Compactness proxy: Sigma_M = log10(Mass / R_e^2)
- Age proxy: D4000
- Formation timescale proxy: [Mg/Fe]

Method:
1. Merge galaxy size data with spectral indices.
2. Calculate Surface Mass Density (Compactness).
3. Regress D4000 against [Mg/Fe] and Compactness.
4. Check sign of Compactness coefficient.
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
    galaxies_path = os.path.join(DATA_DIR, 'sdss_galaxies.csv')
    
    if not os.path.exists(indices_path) or not os.path.exists(galaxies_path):
        raise FileNotFoundError("Required data files not found.")
        
    print(f"Loading indices from {indices_path}...")
    df_idx = pd.read_csv(indices_path)
    
    print(f"Loading galaxies from {galaxies_path}...")
    df_gal = pd.read_csv(galaxies_path)
    
    # Merge
    print("Merging datasets...")
    # Ensure specobjid is string or consistent type if needed, but usually pandas handles int64 fine
    df = pd.merge(df_idx, df_gal[['specobjid', 'petroR50_r', 'petroR90_r', 'expAB_r']], on='specobjid', how='inner')
    print(f"Merged dataset size: {len(df)}")
    return df

def analyze_size_age(df):
    print("Analyzing Size-Age relation...")
    
    # 1. Variables
    # Mass
    df['log_Mass'] = df['log_mass']
    
    # Size (Petrosian half-light radius in arcsec -> need kpc ideally, but at fixed z it's fine. 
    # Actually we should use physical size if we want density.
    # We have redshift. 
    # Approx conversion: 1 arcsec ~ 8 kpc at z=1? No.
    # At z=0.1, 1 arcsec ~ 1.8 kpc.
    # Let's use astropy for robust conversion or just use relative if z range is small.
    # The control tests script used astropy. Let's do a simple approx or use redshift.
    from astropy.cosmology import FlatLambdaCDM
    cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
    kpc_per_arcsec = cosmo.kpc_proper_per_arcmin(df['redshift'].values).value / 60.0
    df['R_e_kpc'] = df['petroR50_r'] * kpc_per_arcsec
    
    # Surface Mass Density (Compactness)
    # Sigma ~ M / R^2
    # log Sigma = log M - 2 log R
    df['log_R_e'] = np.log10(np.maximum(0.1, df['R_e_kpc']))
    df['Compactness'] = df['log_Mass'] - 2 * df['log_R_e']
    
    # [Mg/Fe] proxy
    df['Fe_avg'] = (df['fe5270'] + df['fe5335']) / 2
    df['MgFe_ratio'] = df['mgb'] / df['Fe_avg']
    df['log_MgFe'] = np.log10(np.maximum(0.1, df['MgFe_ratio']))
    
    # Age
    df['log_Age'] = np.log10(df['d4000'])
    
    # Sigma (Velocity Dispersion) for reference
    df['log_sigma'] = np.log10(df['veldisp'])
    
    # 2. Quality Cuts
    # Early-type like selection to ensure D4000 is meaningful age indicator (not just current SFR)
    mask = (
        (df['d4000'] > 1.3) & (df['d4000'] < 2.5) &
        (df['mgb'] > 0) & (df['fe5270'] > 0) & 
        (df['redshift'] > 0.02) & (df['redshift'] < 0.2) &
        (df['petroR50_r'] > 0.5) & # Resolved
        (df['expAB_r'] > 0.5) # Not edge-on disks ideally, but we want spheroidals mostly
    )
    df_clean = df[mask].copy()
    print(f"  Selected {len(df_clean)} galaxies")
    
    # 3. Regression
    # We want to explain Age (D4000).
    # Predictors: [Mg/Fe] (Formation Timescale), Compactness
    # Also control for Mass (standard downsizing)
    
    X = df_clean[['log_MgFe', 'Compactness', 'log_Mass']].values
    y = df_clean['log_Age'].values
    
    reg = LinearRegression().fit(X, y)
    
    # Check coefficient of Compactness
    coeffs = {
        'MgFe_coeff': float(reg.coef_[0]),
        'Compactness_coeff': float(reg.coef_[1]),
        'Mass_coeff': float(reg.coef_[2]),
        'Intercept': float(reg.intercept_)
    }
    
    print("\nRegression Results (Target: log D4000):")
    print(f"  [Mg/Fe] coeff:    {coeffs['MgFe_coeff']:.4f} (Expected > 0: High Mg/Fe -> Old)")
    print(f"  Compactness coeff:{coeffs['Compactness_coeff']:.4f} (Standard > 0, TEP < 0)")
    print(f"  Mass coeff:       {coeffs['Mass_coeff']:.4f} (Downsizing > 0)")
    
    # 4. Partial Correlation for verification
    # Residualize Age against Mass and MgFe
    reg_control = LinearRegression().fit(df_clean[['log_MgFe', 'log_Mass']].values, y)
    age_resid = y - reg_control.predict(df_clean[['log_MgFe', 'log_Mass']].values)
    
    # Residualize Compactness against Mass and MgFe
    reg_compact = LinearRegression().fit(df_clean[['log_MgFe', 'log_Mass']].values, df_clean['Compactness'].values)
    compact_resid = df_clean['Compactness'].values - reg_compact.predict(df_clean[['log_MgFe', 'log_Mass']].values)
    
    r_part, p_part = stats.pearsonr(compact_resid, age_resid)
    print(f"\nPartial Correlation r(Age, Compactness | Mass, [Mg/Fe]): {r_part:.4f} (p={p_part:.2e})")
    
    return {
        'coeffs': coeffs,
        'r_partial': float(r_part),
        'p_partial': float(p_part),
        'n_sample': int(len(df_clean))
    }, df_clean

def create_figure(df, results):
    print("Generating figure...")
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # 1. Age vs Compactness (Raw)
    ax = axes[0]
    ax.scatter(df['Compactness'], df['log_Age'], alpha=0.1, s=2, c='blue')
    m, b = np.polyfit(df['Compactness'], df['log_Age'], 1)
    x = np.linspace(df['Compactness'].min(), df['Compactness'].max(), 100)
    ax.plot(x, m*x + b, 'r-', label='Raw Fit')
    ax.set_xlabel(r'Compactness ($\Sigma \propto M/R_e^2$)')
    ax.set_ylabel(r'$\log(\mathrm{D4000})$')
    ax.set_title('Raw Age-Size Relation')
    
    # 2. Partial Plot
    # We need to re-calculate residuals for plotting or pass them
    # Quick re-calc
    X_control = df[['log_MgFe', 'log_Mass']].values
    y_age = df['log_Age'].values
    y_compact = df['Compactness'].values
    
    res_age = y_age - LinearRegression().fit(X_control, y_age).predict(X_control)
    res_compact = y_compact - LinearRegression().fit(X_control, y_compact).predict(X_control)
    
    ax = axes[1]
    ax.scatter(res_compact, res_age, alpha=0.1, s=2, c='purple')
    m, b = np.polyfit(res_compact, res_age, 1)
    x = np.linspace(res_compact.min(), res_compact.max(), 100)
    ax.plot(x, m*x + b, 'r-', label=f'Partial Fit (r={results["r_partial"]:.3f})')
    ax.set_xlabel(r'Residual Compactness')
    ax.set_ylabel(r'Residual Age')
    ax.set_title('Controlled Size-Age Relation\n(Fixed Mass & [Mg/Fe])')
    ax.legend()
    
    plt.tight_layout()
    out_path = os.path.join(FIGURES_DIR, 'sdss_test_k_size_age.png')
    plt.savefig(out_path, dpi=150)
    print(f"Figure saved to {out_path}")

def main():
    df = load_data()
    results, df_clean = analyze_size_age(df)
    
    out_path = os.path.join(RESULTS_DIR, 'sdss_test_k_results.json')
    with open(out_path, 'w') as f:
        json.dump(results, f, indent=2)
        
    create_figure(df_clean, results)
    
    print("\nSUMMARY TEST K:")
    print("Standard Prediction: r > 0 (Compact = Older)")
    print("TEP Prediction: r < 0 (Compact = Younger/Slower Time)")
    print(f"Observed Partial r: {results['r_partial']:.4f}")
    
    if results['r_partial'] < 0:
        print("RESULT: CONSISTENT with TEP (Sign Flip Observed)")
    elif results['r_partial'] < 0.05:
         print("RESULT: WEAK/NULL (Standard relation suppressed)")
    else:
        print("RESULT: STANDARD PHYSICS (Compact galaxies are older)")

if __name__ == "__main__":
    main()
