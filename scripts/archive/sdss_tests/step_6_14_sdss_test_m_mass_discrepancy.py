#!/usr/bin/env python3
"""
Step 6.14: SDSS Test M - Dynamical vs Stellar Mass Discrepancy

Hypothesis:
Dynamical mass (M_dyn ~ sigma^2 Re / G) measures the total gravitational potential.
Stellar mass (M_star) is derived from SED/spectral fitting assuming standard evolutionary tracks.
Under TEP, time dilation in deep potentials affects the apparent evolution rate, potentially altering 
the M/L ratio inferred from stellar populations (making them appear older/more massive per unit light?).

Prediction:
If time flows slower in deep potentials, stars evolve slower. 
However, standard models assume standard rates. 
If a galaxy appears "older" or "younger" than it is, M/L changes.
Wait, if TEP makes them appear younger (Test K prediction, though contradicted), M/L would be lower.
If TEP makes them appear older (Test B/D results suggested this? No, Test B said SED masses were LOWER than PCA masses in high sigma).

Let's stick to the hypothesis in the plan:
"Stellar population masses are OVERESTIMATED at high sigma relative to dynamical masses" 
-> r(log(M_star/M_dyn), sigma) > 0.

Wait, if Test B showed SED masses are UNDERESTIMATED relative to PCA (geometrical) masses?
Let's just look for the trend. The deviation from the "Fundamental Plane" or standard M_star/M_dyn relation.

Method:
1. Load galaxies.
2. Calculate R_e in kpc.
3. Calculate M_dyn = 5 * sigma^2 * R_e / G.
4. Calculate Delta = log10(M_star) - log10(M_dyn).
5. Correlate Delta with sigma, controlling for structure (concentration, n).
"""

import pandas as pd
import numpy as np
from scipy import stats
import os
import json
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
from astropy.cosmology import FlatLambdaCDM

# Paths
DATA_DIR = os.path.join(os.path.dirname(__file__), '..', '..', 'data', 'sdss')
RESULTS_DIR = os.path.join(os.path.dirname(__file__), '..', '..', 'results', 'outputs')
FIGURES_DIR = os.path.join(os.path.dirname(__file__), '..', '..', 'results', 'figures')

def load_data():
    path = os.path.join(DATA_DIR, 'sdss_galaxies.csv')
    if not os.path.exists(path):
        raise FileNotFoundError(f"Data file not found: {path}")
    print(f"Loading {path}...")
    df = pd.read_csv(path)
    return df

def analyze_mass_discrepancy(df):
    print("Analyzing Mass Discrepancy...")
    
    # 1. Cosmology for Size
    cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
    # Filter valid redshift first to avoid warnings
    df = df[(df['redshift'] > 0.01) & (df['redshift'] < 0.3)].copy()
    
    kpc_per_arcsec = cosmo.kpc_proper_per_arcmin(df['redshift'].values).value / 60.0
    df['Re_kpc'] = df['petroR50_r'] * kpc_per_arcsec
    
    # 2. Dynamical Mass
    # M_dyn = K * sigma^2 * Re / G
    # G = 4.302e-6 kpc (km/s)^2 / M_sun
    G = 4.302e-6
    K = 5.0 # Standard virial coefficient
    
    # Sigma in km/s
    # veldisp is in km/s.
    # Filter valid sigma
    df = df[(df['veldisp'] > 50) & (df['veldisp'] < 450) & (df['Re_kpc'] > 0.5) & (df['Re_kpc'] < 30)].copy()
    
    df['M_dyn'] = K * (df['veldisp']**2) * df['Re_kpc'] / G
    df['log_M_dyn'] = np.log10(df['M_dyn'])
    
    # 3. Stellar Mass
    # log_mass is already log10(M_sun)
    df = df[(df['log_mass'] > 8) & (df['log_mass'] < 13)].copy()
    
    # 4. Discrepancy
    # IMF mismatch, Dark Matter fraction, etc. are standard explanations.
    # We look for a residual correlation with sigma AFTER controlling for structure.
    df['Delta_M'] = df['log_mass'] - df['log_M_dyn']
    
    df['log_sigma'] = np.log10(df['veldisp'])
    
    # Structure proxy: concentration = R90/R50
    df['conc'] = df['petroR90_r'] / df['petroR50_r']
    
    # 5. Regression / Control
    # M_dyn includes Dark Matter. M_star is just stars.
    # Typically M_star/M_dyn < 1. 
    # The fraction M_star/M_dyn usually increases with Mass/Sigma until huge masses where DM dominates?
    # Or fundamental plane tilt.
    
    # Control for Concentration (morphology) and Size (Re)
    # We want to isolate the sigma dependence.
    
    features = ['conc', 'Re_kpc']
    X = df[features].values
    y = df['Delta_M'].values
    
    reg = LinearRegression().fit(X, y)
    df['Delta_M_resid'] = y - reg.predict(X)
    
    # Correlation
    r_simple, p_simple = stats.pearsonr(df['log_sigma'], df['Delta_M'])
    r_controlled, p_controlled = stats.pearsonr(df['log_sigma'], df['Delta_M_resid'])
    
    print(f"N = {len(df)}")
    print(f"Simple r(Delta M, sigma): {r_simple:.4f} (p={p_simple:.2e})")
    print(f"Controlled r(Delta M resid, sigma): {r_controlled:.4f} (p={p_controlled:.2e})")
    
    return {
        'r_simple': float(r_simple),
        'p_simple': float(p_simple),
        'r_controlled': float(r_controlled),
        'p_controlled': float(p_controlled),
        'n_sample': int(len(df))
    }, df

def create_figure(df, results):
    print("Generating figure...")
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # Plot 1: Raw
    # Scatter plot
    ax = axes[0]
    # Hexbin for density
    hb = ax.hexbin(df['log_sigma'], df['Delta_M'], gridsize=50, cmap='inferno', mincnt=1)
    ax.set_xlabel(r'$\log(\sigma)$')
    ax.set_ylabel(r'$\log(M_*) - \log(M_{dyn})$')
    ax.set_title(f"Raw Mass Discrepancy\nr={results['r_simple']:.3f}")
    cb = plt.colorbar(hb, ax=ax)
    cb.set_label('Count')
    
    # Fit
    m, b = np.polyfit(df['log_sigma'], df['Delta_M'], 1)
    x = np.linspace(df['log_sigma'].min(), df['log_sigma'].max(), 100)
    ax.plot(x, m*x + b, 'c--', lw=2)
    
    # Plot 2: Controlled
    ax = axes[1]
    hb2 = ax.hexbin(df['log_sigma'], df['Delta_M_resid'], gridsize=50, cmap='inferno', mincnt=1)
    ax.set_xlabel(r'$\log(\sigma)$')
    ax.set_ylabel(r'Residual $\Delta M$ (Structure Controlled)')
    ax.set_title(f"Controlled Discrepancy\nr={results['r_controlled']:.3f}")
    cb2 = plt.colorbar(hb2, ax=ax)
    
    # Fit
    m2, b2 = np.polyfit(df['log_sigma'], df['Delta_M_resid'], 1)
    ax.plot(x, m2*x + b2, 'c--', lw=2)
    
    plt.tight_layout()
    out_path = os.path.join(FIGURES_DIR, 'sdss_test_m_mass_discrepancy.png')
    plt.savefig(out_path, dpi=150)
    print(f"Figure saved to {out_path}")

def main():
    df = load_data()
    results, df_clean = analyze_mass_discrepancy(df)
    
    out_path = os.path.join(RESULTS_DIR, 'sdss_test_m_results.json')
    with open(out_path, 'w') as f:
        json.dump(results, f, indent=2)
        
    create_figure(df_clean, results)
    
    print("\nSUMMARY TEST M:")
    print(f"TEP Prediction: r > 0 (Stellar mass over-estimated at high sigma)")
    print(f"Observed: r = {results['r_controlled']:.4f}")
    if results['r_controlled'] > 0.05:
        print("RESULT: CONSISTENT with TEP.")
    elif results['r_controlled'] < -0.05:
        print("RESULT: CONTRADICTED (Stellar mass fraction decreases with sigma).")
    else:
        print("RESULT: NULL/INCONCLUSIVE.")

if __name__ == "__main__":
    main()
