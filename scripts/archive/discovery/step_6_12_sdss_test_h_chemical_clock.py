#!/usr/bin/env python3
"""
Step 6.12: SDSS Test H - Chemical Clock Discrepancy

Hypothesis:
[Mg/Fe] is a "chemical clock" (SN II vs SN Ia timescale).
Spectroscopic age (D4000) is a "cooling clock".
Under TEP, time dilation affects cooling/evolution (D4000) but NOT nucleosynthesis yields ([Mg/Fe] is frozen at formation).
High-sigma galaxies should show DISCREPANT [Mg/Fe] at fixed spectroscopic age.

TEP Prediction:
r(Δ[Mg/Fe], σ) > 0 at fixed Age and Z.
(Enhanced alpha-abundance relative to what is expected for that age, 
because the 'age' clock ran slower than the 'enrichment' clock? 
Actually: If time flows slower, the Ia delay (1 Gyr) takes 'longer' in proper time.
So for a fixed proper time elapsed (Age), fewer Ia's have exploded.
Thus [Mg/Fe] stays high longer. 
So at fixed Age, higher sigma -> Higher [Mg/Fe]. Correct.)

Data:
- sdss_spectral_indices.csv (Mgb, Fe5270, Fe5335, D4000, Hbeta)
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
    path = os.path.join(DATA_DIR, 'sdss_spectral_indices.csv')
    if not os.path.exists(path):
        raise FileNotFoundError(f"Data file not found: {path}")
    
    print(f"Loading {path}...")
    df = pd.read_csv(path)
    return df

def analyze_chemical_clock(df):
    print("Analyzing chemical clock discrepancy...")
    
    # 1. Compute Indices
    # Thomas et al. (2003) [MgFe]' index as total metallicity proxy
    # [MgFe]' = sqrt(Mgb * (0.72*Fe5270 + 0.28*Fe5335))
    df['Fe_idx'] = 0.72 * df['fe5270'] + 0.28 * df['fe5335']
    # Avoid negative sqrt
    df['MgFe_prime'] = np.sqrt(np.maximum(0, df['mgb'] * df['Fe_idx']))
    
    # [Mg/Fe] proxy
    # We use log(Mgb / <Fe>) as a simple proxy, or Mgb - <Fe> in index space
    # Common proxy: Mgb / <Fe> where <Fe> = (Fe5270+Fe5335)/2
    df['Fe_avg'] = (df['fe5270'] + df['fe5335']) / 2
    df['MgFe_ratio'] = df['mgb'] / df['Fe_avg']
    df['log_MgFe'] = np.log10(df['MgFe_ratio'])
    
    # Spectroscopic Age Proxy: D4000 (primary)
    # Hbeta is also useful but more sensitive to recent SF. D4000 tracks old population age better.
    df['log_Age'] = np.log10(df['d4000'])
    
    # Sigma
    df['log_sigma'] = np.log10(df['veldisp'])
    
    # Quality Cuts
    # ETG selection: D4000 > 1.5, bptclass often -1 (inactive) or AGN. 
    # We want valid indices.
    mask = (
        (df['d4000'] > 1.2) & (df['d4000'] < 2.5) &
        (df['mgb'] > 0) & (df['fe5270'] > 0) & (df['fe5335'] > 0) &
        (df['veldisp'] > 50) & (df['veldisp'] < 450) &
        (df['z_err'] < 0.001) &
        (df['MgFe_prime'] > 0)
    )
    df_clean = df[mask].copy()
    print(f"  Selected {len(df_clean)} / {len(df)} galaxies (Quality + ETG-like cuts)")
    
    # 2. Control for Age and Metallicity
    # We want to predict [Mg/Fe] from Age and Z (standard physics)
    # And look for residual correlation with sigma.
    
    X = df_clean[['log_Age', 'MgFe_prime']].values
    y = df_clean['log_MgFe'].values
    
    reg = LinearRegression().fit(X, y)
    df_clean['log_MgFe_pred'] = reg.predict(X)
    df_clean['MgFe_resid'] = df_clean['log_MgFe'] - df_clean['log_MgFe_pred']
    
    # 3. Correlate Residual with Sigma
    r_simple, p_simple = stats.pearsonr(df_clean['log_sigma'], df_clean['log_MgFe'])
    r_resid, p_resid = stats.pearsonr(df_clean['log_sigma'], df_clean['MgFe_resid'])
    
    print(f"\nResults:")
    print(f"  Simple r(log[Mg/Fe], log σ): {r_simple:.4f} (p={p_simple:.2e})")
    print(f"  Controlled r(Δ[Mg/Fe], log σ): {r_resid:.4f} (p={p_resid:.2e})")
    
    return {
        'r_simple': float(r_simple),
        'p_simple': float(p_simple),
        'r_controlled': float(r_resid),
        'p_controlled': float(p_resid),
        'n_sample': int(len(df_clean)),
        'coeffs': {
            'age_coeff': float(reg.coef_[0]),
            'z_coeff': float(reg.coef_[1]),
            'intercept': float(reg.intercept_)
        }
    }, df_clean

def create_figure(df, results):
    print("Generating figure...")
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # Plot 1: [Mg/Fe] vs Sigma (Simple)
    # sns.regplot(data=df, x='log_sigma', y='log_MgFe', ax=axes[0], 
    #             scatter_kws={'alpha': 0.1, 's': 2}, line_kws={'color': 'red'})
    ax = axes[0]
    ax.scatter(df['log_sigma'], df['log_MgFe'], alpha=0.1, s=2, c='blue', label='Data')
    # Fit line
    m, b = np.polyfit(df['log_sigma'], df['log_MgFe'], 1)
    x_range = np.linspace(df['log_sigma'].min(), df['log_sigma'].max(), 100)
    ax.plot(x_range, m*x_range + b, color='red', label=f'Fit (r={results["r_simple"]:.3f})')
    
    ax.set_title(f"Raw Relation\nr={results['r_simple']:.3f}")
    ax.set_xlabel(r'$\log(\sigma)$')
    ax.set_ylabel(r'$\log(\mathrm{Mg}/\mathrm{Fe})$')
    ax.legend()
    
    # Plot 2: Residual [Mg/Fe] vs Sigma
    # sns.regplot(data=df, x='log_sigma', y='MgFe_resid', ax=axes[1],
    #             scatter_kws={'alpha': 0.1, 's': 2}, line_kws={'color': 'green'})
    ax = axes[1]
    ax.scatter(df['log_sigma'], df['MgFe_resid'], alpha=0.1, s=2, c='green', label='Residuals')
    # Fit line
    m_res, b_res = np.polyfit(df['log_sigma'], df['MgFe_resid'], 1)
    ax.plot(x_range, m_res*x_range + b_res, color='red', label=f'Fit (r={results["r_controlled"]:.3f})')
    
    ax.set_title(f"Controlled Residual\n(Fixed Age & Z)\nr={results['r_controlled']:.3f}")
    ax.set_xlabel(r'$\log(\sigma)$')
    ax.set_ylabel(r'$\Delta \log(\mathrm{Mg}/\mathrm{Fe})$')
    ax.legend()
    
    plt.tight_layout()
    fig_path = os.path.join(FIGURES_DIR, 'sdss_test_h_chemical_clock.png')
    plt.savefig(fig_path, dpi=150)
    print(f"Figure saved to {fig_path}")

def main():
    # Load
    df = load_data()
    
    # Analyze
    results, df_clean = analyze_chemical_clock(df)
    
    # Save Results
    out_path = os.path.join(RESULTS_DIR, 'sdss_test_h_results.json')
    with open(out_path, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"Results saved to {out_path}")
    
    # Figure
    create_figure(df_clean, results)
    
    # Summary
    print("\nSUMMARY TEST H:")
    print(f"TEP Prediction: r > 0 (Positive residual correlation)")
    print(f"Observed: r = {results['r_controlled']:.4f}")
    if results['r_controlled'] > 0.05 and results['p_controlled'] < 0.05:
        print("RESULT: CONSISTENT with TEP.")
    else:
        print("RESULT: INCONSISTENT or NULL.")

if __name__ == "__main__":
    main()
