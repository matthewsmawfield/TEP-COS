#!/usr/bin/env python3
"""
Generate publication-quality figures for TEP-COS manuscript.

Figures:
1. Primary dipole detection (δV vs x_CMB)
2. Two-fluid comparison (Stars vs Gas)
3. Stratification summary (mass, inclination, morphology)
4. Lopsidedness correlation
5. Sky map of velocity asymmetry

Author: Matthew Lukin Smawfield
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from pathlib import Path
import json

# Publication style
plt.rcParams.update({
    'font.family': 'serif',
    'font.serif': ['Times New Roman', 'Times', 'DejaVu Serif'],
    'font.size': 10,
    'axes.labelsize': 11,
    'axes.titlesize': 12,
    'legend.fontsize': 9,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'figure.dpi': 150,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'mathtext.fontset': 'stix',
})

# Colors
COLOR_STELLAR = '#2E86AB'  # Blue
COLOR_GAS = '#E94F37'      # Red
COLOR_POSITIVE = '#D64933'
COLOR_NEGATIVE = '#1E88E5'
COLOR_NEUTRAL = '#757575'

def load_data():
    """Load all analysis results."""
    base = Path('/Users/matthewsmawfield/www/TEP-COS/results/outputs')
    
    data = {}
    
    # Stellar per-galaxy
    stellar_path = base / 'step_2_0_per_galaxy_N2000_corrected.csv'
    if stellar_path.exists():
        data['stellar'] = pd.read_csv(stellar_path)

    stellar_summary_path = base / 'step_2_0_cosmic_coriolis_summary_N2000_corrected.json'
    if stellar_summary_path.exists():
        with open(stellar_summary_path) as f:
            data['stellar_summary'] = json.load(f)
    
    # Gas per-galaxy
    gas_path = base / 'step_2_0_per_galaxy_N2000_gas_corrected.csv'
    if gas_path.exists():
        data['gas'] = pd.read_csv(gas_path)

    gas_summary_path = base / 'step_2_0_cosmic_coriolis_summary_N2000_gas_corrected.json'
    if gas_summary_path.exists():
        with open(gas_summary_path) as f:
            data['gas_summary'] = json.load(f)
    
    # Comprehensive discovery
    disc_path = base / 'discovery_comprehensive_corrected_v2/step_2_9_comprehensive_per_galaxy.csv'
    if disc_path.exists():
        data['discovery'] = pd.read_csv(disc_path)
    
    # Metallicity/peculiar
    metal_path = base / 'discovery_corrected_v4/step_2_10_metallicity_peculiar.csv'
    if metal_path.exists():
        data['metallicity'] = pd.read_csv(metal_path)

    metal_summary_path = base / 'discovery_corrected_v4/step_2_10_metallicity_peculiar_summary.json'
    if metal_summary_path.exists():
        with open(metal_summary_path) as f:
            data['metallicity_summary'] = json.load(f)
    
    # Stratification
    strat_path = base / 'tomography_N2000/step_2_8_expanded_stratification.json'
    if strat_path.exists():
        with open(strat_path) as f:
            data['stratification'] = json.load(f)
    
    return data

def figure_1_primary_dipole(data, out_dir):
    """Figure 1: Primary dipole detection."""
    df = data['stellar']
    
    fig, ax = plt.subplots(figsize=(7, 5))
    
    # Scatter plot
    x = df['x_cmb'].values
    y = df['delta_v'].values
    
    # Color by sign
    colors = np.where(y > 0, COLOR_POSITIVE, COLOR_NEGATIVE)
    
    ax.scatter(x, y, c=colors, alpha=0.3, s=15, edgecolors='none')
    
    summary = data.get('stellar_summary', {})
    a = summary.get('fit', {}).get('a', np.nan)
    a_se = summary.get('fit', {}).get('a_se', np.nan)
    b = summary.get('fit', {}).get('b', np.nan)
    p = summary.get('permutation', {}).get('p_value_pair', np.nan)

    x_fit = np.linspace(-1, 1, 100)
    y_fit = a * x_fit + b

    ax.plot(x_fit, y_fit, 'k-', lw=2, label=f'Slope = {a:.2f} ± {a_se:.2f} km/s (p = {p:.2f})')
    
    # Binned means
    bins = np.linspace(-1, 1, 11)
    bin_centers = (bins[:-1] + bins[1:]) / 2
    bin_means = []
    bin_errs = []
    for i in range(len(bins)-1):
        mask = (x >= bins[i]) & (x < bins[i+1])
        if mask.sum() > 10:
            bin_means.append(np.mean(y[mask]))
            bin_errs.append(np.std(y[mask]) / np.sqrt(mask.sum()))
        else:
            bin_means.append(np.nan)
            bin_errs.append(np.nan)
    
    ax.errorbar(bin_centers, bin_means, yerr=bin_errs, fmt='ko', ms=8, 
                capsize=3, capthick=1.5, elinewidth=1.5, zorder=10,
                label='Binned means')
    
    ax.axhline(0, color='gray', ls='--', lw=0.8, alpha=0.5)
    ax.axvline(0, color='gray', ls='--', lw=0.8, alpha=0.5)
    
    ax.set_xlabel(r'CMB Projection ($x_{\rm CMB} = \cos\theta$)')
    ax.set_ylabel(r'Velocity Asymmetry $\delta V$ (km/s)')
    ax.set_title('Cosmic Coriolis: Stellar Velocity Dipole')
    ax.legend(loc='upper left')
    ax.set_xlim(-1.1, 1.1)
    
    # Add annotation
    ax.annotate('Toward CMB\nApex', xy=(0.8, ax.get_ylim()[1]*0.8), 
                ha='center', fontsize=9, style='italic')
    ax.annotate('Away from\nCMB Apex', xy=(-0.8, ax.get_ylim()[1]*0.8), 
                ha='center', fontsize=9, style='italic')
    
    plt.tight_layout()
    fig.savefig(out_dir / 'figure_1_primary_dipole.png')
    fig.savefig(out_dir / 'figure_1_primary_dipole.pdf')
    plt.close()
    print(f"[SUCCESS] Saved Figure 1: Primary Dipole")

def figure_2_two_fluid(data, out_dir):
    """Figure 2: Two-fluid comparison."""
    df_star = data['stellar']
    df_gas = data['gas']

    star_summary = data.get('stellar_summary', {})
    gas_summary = data.get('gas_summary', {})
    p_star = star_summary.get('permutation', {}).get('p_value_pair', np.nan)
    p_gas = gas_summary.get('permutation', {}).get('p_value_pair', np.nan)
    
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5), sharey=True)
    
    for ax, df, label, color, title in [
        (axes[0], df_star, 'Stellar', COLOR_STELLAR, 'Stellar Kinematics'),
        (axes[1], df_gas, 'Gas (Hα)', COLOR_GAS, 'Gas Kinematics')
    ]:
        x = df['x_cmb'].values
        y = df['delta_v'].values
        
        ax.scatter(x, y, c=color, alpha=0.2, s=10, edgecolors='none')
        
        # Fit
        from sklearn.linear_model import HuberRegressor
        mask = np.isfinite(x) & np.isfinite(y)
        X = x[mask].reshape(-1, 1)
        Y = y[mask]
        model = HuberRegressor().fit(X, Y)
        
        x_fit = np.linspace(-1, 1, 100)
        y_fit = model.predict(x_fit.reshape(-1, 1))
        
        ax.plot(x_fit, y_fit, color='black', lw=2)
        
        # Binned means
        bins = np.linspace(-1, 1, 9)
        bin_centers = (bins[:-1] + bins[1:]) / 2
        bin_means = []
        bin_errs = []
        for i in range(len(bins)-1):
            m = (x >= bins[i]) & (x < bins[i+1])
            if m.sum() > 10:
                bin_means.append(np.mean(y[m]))
                bin_errs.append(np.std(y[m]) / np.sqrt(m.sum()))
            else:
                bin_means.append(np.nan)
                bin_errs.append(np.nan)
        
        ax.errorbar(bin_centers, bin_means, yerr=bin_errs, fmt='o', 
                    color='black', ms=7, capsize=3, zorder=10)
        
        ax.axhline(0, color='gray', ls='--', lw=0.8, alpha=0.5)
        ax.set_xlabel(r'$x_{\rm CMB}$')
        ax.set_title(f'{title}\nSlope = {model.coef_[0]:.1f} km/s')
        ax.set_xlim(-1.1, 1.1)
    
    axes[0].set_ylabel(r'$\delta V$ (km/s)')
    
    # Add significance annotations
    axes[0].annotate(f'p = {p_star:.2f}', xy=(0.95, 0.95), xycoords='axes fraction',
                     ha='right', va='top', fontsize=10,
                     bbox=dict(boxstyle='round', facecolor='lightgray', alpha=0.8))
    axes[1].annotate(f'p = {p_gas:.2f}', xy=(0.95, 0.95), xycoords='axes fraction',
                     ha='right', va='top', fontsize=10,
                     bbox=dict(boxstyle='round', facecolor='lightgray', alpha=0.8))
    
    plt.tight_layout()
    fig.savefig(out_dir / 'figure_2_two_fluid.png')
    fig.savefig(out_dir / 'figure_2_two_fluid.pdf')
    plt.close()
    print(f"[SUCCESS] Saved Figure 2: Two-Fluid Comparison")

def figure_3_stratification(data, out_dir):
    """Figure 3: Stratification summary."""
    strat = data['stratification']['stellar']
    
    fig, axes = plt.subplots(1, 3, figsize=(12, 4))
    
    # Mass stratification
    ax = axes[0]
    mass_data = strat['mass']
    labels = [d['label'].split('(')[0].strip() for d in mass_data]
    slopes = [d['slope'] for d in mass_data]
    errs = [d['slope_err'] for d in mass_data]
    colors = [COLOR_POSITIVE if s > 5 else COLOR_NEUTRAL for s in slopes]
    
    bars = ax.barh(range(len(labels)), slopes, xerr=errs, color=colors, 
                   capsize=3, alpha=0.8)
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels)
    ax.axvline(0, color='black', lw=0.8)
    ax.set_xlabel('Slope (km/s)')
    ax.set_title('Mass Stratification\n(Velocity Dispersion)')
    
    # Inclination stratification
    ax = axes[1]
    inc_data = strat['inclination']
    labels = [d['label'].split('(')[0].strip() for d in inc_data]
    slopes = [d['slope'] for d in inc_data]
    errs = [d['slope_err'] for d in inc_data]
    colors = [COLOR_POSITIVE if s > 5 else COLOR_NEUTRAL for s in slopes]
    
    ax.barh(range(len(labels)), slopes, xerr=errs, color=colors, 
            capsize=3, alpha=0.8)
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels)
    ax.axvline(0, color='black', lw=0.8)
    ax.set_xlabel('Slope (km/s)')
    ax.set_title('Inclination Stratification\n(Axial Ratio b/a)')
    
    # Morphology stratification
    ax = axes[2]
    morph_data = strat['morphology']
    labels = [d['label'].split('(')[0].strip() for d in morph_data]
    slopes = [d['slope'] for d in morph_data]
    errs = [d['slope_err'] for d in morph_data]
    colors = [COLOR_POSITIVE if s > 5 else COLOR_NEUTRAL for s in slopes]
    
    ax.barh(range(len(labels)), slopes, xerr=errs, color=colors, 
            capsize=3, alpha=0.8)
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels)
    ax.axvline(0, color='black', lw=0.8)
    ax.set_xlabel('Slope (km/s)')
    ax.set_title('Morphology Stratification\n(Sérsic Index)')
    
    plt.tight_layout()
    fig.savefig(out_dir / 'figure_3_stratification.png')
    fig.savefig(out_dir / 'figure_3_stratification.pdf')
    plt.close()
    print(f"[SUCCESS] Saved Figure 3: Stratification")

def figure_4_lopsidedness(data, out_dir):
    """Figure 4: Lopsidedness correlation."""
    df = data['metallicity']

    summary = data.get('metallicity_summary', {})
    lop = summary.get('lopsidedness_slope', {})
    slope = lop.get('slope', np.nan)
    err = lop.get('err', np.nan)
    
    if 'lopsidedness' not in df.columns:
        print("[SKIP] No lopsidedness data")
        return
    
    fig, ax = plt.subplots(figsize=(6, 5))
    
    valid = df.dropna(subset=['lopsidedness', 'x_cmb'])
    x = valid['x_cmb'].values
    y = valid['lopsidedness'].values
    
    ax.scatter(x, y, c=COLOR_STELLAR, alpha=0.3, s=15, edgecolors='none')
    
    # Fit
    from sklearn.linear_model import HuberRegressor
    mask = np.isfinite(x) & np.isfinite(y)
    X = x[mask].reshape(-1, 1)
    Y = y[mask]
    model = HuberRegressor().fit(X, Y)
    
    x_fit = np.linspace(-1, 1, 100)
    y_fit = model.predict(x_fit.reshape(-1, 1))
    
    ax.plot(x_fit, y_fit, 'k-', lw=2, label=f'Slope = {slope:.2f} ± {err:.2f} km/s (null)')
    
    # Binned means
    bins = np.linspace(-1, 1, 9)
    bin_centers = (bins[:-1] + bins[1:]) / 2
    bin_means = []
    bin_errs = []
    for i in range(len(bins)-1):
        m = (x >= bins[i]) & (x < bins[i+1])
        if m.sum() > 10:
            bin_means.append(np.mean(y[m]))
            bin_errs.append(np.std(y[m]) / np.sqrt(m.sum()))
        else:
            bin_means.append(np.nan)
            bin_errs.append(np.nan)
    
    ax.errorbar(bin_centers, bin_means, yerr=bin_errs, fmt='ko', ms=8, 
                capsize=3, zorder=10)
    
    ax.set_xlabel(r'CMB Projection ($x_{\rm CMB}$)')
    ax.set_ylabel('Kinematic Lopsidedness (km/s)')
    ax.set_title('Lopsidedness-CMB Correlation')
    ax.legend(loc='upper left')
    
    # Significance box
    ax.annotate('Null', xy=(0.95, 0.95), xycoords='axes fraction',
                ha='right', va='top', fontsize=11,
                bbox=dict(boxstyle='round', facecolor='lightgray', alpha=0.8))
    
    plt.tight_layout()
    fig.savefig(out_dir / 'figure_4_lopsidedness.png')
    fig.savefig(out_dir / 'figure_4_lopsidedness.pdf')
    plt.close()
    print(f"[SUCCESS] Saved Figure 4: Lopsidedness")

def figure_5_sky_map(data, out_dir):
    """Figure 5: Sky map of velocity asymmetry."""
    df = data['stellar']
    
    fig = plt.figure(figsize=(10, 5))
    ax = fig.add_subplot(111, projection='mollweide')
    
    # Convert RA/Dec to radians for Mollweide
    ra = np.radians(df['ra_deg'].values - 180)  # Center at RA=180
    dec = np.radians(df['dec_deg'].values)
    dv = df['delta_v'].values
    
    # Color normalization centered at 0
    vmax = np.percentile(np.abs(dv[np.isfinite(dv)]), 95)
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)
    
    sc = ax.scatter(ra, dec, c=dv, cmap='RdBu_r', norm=norm, 
                    s=8, alpha=0.6, edgecolors='none')
    
    # Mark CMB apex
    cmb_ra = np.radians(168 - 180)
    cmb_dec = np.radians(-7)
    ax.scatter([cmb_ra], [cmb_dec], marker='*', s=200, c='gold', 
               edgecolors='black', linewidths=1, zorder=10, label='CMB Apex')
    
    # Mark CMB anti-apex
    anti_ra = np.radians(168 + 180 - 180)
    anti_dec = np.radians(7)
    ax.scatter([anti_ra], [anti_dec], marker='*', s=200, c='white', 
               edgecolors='black', linewidths=1, zorder=10, label='CMB Anti-apex')
    
    ax.set_title('Sky Distribution of Velocity Asymmetry', fontsize=12, pad=20)
    ax.grid(True, alpha=0.3)
    
    cbar = plt.colorbar(sc, ax=ax, orientation='horizontal', 
                        fraction=0.05, pad=0.1, aspect=40)
    cbar.set_label(r'$\delta V$ (km/s)')
    
    ax.legend(loc='lower right', fontsize=9)
    
    plt.tight_layout()
    fig.savefig(out_dir / 'figure_5_sky_map.png')
    fig.savefig(out_dir / 'figure_5_sky_map.pdf')
    plt.close()
    print(f"[SUCCESS] Saved Figure 5: Sky Map")

def main():
    out_dir = Path('/Users/matthewsmawfield/www/TEP-COS/results/figures/manuscript')
    out_dir.mkdir(parents=True, exist_ok=True)
    
    print("[PROCESS] Loading data...")
    data = load_data()
    
    print("[PROCESS] Generating figures...")
    
    if 'stellar' in data:
        figure_1_primary_dipole(data, out_dir)
    
    if 'stellar' in data and 'gas' in data:
        figure_2_two_fluid(data, out_dir)
    
    if 'stratification' in data:
        figure_3_stratification(data, out_dir)
    
    if 'metallicity' in data:
        figure_4_lopsidedness(data, out_dir)
    
    if 'stellar' in data:
        figure_5_sky_map(data, out_dir)
    
    print(f"\n[SUCCESS] All figures saved to {out_dir}")

if __name__ == '__main__':
    main()
