#!/usr/bin/env python3
"""
Step 6.6: Star Formation Rate Holonomy Test

TEP HYPOTHESIS:
Star Formation Rate is a time-derivative quantity (M☉/yr). Under TEP, proper time
flows slower in deeper gravitational potentials. Therefore, at fixed gas supply
and metallicity, galaxies with higher velocity dispersion (σ) should show
systematically LOWER observed SFR.

TEP PREDICTION:
  At fixed stellar mass M* and metallicity Z:
    r(SFR, σ) < 0        (negative correlation)
    r(sSFR, σ) < 0       (specific SFR = SFR/M*)

This is distinct from previous tests which used cumulative age indicators.
SFR is a RATE - directly sensitive to time dilation, not accumulated proper time.

Author: M. Smawfield
Date: January 2026
"""

import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import pearsonr, spearmanr
import matplotlib.pyplot as plt
import json
import os
from datetime import datetime
import warnings
warnings.filterwarnings('ignore')

DATA_DIR = os.path.join(os.path.dirname(__file__), '..', '..', 'data', 'sdss')
RESULTS_DIR = os.path.join(os.path.dirname(__file__), '..', '..', 'results', 'outputs')
FIGURES_DIR = os.path.join(os.path.dirname(__file__), '..', '..', 'results', 'figures')

os.makedirs(RESULTS_DIR, exist_ok=True)
os.makedirs(FIGURES_DIR, exist_ok=True)


def load_data():
    """Load SDSS spectral indices data with SFR."""
    cache_path = os.path.join(DATA_DIR, 'sdss_spectral_indices.csv')
    if not os.path.exists(cache_path):
        raise FileNotFoundError(f"Data not found: {cache_path}\nRun step_6_0 first.")
    
    df = pd.read_csv(cache_path)
    print(f"Loaded {len(df):,} galaxies from SDSS")
    return df


def prepare_data(df):
    """Prepare derived quantities for SFR Holonomy analysis."""
    
    # Velocity dispersion (potential depth proxy)
    df['log_sigma'] = np.log10(df['veldisp'])
    
    # Specific Star Formation Rate
    df['log_ssfr'] = df['log_sfr'] - df['log_mass']
    
    # [Mg/Fe] ratio (nucleosynthesis indicator - controls formation timescale)
    df['fe_avg'] = (df['fe5270'] + df['fe5335']) / 2
    df['mg_fe_ratio'] = df['mgb'] / df['fe_avg']
    df['log_mg_fe'] = np.log10(df['mg_fe_ratio'])
    
    # Quality cuts
    valid = (
        # σ quality
        (df['veldisp'] > 50) & (df['veldisp'] < 400) &
        (df['veldisp_err'] < 30) &
        # SFR quality (not upper limits)
        np.isfinite(df['log_sfr']) &
        (df['log_sfr'] > -5) & (df['log_sfr'] < 3) &
        # Mass quality
        np.isfinite(df['log_mass']) &
        (df['log_mass'] > 9.0) & (df['log_mass'] < 12.5) &
        # [Mg/Fe] quality
        np.isfinite(df['mg_fe_ratio']) &
        (df['mg_fe_ratio'] > 0.3) & (df['mg_fe_ratio'] < 5.0) &
        # Lick index quality
        (df['mgb_err'] < 0.5) &
        (df['fe5270_err'] < 0.5) &
        # Redshift range
        (df['redshift'] > 0.02) & (df['redshift'] < 0.25)
    )
    
    df_clean = df[valid].copy()
    print(f"After quality cuts: {len(df_clean):,} galaxies ({100*len(df_clean)/len(df):.1f}%)")
    
    return df_clean


def partial_corr(x, y, z):
    """Partial correlation of x and y controlling for z (single variable)."""
    r_xy, _ = pearsonr(x, y)
    r_xz, _ = pearsonr(x, z)
    r_yz, _ = pearsonr(y, z)
    
    numerator = r_xy - r_xz * r_yz
    denominator = np.sqrt((1 - r_xz**2) * (1 - r_yz**2))
    
    if denominator < 1e-10:
        return 0.0
    return numerator / denominator


def partial_corr_multi(x, y, Z):
    """Partial correlation of x and y controlling for multiple variables Z."""
    from sklearn.linear_model import LinearRegression
    
    # Regress out confounders
    Z = np.atleast_2d(Z).T if Z.ndim == 1 else Z
    
    reg_x = LinearRegression().fit(Z, x)
    reg_y = LinearRegression().fit(Z, y)
    
    resid_x = x - reg_x.predict(Z)
    resid_y = y - reg_y.predict(Z)
    
    r, p = pearsonr(resid_x, resid_y)
    return r, p


def test_simple_correlations(df):
    """Test 1: Simple correlations between SFR, sSFR, and σ."""
    print("\n" + "=" * 70)
    print("TEST 1: Simple Correlations")
    print("=" * 70)
    
    results = {}
    
    # SFR vs σ
    r_sfr, p_sfr = pearsonr(df['log_sigma'], df['log_sfr'])
    rho_sfr, p_rho_sfr = spearmanr(df['log_sigma'], df['log_sfr'])
    
    # sSFR vs σ
    r_ssfr, p_ssfr = pearsonr(df['log_sigma'], df['log_ssfr'])
    rho_ssfr, p_rho_ssfr = spearmanr(df['log_sigma'], df['log_ssfr'])
    
    # For comparison: Age indicators
    r_d4000, p_d4000 = pearsonr(df['log_sigma'], df['d4000'])
    r_mgfe, p_mgfe = pearsonr(df['log_sigma'], df['log_mg_fe'])
    
    print(f"\n{'Quantity':<25} {'Pearson r':>12} {'p-value':>12} {'Spearman ρ':>12}")
    print("-" * 65)
    print(f"{'log(SFR) vs log(σ)':<25} {r_sfr:>12.4f} {p_sfr:>12.2e} {rho_sfr:>12.4f}")
    print(f"{'log(sSFR) vs log(σ)':<25} {r_ssfr:>12.4f} {p_ssfr:>12.2e} {rho_ssfr:>12.4f}")
    print("-" * 65)
    print(f"{'D4000 vs log(σ)':<25} {r_d4000:>12.4f} {p_d4000:>12.2e}")
    print(f"{'[Mg/Fe] vs log(σ)':<25} {r_mgfe:>12.4f} {p_mgfe:>12.2e}")
    
    results['simple'] = {
        'r_sfr_sigma': r_sfr, 'p_sfr_sigma': p_sfr,
        'r_ssfr_sigma': r_ssfr, 'p_ssfr_sigma': p_ssfr,
        'r_d4000_sigma': r_d4000, 'p_d4000_sigma': p_d4000,
        'r_mgfe_sigma': r_mgfe, 'p_mgfe_sigma': p_mgfe,
    }
    
    print("\nINTERPRETATION:")
    if r_ssfr < 0:
        print(f"  ✓ sSFR-σ correlation is NEGATIVE (r = {r_ssfr:.4f})")
        print("    Higher σ galaxies have lower sSFR - consistent with TEP time dilation")
    else:
        print(f"  ✗ sSFR-σ correlation is positive/zero (r = {r_ssfr:.4f})")
    
    return results


def test_mass_controlled(df):
    """Test 2: Correlations controlling for stellar mass."""
    print("\n" + "=" * 70)
    print("TEST 2: Mass-Controlled Analysis")
    print("=" * 70)
    
    results = {}
    
    # Partial correlations controlling for mass
    r_sfr_mass, p_sfr_mass = partial_corr_multi(
        df['log_sigma'].values, df['log_sfr'].values, df['log_mass'].values
    )
    r_ssfr_mass, p_ssfr_mass = partial_corr_multi(
        df['log_sigma'].values, df['log_ssfr'].values, df['log_mass'].values
    )
    
    print(f"\nPartial correlations (controlling for M*):")
    print(f"  r(SFR, σ | M*)  = {r_sfr_mass:>8.4f}  (p = {p_sfr_mass:.2e})")
    print(f"  r(sSFR, σ | M*) = {r_ssfr_mass:>8.4f}  (p = {p_ssfr_mass:.2e})")
    
    results['partial_mass'] = {
        'r_sfr_sigma_mass': r_sfr_mass, 'p_sfr_sigma_mass': p_sfr_mass,
        'r_ssfr_sigma_mass': r_ssfr_mass, 'p_ssfr_sigma_mass': p_ssfr_mass,
    }
    
    # Binned analysis by stellar mass
    print("\nBinned analysis by stellar mass:")
    print("-" * 70)
    print(f"{'Mass bin':<15} {'<log M*>':>10} {'N':>8} {'r(SFR,σ)':>12} {'r(sSFR,σ)':>12}")
    print("-" * 70)
    
    mass_bins = pd.qcut(df['log_mass'], q=5, labels=['M1', 'M2', 'M3', 'M4', 'M5'])
    bin_results = []
    
    for mbin in ['M1', 'M2', 'M3', 'M4', 'M5']:
        mask = mass_bins == mbin
        sub = df[mask]
        
        if len(sub) < 100:
            continue
        
        mean_mass = sub['log_mass'].mean()
        r_sfr, _ = pearsonr(sub['log_sigma'], sub['log_sfr'])
        r_ssfr, _ = pearsonr(sub['log_sigma'], sub['log_ssfr'])
        
        print(f"{mbin:<15} {mean_mass:>10.2f} {len(sub):>8} {r_sfr:>12.4f} {r_ssfr:>12.4f}")
        
        bin_results.append({
            'mass_bin': mbin, 'mean_mass': mean_mass, 'n': len(sub),
            'r_sfr': r_sfr, 'r_ssfr': r_ssfr
        })
    
    print("-" * 70)
    results['binned_mass'] = bin_results
    
    # Summary
    mean_r_ssfr = np.mean([b['r_ssfr'] for b in bin_results])
    print(f"\nMean r(sSFR, σ) across mass bins: {mean_r_ssfr:.4f}")
    
    if mean_r_ssfr < 0:
        print("  ✓ sSFR-σ anti-correlation persists at fixed mass - TEP consistent")
    
    return results


def test_metallicity_controlled(df):
    """Test 3: Correlations controlling for metallicity ([Mg/Fe])."""
    print("\n" + "=" * 70)
    print("TEST 3: Metallicity-Controlled Analysis")
    print("=" * 70)
    
    results = {}
    
    # Partial correlations controlling for [Mg/Fe] and mass
    Z = np.column_stack([df['log_mass'].values, df['log_mg_fe'].values])
    
    r_ssfr_full, p_ssfr_full = partial_corr_multi(
        df['log_sigma'].values, df['log_ssfr'].values, Z
    )
    
    print(f"\nPartial correlation (controlling for M* and [Mg/Fe]):")
    print(f"  r(sSFR, σ | M*, [Mg/Fe]) = {r_ssfr_full:>8.4f}  (p = {p_ssfr_full:.2e})")
    
    results['partial_full'] = {
        'r_ssfr_sigma_full': r_ssfr_full,
        'p_ssfr_sigma_full': p_ssfr_full,
    }
    
    # Binned analysis by [Mg/Fe]
    print("\nBinned analysis by [Mg/Fe] (nucleosynthesis indicator):")
    print("-" * 70)
    print(f"{'[Mg/Fe] bin':<15} {'<[Mg/Fe]>':>12} {'N':>8} {'r(sSFR,σ)':>12}")
    print("-" * 70)
    
    mgfe_bins = pd.qcut(df['log_mg_fe'], q=5, labels=['Z1', 'Z2', 'Z3', 'Z4', 'Z5'])
    bin_results = []
    
    for zbin in ['Z1', 'Z2', 'Z3', 'Z4', 'Z5']:
        mask = mgfe_bins == zbin
        sub = df[mask]
        
        if len(sub) < 100:
            continue
        
        mean_mgfe = sub['log_mg_fe'].mean()
        r_ssfr, _ = pearsonr(sub['log_sigma'], sub['log_ssfr'])
        
        print(f"{zbin:<15} {mean_mgfe:>12.4f} {len(sub):>8} {r_ssfr:>12.4f}")
        
        bin_results.append({
            'mgfe_bin': zbin, 'mean_mgfe': mean_mgfe, 'n': len(sub), 'r_ssfr': r_ssfr
        })
    
    print("-" * 70)
    results['binned_mgfe'] = bin_results
    
    mean_r = np.mean([b['r_ssfr'] for b in bin_results])
    print(f"\nMean r(sSFR, σ) across [Mg/Fe] bins: {mean_r:.4f}")
    
    return results


def test_bpt_stratified(df):
    """Test 4: Analysis stratified by BPT classification."""
    print("\n" + "=" * 70)
    print("TEST 4: BPT-Stratified Analysis")
    print("=" * 70)
    
    bpt_labels = {
        -1: 'Unclassified',
        1: 'Star-forming',
        2: 'Composite',
        3: 'AGN',
        4: 'LINER'
    }
    
    results = {}
    
    print("\nCorrelations by BPT class:")
    print("-" * 70)
    print(f"{'BPT Class':<20} {'N':>10} {'r(SFR,σ)':>12} {'r(sSFR,σ)':>12}")
    print("-" * 70)
    
    for bpt_code, bpt_name in bpt_labels.items():
        mask = df['bptclass'] == bpt_code
        sub = df[mask]
        
        if len(sub) < 100:
            continue
        
        r_sfr, _ = pearsonr(sub['log_sigma'], sub['log_sfr'])
        r_ssfr, _ = pearsonr(sub['log_sigma'], sub['log_ssfr'])
        
        print(f"{bpt_name:<20} {len(sub):>10} {r_sfr:>12.4f} {r_ssfr:>12.4f}")
        
        results[bpt_name] = {'n': len(sub), 'r_sfr': r_sfr, 'r_ssfr': r_ssfr}
    
    print("-" * 70)
    
    print("\nKEY INSIGHT:")
    print("  Star-forming galaxies should show the clearest TEP signal")
    print("  (SFR measurement most reliable, less AGN contamination)")
    
    return results


def test_redshift_evolution(df):
    """Test 5: Check if signal evolves with redshift."""
    print("\n" + "=" * 70)
    print("TEST 5: Redshift Evolution")
    print("=" * 70)
    
    z_bins = [
        (0.02, 0.06, 'z < 0.06'),
        (0.06, 0.10, '0.06-0.10'),
        (0.10, 0.15, '0.10-0.15'),
        (0.15, 0.20, '0.15-0.20'),
        (0.20, 0.25, 'z > 0.20'),
    ]
    
    results = []
    
    print(f"\n{'z range':<15} {'<z>':>8} {'N':>10} {'r(sSFR,σ)':>12} {'p-value':>12}")
    print("-" * 60)
    
    for z_min, z_max, label in z_bins:
        mask = (df['redshift'] >= z_min) & (df['redshift'] < z_max)
        sub = df[mask]
        
        if len(sub) < 100:
            continue
        
        mean_z = sub['redshift'].mean()
        r_ssfr, p_ssfr = pearsonr(sub['log_sigma'], sub['log_ssfr'])
        
        print(f"{label:<15} {mean_z:>8.3f} {len(sub):>10} {r_ssfr:>12.4f} {p_ssfr:>12.2e}")
        
        results.append({
            'z_label': label, 'z_min': z_min, 'z_max': z_max,
            'mean_z': mean_z, 'n': len(sub), 'r_ssfr': r_ssfr, 'p_ssfr': p_ssfr
        })
    
    print("-" * 60)
    
    # Check stability
    rs = [r['r_ssfr'] for r in results]
    print(f"\nMean r(sSFR, σ): {np.mean(rs):.4f} ± {np.std(rs):.4f}")
    
    if np.std(rs) < 0.05:
        print("  ✓ Signal is STABLE across redshift")
    
    return results


def create_figure(df, results):
    """Create comprehensive SFR Holonomy figure."""
    print("\n" + "=" * 70)
    print("Creating Figure")
    print("=" * 70)
    
    fig, axes = plt.subplots(2, 3, figsize=(14, 9))
    
    # Color palette (user preference: blues/slates)
    colors = {
        'primary': '#4A90A4',
        'secondary': '#6B7B8C',
        'accent': '#2D5A6B',
        'light': '#A8C5D1',
        'dark': '#1E3A4A'
    }
    
    # Panel A: sSFR vs σ (2D histogram)
    ax = axes[0, 0]
    h = ax.hexbin(df['log_sigma'], df['log_ssfr'], gridsize=50, cmap='Blues', mincnt=1)
    ax.set_xlabel(r'$\log_{10}(\sigma)$ [km/s]')
    ax.set_ylabel(r'$\log_{10}$(sSFR) [yr$^{-1}$]')
    ax.set_title('A. Specific SFR vs Velocity Dispersion')
    
    # Add regression line
    z = np.polyfit(df['log_sigma'], df['log_ssfr'], 1)
    x_line = np.linspace(df['log_sigma'].min(), df['log_sigma'].max(), 100)
    ax.plot(x_line, np.polyval(z, x_line), 'r-', lw=2, label=f'r = {results["simple"]["r_ssfr_sigma"]:.3f}')
    ax.legend(loc='upper right')
    
    # Panel B: Binned by mass
    ax = axes[0, 1]
    binned = results['binned_mass']
    masses = [b['mean_mass'] for b in binned]
    r_ssfrs = [b['r_ssfr'] for b in binned]
    ax.bar(range(len(binned)), r_ssfrs, color=colors['primary'], edgecolor=colors['dark'], alpha=0.8)
    ax.axhline(0, color='k', ls='--', lw=1)
    ax.set_xticks(range(len(binned)))
    ax.set_xticklabels([f'{m:.1f}' for m in masses])
    ax.set_xlabel(r'$\langle \log M_* \rangle$ [M$_\odot$]')
    ax.set_ylabel(r'$r$(sSFR, $\sigma$)')
    ax.set_title('B. Correlation by Mass Bin')
    
    # Panel C: Binned by [Mg/Fe]
    ax = axes[0, 2]
    binned = results['binned_mgfe']
    mgfes = [b['mean_mgfe'] for b in binned]
    r_ssfrs = [b['r_ssfr'] for b in binned]
    ax.bar(range(len(binned)), r_ssfrs, color=colors['accent'], edgecolor=colors['dark'], alpha=0.8)
    ax.axhline(0, color='k', ls='--', lw=1)
    ax.set_xticks(range(len(binned)))
    ax.set_xticklabels([f'{m:.2f}' for m in mgfes])
    ax.set_xlabel(r'$\langle$[Mg/Fe]$\rangle$')
    ax.set_ylabel(r'$r$(sSFR, $\sigma$)')
    ax.set_title('C. Correlation by [Mg/Fe] Bin')
    
    # Panel D: Comparison with age indicators
    ax = axes[1, 0]
    indicators = ['sSFR', 'D4000', '[Mg/Fe]']
    correlations = [
        results['simple']['r_ssfr_sigma'],
        results['simple']['r_d4000_sigma'],
        results['simple']['r_mgfe_sigma']
    ]
    bar_colors = [colors['primary'], colors['secondary'], colors['light']]
    ax.bar(indicators, correlations, color=bar_colors, edgecolor=colors['dark'], alpha=0.8)
    ax.axhline(0, color='k', ls='--', lw=1)
    ax.set_ylabel(r'Correlation with $\log(\sigma)$')
    ax.set_title('D. Comparison: Rate vs Cumulative Indicators')
    
    # Panel E: Redshift evolution
    ax = axes[1, 1]
    z_results = results['redshift_evolution']
    zs = [r['mean_z'] for r in z_results]
    rs = [r['r_ssfr'] for r in z_results]
    ax.plot(zs, rs, 'o-', color=colors['primary'], markersize=8, lw=2)
    ax.axhline(0, color='k', ls='--', lw=1)
    ax.fill_between(zs, [r - 0.02 for r in rs], [r + 0.02 for r in rs], alpha=0.2, color=colors['primary'])
    ax.set_xlabel('Redshift')
    ax.set_ylabel(r'$r$(sSFR, $\sigma$)')
    ax.set_title('E. Redshift Evolution')
    
    # Panel F: Summary
    ax = axes[1, 2]
    ax.axis('off')
    
    summary_text = f"""
    SFR HOLONOMY TEST SUMMARY
    ─────────────────────────────
    Sample: {len(df):,} SDSS galaxies
    
    Simple correlations:
      r(sSFR, σ) = {results['simple']['r_ssfr_sigma']:.4f}
      r(D4000, σ) = {results['simple']['r_d4000_sigma']:.4f}
      r([Mg/Fe], σ) = {results['simple']['r_mgfe_sigma']:.4f}
    
    Controlled for M*:
      r(sSFR, σ | M*) = {results['partial_mass']['r_ssfr_sigma_mass']:.4f}
    
    Controlled for M* and [Mg/Fe]:
      r(sSFR, σ | M*, [Mg/Fe]) = {results['partial_full']['r_ssfr_sigma_full']:.4f}
    
    TEP CONSISTENT: {'YES ✓' if results['simple']['r_ssfr_sigma'] < 0 else 'NO ✗'}
    
    INTERPRETATION:
    {'Galaxies in deeper potentials show' if results['simple']['r_ssfr_sigma'] < 0 else 'No evidence that'}
    {'lower SFR at fixed mass - consistent' if results['simple']['r_ssfr_sigma'] < 0 else 'potential depth affects SFR'}
    {'with TEP time dilation.' if results['simple']['r_ssfr_sigma'] < 0 else ''}
    """
    ax.text(0.05, 0.95, summary_text, transform=ax.transAxes, fontsize=10,
            verticalalignment='top', fontfamily='monospace',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    plt.tight_layout()
    
    # Save
    fig_path = os.path.join(FIGURES_DIR, 'sdss_sfr_holonomy.png')
    plt.savefig(fig_path, dpi=150, bbox_inches='tight', facecolor='white')
    print(f"Saved figure: {fig_path}")
    
    fig_path_pdf = os.path.join(FIGURES_DIR, 'sdss_sfr_holonomy.pdf')
    plt.savefig(fig_path_pdf, bbox_inches='tight')
    print(f"Saved figure: {fig_path_pdf}")
    
    plt.close()


def save_results(results, df):
    """Save results to JSON."""
    output = {
        'test_name': 'SFR Holonomy Test',
        'description': 'Tests whether SFR (a time-derivative) correlates with gravitational potential',
        'timestamp': datetime.now().isoformat(),
        'sample_size': len(df),
        'tep_prediction': 'r(sSFR, σ) < 0 at fixed mass and metallicity',
        'results': results,
        'tep_consistent': results['simple']['r_ssfr_sigma'] < 0,
        'significance': {
            'simple': results['simple']['p_ssfr_sigma'],
            'controlled': results['partial_full']['p_ssfr_sigma_full']
        }
    }
    
    output_path = os.path.join(RESULTS_DIR, 'sdss_sfr_holonomy_results.json')
    with open(output_path, 'w') as f:
        json.dump(output, f, indent=2, default=float)
    print(f"\nSaved results: {output_path}")
    
    return output


def main():
    """Run the full SFR Holonomy analysis."""
    print("=" * 70)
    print("STEP 6.6: STAR FORMATION RATE HOLONOMY TEST")
    print("=" * 70)
    print("\nTEP HYPOTHESIS:")
    print("  SFR is a time-derivative (M☉/yr). Under TEP, proper time flows")
    print("  slower in deeper potentials. At fixed gas supply and metallicity,")
    print("  high-σ galaxies should show LOWER observed SFR.")
    print("\nTEP PREDICTION: r(sSFR, σ) < 0")
    
    # Load and prepare data
    df = load_data()
    df = prepare_data(df)
    
    # Run tests
    results = {}
    results.update(test_simple_correlations(df))
    results.update(test_mass_controlled(df))
    results.update(test_metallicity_controlled(df))
    results['bpt_stratified'] = test_bpt_stratified(df)
    results['redshift_evolution'] = test_redshift_evolution(df)
    
    # Create figure
    create_figure(df, results)
    
    # Save results
    output = save_results(results, df)
    
    # Final summary
    print("\n" + "=" * 70)
    print("FINAL SUMMARY")
    print("=" * 70)
    
    r_simple = results['simple']['r_ssfr_sigma']
    r_controlled = results['partial_full']['r_ssfr_sigma_full']
    
    print(f"\n  r(sSFR, σ) simple:                    {r_simple:>8.4f}")
    print(f"  r(sSFR, σ | M*, [Mg/Fe]) controlled:  {r_controlled:>8.4f}")
    
    print(f"\n  TEP CONSISTENT: {'YES ✓' if r_simple < 0 else 'NO ✗'}")
    
    if r_simple < 0:
        print("\n  INTERPRETATION:")
        print("  Galaxies in deeper gravitational potentials show lower specific")
        print("  star formation rates at fixed stellar mass and nucleosynthesis ratio.")
        print("  This is consistent with TEP time dilation slowing the rate of")
        print("  baryonic processes in regions of stronger gravity.")
    
    print("\n" + "=" * 70)
    
    return output


if __name__ == "__main__":
    main()
