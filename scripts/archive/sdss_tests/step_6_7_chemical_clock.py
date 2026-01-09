#!/usr/bin/env python3
"""
Step 6.7: Chemical Clock vs Spectroscopic Clock Test

SMOKING GUN TEST #3: Two independent time standards in the same objects

TEP HYPOTHESIS:
- [Mg/Fe] ratio is a NUCLEOSYNTHESIS clock (Type Ia delay timescale)
  → Measures coordinate time (unaffected by local time dilation)
- D4000 and Hβ are SPECTROSCOPIC age indicators (stellar evolution)
  → Measure proper time (affected by local time dilation)

Under standard physics: At fixed [Mg/Fe], D4000 should track age consistently.
Under TEP: At fixed [Mg/Fe], high-σ galaxies experienced LESS proper time,
           so they should appear YOUNGER (lower D4000, higher Hβ).

TEP PREDICTION:
  At fixed [Mg/Fe] and M*:
    r(D4000, σ) < 0   (high-σ appears younger)
    r(Hβ, σ) > 0      (high-σ has stronger Hβ = younger)

This is a SIGN TEST - standard physics predicts the OPPOSITE or null.

Author: M. Smawfield
Date: January 2026
"""

import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import pearsonr, spearmanr
from sklearn.linear_model import LinearRegression
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
    """Load SDSS spectral indices data."""
    df = pd.read_csv(os.path.join(DATA_DIR, 'sdss_spectral_indices.csv'))
    print(f"Loaded {len(df):,} galaxies")
    return df


def prepare_data(df):
    """Prepare derived quantities."""
    df['log_sigma'] = np.log10(df['veldisp'])
    df['fe_avg'] = (df['fe5270'] + df['fe5335']) / 2
    df['mg_fe_ratio'] = df['mgb'] / df['fe_avg']
    df['log_mg_fe'] = np.log10(df['mg_fe_ratio'])
    
    # Quality cuts
    valid = (
        (df['veldisp'] > 50) & (df['veldisp'] < 400) &
        (df['veldisp_err'] < 30) &
        np.isfinite(df['d4000']) & (df['d4000'] > 1.0) & (df['d4000'] < 2.5) &
        np.isfinite(df['hbeta']) & (df['hbeta'] > 0) & (df['hbeta'] < 6) &
        np.isfinite(df['log_mass']) & (df['log_mass'] > 9.0) & (df['log_mass'] < 12.5) &
        np.isfinite(df['mg_fe_ratio']) & (df['mg_fe_ratio'] > 0.3) & (df['mg_fe_ratio'] < 5.0) &
        (df['mgb_err'] < 0.5) & (df['fe5270_err'] < 0.5) &
        (df['d4000_err'] < 0.1) & (df['hbeta_err'] < 0.5) &
        (df['redshift'] > 0.02) & (df['redshift'] < 0.25)
    )
    
    df_clean = df[valid].copy()
    print(f"After quality cuts: {len(df_clean):,} galaxies ({100*len(df_clean)/len(df):.1f}%)")
    return df_clean


def partial_corr(x, y, Z):
    """Partial correlation controlling for multiple variables."""
    Z = np.atleast_2d(Z).T if Z.ndim == 1 else Z
    reg_x = LinearRegression().fit(Z, x)
    reg_y = LinearRegression().fit(Z, y)
    resid_x = x - reg_x.predict(Z)
    resid_y = y - reg_y.predict(Z)
    r, p = pearsonr(resid_x, resid_y)
    return r, p, resid_x, resid_y


def test_simple_correlations(df):
    """Test 1: Simple correlations."""
    print("\n" + "=" * 70)
    print("TEST 1: Simple Correlations (Age Indicators vs σ)")
    print("=" * 70)
    
    r_d4000, p_d4000 = pearsonr(df['log_sigma'], df['d4000'])
    r_hbeta, p_hbeta = pearsonr(df['log_sigma'], df['hbeta'])
    r_mgfe, p_mgfe = pearsonr(df['log_sigma'], df['log_mg_fe'])
    
    print(f"\nr(D4000, σ):   {r_d4000:+.4f}  (p = {p_d4000:.2e})")
    print(f"r(Hβ, σ):      {r_hbeta:+.4f}  (p = {p_hbeta:.2e})")
    print(f"r([Mg/Fe], σ): {r_mgfe:+.4f}  (p = {p_mgfe:.2e})")
    
    print("\nInterpretation:")
    print("  Standard: High-σ = older = higher D4000, lower Hβ, higher [Mg/Fe]")
    print("  All correlations follow standard expectation (so far)")
    
    return {
        'r_d4000_sigma': r_d4000,
        'r_hbeta_sigma': r_hbeta,
        'r_mgfe_sigma': r_mgfe
    }


def test_chemical_clock_controlled(df):
    """
    Test 2: THE SMOKING GUN
    
    At FIXED [Mg/Fe] (nucleosynthesis clock), what is the correlation
    between D4000 (spectroscopic age) and σ (potential depth)?
    
    Standard: r ≈ 0 or slightly positive (no residual after age control)
    TEP: r < 0 (high-σ appears younger than [Mg/Fe] predicts)
    """
    print("\n" + "=" * 70)
    print("TEST 2: Chemical Clock Discrepancy (THE SMOKING GUN)")
    print("=" * 70)
    
    print("\nKey insight: [Mg/Fe] is set by Type Ia delay timescale")
    print("             D4000/Hβ are set by stellar evolution")
    print("             If time flows differently, these clocks disagree")
    
    # Control for [Mg/Fe] only
    Z_mgfe = df['log_mg_fe'].values
    r_d4000_mgfe, p_d4000_mgfe, _, _ = partial_corr(
        df['log_sigma'].values, df['d4000'].values, Z_mgfe
    )
    r_hbeta_mgfe, p_hbeta_mgfe, _, _ = partial_corr(
        df['log_sigma'].values, df['hbeta'].values, Z_mgfe
    )
    
    print(f"\nControlling for [Mg/Fe] only:")
    print(f"  r(D4000, σ | [Mg/Fe]) = {r_d4000_mgfe:+.4f}  (p = {p_d4000_mgfe:.2e})")
    print(f"  r(Hβ, σ | [Mg/Fe])    = {r_hbeta_mgfe:+.4f}  (p = {p_hbeta_mgfe:.2e})")
    
    # Control for [Mg/Fe] AND M*
    Z_full = np.column_stack([df['log_mg_fe'].values, df['log_mass'].values])
    r_d4000_full, p_d4000_full, resid_sigma, resid_d4000 = partial_corr(
        df['log_sigma'].values, df['d4000'].values, Z_full
    )
    r_hbeta_full, p_hbeta_full, _, resid_hbeta = partial_corr(
        df['log_sigma'].values, df['hbeta'].values, Z_full
    )
    
    print(f"\nControlling for [Mg/Fe] AND M*:")
    print(f"  r(D4000, σ | [Mg/Fe], M*) = {r_d4000_full:+.4f}  (p = {p_d4000_full:.2e})")
    print(f"  r(Hβ, σ | [Mg/Fe], M*)    = {r_hbeta_full:+.4f}  (p = {p_hbeta_full:.2e})")
    
    # Also control for redshift
    Z_fullz = np.column_stack([df['log_mg_fe'].values, df['log_mass'].values, df['redshift'].values])
    r_d4000_z, p_d4000_z, _, _ = partial_corr(
        df['log_sigma'].values, df['d4000'].values, Z_fullz
    )
    r_hbeta_z, p_hbeta_z, _, _ = partial_corr(
        df['log_sigma'].values, df['hbeta'].values, Z_fullz
    )
    
    print(f"\nControlling for [Mg/Fe], M*, AND z:")
    print(f"  r(D4000, σ | [Mg/Fe], M*, z) = {r_d4000_z:+.4f}  (p = {p_d4000_z:.2e})")
    print(f"  r(Hβ, σ | [Mg/Fe], M*, z)    = {r_hbeta_z:+.4f}  (p = {p_hbeta_z:.2e})")
    
    # VERDICT
    print("\n" + "-" * 50)
    print("SMOKING GUN VERDICT:")
    if r_d4000_full < -0.05:
        print(f"  ✓ D4000-σ correlation is NEGATIVE ({r_d4000_full:+.4f})")
        print("    At fixed [Mg/Fe], high-σ galaxies appear YOUNGER")
        print("    This is OPPOSITE to standard expectation!")
        tep_consistent = True
    elif r_d4000_full > 0.05:
        print(f"  ✗ D4000-σ correlation is POSITIVE ({r_d4000_full:+.4f})")
        print("    Standard physics: high-σ = older even at fixed [Mg/Fe]")
        tep_consistent = False
    else:
        print(f"  ~ D4000-σ correlation is WEAK ({r_d4000_full:+.4f})")
        print("    Inconclusive - need more controls or different approach")
        tep_consistent = None
    
    return {
        'r_d4000_mgfe': r_d4000_mgfe,
        'r_hbeta_mgfe': r_hbeta_mgfe,
        'r_d4000_full': r_d4000_full,
        'r_hbeta_full': r_hbeta_full,
        'r_d4000_z': r_d4000_z,
        'r_hbeta_z': r_hbeta_z,
        'tep_consistent': tep_consistent,
        'resid_sigma': resid_sigma,
        'resid_d4000': resid_d4000,
        'resid_hbeta': resid_hbeta
    }


def test_binned_analysis(df):
    """Test 3: Binned analysis within [Mg/Fe] bins."""
    print("\n" + "=" * 70)
    print("TEST 3: Binned Analysis (within [Mg/Fe] bins)")
    print("=" * 70)
    
    df['mgfe_bin'] = pd.qcut(df['log_mg_fe'], q=5, labels=['Z1','Z2','Z3','Z4','Z5'])
    
    print(f"\n{'[Mg/Fe] bin':<12} {'<[Mg/Fe]>':>10} {'N':>8} {'r(D4000,σ)':>12} {'r(Hβ,σ)':>10}")
    print("-" * 60)
    
    results = []
    for zbin in ['Z1', 'Z2', 'Z3', 'Z4', 'Z5']:
        subset = df[df['mgfe_bin'] == zbin]
        r_d4000, _ = pearsonr(subset['log_sigma'], subset['d4000'])
        r_hbeta, _ = pearsonr(subset['log_sigma'], subset['hbeta'])
        mean_mgfe = subset['log_mg_fe'].mean()
        
        print(f"{zbin:<12} {mean_mgfe:>10.4f} {len(subset):>8} {r_d4000:>+12.4f} {r_hbeta:>+10.4f}")
        results.append({
            'bin': zbin, 
            'mean_mgfe': mean_mgfe, 
            'n': len(subset),
            'r_d4000': r_d4000, 
            'r_hbeta': r_hbeta
        })
    
    # Mean across bins
    mean_r_d4000 = np.mean([r['r_d4000'] for r in results])
    mean_r_hbeta = np.mean([r['r_hbeta'] for r in results])
    print("-" * 60)
    print(f"{'Mean':<12} {'':>10} {'':>8} {mean_r_d4000:>+12.4f} {mean_r_hbeta:>+10.4f}")
    
    return results


def test_matched_pairs(df):
    """Test 4: Matched-pair analysis."""
    print("\n" + "=" * 70)
    print("TEST 4: Matched-Pair Analysis")
    print("=" * 70)
    
    # Create matched cells
    df['mgfe_cell'] = pd.qcut(df['log_mg_fe'], q=10, labels=False, duplicates='drop')
    df['mass_cell'] = pd.qcut(df['log_mass'], q=10, labels=False, duplicates='drop')
    df['z_cell'] = pd.qcut(df['redshift'], q=5, labels=False, duplicates='drop')
    df['match_cell'] = df['mgfe_cell'].astype(str) + '_' + df['mass_cell'].astype(str) + '_' + df['z_cell'].astype(str)
    
    cell_results = []
    for cell, group in df.groupby('match_cell'):
        if len(group) >= 30:
            r, p = pearsonr(group['log_sigma'], group['d4000'])
            cell_results.append({'cell': cell, 'n': len(group), 'r': r, 'p': p})
    
    print(f"\nMatched-pair cells: {len(cell_results)}")
    
    rs = [c['r'] for c in cell_results]
    ns = [c['n'] for c in cell_results]
    weighted_r = np.average(rs, weights=ns)
    
    positive = sum(1 for c in cell_results if c['r'] > 0)
    negative = sum(1 for c in cell_results if c['r'] < 0)
    
    print(f"Weighted mean r(D4000, σ): {weighted_r:+.4f}")
    print(f"Cells with r > 0: {positive} ({100*positive/len(cell_results):.1f}%)")
    print(f"Cells with r < 0: {negative} ({100*negative/len(cell_results):.1f}%)")
    
    # Sign test
    n_total = positive + negative
    p_sign = 2 * min(
        stats.binom.cdf(min(positive, negative), n_total, 0.5),
        1 - stats.binom.cdf(max(positive, negative) - 1, n_total, 0.5)
    )
    print(f"Sign test p-value: {p_sign:.4f}")
    
    if weighted_r < 0 and negative > positive:
        print("\n✓ CONSISTENT: Majority of cells show negative D4000-σ correlation")
    else:
        print("\n✗ INCONSISTENT or ambiguous")
    
    return {
        'n_cells': len(cell_results),
        'weighted_r': weighted_r,
        'n_positive': positive,
        'n_negative': negative,
        'p_sign': p_sign
    }


def create_figure(df, results):
    """Create publication figure."""
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # Panel A: D4000 vs σ colored by [Mg/Fe]
    ax = axes[0, 0]
    scatter = ax.scatter(df['log_sigma'], df['d4000'], 
                        c=df['log_mg_fe'], cmap='coolwarm', 
                        s=1, alpha=0.3, rasterized=True)
    ax.set_xlabel('log(σ) [km/s]')
    ax.set_ylabel('D4000')
    ax.set_title('A. D4000 vs σ (colored by [Mg/Fe])')
    plt.colorbar(scatter, ax=ax, label='log([Mg/Fe])')
    
    # Panel B: Residual D4000 vs residual σ
    ax = axes[0, 1]
    if 'resid_sigma' in results['controlled']:
        resid_sigma = results['controlled']['resid_sigma']
        resid_d4000 = results['controlled']['resid_d4000']
        ax.scatter(resid_sigma, resid_d4000, s=1, alpha=0.2, c='steelblue', rasterized=True)
        
        # Trend line
        slope, intercept = np.polyfit(resid_sigma, resid_d4000, 1)
        x_line = np.linspace(resid_sigma.min(), resid_sigma.max(), 100)
        ax.plot(x_line, slope * x_line + intercept, 'r-', lw=2)
        
        r = results['controlled']['r_d4000_full']
        ax.text(0.05, 0.95, f'r = {r:+.3f}', transform=ax.transAxes, 
                fontsize=12, fontweight='bold', va='top')
    ax.set_xlabel('σ residual (after [Mg/Fe], M*)')
    ax.set_ylabel('D4000 residual')
    ax.set_title('B. Residual Correlation (THE TEST)')
    ax.axhline(0, c='gray', ls='--', alpha=0.5)
    ax.axvline(0, c='gray', ls='--', alpha=0.5)
    
    # Panel C: Binned analysis
    ax = axes[1, 0]
    binned = results['binned']
    x = range(len(binned))
    r_d4000 = [b['r_d4000'] for b in binned]
    r_hbeta = [b['r_hbeta'] for b in binned]
    width = 0.35
    ax.bar([i - width/2 for i in x], r_d4000, width, label='r(D4000, σ)', color='steelblue')
    ax.bar([i + width/2 for i in x], r_hbeta, width, label='r(Hβ, σ)', color='coral')
    ax.axhline(0, c='black', ls='-', lw=0.5)
    ax.set_xticks(x)
    ax.set_xticklabels([f"Z{i+1}" for i in x])
    ax.set_xlabel('[Mg/Fe] quintile')
    ax.set_ylabel('Correlation with log(σ)')
    ax.set_title('C. Within [Mg/Fe] Bins')
    ax.legend()
    
    # Panel D: Matched-pair distribution
    ax = axes[1, 1]
    matched = results['matched']
    ax.text(0.5, 0.9, f"Matched-pair cells: {matched['n_cells']}", 
            transform=ax.transAxes, ha='center', fontsize=12)
    ax.text(0.5, 0.7, f"Weighted r(D4000, σ) = {matched['weighted_r']:+.4f}", 
            transform=ax.transAxes, ha='center', fontsize=14, fontweight='bold')
    ax.text(0.5, 0.5, f"Cells with r < 0: {matched['n_negative']} ({100*matched['n_negative']/matched['n_cells']:.1f}%)", 
            transform=ax.transAxes, ha='center', fontsize=12)
    ax.text(0.5, 0.3, f"Sign test p = {matched['p_sign']:.4f}", 
            transform=ax.transAxes, ha='center', fontsize=12)
    
    verdict = "TEP CONSISTENT" if matched['weighted_r'] < 0 else "STANDARD PHYSICS"
    color = 'green' if matched['weighted_r'] < 0 else 'red'
    ax.text(0.5, 0.1, verdict, transform=ax.transAxes, ha='center', 
            fontsize=14, fontweight='bold', color=color)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis('off')
    ax.set_title('D. Matched-Pair Summary')
    
    plt.suptitle('Chemical Clock vs Spectroscopic Clock Test\n'
                 'TEP Prediction: At fixed [Mg/Fe], high-σ galaxies appear YOUNGER (lower D4000)',
                 fontsize=14, fontweight='bold')
    plt.tight_layout()
    
    for ext in ['png', 'pdf']:
        path = os.path.join(FIGURES_DIR, f'sdss_chemical_clock.{ext}')
        plt.savefig(path, dpi=150 if ext == 'png' else None, bbox_inches='tight')
        print(f"Saved: {path}")
    plt.close()


def main():
    """Run Chemical Clock test."""
    print("=" * 70)
    print("STEP 6.7: CHEMICAL CLOCK VS SPECTROSCOPIC CLOCK TEST")
    print("=" * 70)
    
    print("\nSMOKING GUN RATIONALE:")
    print("  [Mg/Fe] = nucleosynthesis clock (coordinate time)")
    print("  D4000/Hβ = spectroscopic clock (proper time)")
    print("  If clocks disagree tracking σ, time itself is modified!")
    
    df = load_data()
    df = prepare_data(df)
    
    results = {}
    results['simple'] = test_simple_correlations(df)
    results['controlled'] = test_chemical_clock_controlled(df)
    results['binned'] = test_binned_analysis(df)
    results['matched'] = test_matched_pairs(df)
    
    # Create figure
    print("\n" + "=" * 70)
    print("Creating Figure")
    print("=" * 70)
    create_figure(df, results)
    
    # Final summary
    print("\n" + "=" * 70)
    print("FINAL SUMMARY")
    print("=" * 70)
    
    r_controlled = results['controlled']['r_d4000_full']
    r_matched = results['matched']['weighted_r']
    
    print(f"\n  r(D4000, σ | [Mg/Fe], M*) = {r_controlled:+.4f}")
    print(f"  Matched-pair r(D4000, σ)  = {r_matched:+.4f}")
    
    if r_controlled < 0 and r_matched < 0:
        print("\n  SMOKING GUN VERDICT: TEP CONSISTENT ✓")
        print("  At fixed nucleosynthesis age, high-σ galaxies appear younger")
        print("  This is the OPPOSITE of what standard physics predicts!")
        verdict = "TEP_CONSISTENT"
    elif r_controlled > 0 and r_matched > 0:
        print("\n  SMOKING GUN VERDICT: STANDARD PHYSICS")
        print("  High-σ galaxies appear older even at fixed [Mg/Fe]")
        verdict = "STANDARD"
    else:
        print("\n  SMOKING GUN VERDICT: AMBIGUOUS")
        print("  Results are mixed - need further investigation")
        verdict = "AMBIGUOUS"
    
    # Save results
    output = {
        'test_name': 'Chemical Clock vs Spectroscopic Clock',
        'timestamp': datetime.now().isoformat(),
        'sample_size': len(df),
        'results': {
            'simple': results['simple'],
            'controlled': {k: v for k, v in results['controlled'].items() 
                          if not isinstance(v, np.ndarray)},
            'binned': results['binned'],
            'matched': results['matched']
        },
        'verdict': verdict
    }
    
    output_path = os.path.join(RESULTS_DIR, 'sdss_chemical_clock_results.json')
    with open(output_path, 'w') as f:
        json.dump(output, f, indent=2, default=float)
    print(f"\nSaved: {output_path}")


if __name__ == "__main__":
    main()
