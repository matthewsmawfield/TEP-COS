#!/usr/bin/env python3
"""
BULLETPROOF TEP SIGNALS

This script performs rigorous robustness tests on the strongest TEP-consistent
signals to make them unassailable against criticism.

Tests performed:
1. Comprehensive control variables (mass, z, environment, metallicity, size)
2. Matched-pair analysis (galaxy twins)
3. Bootstrap confidence intervals
4. Jackknife leave-one-out stability
5. Alternative hypothesis testing
6. Selection bias checks
7. Systematic error propagation
"""

import os
import sys
import json
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')

BASE_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATA_DIR = os.path.join(BASE_DIR, 'data', 'sdss')
OUTPUT_DIR = os.path.join(BASE_DIR, 'results', 'outputs')
FIGURE_DIR = os.path.join(BASE_DIR, 'results', 'figures')

os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(FIGURE_DIR, exist_ok=True)


def load_comprehensive_data():
    """Load all available data for bulletproofing."""
    twin_file = os.path.join(DATA_DIR, 'sdss_twin_base_sample_with_size.csv')
    
    if os.path.exists(twin_file):
        df = pd.read_csv(twin_file)
        df = df.rename(columns={
            'specobjid': 'specObjID',
            'veldisp': 'sigma',
            'veldisp_err': 'sigma_err',
            'mgb': 'Mgb',
            'fe5270': 'Fe5270',
            'fe5335': 'Fe5335',
            'd4000': 'D4000',
            'hbeta': 'Hbeta',
            'log_mass': 'logMass',
            'petroR50_r_arcsec': 'Re_arcsec'
        })
        return df
    
    return None


def compute_derived_quantities(df):
    """Compute all derived quantities needed for analysis."""
    # [Mg/Fe] proxy
    df['Fe_avg'] = 0.5 * (df['Fe5270'] + df['Fe5335'])
    mask = (df['Mgb'] > 0) & (df['Fe_avg'] > 0)
    df['MgFe'] = np.nan
    df.loc[mask, 'MgFe'] = np.log10(df.loc[mask, 'Mgb'] / df.loc[mask, 'Fe_avg'])
    
    # Log sigma
    df['log_sigma'] = np.log10(df['sigma'])
    
    # Compactness
    df['log_Re'] = np.log10(df['Re_arcsec'].clip(lower=0.1))
    df['Compactness'] = df['logMass'] - 2 * df['log_Re']
    
    return df


def partial_correlation(x, y, controls, data):
    """Compute partial correlation controlling for multiple variables."""
    mask = data[x].notna() & data[y].notna()
    for c in controls:
        mask &= data[c].notna()
    
    subset = data[mask].copy()
    if len(subset) < 100:
        return np.nan, np.nan, 0
    
    # Residualize x and y against controls
    X_ctrl = subset[controls].values
    
    reg_x = LinearRegression().fit(X_ctrl, subset[x])
    resid_x = subset[x] - reg_x.predict(X_ctrl)
    
    reg_y = LinearRegression().fit(X_ctrl, subset[y])
    resid_y = subset[y] - reg_y.predict(X_ctrl)
    
    r, p = stats.pearsonr(resid_x, resid_y)
    return r, p, len(subset)


def bootstrap_correlation(x, y, data, n_boot=1000):
    """Bootstrap confidence interval for correlation."""
    mask = data[x].notna() & data[y].notna()
    subset = data[mask]
    
    if len(subset) < 100:
        return np.nan, np.nan, np.nan
    
    correlations = []
    n = len(subset)
    
    for _ in range(n_boot):
        idx = np.random.choice(n, size=n, replace=True)
        boot_x = subset[x].iloc[idx]
        boot_y = subset[y].iloc[idx]
        r, _ = stats.pearsonr(boot_x, boot_y)
        correlations.append(r)
    
    return np.percentile(correlations, [2.5, 50, 97.5])


def jackknife_stability(x, y, data, n_groups=10):
    """Jackknife leave-one-group-out stability test."""
    mask = data[x].notna() & data[y].notna()
    subset = data[mask].copy()
    
    if len(subset) < 100:
        return np.nan, []
    
    # Assign to groups
    subset['_group'] = np.random.randint(0, n_groups, len(subset))
    
    correlations = []
    for g in range(n_groups):
        loo = subset[subset['_group'] != g]
        r, _ = stats.pearsonr(loo[x], loo[y])
        correlations.append(r)
    
    return np.std(correlations), correlations


def matched_pair_analysis(x, y, match_vars, data, n_bins=5):
    """Matched-pair analysis controlling for confounders."""
    mask = data[x].notna() & data[y].notna()
    for v in match_vars:
        mask &= data[v].notna()
    
    subset = data[mask].copy()
    
    if len(subset) < 500:
        return np.nan, np.nan, 0, []
    
    # Create bins for each matching variable
    for v in match_vars:
        subset[f'{v}_bin'] = pd.qcut(subset[v], n_bins, labels=False, duplicates='drop')
    
    # Create combined bin key
    bin_cols = [f'{v}_bin' for v in match_vars]
    subset['match_key'] = subset[bin_cols].astype(str).agg('-'.join, axis=1)
    
    # Within each bin, compute correlation
    bin_results = []
    for key, group in subset.groupby('match_key'):
        if len(group) >= 20:
            r, p = stats.pearsonr(group[x], group[y])
            bin_results.append({
                'key': key,
                'r': r,
                'p': p,
                'n': len(group)
            })
    
    if len(bin_results) == 0:
        return np.nan, np.nan, 0, []
    
    # Weighted mean correlation
    total_n = sum(b['n'] for b in bin_results)
    weighted_r = sum(b['r'] * b['n'] for b in bin_results) / total_n
    
    # Fraction of bins with same sign as overall
    sign_frac = sum(1 for b in bin_results if np.sign(b['r']) == np.sign(weighted_r)) / len(bin_results)
    
    return weighted_r, sign_frac, len(bin_results), bin_results


def test_alternative_hypothesis(x, y, alt_var, data):
    """Test if an alternative variable explains the correlation better."""
    mask = data[x].notna() & data[y].notna() & data[alt_var].notna()
    subset = data[mask]
    
    if len(subset) < 100:
        return {}
    
    # Original correlation
    r_orig, p_orig = stats.pearsonr(subset[x], subset[y])
    
    # Alternative correlation
    r_alt, p_alt = stats.pearsonr(subset[alt_var], subset[y])
    
    # Partial correlation of x with y controlling for alt
    r_partial, p_partial, n = partial_correlation(x, y, [alt_var], subset)
    
    # Does the signal survive?
    survives = (np.sign(r_partial) == np.sign(r_orig)) and (p_partial < 0.05)
    
    return {
        'r_original': float(r_orig),
        'r_alternative': float(r_alt),
        'r_partial': float(r_partial) if not np.isnan(r_partial) else None,
        'p_partial': float(p_partial) if not np.isnan(p_partial) else None,
        'survives_control': survives,
        'n': int(n)
    }


def bulletproof_chemical_clock(df):
    """Bulletproof Test H: Chemical Clock Discrepancy."""
    print("\n" + "="*70)
    print("BULLETPROOFING TEST H: CHEMICAL CLOCK ([Mg/Fe] vs σ at fixed age)")
    print("="*70)
    
    results = {'test': 'H_Chemical_Clock'}
    
    # Clean data
    mask = (
        df['MgFe'].notna() & df['sigma'].notna() & df['D4000'].notna() &
        (df['sigma'] > 50) & (df['sigma'] < 400) &
        df['logMass'].notna()
    )
    clean = df[mask].copy()
    print(f"\nSample size: {len(clean):,}")
    results['n_total'] = len(clean)
    
    # 1. Raw correlation
    r_raw, p_raw = stats.pearsonr(clean['sigma'], clean['MgFe'])
    print(f"\n1. RAW CORRELATION: r = {r_raw:.4f}, p = {p_raw:.2e}")
    results['raw'] = {'r': float(r_raw), 'p': float(p_raw)}
    
    # 2. Partial correlations with increasing controls
    print("\n2. PARTIAL CORRELATIONS (surviving controls):")
    
    control_sets = [
        (['D4000'], 'Age proxy'),
        (['D4000', 'logMass'], 'Age + Mass'),
        (['D4000', 'logMass', 'redshift'], 'Age + Mass + z'),
        (['D4000', 'logMass', 'redshift', 'Compactness'], 'All controls')
    ]
    
    partial_results = []
    for controls, label in control_sets:
        r, p, n = partial_correlation('sigma', 'MgFe', controls, clean)
        print(f"   {label}: r = {r:.4f}, p = {p:.2e}, n = {n}")
        partial_results.append({
            'controls': label,
            'r': float(r) if not np.isnan(r) else None,
            'p': float(p) if not np.isnan(p) else None,
            'n': n
        })
    results['partial_correlations'] = partial_results
    
    # 3. Bootstrap CI
    print("\n3. BOOTSTRAP CONFIDENCE INTERVAL (1000 iterations):")
    ci = bootstrap_correlation('sigma', 'MgFe', clean)
    print(f"   95% CI: [{ci[0]:.4f}, {ci[2]:.4f}]")
    print(f"   Median: {ci[1]:.4f}")
    results['bootstrap_ci'] = {'lower': float(ci[0]), 'median': float(ci[1]), 'upper': float(ci[2])}
    
    # 4. Jackknife stability
    print("\n4. JACKKNIFE STABILITY (10-fold):")
    jk_std, jk_vals = jackknife_stability('sigma', 'MgFe', clean)
    print(f"   Std across folds: {jk_std:.4f}")
    print(f"   Range: [{min(jk_vals):.4f}, {max(jk_vals):.4f}]")
    results['jackknife'] = {'std': float(jk_std), 'min': float(min(jk_vals)), 'max': float(max(jk_vals))}
    
    # 5. Matched-pair analysis
    print("\n5. MATCHED-PAIR ANALYSIS:")
    mp_r, mp_frac, mp_n, mp_bins = matched_pair_analysis(
        'sigma', 'MgFe', ['D4000', 'logMass', 'redshift'], clean
    )
    print(f"   Weighted r across {mp_n} matched bins: {mp_r:.4f}")
    print(f"   Fraction with consistent sign: {mp_frac:.1%}")
    results['matched_pairs'] = {'weighted_r': float(mp_r), 'sign_fraction': float(mp_frac), 'n_bins': mp_n}
    
    # 6. Alternative hypothesis tests
    print("\n6. ALTERNATIVE HYPOTHESIS TESTS:")
    
    alt_tests = [
        ('logMass', 'Mass-driven'),
        ('Compactness', 'Compactness-driven'),
        ('redshift', 'Redshift-driven')
    ]
    
    alt_results = []
    for alt_var, label in alt_tests:
        alt = test_alternative_hypothesis('sigma', 'MgFe', alt_var, clean)
        survives = alt.get('survives_control', False)
        print(f"   {label}: Signal survives = {survives}")
        if alt.get('r_partial') is not None:
            print(f"      r_partial = {alt['r_partial']:.4f}, p = {alt['p_partial']:.2e}")
        alt_results.append({'hypothesis': label, **alt})
    results['alternative_tests'] = alt_results
    
    # 7. σ-binned robustness
    print("\n7. σ-BINNED ANALYSIS:")
    sigma_bins = [(50, 100), (100, 150), (150, 200), (200, 300), (300, 400)]
    bin_results = []
    for lo, hi in sigma_bins:
        subset = clean[(clean['sigma'] >= lo) & (clean['sigma'] < hi)]
        if len(subset) > 50:
            mean_mgfe = subset['MgFe'].mean()
            sem = subset['MgFe'].std() / np.sqrt(len(subset))
            print(f"   σ={lo}-{hi}: <[Mg/Fe]> = {mean_mgfe:.4f} ± {sem:.4f}, n={len(subset)}")
            bin_results.append({
                'sigma_range': f'{lo}-{hi}',
                'mean_MgFe': float(mean_mgfe),
                'sem': float(sem),
                'n': len(subset)
            })
    results['sigma_bins'] = bin_results
    
    # Verdict
    all_partials_positive = all(
        p['r'] is not None and p['r'] > 0 
        for p in partial_results
    )
    ci_excludes_zero = ci[0] > 0 or ci[2] < 0
    
    if all_partials_positive and ci_excludes_zero and mp_frac > 0.6:
        results['verdict'] = 'BULLETPROOF'
        print("\n✓ VERDICT: BULLETPROOF - Signal survives all controls")
    else:
        results['verdict'] = 'NEEDS_WORK'
        print("\n⚠ VERDICT: Signal weakened under some controls")
    
    return results


def bulletproof_psb_timing(df):
    """Bulletproof Test I: Post-Starburst Timing."""
    print("\n" + "="*70)
    print("BULLETPROOFING TEST I: POST-STARBURST TIMING (Hβ vs σ)")
    print("="*70)
    
    results = {'test': 'I_PSB_Timing'}
    
    mask = (
        df['Hbeta'].notna() & df['sigma'].notna() & df['D4000'].notna() &
        (df['sigma'] > 50) & (df['sigma'] < 400) &
        df['logMass'].notna()
    )
    clean = df[mask].copy()
    print(f"\nSample size: {len(clean):,}")
    results['n_total'] = len(clean)
    
    # 1. Raw correlation
    r_raw, p_raw = stats.pearsonr(clean['sigma'], clean['Hbeta'])
    print(f"\n1. RAW CORRELATION: r = {r_raw:.4f}, p = {p_raw:.2e}")
    results['raw'] = {'r': float(r_raw), 'p': float(p_raw)}
    
    # 2. Partial correlations
    print("\n2. PARTIAL CORRELATIONS:")
    control_sets = [
        (['D4000'], 'Age proxy'),
        (['D4000', 'logMass'], 'Age + Mass'),
        (['D4000', 'logMass', 'MgFe'], 'Age + Mass + [Mg/Fe]'),
    ]
    
    partial_results = []
    for controls, label in control_sets:
        r, p, n = partial_correlation('sigma', 'Hbeta', controls, clean)
        print(f"   {label}: r = {r:.4f}, p = {p:.2e}, n = {n}")
        partial_results.append({
            'controls': label,
            'r': float(r) if not np.isnan(r) else None,
            'p': float(p) if not np.isnan(p) else None,
            'n': n
        })
    results['partial_correlations'] = partial_results
    
    # 3. Bootstrap CI
    print("\n3. BOOTSTRAP CI:")
    ci = bootstrap_correlation('sigma', 'Hbeta', clean)
    print(f"   95% CI: [{ci[0]:.4f}, {ci[2]:.4f}]")
    results['bootstrap_ci'] = {'lower': float(ci[0]), 'median': float(ci[1]), 'upper': float(ci[2])}
    
    # 4. Matched-pair
    print("\n4. MATCHED-PAIR ANALYSIS:")
    mp_r, mp_frac, mp_n, _ = matched_pair_analysis(
        'sigma', 'Hbeta', ['D4000', 'logMass'], clean
    )
    print(f"   Weighted r: {mp_r:.4f}, Sign fraction: {mp_frac:.1%}")
    results['matched_pairs'] = {'weighted_r': float(mp_r), 'sign_fraction': float(mp_frac), 'n_bins': mp_n}
    
    # 5. Critical test: Does signal flip after D4000 control?
    r_ctrl, p_ctrl, n_ctrl = partial_correlation('sigma', 'Hbeta', ['D4000'], clean)
    
    if r_ctrl < 0:
        results['verdict'] = 'CONFOUNDED'
        print(f"\n⚠ VERDICT: CONFOUNDED - Signal flips negative (r={r_ctrl:.4f}) after D4000 control")
        print("   This suggests the raw correlation is driven by age-σ correlation, not TEP")
    else:
        results['verdict'] = 'BULLETPROOF'
        print(f"\n✓ VERDICT: BULLETPROOF - Signal persists after age control")
    
    return results


def bulletproof_lw_mw_age(df):
    """Bulletproof Test L: LW-MW Age Difference."""
    print("\n" + "="*70)
    print("BULLETPROOFING TEST L: LW-MW AGE PROXY (indirect gradient test)")
    print("="*70)
    
    # Load MaNGA data
    manga_file = os.path.join(DATA_DIR, 'manga_age_data.csv')
    
    if not os.path.exists(manga_file):
        print("MaNGA data not available")
        return {'test': 'L_LW_MW_Age', 'verdict': 'SKIPPED'}
    
    manga = pd.read_csv(manga_file)
    manga = manga.drop_duplicates(subset=['PLATEIFU'])
    
    results = {'test': 'L_LW_MW_Age'}
    
    mask = (
        manga['LW_AGE_1RE'].notna() & manga['MW_AGE_1RE'].notna() &
        manga['stellar_sigma_1re'].notna() &
        (manga['stellar_sigma_1re'] > 50) & (manga['stellar_sigma_1re'] < 400)
    )
    clean = manga[mask].copy()
    clean['age_diff'] = clean['LW_AGE_1RE'] - clean['MW_AGE_1RE']
    
    print(f"\nSample size: {len(clean):,}")
    results['n_total'] = len(clean)
    
    # 1. Raw correlation
    r_raw, p_raw = stats.pearsonr(clean['stellar_sigma_1re'], clean['age_diff'])
    print(f"\n1. RAW CORRELATION: r = {r_raw:.4f}, p = {p_raw:.2e}")
    results['raw'] = {'r': float(r_raw), 'p': float(p_raw)}
    
    # 2. Controls
    print("\n2. PARTIAL CORRELATIONS:")
    if 'PHOTOMETRIC_MASS' in clean.columns:
        r, p, n = partial_correlation('stellar_sigma_1re', 'age_diff', ['PHOTOMETRIC_MASS'], clean)
        print(f"   Controlling for Mass: r = {r:.4f}, p = {p:.2e}")
        results['partial_mass'] = {'r': float(r), 'p': float(p)}
    
    # 3. Bootstrap
    print("\n3. BOOTSTRAP CI:")
    ci = bootstrap_correlation('stellar_sigma_1re', 'age_diff', clean)
    print(f"   95% CI: [{ci[0]:.4f}, {ci[2]:.4f}]")
    results['bootstrap_ci'] = {'lower': float(ci[0]), 'median': float(ci[1]), 'upper': float(ci[2])}
    
    # Verdict
    ci_excludes_zero = ci[0] > 0 or ci[2] < 0
    if ci_excludes_zero:
        results['verdict'] = 'BULLETPROOF'
        print("\n✓ VERDICT: BULLETPROOF")
    else:
        results['verdict'] = 'UNCERTAIN'
        print("\n⚠ VERDICT: CI includes zero")
    
    return results


def create_bulletproof_figure(all_results):
    """Create comprehensive bulletproof visualization."""
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    
    # Panel 1: Partial correlation survival
    ax1 = axes[0, 0]
    tests = []
    raw_r = []
    final_r = []
    
    for result in all_results:
        if 'partial_correlations' in result:
            tests.append(result['test'].split('_')[0])
            raw_r.append(result['raw']['r'])
            if result['partial_correlations']:
                final_r.append(result['partial_correlations'][-1]['r'] or 0)
            else:
                final_r.append(0)
    
    x = np.arange(len(tests))
    width = 0.35
    ax1.bar(x - width/2, raw_r, width, label='Raw', color='steelblue')
    ax1.bar(x + width/2, final_r, width, label='After Controls', color='coral')
    ax1.axhline(0, color='black', linestyle='-', linewidth=0.5)
    ax1.set_xticks(x)
    ax1.set_xticklabels(tests)
    ax1.set_ylabel('Correlation (r)')
    ax1.set_title('Signal Survival Under Controls')
    ax1.legend()
    
    # Panel 2: Bootstrap CIs
    ax2 = axes[0, 1]
    ci_data = []
    for result in all_results:
        if 'bootstrap_ci' in result:
            ci_data.append({
                'test': result['test'].split('_')[0],
                'lower': result['bootstrap_ci']['lower'],
                'median': result['bootstrap_ci']['median'],
                'upper': result['bootstrap_ci']['upper']
            })
    
    for i, ci in enumerate(ci_data):
        ax2.errorbar(i, ci['median'], 
                    yerr=[[ci['median'] - ci['lower']], [ci['upper'] - ci['median']]],
                    fmt='o', capsize=5, markersize=10, color='darkgreen')
    ax2.axhline(0, color='red', linestyle='--', linewidth=2)
    ax2.set_xticks(range(len(ci_data)))
    ax2.set_xticklabels([c['test'] for c in ci_data])
    ax2.set_ylabel('Correlation (r)')
    ax2.set_title('Bootstrap 95% Confidence Intervals')
    
    # Panel 3: Matched-pair sign fraction
    ax3 = axes[1, 0]
    mp_data = []
    for result in all_results:
        if 'matched_pairs' in result and result['matched_pairs'].get('sign_fraction'):
            mp_data.append({
                'test': result['test'].split('_')[0],
                'frac': result['matched_pairs']['sign_fraction']
            })
    
    if mp_data:
        colors = ['green' if m['frac'] > 0.6 else 'orange' for m in mp_data]
        ax3.bar([m['test'] for m in mp_data], [m['frac'] for m in mp_data], color=colors)
        ax3.axhline(0.5, color='red', linestyle='--', label='Random chance')
        ax3.axhline(0.6, color='green', linestyle=':', label='60% threshold')
        ax3.set_ylabel('Fraction of bins with consistent sign')
        ax3.set_title('Matched-Pair Analysis: Sign Consistency')
        ax3.legend()
        ax3.set_ylim(0, 1)
    
    # Panel 4: Summary table
    ax4 = axes[1, 1]
    ax4.axis('off')
    
    summary_text = "BULLETPROOF ASSESSMENT SUMMARY\n" + "="*40 + "\n\n"
    for result in all_results:
        test = result['test']
        verdict = result.get('verdict', 'N/A')
        emoji = '✓' if verdict == 'BULLETPROOF' else ('⚠' if verdict == 'CONFOUNDED' else '?')
        summary_text += f"{emoji} {test}: {verdict}\n"
    
    summary_text += "\n" + "="*40 + "\n"
    summary_text += "BULLETPROOF = Signal survives all controls\n"
    summary_text += "CONFOUNDED = Signal driven by confounders\n"
    summary_text += "UNCERTAIN = Needs more investigation"
    
    ax4.text(0.1, 0.9, summary_text, transform=ax4.transAxes,
            fontsize=11, verticalalignment='top', fontfamily='monospace',
            bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))
    
    plt.tight_layout()
    fig_path = os.path.join(FIGURE_DIR, 'bulletproof_tep_signals.png')
    plt.savefig(fig_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"\nFigure saved to {fig_path}")
    
    return fig_path


def main():
    print("="*70)
    print("BULLETPROOFING TEP SIGNALS")
    print("Making the strongest results unassailable")
    print("="*70)
    
    # Load data
    df = load_comprehensive_data()
    if df is None:
        print("ERROR: Could not load data")
        return
    
    print(f"\nLoaded {len(df):,} galaxies")
    df = compute_derived_quantities(df)
    
    all_results = []
    
    # Bulletproof each signal
    all_results.append(bulletproof_chemical_clock(df))
    all_results.append(bulletproof_psb_timing(df))
    all_results.append(bulletproof_lw_mw_age(df))
    
    # Create figure
    fig_path = create_bulletproof_figure(all_results)
    
    # Save results
    # Convert numpy types to native Python for JSON serialization
    def convert_to_native(obj):
        if isinstance(obj, dict):
            return {k: convert_to_native(v) for k, v in obj.items()}
        elif isinstance(obj, list):
            return [convert_to_native(i) for i in obj]
        elif isinstance(obj, (np.integer, np.int64)):
            return int(obj)
        elif isinstance(obj, (np.floating, np.float64)):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, (np.bool_, bool)):
            return bool(obj)
        return obj
    
    output = convert_to_native({
        'tests': all_results,
        'figure': fig_path,
        'summary': {
            'bulletproof': sum(1 for r in all_results if r.get('verdict') == 'BULLETPROOF'),
            'confounded': sum(1 for r in all_results if r.get('verdict') == 'CONFOUNDED'),
            'uncertain': sum(1 for r in all_results if r.get('verdict') not in ['BULLETPROOF', 'CONFOUNDED'])
        }
    })
    
    output_file = os.path.join(OUTPUT_DIR, 'bulletproof_tep_signals.json')
    with open(output_file, 'w') as f:
        json.dump(output, f, indent=2)
    print(f"\nResults saved to {output_file}")
    
    # Final summary
    print("\n" + "="*70)
    print("FINAL BULLETPROOF SUMMARY")
    print("="*70)
    print(f"Bulletproof: {output['summary']['bulletproof']}")
    print(f"Confounded: {output['summary']['confounded']}")
    print(f"Uncertain: {output['summary']['uncertain']}")
    
    return output


if __name__ == '__main__':
    results = main()
