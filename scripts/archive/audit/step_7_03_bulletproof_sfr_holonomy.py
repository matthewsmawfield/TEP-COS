#!/usr/bin/env python3
"""
BULLETPROOF SFR HOLONOMY

The SFR vs σ correlation is one of the strongest TEP-consistent signals.
Raw r = -0.59, after mass control r = -0.49, after full controls r = -0.42.

Key question: Is this TEP or standard downsizing?
"""

import os
import json
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.linear_model import LinearRegression
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')

BASE_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATA_DIR = os.path.join(BASE_DIR, 'data', 'sdss')
OUTPUT_DIR = os.path.join(BASE_DIR, 'results', 'outputs')
FIGURE_DIR = os.path.join(BASE_DIR, 'results', 'figures')


def partial_correlation(x, y, controls, data):
    """Compute partial correlation."""
    mask = data[x].notna() & data[y].notna()
    for c in controls:
        mask &= data[c].notna()
    
    subset = data[mask].copy()
    if len(subset) < 100:
        return np.nan, np.nan, 0
    
    X_ctrl = subset[controls].values
    
    reg_x = LinearRegression().fit(X_ctrl, subset[x])
    resid_x = subset[x] - reg_x.predict(X_ctrl)
    
    reg_y = LinearRegression().fit(X_ctrl, subset[y])
    resid_y = subset[y] - reg_y.predict(X_ctrl)
    
    r, p = stats.pearsonr(resid_x, resid_y)
    return r, p, len(subset)


def bootstrap_correlation(x, y, data, n_boot=1000):
    """Bootstrap CI for correlation."""
    mask = data[x].notna() & data[y].notna()
    subset = data[mask]
    
    correlations = []
    n = len(subset)
    
    for _ in range(n_boot):
        idx = np.random.choice(n, size=n, replace=True)
        r, _ = stats.pearsonr(subset[x].iloc[idx], subset[y].iloc[idx])
        correlations.append(r)
    
    return np.percentile(correlations, [2.5, 50, 97.5])


def matched_pair_analysis(x, y, match_vars, data, n_bins=5):
    """Matched-pair analysis."""
    mask = data[x].notna() & data[y].notna()
    for v in match_vars:
        mask &= data[v].notna()
    
    subset = data[mask].copy()
    
    for v in match_vars:
        subset[f'{v}_bin'] = pd.qcut(subset[v], n_bins, labels=False, duplicates='drop')
    
    bin_cols = [f'{v}_bin' for v in match_vars]
    subset['match_key'] = subset[bin_cols].astype(str).agg('-'.join, axis=1)
    
    bin_results = []
    for key, group in subset.groupby('match_key'):
        if len(group) >= 20:
            r, p = stats.pearsonr(group[x], group[y])
            bin_results.append({'key': key, 'r': r, 'n': len(group)})
    
    if not bin_results:
        return np.nan, np.nan, 0, []
    
    total_n = sum(b['n'] for b in bin_results)
    weighted_r = sum(b['r'] * b['n'] for b in bin_results) / total_n
    
    # Count consistent sign
    n_consistent = sum(1 for b in bin_results if np.sign(b['r']) == np.sign(weighted_r))
    sign_frac = n_consistent / len(bin_results)
    
    return weighted_r, sign_frac, len(bin_results), bin_results


def main():
    print("="*70)
    print("BULLETPROOFING SFR HOLONOMY (sSFR vs σ)")
    print("="*70)
    
    results = {'test': 'SFR_Holonomy'}
    
    # Load data
    twin_file = os.path.join(DATA_DIR, 'sdss_twin_base_sample_with_size.csv')
    if not os.path.exists(twin_file):
        print("Data not found")
        return {'verdict': 'SKIPPED'}
    
    df = pd.read_csv(twin_file)
    df = df.rename(columns={
        'veldisp': 'sigma',
        'log_mass': 'logMass',
        'd4000': 'D4000',
        'hbeta': 'Hbeta',
        'mgb': 'Mgb',
        'fe5270': 'Fe5270',
        'fe5335': 'Fe5335'
    })
    
    # Compute MgFe
    df['Fe_avg'] = 0.5 * (df['Fe5270'] + df['Fe5335'])
    mask = (df['Mgb'] > 0) & (df['Fe_avg'] > 0)
    df['MgFe'] = np.nan
    df.loc[mask, 'MgFe'] = np.log10(df.loc[mask, 'Mgb'] / df.loc[mask, 'Fe_avg'])
    
    # Use D4000 as SFR proxy (inverse relationship)
    # Low D4000 = high sSFR
    # Create sSFR proxy: -D4000 or use Hbeta as young star indicator
    df['sSFR_proxy'] = -df['D4000']  # Higher = more star forming
    
    # Clean
    mask = (
        df['sigma'].notna() & df['sSFR_proxy'].notna() &
        (df['sigma'] > 50) & (df['sigma'] < 400) &
        df['logMass'].notna()
    )
    clean = df[mask].copy()
    
    print(f"\nSample size: {len(clean):,}")
    results['n_total'] = len(clean)
    
    # 1. Raw correlation
    r_raw, p_raw = stats.pearsonr(clean['sigma'], clean['sSFR_proxy'])
    print(f"\n1. RAW CORRELATION (σ vs -D4000):")
    print(f"   r = {r_raw:.4f}, p = {p_raw:.2e}")
    results['raw'] = {'r': float(r_raw), 'p': float(p_raw)}
    
    # 2. Partial correlations
    print("\n2. PARTIAL CORRELATIONS:")
    
    control_sets = [
        (['logMass'], 'Mass'),
        (['logMass', 'MgFe'], 'Mass + [Mg/Fe]'),
        (['logMass', 'MgFe', 'redshift'], 'Mass + [Mg/Fe] + z'),
    ]
    
    partial_results = []
    for controls, label in control_sets:
        available = [c for c in controls if c in clean.columns]
        if available:
            r, p, n = partial_correlation('sigma', 'sSFR_proxy', available, clean)
            print(f"   {label}: r = {r:.4f}, p = {p:.2e}")
            partial_results.append({
                'controls': label,
                'r': float(r) if not np.isnan(r) else None,
                'p': float(p) if not np.isnan(p) else None
            })
    results['partial_correlations'] = partial_results
    
    # 3. Bootstrap CI
    print("\n3. BOOTSTRAP CI (1000 iterations):")
    ci = bootstrap_correlation('sigma', 'sSFR_proxy', clean)
    print(f"   95% CI: [{ci[0]:.4f}, {ci[2]:.4f}]")
    results['bootstrap_ci'] = {
        'lower': float(ci[0]), 
        'median': float(ci[1]), 
        'upper': float(ci[2])
    }
    
    # 4. Matched-pair analysis
    print("\n4. MATCHED-PAIR ANALYSIS:")
    mp_r, mp_frac, mp_n, mp_bins = matched_pair_analysis(
        'sigma', 'sSFR_proxy', ['logMass', 'MgFe'], clean
    )
    print(f"   Weighted r across {mp_n} bins: {mp_r:.4f}")
    print(f"   Sign consistency: {mp_frac:.1%}")
    results['matched_pairs'] = {
        'weighted_r': float(mp_r),
        'sign_fraction': float(mp_frac),
        'n_bins': mp_n
    }
    
    # 5. Binned analysis by σ
    print("\n5. BINNED ANALYSIS BY σ:")
    sigma_bins = [(50, 100), (100, 150), (150, 200), (200, 300), (300, 400)]
    binned = []
    for lo, hi in sigma_bins:
        subset = clean[(clean['sigma'] >= lo) & (clean['sigma'] < hi)]
        if len(subset) > 50:
            mean_ssfr = subset['sSFR_proxy'].mean()
            sem = subset['sSFR_proxy'].std() / np.sqrt(len(subset))
            print(f"   σ={lo}-{hi}: <-D4000> = {mean_ssfr:.4f} ± {sem:.4f}, n={len(subset)}")
            binned.append({
                'range': f'{lo}-{hi}',
                'mean': float(mean_ssfr),
                'sem': float(sem),
                'n': len(subset)
            })
    results['binned'] = binned
    
    # 6. Critical test: Does signal survive after controlling for [Mg/Fe]?
    # This tests whether it's truly TEP (time-related) vs standard downsizing
    r_mgfe, p_mgfe, n = partial_correlation('sigma', 'sSFR_proxy', ['MgFe'], clean)
    print(f"\n6. CRITICAL TEST - After [Mg/Fe] control only:")
    print(f"   r = {r_mgfe:.4f}, p = {p_mgfe:.2e}")
    results['mgfe_control'] = {'r': float(r_mgfe), 'p': float(p_mgfe)}
    
    # 7. Alternative test: Hbeta as young star proxy
    if 'Hbeta' in clean.columns and clean['Hbeta'].notna().sum() > 1000:
        print("\n7. ALTERNATIVE: Hβ (young star proxy) vs σ:")
        hb_mask = clean['Hbeta'].notna() & np.isfinite(clean['Hbeta'])
        hb_clean = clean[hb_mask]
        r_hb, p_hb = stats.pearsonr(hb_clean['sigma'], hb_clean['Hbeta'])
        print(f"   Raw r(Hβ, σ) = {r_hb:.4f}")
        
        r_hb_ctrl, p_hb_ctrl, _ = partial_correlation('sigma', 'Hbeta', ['D4000', 'logMass'], hb_clean)
        print(f"   After D4000+Mass control: r = {r_hb_ctrl:.4f}")
        results['hbeta_test'] = {
            'raw_r': float(r_hb),
            'controlled_r': float(r_hb_ctrl) if not np.isnan(r_hb_ctrl) else None
        }
    
    # Verdict
    print("\n" + "="*70)
    print("VERDICT ASSESSMENT:")
    
    ci_excludes_zero = (ci[0] > 0 and ci[2] > 0) or (ci[0] < 0 and ci[2] < 0)
    survives_controls = all(
        p['r'] is not None and np.sign(p['r']) == np.sign(r_raw)
        for p in partial_results
    )
    
    print(f"   ✓ CI excludes zero: {ci_excludes_zero}")
    print(f"   ✓ Sign survives all controls: {survives_controls}")
    print(f"   ✓ Matched-pair consistency > 60%: {mp_frac > 0.6}")
    
    # Note: This test is DEGENERATE with standard astrophysics
    # Both TEP and downsizing predict sSFR decreases with σ
    # The signal is REAL but NOT DISCRIMINATING
    
    if ci_excludes_zero and survives_controls and mp_frac > 0.6:
        results['verdict'] = 'BULLETPROOF_BUT_DEGENERATE'
        print("\n✓ VERDICT: BULLETPROOF but DEGENERATE")
        print("   Signal is robust but indistinguishable from standard downsizing.")
    elif ci_excludes_zero and survives_controls:
        results['verdict'] = 'ROBUST_BUT_DEGENERATE'
        print("\n✓ VERDICT: ROBUST but DEGENERATE")
    else:
        results['verdict'] = 'UNCERTAIN'
        print("\n⚠ VERDICT: UNCERTAIN")
    
    results['interpretation'] = (
        "The sSFR-σ correlation is statistically robust but does not discriminate TEP "
        "from standard downsizing. Both predict lower sSFR at higher σ. "
        "This test confirms consistency but cannot provide unique TEP evidence."
    )
    
    # Save
    def convert(obj):
        if isinstance(obj, (np.integer,)): return int(obj)
        if isinstance(obj, (np.floating,)): return float(obj)
        if isinstance(obj, np.ndarray): return obj.tolist()
        if isinstance(obj, (np.bool_, bool)): return bool(obj)
        if isinstance(obj, dict): return {k: convert(v) for k, v in obj.items()}
        if isinstance(obj, list): return [convert(i) for i in obj]
        return obj
    
    output_file = os.path.join(OUTPUT_DIR, 'bulletproof_sfr_holonomy.json')
    with open(output_file, 'w') as f:
        json.dump(convert(results), f, indent=2)
    print(f"\nResults saved to {output_file}")
    
    return results


if __name__ == '__main__':
    results = main()
