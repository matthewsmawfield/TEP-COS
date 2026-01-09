#!/usr/bin/env python3
"""
BULLETPROOF ADDITIONAL TEP SIGNALS

Tests M (Mass Discrepancy), DQ (Satellite Abundance), DH (Dust-to-Gas)
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

os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(FIGURE_DIR, exist_ok=True)


def partial_correlation(x, y, controls, data):
    """Compute partial correlation controlling for variables."""
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


def permutation_test(x, y, data, n_perm=5000):
    """Permutation test for significance."""
    mask = data[x].notna() & data[y].notna()
    subset = data[mask]
    
    if len(subset) < 100:
        return np.nan, np.nan
    
    r_obs, _ = stats.pearsonr(subset[x], subset[y])
    
    null_r = []
    y_values = subset[y].values.copy()
    
    for _ in range(n_perm):
        np.random.shuffle(y_values)
        r_null, _ = stats.pearsonr(subset[x], y_values)
        null_r.append(r_null)
    
    p_perm = np.mean(np.abs(null_r) >= np.abs(r_obs))
    return r_obs, p_perm


def bulletproof_mass_discrepancy():
    """Bulletproof Test M: Mass Discrepancy."""
    print("\n" + "="*70)
    print("BULLETPROOFING TEST M: MASS DISCREPANCY")
    print("="*70)
    
    results = {'test': 'M_Mass_Discrepancy'}
    
    # Load data
    cache_file = os.path.join(DATA_DIR, 'sdss_mass_comparison.csv')
    if not os.path.exists(cache_file):
        print("Mass comparison data not found")
        return {'test': 'M_Mass_Discrepancy', 'verdict': 'SKIPPED'}
    
    df = pd.read_csv(cache_file)
    
    # Clean data
    mask = (
        df['sigma_stars'].notna() &
        df['logM_FSPS'].notna() &
        df['logM_PCA'].notna() &
        (df['sigma_stars'] > 50) & (df['sigma_stars'] < 400)
    )
    clean = df[mask].copy()
    clean['delta_M'] = clean['logM_FSPS'] - clean['logM_PCA']
    
    print(f"\nSample size: {len(clean):,}")
    results['n_total'] = len(clean)
    
    # 1. Raw correlation
    r_raw, p_raw = stats.pearsonr(clean['sigma_stars'], clean['delta_M'])
    print(f"\n1. RAW CORRELATION: r = {r_raw:.4f}, p = {p_raw:.2e}")
    results['raw'] = {'r': float(r_raw), 'p': float(p_raw)}
    
    # 2. Bootstrap CI
    print("\n2. BOOTSTRAP CI (1000 iterations):")
    ci = bootstrap_correlation('sigma_stars', 'delta_M', clean)
    print(f"   95% CI: [{ci[0]:.4f}, {ci[2]:.4f}]")
    results['bootstrap_ci'] = {'lower': float(ci[0]), 'median': float(ci[1]), 'upper': float(ci[2])}
    
    # 3. Permutation test
    print("\n3. PERMUTATION TEST (5000 iterations):")
    r_obs, p_perm = permutation_test('sigma_stars', 'delta_M', clean)
    print(f"   Permutation p = {p_perm:.6f}")
    results['permutation_p'] = float(p_perm)
    
    # 4. Partial correlations
    print("\n4. PARTIAL CORRELATIONS:")
    control_sets = [
        (['redshift'], 'z'),
        (['redshift', 'g_minus_r'], 'z + color'),
    ]
    
    partial_results = []
    for controls, label in control_sets:
        available = [c for c in controls if c in clean.columns and clean[c].notna().sum() > 100]
        if available:
            r, p, n = partial_correlation('sigma_stars', 'delta_M', available, clean)
            print(f"   Controlling for {label}: r = {r:.4f}, p = {p:.2e}")
            partial_results.append({'controls': label, 'r': float(r), 'p': float(p)})
    results['partial_correlations'] = partial_results
    
    # 5. Binned analysis
    print("\n5. BINNED ANALYSIS:")
    sigma_bins = [(50, 100), (100, 150), (150, 200), (200, 300), (300, 400)]
    binned = []
    for lo, hi in sigma_bins:
        subset = clean[(clean['sigma_stars'] >= lo) & (clean['sigma_stars'] < hi)]
        if len(subset) > 20:
            mean = subset['delta_M'].mean()
            sem = subset['delta_M'].std() / np.sqrt(len(subset))
            print(f"   σ={lo}-{hi}: ΔM = {mean:.4f} ± {sem:.4f}, n={len(subset)}")
            binned.append({'range': f'{lo}-{hi}', 'mean': float(mean), 'sem': float(sem), 'n': len(subset)})
    results['binned'] = binned
    
    # Verdict
    ci_excludes_zero = (ci[0] > 0 and ci[2] > 0) or (ci[0] < 0 and ci[2] < 0)
    
    if ci_excludes_zero and p_perm < 0.001:
        results['verdict'] = 'BULLETPROOF'
        print("\n✓ VERDICT: BULLETPROOF")
    elif p_perm < 0.01:
        results['verdict'] = 'ROBUST'
        print("\n✓ VERDICT: ROBUST")
    else:
        results['verdict'] = 'UNCERTAIN'
        print("\n⚠ VERDICT: UNCERTAIN")
    
    return results


def bulletproof_satellite_abundance():
    """Bulletproof Test DQ: Satellite Abundance."""
    print("\n" + "="*70)
    print("BULLETPROOFING TEST DQ: SATELLITE ABUNDANCE")
    print("="*70)
    
    results = {'test': 'DQ_Satellite_Abundance'}
    
    # Load satellite data
    cache_file = os.path.join(DATA_DIR, 'sdss_cmg_survival.csv')
    if not os.path.exists(cache_file):
        print("Satellite abundance data not found - checking alternatives")
        # Try alternative files
        alt_files = ['sdss_satellite_abundance.csv', 'sdss_env_data.csv']
        for alt in alt_files:
            alt_path = os.path.join(DATA_DIR, alt)
            if os.path.exists(alt_path):
                cache_file = alt_path
                break
        else:
            return {'test': 'DQ_Satellite_Abundance', 'verdict': 'SKIPPED', 
                    'reason': 'No satellite data available'}
    
    df = pd.read_csv(cache_file)
    print(f"Loaded data with columns: {list(df.columns)[:10]}")
    
    # Look for relevant columns
    sigma_col = None
    sat_col = None
    mass_col = None
    
    for col in df.columns:
        col_lower = col.lower()
        if 'sigma' in col_lower or 'veldisp' in col_lower:
            sigma_col = col
        if 'sat' in col_lower or 'n_' in col_lower or 'neighbor' in col_lower:
            sat_col = col
        if 'mass' in col_lower:
            mass_col = col
    
    if sigma_col is None:
        print("No sigma column found")
        return {'test': 'DQ_Satellite_Abundance', 'verdict': 'SKIPPED', 
                'reason': 'No sigma column'}
    
    print(f"Using sigma: {sigma_col}, mass: {mass_col}")
    
    # If no satellite column, try compactness-based analysis
    if sat_col is None and mass_col and 'petroR50' in str(df.columns):
        print("No satellite column - using compactness proxy")
        # Use existing CMG survival data differently
        return {'test': 'DQ_Satellite_Abundance', 'verdict': 'SKIPPED',
                'reason': 'No direct satellite count available'}
    
    results['columns_found'] = {
        'sigma': sigma_col,
        'satellites': sat_col,
        'mass': mass_col
    }
    
    if sat_col is None:
        return {'test': 'DQ_Satellite_Abundance', 'verdict': 'SKIPPED',
                'reason': 'No satellite column found'}
    
    # Clean data
    mask = df[sigma_col].notna() & df[sat_col].notna()
    if mass_col:
        mask &= df[mass_col].notna()
    clean = df[mask].copy()
    
    print(f"\nSample size: {len(clean):,}")
    results['n_total'] = len(clean)
    
    if len(clean) < 100:
        return {'test': 'DQ_Satellite_Abundance', 'verdict': 'INSUFFICIENT_DATA'}
    
    # Analysis
    r_raw, p_raw = stats.pearsonr(clean[sigma_col], clean[sat_col])
    print(f"\n1. RAW CORRELATION: r = {r_raw:.4f}, p = {p_raw:.2e}")
    results['raw'] = {'r': float(r_raw), 'p': float(p_raw)}
    
    # Bootstrap
    ci = bootstrap_correlation(sigma_col, sat_col, clean)
    print(f"\n2. BOOTSTRAP CI: [{ci[0]:.4f}, {ci[2]:.4f}]")
    results['bootstrap_ci'] = {'lower': float(ci[0]), 'median': float(ci[1]), 'upper': float(ci[2])}
    
    # Partial correlation controlling for mass
    if mass_col:
        r_partial, p_partial, n = partial_correlation(sigma_col, sat_col, [mass_col], clean)
        print(f"\n3. PARTIAL (|Mass): r = {r_partial:.4f}, p = {p_partial:.2e}")
        results['partial_mass'] = {'r': float(r_partial), 'p': float(p_partial)}
    
    # Verdict
    ci_excludes_zero = (ci[0] > 0 and ci[2] > 0) or (ci[0] < 0 and ci[2] < 0)
    
    if ci_excludes_zero:
        results['verdict'] = 'BULLETPROOF'
        print("\n✓ VERDICT: BULLETPROOF")
    else:
        results['verdict'] = 'UNCERTAIN'
        print("\n⚠ VERDICT: UNCERTAIN")
    
    return results


def bulletproof_dust_gas():
    """Bulletproof Test DH: Dust-to-Gas Ratio."""
    print("\n" + "="*70)
    print("BULLETPROOFING TEST DH: DUST-TO-GAS RATIO")
    print("="*70)
    
    results = {'test': 'DH_Dust_Gas'}
    
    cache_file = os.path.join(DATA_DIR, 'sdss_dust_gas.csv')
    if not os.path.exists(cache_file):
        print("Dust-gas data not found")
        return {'test': 'DH_Dust_Gas', 'verdict': 'SKIPPED'}
    
    df = pd.read_csv(cache_file)
    print(f"Columns: {list(df.columns)}")
    
    # Find relevant columns
    sigma_col = None
    dust_col = None
    
    for col in df.columns:
        col_lower = col.lower()
        if 'sigma' in col_lower or 'veldisp' in col_lower:
            sigma_col = col
        if 'dust' in col_lower or 'ebv' in col_lower or 'av' in col_lower:
            dust_col = col
    
    if sigma_col is None or dust_col is None:
        print(f"Missing columns: sigma={sigma_col}, dust={dust_col}")
        return {'test': 'DH_Dust_Gas', 'verdict': 'SKIPPED', 
                'reason': f'Missing sigma or dust column'}
    
    print(f"Using sigma: {sigma_col}, dust: {dust_col}")
    
    # Clean
    mask = df[sigma_col].notna() & df[dust_col].notna()
    clean = df[mask].copy()
    
    # Quality cuts
    if sigma_col in clean.columns:
        mask = (clean[sigma_col] > 30) & (clean[sigma_col] < 400)
        clean = clean[mask]
    
    print(f"\nSample size: {len(clean):,}")
    results['n_total'] = len(clean)
    
    if len(clean) < 100:
        return {'test': 'DH_Dust_Gas', 'verdict': 'INSUFFICIENT_DATA'}
    
    # Raw correlation
    r_raw, p_raw = stats.pearsonr(clean[sigma_col], clean[dust_col])
    print(f"\n1. RAW CORRELATION: r = {r_raw:.4f}, p = {p_raw:.2e}")
    results['raw'] = {'r': float(r_raw), 'p': float(p_raw)}
    
    # Bootstrap
    ci = bootstrap_correlation(sigma_col, dust_col, clean)
    print(f"\n2. BOOTSTRAP CI: [{ci[0]:.4f}, {ci[2]:.4f}]")
    results['bootstrap_ci'] = {'lower': float(ci[0]), 'median': float(ci[1]), 'upper': float(ci[2])}
    
    # Permutation
    r_obs, p_perm = permutation_test(sigma_col, dust_col, clean)
    print(f"\n3. PERMUTATION p = {p_perm:.6f}")
    results['permutation_p'] = float(p_perm)
    
    # Verdict
    ci_excludes_zero = (ci[0] > 0 and ci[2] > 0) or (ci[0] < 0 and ci[2] < 0)
    
    if ci_excludes_zero and p_perm < 0.001:
        results['verdict'] = 'BULLETPROOF'
        print("\n✓ VERDICT: BULLETPROOF")
    elif p_perm < 0.01:
        results['verdict'] = 'ROBUST'
        print("\n✓ VERDICT: ROBUST")
    else:
        results['verdict'] = 'UNCERTAIN'
        print("\n⚠ VERDICT: UNCERTAIN")
    
    return results


def convert_to_native(obj):
    """Convert numpy types for JSON serialization."""
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


def main():
    print("="*70)
    print("BULLETPROOFING ADDITIONAL TEP SIGNALS")
    print("="*70)
    
    all_results = []
    
    # Test M
    all_results.append(bulletproof_mass_discrepancy())
    
    # Test DQ
    all_results.append(bulletproof_satellite_abundance())
    
    # Test DH
    all_results.append(bulletproof_dust_gas())
    
    # Summary
    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)
    
    for r in all_results:
        test = r.get('test', 'Unknown')
        verdict = r.get('verdict', 'N/A')
        print(f"  {test}: {verdict}")
    
    # Save
    output = convert_to_native({
        'tests': all_results,
        'summary': {
            'bulletproof': sum(1 for r in all_results if r.get('verdict') == 'BULLETPROOF'),
            'robust': sum(1 for r in all_results if r.get('verdict') == 'ROBUST'),
            'uncertain': sum(1 for r in all_results if r.get('verdict') == 'UNCERTAIN'),
            'skipped': sum(1 for r in all_results if r.get('verdict') == 'SKIPPED')
        }
    })
    
    output_file = os.path.join(OUTPUT_DIR, 'bulletproof_additional_signals.json')
    with open(output_file, 'w') as f:
        json.dump(output, f, indent=2)
    print(f"\nResults saved to {output_file}")
    
    return output


if __name__ == '__main__':
    results = main()
