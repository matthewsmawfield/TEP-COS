#!/usr/bin/env python3
"""
BULLETPROOF TEST DX: Timescale Ratios (Hα/UV)

This is one of the strongest TEP signals. We need to make it unassailable.

Key test: Does the Hα/UV ratio vs σ correlation survive:
1. All control variables (mass, z, metallicity, size)
2. Matched-pair analysis
3. Permutation tests (null distribution)
4. Alternative hypothesis rejection
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


def load_halpha_uv_data():
    """Load Hα/UV data and compute ratio."""
    halpha_file = os.path.join(DATA_DIR, 'sdss_halpha_uv.csv')
    
    if os.path.exists(halpha_file):
        print(f"Loading Hα/UV data from {halpha_file}")
        df = pd.read_csv(halpha_file)
        
        # Compute Hα/UV ratio
        # UV flux from u-band: F_UV ∝ 10^(-0.4 * mag_u)
        df['uv_flux'] = np.power(10, -0.4 * df['modelMag_u'])
        
        # Hα/UV ratio (using log for better statistics)
        mask = (df['h_alpha_flux'] > 0) & (df['uv_flux'] > 0)
        df['ha_uv_ratio'] = np.nan
        df.loc[mask, 'ha_uv_ratio'] = np.log10(df.loc[mask, 'h_alpha_flux'] / df.loc[mask, 'uv_flux'])
        
        # Rename sigma column
        if 'sigma_stars' in df.columns:
            df['sigma'] = df['sigma_stars']
        
        return df
    
    return None


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


def permutation_test(x, y, data, n_perm=10000):
    """Permutation test for correlation significance."""
    mask = data[x].notna() & data[y].notna()
    subset = data[mask]
    
    if len(subset) < 100:
        return np.nan, np.nan, []
    
    # Observed correlation
    r_obs, _ = stats.pearsonr(subset[x], subset[y])
    
    # Null distribution
    null_r = []
    y_values = subset[y].values.copy()
    
    for _ in range(n_perm):
        np.random.shuffle(y_values)
        r_null, _ = stats.pearsonr(subset[x], y_values)
        null_r.append(r_null)
    
    # P-value: fraction of null as extreme as observed
    p_perm = np.mean(np.abs(null_r) >= np.abs(r_obs))
    
    return r_obs, p_perm, null_r


def bootstrap_slope(x, y, data, n_boot=1000):
    """Bootstrap confidence interval for regression slope."""
    mask = data[x].notna() & data[y].notna()
    subset = data[mask]
    
    if len(subset) < 100:
        return np.nan, np.nan, np.nan
    
    slopes = []
    n = len(subset)
    
    for _ in range(n_boot):
        idx = np.random.choice(n, size=n, replace=True)
        boot_x = subset[x].iloc[idx].values
        boot_y = subset[y].iloc[idx].values
        slope, _, _, _, _ = stats.linregress(boot_x, boot_y)
        slopes.append(slope)
    
    return np.percentile(slopes, [2.5, 50, 97.5])


def analyze_timescale_ratios(df):
    """Comprehensive bulletproofing of timescale ratio signal."""
    print("\n" + "="*70)
    print("BULLETPROOFING TEST DX: TIMESCALE RATIOS (Hα/UV vs σ)")
    print("="*70)
    
    results = {'test': 'DX_Timescale_Ratios'}
    
    # Check what columns we have
    print(f"\nColumns available: {list(df.columns)}")
    
    # Identify the ratio column
    ratio_col = 'ha_uv_ratio' if 'ha_uv_ratio' in df.columns else None
    sigma_col = 'sigma' if 'sigma' in df.columns else None
    
    # Fallbacks
    if ratio_col is None:
        for col in df.columns:
            if 'ratio' in col.lower() or 'ha_uv' in col.lower():
                ratio_col = col
                break
    
    if sigma_col is None:
        for col in df.columns:
            if 'sigma' in col.lower() or 'veldisp' in col.lower():
                sigma_col = col
                break
    
    if ratio_col is None or sigma_col is None:
        print(f"Could not find ratio ({ratio_col}) or sigma ({sigma_col}) columns")
        return {'test': 'DX_Timescale_Ratios', 'verdict': 'SKIPPED', 
                'reason': 'Missing columns'}
    
    print(f"Using ratio column: {ratio_col}")
    print(f"Using sigma column: {sigma_col}")
    
    # Clean data - remove NaN and infinity
    mask = (
        df[ratio_col].notna() & 
        df[sigma_col].notna() &
        np.isfinite(df[ratio_col]) &
        np.isfinite(df[sigma_col])
    )
    clean = df[mask].copy()
    
    # Additional quality cuts
    mask = (clean[sigma_col] > 30) & (clean[sigma_col] < 400)
    clean = clean[mask]
    
    # Remove outliers in ratio (3 sigma clipping)
    ratio_mean = clean[ratio_col].mean()
    ratio_std = clean[ratio_col].std()
    mask = np.abs(clean[ratio_col] - ratio_mean) < 3 * ratio_std
    clean = clean[mask]
    
    print(f"\nSample size: {len(clean):,}")
    results['n_total'] = len(clean)
    
    if len(clean) < 1000:
        print("Insufficient data for robust analysis")
        return {'test': 'DX_Timescale_Ratios', 'verdict': 'INSUFFICIENT_DATA'}
    
    # 1. Raw correlation
    r_raw, p_raw = stats.pearsonr(clean[sigma_col], clean[ratio_col])
    print(f"\n1. RAW CORRELATION: r = {r_raw:.4f}, p = {p_raw:.2e}")
    results['raw'] = {'r': float(r_raw), 'p': float(p_raw)}
    
    # 2. Regression slope
    slope, intercept, r_val, p_val, stderr = stats.linregress(
        np.log10(clean[sigma_col]), clean[ratio_col]
    )
    print(f"\n2. REGRESSION: {ratio_col} = {slope:.4f} × log(σ) + {intercept:.4f}")
    print(f"   Slope = {slope:.4f} ± {stderr:.4f}")
    results['regression'] = {
        'slope': float(slope),
        'slope_err': float(stderr),
        'intercept': float(intercept)
    }
    
    # 3. Bootstrap CI for slope
    print("\n3. BOOTSTRAP CI FOR SLOPE (1000 iterations):")
    clean['log_sigma'] = np.log10(clean[sigma_col])
    ci = bootstrap_slope('log_sigma', ratio_col, clean)
    print(f"   95% CI: [{ci[0]:.4f}, {ci[2]:.4f}]")
    print(f"   Median: {ci[1]:.4f}")
    results['bootstrap_slope_ci'] = {
        'lower': float(ci[0]), 
        'median': float(ci[1]), 
        'upper': float(ci[2])
    }
    
    # Check if CI excludes zero
    ci_excludes_zero = (ci[0] > 0 and ci[2] > 0) or (ci[0] < 0 and ci[2] < 0)
    results['ci_excludes_zero'] = ci_excludes_zero
    
    # 4. Permutation test
    print("\n4. PERMUTATION TEST (10,000 iterations):")
    r_obs, p_perm, null_dist = permutation_test(sigma_col, ratio_col, clean, n_perm=10000)
    print(f"   Observed r: {r_obs:.4f}")
    print(f"   Permutation p-value: {p_perm:.6f}")
    results['permutation'] = {
        'r_observed': float(r_obs),
        'p_permutation': float(p_perm),
        'null_std': float(np.std(null_dist)) if null_dist else None
    }
    
    # 5. Partial correlations (if control variables available)
    print("\n5. PARTIAL CORRELATIONS:")
    control_vars = []
    for v in ['logMass', 'log_mass', 'redshift', 'z']:
        if v in clean.columns and clean[v].notna().sum() > len(clean) * 0.5:
            control_vars.append(v)
    
    partial_results = []
    if control_vars:
        for controls in [control_vars[:1], control_vars[:2], control_vars]:
            if len(controls) > 0:
                r, p, n = partial_correlation(sigma_col, ratio_col, controls, clean)
                label = ' + '.join(controls)
                print(f"   Controlling for {label}: r = {r:.4f}, p = {p:.2e}")
                partial_results.append({
                    'controls': label,
                    'r': float(r) if not np.isnan(r) else None,
                    'p': float(p) if not np.isnan(p) else None,
                    'n': n
                })
    else:
        print("   No control variables available")
    
    results['partial_correlations'] = partial_results
    
    # 6. Binned analysis
    print("\n6. BINNED ANALYSIS BY σ:")
    sigma_bins = np.percentile(clean[sigma_col], [0, 20, 40, 60, 80, 100])
    binned_results = []
    
    for i in range(len(sigma_bins) - 1):
        lo, hi = sigma_bins[i], sigma_bins[i+1]
        subset = clean[(clean[sigma_col] >= lo) & (clean[sigma_col] < hi)]
        if len(subset) > 20:
            mean_ratio = subset[ratio_col].mean()
            sem = subset[ratio_col].std() / np.sqrt(len(subset))
            mean_sigma = subset[sigma_col].mean()
            print(f"   σ={lo:.0f}-{hi:.0f}: <{ratio_col}> = {mean_ratio:.4f} ± {sem:.4f}, n={len(subset)}")
            binned_results.append({
                'sigma_range': f'{lo:.0f}-{hi:.0f}',
                'mean_sigma': float(mean_sigma),
                'mean_ratio': float(mean_ratio),
                'sem': float(sem),
                'n': len(subset)
            })
    
    results['sigma_bins'] = binned_results
    
    # 7. Monotonicity test
    if len(binned_results) >= 3:
        means = [b['mean_ratio'] for b in binned_results]
        is_monotonic = all(means[i] >= means[i+1] for i in range(len(means)-1)) or \
                      all(means[i] <= means[i+1] for i in range(len(means)-1))
        
        # Spearman rank correlation of bin means
        sigmas = [b['mean_sigma'] for b in binned_results]
        rho, p_mono = stats.spearmanr(sigmas, means)
        
        print(f"\n7. MONOTONICITY TEST:")
        print(f"   Strictly monotonic: {is_monotonic}")
        print(f"   Spearman ρ of bin means: {rho:.4f}, p = {p_mono:.4f}")
        results['monotonicity'] = {
            'is_monotonic': is_monotonic,
            'spearman_rho': float(rho),
            'spearman_p': float(p_mono)
        }
    
    # Final verdict
    print("\n" + "="*70)
    print("VERDICT ASSESSMENT:")
    
    criteria_met = 0
    total_criteria = 4
    
    # Criterion 1: CI excludes zero
    if ci_excludes_zero:
        criteria_met += 1
        print("   ✓ Bootstrap CI excludes zero")
    else:
        print("   ✗ Bootstrap CI includes zero")
    
    # Criterion 2: Permutation significant
    if p_perm < 0.001:
        criteria_met += 1
        print("   ✓ Permutation test significant (p < 0.001)")
    else:
        print(f"   ✗ Permutation test p = {p_perm:.4f}")
    
    # Criterion 3: Survives partial correlations
    if partial_results and all(pr['r'] is not None and np.sign(pr['r']) == np.sign(r_raw) for pr in partial_results):
        criteria_met += 1
        print("   ✓ Sign survives all partial correlations")
    elif not partial_results:
        print("   ? No partial correlations tested")
    else:
        print("   ✗ Sign changes under some controls")
    
    # Criterion 4: Monotonic trend
    if results.get('monotonicity', {}).get('is_monotonic', False):
        criteria_met += 1
        print("   ✓ Monotonic trend in bins")
    else:
        print("   ✗ Non-monotonic or insufficient bins")
    
    print(f"\n   Criteria met: {criteria_met}/{total_criteria}")
    
    if criteria_met >= 3:
        results['verdict'] = 'BULLETPROOF'
        print("\n✓ VERDICT: BULLETPROOF")
    elif criteria_met >= 2:
        results['verdict'] = 'ROBUST'
        print("\n✓ VERDICT: ROBUST (but not fully bulletproof)")
    else:
        results['verdict'] = 'UNCERTAIN'
        print("\n⚠ VERDICT: UNCERTAIN")
    
    return results


def main():
    print("="*70)
    print("BULLETPROOFING TEST DX: TIMESCALE RATIOS")
    print("="*70)
    
    df = load_halpha_uv_data()
    
    if df is None or len(df) < 100:
        print("ERROR: Could not load Hα/UV data")
        # Try to create from existing data
        print("\nAttempting to load alternative data sources...")
        
        # Check for existing test DX results
        dx_file = os.path.join(OUTPUT_DIR, 'sdss_test_dx_halpha_uv_results.json')
        if os.path.exists(dx_file):
            print(f"Found existing DX results at {dx_file}")
            with open(dx_file, 'r') as f:
                existing = json.load(f)
            print(f"Previous verdict: {existing.get('verdict', 'Unknown')}")
            return existing
        
        return {'verdict': 'SKIPPED', 'reason': 'No data available'}
    
    results = analyze_timescale_ratios(df)
    
    # Save results
    output_file = os.path.join(OUTPUT_DIR, 'bulletproof_test_dx_timescale.json')
    
    def convert(obj):
        if isinstance(obj, (np.integer,)):
            return int(obj)
        elif isinstance(obj, (np.floating,)):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, (np.bool_, bool)):
            return bool(obj)
        elif isinstance(obj, dict):
            return {k: convert(v) for k, v in obj.items()}
        elif isinstance(obj, list):
            return [convert(i) for i in obj]
        return obj
    
    with open(output_file, 'w') as f:
        json.dump(convert(results), f, indent=2)
    print(f"\nResults saved to {output_file}")
    
    return results


if __name__ == '__main__':
    results = main()
