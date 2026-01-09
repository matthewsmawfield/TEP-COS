#!/usr/bin/env python3
"""
Test T: Fundamental Plane Residual Analysis

Hypothesis:
The Fundamental Plane relates R_e, σ, and surface brightness. Under TEP, 
if time dilation affects stellar evolution in deep potentials, the M/L 
ratio could be biased, leading to FP residuals that correlate with σ.

TEP Prediction:
  FP residuals should correlate with σ if M/L is systematically biased.
"""

import os
import json
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt

BASE_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATA_DIR = os.path.join(BASE_DIR, 'data', 'sdss')
OUTPUT_DIR = os.path.join(BASE_DIR, 'results', 'outputs')
FIGURE_DIR = os.path.join(BASE_DIR, 'results', 'figures')

os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(FIGURE_DIR, exist_ok=True)


def load_fj_data():
    """Load Faber-Jackson data."""
    fj_file = os.path.join(DATA_DIR, 'sdss_faber_jackson.csv')
    twin_file = os.path.join(DATA_DIR, 'sdss_twin_base_sample_with_size.csv')
    
    if os.path.exists(twin_file):
        print(f"Loading comprehensive data from {twin_file}")
        df = pd.read_csv(twin_file)
        df = df.rename(columns={
            'specobjid': 'specObjID',
            'veldisp': 'sigma',
            'log_mass': 'logMass',
            'petroR50_r_arcsec': 'Re_arcsec',
            'd4000': 'D4000',
            'mgb': 'Mgb',
            'fe5270': 'Fe5270',
            'fe5335': 'Fe5335'
        })
        return df
    
    if os.path.exists(fj_file):
        print(f"Loading Faber-Jackson data from {fj_file}")
        return pd.read_csv(fj_file)
    
    return None


def analyze_fp_residuals(df):
    """Analyze Fundamental Plane residuals."""
    results = {}
    
    # Clean data
    required_cols = ['sigma', 'logMass']
    for col in required_cols:
        if col not in df.columns:
            print(f"Missing column: {col}")
            return {'verdict': 'Skipped', 'interpretation': f'Missing {col} column'}
    
    mask = (
        df['sigma'].notna() & df['logMass'].notna() &
        (df['sigma'] > 50) & (df['sigma'] < 400) &
        (df['logMass'] > 9) & (df['logMass'] < 12.5)
    )
    clean = df[mask].copy()
    
    print(f"\nAnalyzing {len(clean)} galaxies")
    results['n_galaxies'] = len(clean)
    
    # 1. Faber-Jackson relation: L ∝ σ^n (or M ∝ σ^n)
    # Standard: n ≈ 4
    clean['log_sigma'] = np.log10(clean['sigma'])
    
    slope_fj, intercept_fj, r_fj, p_fj, stderr_fj = stats.linregress(
        clean['log_sigma'], clean['logMass']
    )
    
    print(f"\n1. Faber-Jackson Relation:")
    print(f"   log(M*) = {slope_fj:.3f} × log(σ) + {intercept_fj:.3f}")
    print(f"   r = {r_fj:.4f}, stderr = {stderr_fj:.4f}")
    print(f"   Standard FJ slope ≈ 4, observed: {slope_fj:.3f}")
    
    results['faber_jackson'] = {
        'slope': float(slope_fj),
        'slope_err': float(stderr_fj),
        'intercept': float(intercept_fj),
        'r': float(r_fj),
        'deviation_from_4': float(slope_fj - 4)
    }
    
    # 2. FJ residuals vs σ
    clean['M_predicted_FJ'] = slope_fj * clean['log_sigma'] + intercept_fj
    clean['FJ_residual'] = clean['logMass'] - clean['M_predicted_FJ']
    
    # By construction, residuals should have zero mean. 
    # But do they correlate with σ indicating non-linearity?
    r_resid, p_resid = stats.pearsonr(clean['sigma'], clean['FJ_residual'])
    print(f"\n2. FJ Residual vs σ (testing non-linearity):")
    print(f"   r = {r_resid:.6f}, p = {p_resid:.2e}")
    
    results['fj_residual_vs_sigma'] = {
        'r': float(r_resid),
        'p': float(p_resid)
    }
    
    # 3. If we have size, compute FP residual
    if 'Re_arcsec' in clean.columns:
        mask_size = clean['Re_arcsec'].notna() & (clean['Re_arcsec'] > 0.5)
        subset = clean[mask_size].copy()
        
        if len(subset) > 1000:
            # FP: log(R_e) = a × log(σ) + b × <μ> + c
            # Or equivalently: log(M/L) = f(σ, R_e)
            # Simplified: look at M vs σ^2 × R_e (virial theorem)
            
            subset['log_Re'] = np.log10(subset['Re_arcsec'])
            subset['log_virial'] = 2 * subset['log_sigma'] + subset['log_Re']
            
            # Mass vs virial estimator
            slope_vir, intercept_vir, r_vir, p_vir, stderr_vir = stats.linregress(
                subset['log_virial'], subset['logMass']
            )
            
            print(f"\n3. Mass vs Virial Estimator (σ² × R_e):")
            print(f"   slope = {slope_vir:.3f} ± {stderr_vir:.3f}")
            print(f"   r = {r_vir:.4f}")
            print(f"   Expected slope = 1.0 if M/L is constant")
            
            results['virial_relation'] = {
                'slope': float(slope_vir),
                'slope_err': float(stderr_vir),
                'r': float(r_vir),
                'deviation_from_unity': float(slope_vir - 1)
            }
            
            # Virial residuals vs σ
            subset['M_pred_vir'] = slope_vir * subset['log_virial'] + intercept_vir
            subset['vir_residual'] = subset['logMass'] - subset['M_pred_vir']
            
            r_vir_resid, p_vir_resid = stats.pearsonr(subset['sigma'], subset['vir_residual'])
            print(f"\n4. Virial Residual vs σ:")
            print(f"   r = {r_vir_resid:.4f}, p = {p_vir_resid:.2e}")
            
            results['virial_residual_vs_sigma'] = {
                'r': float(r_vir_resid),
                'p': float(p_vir_resid)
            }
    
    # 4. Binned FJ slope analysis
    print("\n5. FJ Slope by σ range (testing TEP tilt):")
    sigma_ranges = [(50, 120), (120, 180), (180, 250), (250, 400)]
    slope_results = []
    
    for lo, hi in sigma_ranges:
        subset = clean[(clean['sigma'] >= lo) & (clean['sigma'] < hi)]
        if len(subset) > 100:
            s, i, r, p, se = stats.linregress(subset['log_sigma'], subset['logMass'])
            print(f"   σ = {lo}-{hi}: slope = {s:.3f} ± {se:.3f}, n = {len(subset)}")
            slope_results.append({
                'sigma_range': f"{lo}-{hi}",
                'slope': float(s),
                'slope_err': float(se),
                'n': len(subset)
            })
    
    results['slope_by_sigma_range'] = slope_results
    
    # Check if slope varies systematically with σ
    if len(slope_results) >= 3:
        slopes = [x['slope'] for x in slope_results]
        slope_trend = slopes[-1] - slopes[0]  # high σ slope - low σ slope
        results['slope_trend'] = float(slope_trend)
        print(f"\n   Slope trend (high - low σ): {slope_trend:.3f}")
    
    # Verdict
    if abs(r_resid) > 0.05 and p_resid < 0.01:
        results['verdict'] = 'Signal'
        results['interpretation'] = ('FJ residuals correlate with σ, indicating '
                                    'non-linear M/L variation with potential depth.')
    else:
        results['verdict'] = 'Null'
        results['interpretation'] = 'FJ residuals do not significantly correlate with σ.'
    
    print(f"\n=== VERDICT: {results['verdict']} ===")
    print(results['interpretation'])
    
    return results, clean


def create_figure(df, results):
    """Create visualization."""
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    clean = df[(df['sigma'] > 50) & (df['sigma'] < 400)].copy()
    clean['log_sigma'] = np.log10(clean['sigma'])
    
    # 1. Faber-Jackson relation
    ax1 = axes[0, 0]
    ax1.scatter(clean['log_sigma'], clean['logMass'], alpha=0.1, s=3)
    
    x_fit = np.linspace(clean['log_sigma'].min(), clean['log_sigma'].max(), 100)
    slope = results['faber_jackson']['slope']
    intercept = results['faber_jackson']['intercept']
    ax1.plot(x_fit, slope * x_fit + intercept, 'r-', linewidth=2,
            label=f'Slope = {slope:.2f}')
    ax1.plot(x_fit, 4 * x_fit + (clean['logMass'].mean() - 4 * clean['log_sigma'].mean()),
            'g--', linewidth=2, label='Standard (slope=4)')
    
    ax1.set_xlabel('log(σ)')
    ax1.set_ylabel('log(M*/M☉)')
    ax1.set_title('Faber-Jackson Relation')
    ax1.legend()
    
    # 2. Slope by σ range
    ax2 = axes[0, 1]
    if results.get('slope_by_sigma_range'):
        x_centers = []
        slopes = []
        slope_errs = []
        for s in results['slope_by_sigma_range']:
            lo, hi = map(int, s['sigma_range'].split('-'))
            x_centers.append((lo + hi) / 2)
            slopes.append(s['slope'])
            slope_errs.append(s['slope_err'])
        
        ax2.errorbar(x_centers, slopes, yerr=slope_errs, fmt='o-', capsize=5, markersize=8)
        ax2.axhline(results['faber_jackson']['slope'], color='red', linestyle='--',
                   label=f"Global slope = {results['faber_jackson']['slope']:.2f}")
        ax2.set_xlabel('σ (km/s)')
        ax2.set_ylabel('Local FJ Slope')
        ax2.set_title('FJ Slope vs σ Range (TEP Tilt Test)')
        ax2.legend()
    
    # 3. Residuals vs σ
    ax3 = axes[1, 0]
    if 'FJ_residual' in clean.columns:
        ax3.scatter(clean['sigma'], clean['FJ_residual'], alpha=0.1, s=3)
        ax3.axhline(0, color='red', linestyle='--')
        ax3.set_xlabel('σ (km/s)')
        ax3.set_ylabel('FJ Residual (dex)')
        ax3.set_title(f"FJ Residual vs σ (r = {results['fj_residual_vs_sigma']['r']:.4f})")
    
    # 4. Summary
    ax4 = axes[1, 1]
    ax4.axis('off')
    
    summary = f"""
TEST T: Fundamental Plane Residual Analysis

HYPOTHESIS:
TEP could cause M/L biases that appear
as FP residuals correlating with σ.

RESULTS:
• N = {results['n_galaxies']:,} galaxies
• FJ slope: {results['faber_jackson']['slope']:.3f} ± {results['faber_jackson']['slope_err']:.3f}
  (standard ≈ 4, deviation: {results['faber_jackson']['deviation_from_4']:.3f})
• Residual vs σ: r = {results['fj_residual_vs_sigma']['r']:.6f}
  p = {results['fj_residual_vs_sigma']['p']:.2e}

VERDICT: {results['verdict']}
"""
    ax4.text(0.05, 0.95, summary, transform=ax4.transAxes,
            fontsize=10, verticalalignment='top', fontfamily='monospace',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    plt.tight_layout()
    fig_path = os.path.join(FIGURE_DIR, 'test_t_fp_residual.png')
    plt.savefig(fig_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"\nFigure saved to {fig_path}")
    
    return fig_path


def main():
    print("=" * 60)
    print("TEST T: Fundamental Plane Residual Analysis")
    print("=" * 60)
    
    df = load_fj_data()
    
    if df is None or len(df) < 100:
        print("ERROR: Insufficient data")
        return {'verdict': 'Skipped', 'interpretation': 'No FP data available.'}
    
    print(f"Loaded {len(df)} galaxies")
    
    results, clean_df = analyze_fp_residuals(df)
    fig_path = create_figure(clean_df, results)
    results['figure_path'] = fig_path
    
    output_file = os.path.join(OUTPUT_DIR, 'sdss_test_t_fp_residual_results.json')
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\nResults saved to {output_file}")
    
    return results


if __name__ == '__main__':
    results = main()
