#!/usr/bin/env python3
"""
Test M: Dynamical vs Stellar Population Mass Discrepancy

Hypothesis:
Dynamical masses (from σ) probe gravity directly. Stellar population masses 
(from SED fitting) assume standard stellar evolution. Under TEP, if time flows 
slower in deeper potentials, spectroscopic ages would be biased, potentially
causing mass discrepancies that correlate with potential depth.

TEP Prediction:
  Δlog(M) = log(M_dyn) - log(M_spec) should correlate with σ
  At high σ, stellar populations appear younger → lower inferred M/L → 
  systematic underestimate of M_spec relative to M_dyn.
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


def load_mass_data():
    """Load mass comparison data."""
    cache_file = os.path.join(DATA_DIR, 'sdss_mass_comparison.csv')
    
    if os.path.exists(cache_file):
        print(f"Loading mass comparison data from {cache_file}")
        df = pd.read_csv(cache_file)
        return df
    
    print("Mass comparison data not found")
    return None


def analyze_mass_discrepancy(df):
    """Analyze mass discrepancy vs velocity dispersion."""
    results = {}
    
    # Clean data
    mask = (
        df['sigma_stars'].notna() &
        df['logM_FSPS'].notna() &
        df['logM_PCA'].notna() &
        (df['sigma_stars'] > 50) & (df['sigma_stars'] < 400)
    )
    clean = df[mask].copy()
    
    print(f"\nAnalyzing {len(clean)} galaxies")
    results['n_galaxies'] = len(clean)
    
    # Compute dynamical mass proxy: M_dyn ∝ σ² R_e / G
    # For FP analysis: log(M_dyn) ≈ 2×log(σ) + constant (at fixed R_e)
    # We'll use σ as the proxy and look at discrepancy between mass estimates
    
    # 1. FSPS vs PCA mass discrepancy
    clean['delta_M_FSPS_PCA'] = clean['logM_FSPS'] - clean['logM_PCA']
    
    r_fsps_pca, p_fsps_pca = stats.pearsonr(clean['sigma_stars'], clean['delta_M_FSPS_PCA'])
    print(f"\n1. FSPS - PCA mass difference vs σ:")
    print(f"   r = {r_fsps_pca:.4f}, p = {p_fsps_pca:.2e}")
    
    results['fsps_pca_discrepancy'] = {
        'r': float(r_fsps_pca),
        'p': float(p_fsps_pca),
        'mean_diff': float(clean['delta_M_FSPS_PCA'].mean())
    }
    
    # 2. FSPS vs Portsmouth mass discrepancy
    if 'logM_Port' in clean.columns:
        clean['delta_M_FSPS_Port'] = clean['logM_FSPS'] - clean['logM_Port']
        mask_port = clean['logM_Port'].notna()
        
        r_fsps_port, p_fsps_port = stats.pearsonr(
            clean.loc[mask_port, 'sigma_stars'], 
            clean.loc[mask_port, 'delta_M_FSPS_Port']
        )
        print(f"\n2. FSPS - Portsmouth mass difference vs σ:")
        print(f"   r = {r_fsps_port:.4f}, p = {p_fsps_port:.2e}")
        
        results['fsps_port_discrepancy'] = {
            'r': float(r_fsps_port),
            'p': float(p_fsps_port),
            'mean_diff': float(clean.loc[mask_port, 'delta_M_FSPS_Port'].mean())
        }
    
    # 3. Virial mass estimator: M_vir ∝ σ² (simplified)
    # Compare with spectroscopic mass
    clean['log_sigma2'] = 2 * np.log10(clean['sigma_stars'])
    clean['M_vir_proxy'] = clean['log_sigma2'] + 5.0  # arbitrary normalization
    clean['delta_vir_FSPS'] = clean['M_vir_proxy'] - clean['logM_FSPS']
    
    r_vir, p_vir = stats.pearsonr(clean['sigma_stars'], clean['delta_vir_FSPS'])
    print(f"\n3. Virial proxy - FSPS mass vs σ:")
    print(f"   r = {r_vir:.4f}, p = {p_vir:.2e}")
    
    # Slope of M_FSPS vs σ²
    slope, intercept, r_val, p_slope, stderr = stats.linregress(
        clean['log_sigma2'], clean['logM_FSPS']
    )
    print(f"\n4. log(M_FSPS) vs log(σ²):")
    print(f"   slope = {slope:.4f} ± {stderr:.4f}")
    print(f"   Expected for virial: slope ≈ 1.0")
    print(f"   Deviation from virial: {slope - 1:.4f}")
    
    results['virial_slope'] = {
        'slope': float(slope),
        'slope_err': float(stderr),
        'deviation_from_unity': float(slope - 1)
    }
    
    # 5. Binned analysis
    print("\n5. Mass discrepancy by σ bin:")
    sigma_bins = [(50, 100), (100, 150), (150, 200), (200, 300), (300, 400)]
    binned = []
    
    for lo, hi in sigma_bins:
        subset = clean[(clean['sigma_stars'] >= lo) & (clean['sigma_stars'] < hi)]
        if len(subset) > 20:
            mean_diff = subset['delta_M_FSPS_PCA'].mean()
            sem = subset['delta_M_FSPS_PCA'].std() / np.sqrt(len(subset))
            print(f"   σ = {lo}-{hi}: Δlog(M) = {mean_diff:.4f} ± {sem:.4f}, n = {len(subset)}")
            binned.append({
                'sigma_range': f"{lo}-{hi}",
                'mean_diff': float(mean_diff),
                'stderr': float(sem),
                'n': len(subset)
            })
    
    results['binned_by_sigma'] = binned
    
    # Verdict
    # If high-σ galaxies have systematically different mass estimates, 
    # this could indicate TEP-related age biases
    if abs(r_fsps_pca) > 0.1 and p_fsps_pca < 0.05:
        if r_fsps_pca > 0:
            results['verdict'] = 'Signal'
            results['interpretation'] = ('Mass discrepancy increases with σ. '
                                        'Different methods disagree more at high σ, '
                                        'potentially due to age-related systematics.')
        else:
            results['verdict'] = 'Signal'
            results['interpretation'] = ('Mass discrepancy decreases with σ. '
                                        'High-σ galaxies show more agreement between methods.')
    else:
        results['verdict'] = 'Null'
        results['interpretation'] = 'No significant σ-dependence in mass discrepancy.'
    
    print(f"\n=== VERDICT: {results['verdict']} ===")
    print(results['interpretation'])
    
    return results, clean


def create_figure(df, results):
    """Create visualization."""
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    clean = df[df['delta_M_FSPS_PCA'].notna()].copy()
    
    # 1. Mass discrepancy vs σ
    ax1 = axes[0, 0]
    ax1.scatter(clean['sigma_stars'], clean['delta_M_FSPS_PCA'], alpha=0.2, s=5)
    ax1.axhline(0, color='red', linestyle='--')
    ax1.set_xlabel('σ (km/s)')
    ax1.set_ylabel('log(M_FSPS) - log(M_PCA)')
    ax1.set_title(f"Mass Discrepancy vs σ (r = {results['fsps_pca_discrepancy']['r']:.3f})")
    
    # 2. M_FSPS vs σ²
    ax2 = axes[0, 1]
    ax2.scatter(clean['log_sigma2'], clean['logM_FSPS'], alpha=0.2, s=5)
    x_fit = np.linspace(clean['log_sigma2'].min(), clean['log_sigma2'].max(), 100)
    slope = results['virial_slope']['slope']
    intercept = clean['logM_FSPS'].mean() - slope * clean['log_sigma2'].mean()
    ax2.plot(x_fit, slope * x_fit + intercept, 'r-', linewidth=2, 
            label=f'Slope = {slope:.2f}')
    ax2.plot(x_fit, 1.0 * x_fit + (clean['logM_FSPS'].mean() - clean['log_sigma2'].mean()), 
            'g--', linewidth=2, label='Virial (slope=1)')
    ax2.set_xlabel('log(σ²)')
    ax2.set_ylabel('log(M_FSPS)')
    ax2.set_title('Mass-σ Relation')
    ax2.legend()
    
    # 3. Binned means
    ax3 = axes[1, 0]
    if results.get('binned_by_sigma'):
        x_centers = []
        y_means = []
        y_errs = []
        for b in results['binned_by_sigma']:
            lo, hi = map(int, b['sigma_range'].split('-'))
            x_centers.append((lo + hi) / 2)
            y_means.append(b['mean_diff'])
            y_errs.append(b['stderr'])
        
        ax3.errorbar(x_centers, y_means, yerr=y_errs, fmt='o-', capsize=5, markersize=8)
        ax3.axhline(0, color='red', linestyle='--')
        ax3.set_xlabel('σ (km/s)')
        ax3.set_ylabel('Mean Δlog(M)')
        ax3.set_title('Binned Mass Discrepancy')
    
    # 4. Summary
    ax4 = axes[1, 1]
    ax4.axis('off')
    
    summary = f"""
TEST M: Dynamical vs Stellar Mass Discrepancy

HYPOTHESIS:
If TEP causes age biases in SED fitting,
mass estimates should disagree at high σ.

RESULTS:
• N = {results['n_galaxies']:,} galaxies
• FSPS-PCA discrepancy vs σ: r = {results['fsps_pca_discrepancy']['r']:.4f}
  p = {results['fsps_pca_discrepancy']['p']:.2e}
• Virial slope: {results['virial_slope']['slope']:.3f} ± {results['virial_slope']['slope_err']:.3f}
  (deviation from unity: {results['virial_slope']['deviation_from_unity']:.3f})

VERDICT: {results['verdict']}
"""
    ax4.text(0.05, 0.95, summary, transform=ax4.transAxes,
            fontsize=10, verticalalignment='top', fontfamily='monospace',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    plt.tight_layout()
    fig_path = os.path.join(FIGURE_DIR, 'test_m_mass_discrepancy.png')
    plt.savefig(fig_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"\nFigure saved to {fig_path}")
    
    return fig_path


def main():
    print("=" * 60)
    print("TEST M: Dynamical vs Stellar Population Mass Discrepancy")
    print("=" * 60)
    
    df = load_mass_data()
    
    if df is None or len(df) < 100:
        print("ERROR: Insufficient data")
        return {'verdict': 'Skipped', 'interpretation': 'No mass comparison data available.'}
    
    print(f"Loaded {len(df)} galaxies")
    
    results, clean_df = analyze_mass_discrepancy(df)
    fig_path = create_figure(clean_df, results)
    results['figure_path'] = fig_path
    
    output_file = os.path.join(OUTPUT_DIR, 'sdss_test_m_mass_discrepancy_results.json')
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\nResults saved to {output_file}")
    
    return results


if __name__ == '__main__':
    results = main()
