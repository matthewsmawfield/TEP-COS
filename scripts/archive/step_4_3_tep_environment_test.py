#!/usr/bin/env python3
"""
Step 4.3: TEP Environment Test

Tests for TEP signatures in SDSS galaxy data by looking for:
1. Anomalous velocity dispersion residuals correlated with environment
2. Systematic age/metallicity offsets in dense vs sparse environments
3. Redshift-dependent patterns that could indicate time-flow variations

TEP Prediction: If time flows slower in deep gravitational wells,
galaxies in clusters should appear "younger" than their redshift suggests
(less evolved for their cosmic epoch).

Author: M. Smawfield
Date: January 2026
"""

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy import stats
from astropy.cosmology import FlatLambdaCDM
import json
import os
from datetime import datetime
import warnings
warnings.filterwarnings('ignore')

cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

DATA_DIR = os.path.join(os.path.dirname(__file__), '..', '..', 'data', 'sdss')
RESULTS_DIR = os.path.join(os.path.dirname(__file__), '..', '..', 'results', 'outputs', 'topology')
FIGURES_DIR = os.path.join(os.path.dirname(__file__), '..', '..', 'results', 'figures', 'topology')

os.makedirs(RESULTS_DIR, exist_ok=True)
os.makedirs(FIGURES_DIR, exist_ok=True)


def load_data():
    """Load SDSS galaxy catalog."""
    print("Loading SDSS data...")
    df = pd.read_csv(os.path.join(DATA_DIR, 'sdss_galaxies.csv'))
    print(f"  Loaded {len(df):,} galaxies")
    return df


def compute_local_density(df, n_neighbors=10, max_search_deg=5.0):
    """
    Compute local galaxy density for each galaxy.
    
    Uses projected distance to Nth nearest neighbor as density proxy.
    """
    print(f"\nComputing local density (N={n_neighbors} neighbors)...")
    
    # Build KD-tree in RA/Dec space
    coords = np.column_stack([df['ra'].values, df['dec'].values])
    tree = cKDTree(coords)
    
    # Query for N+1 neighbors (includes self)
    distances, _ = tree.query(coords, k=n_neighbors + 1)
    
    # Distance to Nth neighbor (in degrees)
    nth_dist = distances[:, -1]
    
    # Convert to density proxy (inverse of area)
    # Higher density = smaller distance to Nth neighbor
    df['local_density'] = 1.0 / (nth_dist**2 + 0.01)  # Add small offset to avoid div by zero
    
    # Log density for better distribution
    df['log_density'] = np.log10(df['local_density'])
    
    # Classify into environment bins
    density_percentiles = np.percentile(df['log_density'], [20, 40, 60, 80])
    df['env_class'] = pd.cut(df['log_density'], 
                             bins=[-np.inf] + list(density_percentiles) + [np.inf],
                             labels=['void', 'sparse', 'average', 'dense', 'cluster'])
    
    print(f"  Environment distribution:")
    for env in ['void', 'sparse', 'average', 'dense', 'cluster']:
        count = (df['env_class'] == env).sum()
        print(f"    {env}: {count:,} ({100*count/len(df):.1f}%)")
    
    return df


def compute_scaling_residuals(df):
    """
    Compute residuals from standard scaling relations.
    
    Key relations:
    1. Faber-Jackson: L ∝ σ^4 (or M ∝ σ^4)
    2. Mass-SFR relation (main sequence)
    3. Mass-concentration relation
    
    Residuals from these relations may reveal TEP effects.
    """
    print("\nComputing scaling relation residuals...")
    
    # 1. Faber-Jackson residual
    # Expected: log(M) = a + b*log(σ)
    # Fit the relation
    mask = np.isfinite(df['log_mass']) & np.isfinite(df['log_sigma'])
    slope, intercept, _, _, _ = stats.linregress(df.loc[mask, 'log_sigma'], 
                                                   df.loc[mask, 'log_mass'])
    df['fj_expected'] = intercept + slope * df['log_sigma']
    df['fj_residual'] = df['log_mass'] - df['fj_expected']
    print(f"  Faber-Jackson: M ∝ σ^{slope:.2f}")
    
    # 2. Mass-SFR residual (specific SFR)
    # Expected: log(SFR) = a + b*log(M)
    mask = np.isfinite(df['log_mass']) & np.isfinite(df['log_sfr'])
    slope_sfr, intercept_sfr, _, _, _ = stats.linregress(df.loc[mask, 'log_mass'],
                                                          df.loc[mask, 'log_sfr'])
    df['sfr_expected'] = intercept_sfr + slope_sfr * df['log_mass']
    df['sfr_residual'] = df['log_sfr'] - df['sfr_expected']
    print(f"  Mass-SFR: SFR ∝ M^{slope_sfr:.2f}")
    
    # 3. Mass-concentration residual
    mask = np.isfinite(df['log_mass']) & np.isfinite(df['concentration'])
    slope_c, intercept_c, _, _, _ = stats.linregress(df.loc[mask, 'log_mass'],
                                                      df.loc[mask, 'concentration'])
    df['conc_expected'] = intercept_c + slope_c * df['log_mass']
    df['conc_residual'] = df['concentration'] - df['conc_expected']
    print(f"  Mass-Concentration: C ∝ M^{slope_c:.2f}")
    
    return df


def test_environment_dependence(df):
    """
    Test if scaling relation residuals depend on environment.
    
    TEP prediction: Galaxies in dense environments (deep potential wells)
    should show systematic offsets due to time dilation effects.
    """
    print("\n" + "=" * 70)
    print("TESTING ENVIRONMENT DEPENDENCE OF SCALING RESIDUALS")
    print("=" * 70)
    
    results = {}
    
    for residual_name in ['fj_residual', 'sfr_residual', 'conc_residual']:
        print(f"\n{residual_name.upper()}:")
        
        env_means = {}
        env_stds = {}
        
        for env in ['void', 'sparse', 'average', 'dense', 'cluster']:
            mask = (df['env_class'] == env) & np.isfinite(df[residual_name])
            values = df.loc[mask, residual_name]
            env_means[env] = float(values.mean())
            env_stds[env] = float(values.std() / np.sqrt(len(values)))
            print(f"  {env:8s}: {env_means[env]:+.4f} ± {env_stds[env]:.4f}")
        
        # Test void vs cluster difference
        void_vals = df.loc[(df['env_class'] == 'void') & np.isfinite(df[residual_name]), residual_name]
        cluster_vals = df.loc[(df['env_class'] == 'cluster') & np.isfinite(df[residual_name]), residual_name]
        
        t_stat, p_value = stats.ttest_ind(void_vals, cluster_vals)
        
        # Effect size (Cohen's d)
        pooled_std = np.sqrt((void_vals.std()**2 + cluster_vals.std()**2) / 2)
        cohens_d = (void_vals.mean() - cluster_vals.mean()) / pooled_std
        
        print(f"  Void vs Cluster: t={t_stat:.2f}, p={p_value:.2e}, Cohen's d={cohens_d:.3f}")
        
        results[residual_name] = {
            'env_means': env_means,
            'env_stds': env_stds,
            't_statistic': float(t_stat),
            'p_value': float(p_value),
            'cohens_d': float(cohens_d),
        }
    
    return results


def test_redshift_evolution(df):
    """
    Test if environment effects change with redshift.
    
    TEP prediction: If time-flow varies, the environment effect
    should accumulate over cosmic time, becoming stronger at higher z.
    """
    print("\n" + "=" * 70)
    print("TESTING REDSHIFT EVOLUTION OF ENVIRONMENT EFFECTS")
    print("=" * 70)
    
    # Split into redshift bins
    z_bins = [(0.01, 0.08), (0.08, 0.15), (0.15, 0.25), (0.25, 0.40), (0.40, 0.60)]
    
    results = []
    
    for z_min, z_max in z_bins:
        mask = (df['redshift'] >= z_min) & (df['redshift'] < z_max)
        df_bin = df[mask]
        
        if len(df_bin) < 1000:
            continue
        
        # Compute void-cluster difference in FJ residual
        void_mask = (df_bin['env_class'] == 'void') & np.isfinite(df_bin['fj_residual'])
        cluster_mask = (df_bin['env_class'] == 'cluster') & np.isfinite(df_bin['fj_residual'])
        
        if void_mask.sum() < 100 or cluster_mask.sum() < 100:
            continue
        
        void_mean = df_bin.loc[void_mask, 'fj_residual'].mean()
        cluster_mean = df_bin.loc[cluster_mask, 'fj_residual'].mean()
        diff = void_mean - cluster_mean
        
        # Error on difference
        void_sem = df_bin.loc[void_mask, 'fj_residual'].std() / np.sqrt(void_mask.sum())
        cluster_sem = df_bin.loc[cluster_mask, 'fj_residual'].std() / np.sqrt(cluster_mask.sum())
        diff_err = np.sqrt(void_sem**2 + cluster_sem**2)
        
        z_mid = (z_min + z_max) / 2
        t_lookback = cosmo.lookback_time(z_mid).value
        
        print(f"  z={z_min:.2f}-{z_max:.2f} (t={t_lookback:.1f} Gyr): "
              f"Δ(FJ) = {diff:+.4f} ± {diff_err:.4f}")
        
        results.append({
            'z_min': z_min,
            'z_max': z_max,
            'z_mid': z_mid,
            't_lookback': t_lookback,
            'fj_diff': diff,
            'fj_diff_err': diff_err,
            'n_void': int(void_mask.sum()),
            'n_cluster': int(cluster_mask.sum()),
        })
    
    # Test for trend with redshift
    if len(results) >= 3:
        z_vals = [r['z_mid'] for r in results]
        diff_vals = [r['fj_diff'] for r in results]
        
        slope, intercept, r_value, p_value, std_err = stats.linregress(z_vals, diff_vals)
        
        print(f"\n  Trend with redshift: slope = {slope:.4f} ± {std_err:.4f}")
        print(f"  Correlation: r = {r_value:.3f}, p = {p_value:.3f}")
        
        return {
            'bins': results,
            'trend_slope': float(slope),
            'trend_slope_err': float(std_err),
            'trend_r': float(r_value),
            'trend_p': float(p_value),
        }
    
    return {'bins': results}


def test_velocity_anomaly(df):
    """
    Test for anomalous velocity dispersion patterns.
    
    TEP prediction: If time flows differently, the relationship between
    velocity dispersion and other properties may show environment-dependent
    anomalies.
    """
    print("\n" + "=" * 70)
    print("TESTING VELOCITY DISPERSION ANOMALIES")
    print("=" * 70)
    
    # Compare velocity dispersion at fixed mass across environments
    mass_bins = [(9.5, 10.0), (10.0, 10.5), (10.5, 11.0), (11.0, 11.5)]
    
    results = []
    
    for m_min, m_max in mass_bins:
        mass_mask = (df['log_mass'] >= m_min) & (df['log_mass'] < m_max)
        
        void_mask = mass_mask & (df['env_class'] == 'void')
        cluster_mask = mass_mask & (df['env_class'] == 'cluster')
        
        if void_mask.sum() < 50 or cluster_mask.sum() < 50:
            continue
        
        void_sigma = df.loc[void_mask, 'log_sigma'].mean()
        cluster_sigma = df.loc[cluster_mask, 'log_sigma'].mean()
        
        void_sem = df.loc[void_mask, 'log_sigma'].std() / np.sqrt(void_mask.sum())
        cluster_sem = df.loc[cluster_mask, 'log_sigma'].std() / np.sqrt(cluster_mask.sum())
        
        diff = cluster_sigma - void_sigma
        diff_err = np.sqrt(void_sem**2 + cluster_sem**2)
        
        # Significance
        z_score = diff / diff_err if diff_err > 0 else 0
        
        print(f"  M = {m_min:.1f}-{m_max:.1f}: "
              f"σ_cluster - σ_void = {diff:+.4f} ± {diff_err:.4f} ({z_score:.1f}σ)")
        
        results.append({
            'mass_min': m_min,
            'mass_max': m_max,
            'sigma_diff': float(diff),
            'sigma_diff_err': float(diff_err),
            'z_score': float(z_score),
            'n_void': int(void_mask.sum()),
            'n_cluster': int(cluster_mask.sum()),
        })
    
    return results


def create_visualization(df, env_results, z_results, sigma_results, output_path):
    """Create visualization of results."""
    import matplotlib.pyplot as plt
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    
    # 1. Environment distribution
    ax = axes[0, 0]
    env_counts = df['env_class'].value_counts()
    colors = ['#2166ac', '#67a9cf', '#d1e5f0', '#fddbc7', '#b2182b']
    ax.bar(range(5), [env_counts.get(e, 0) for e in ['void', 'sparse', 'average', 'dense', 'cluster']],
           color=colors)
    ax.set_xticks(range(5))
    ax.set_xticklabels(['Void', 'Sparse', 'Average', 'Dense', 'Cluster'])
    ax.set_ylabel('Galaxy Count')
    ax.set_title('Environment Classification')
    
    # 2. Scaling residuals by environment
    ax = axes[0, 1]
    envs = ['void', 'sparse', 'average', 'dense', 'cluster']
    x = np.arange(len(envs))
    width = 0.25
    
    for i, (resid, label) in enumerate([('fj_residual', 'Faber-Jackson'),
                                         ('sfr_residual', 'Mass-SFR'),
                                         ('conc_residual', 'Mass-Conc')]):
        means = [env_results[resid]['env_means'][e] for e in envs]
        errs = [env_results[resid]['env_stds'][e] for e in envs]
        ax.bar(x + i*width, means, width, yerr=errs, label=label, alpha=0.8)
    
    ax.axhline(0, color='black', linestyle='--', alpha=0.5)
    ax.set_xticks(x + width)
    ax.set_xticklabels(['Void', 'Sparse', 'Avg', 'Dense', 'Cluster'])
    ax.set_ylabel('Residual')
    ax.set_title('Scaling Relation Residuals by Environment')
    ax.legend()
    
    # 3. Redshift evolution of environment effect
    ax = axes[1, 0]
    if z_results and 'bins' in z_results and len(z_results['bins']) > 0:
        z_vals = [r['z_mid'] for r in z_results['bins']]
        diff_vals = [r['fj_diff'] for r in z_results['bins']]
        diff_errs = [r['fj_diff_err'] for r in z_results['bins']]
        
        ax.errorbar(z_vals, diff_vals, yerr=diff_errs, fmt='o-', capsize=5, 
                   markersize=8, linewidth=2)
        ax.axhline(0, color='gray', linestyle='--', alpha=0.5)
        
        if 'trend_slope' in z_results:
            z_fit = np.linspace(min(z_vals), max(z_vals), 100)
            y_fit = z_results['trend_slope'] * z_fit + (np.mean(diff_vals) - z_results['trend_slope'] * np.mean(z_vals))
            ax.plot(z_fit, y_fit, 'r--', alpha=0.7, 
                   label=f'Trend: {z_results["trend_slope"]:.4f}±{z_results["trend_slope_err"]:.4f}')
            ax.legend()
    
    ax.set_xlabel('Redshift')
    ax.set_ylabel('Void - Cluster FJ Residual')
    ax.set_title('Environment Effect vs Redshift')
    
    # 4. Summary text
    ax = axes[1, 1]
    ax.axis('off')
    
    summary = """
TEP ENVIRONMENT TEST SUMMARY

HYPOTHESIS: If time flows slower in deep gravitational
wells (clusters), galaxies there should appear "younger"
(less evolved) than field galaxies at the same redshift.

KEY RESULTS:

Faber-Jackson Residual (Void vs Cluster):
"""
    
    if 'fj_residual' in env_results:
        fj = env_results['fj_residual']
        summary += f"""  Difference: {fj['env_means']['void'] - fj['env_means']['cluster']:+.4f}
  t-statistic: {fj['t_statistic']:.2f}
  p-value: {fj['p_value']:.2e}
  Cohen's d: {fj['cohens_d']:.3f}
"""
    
    if z_results and 'trend_slope' in z_results:
        summary += f"""
Redshift Evolution:
  Slope: {z_results['trend_slope']:.4f} ± {z_results['trend_slope_err']:.4f}
  Correlation: r = {z_results['trend_r']:.3f}
  p-value: {z_results['trend_p']:.3f}
"""
    
    # Interpretation
    if 'fj_residual' in env_results:
        p = env_results['fj_residual']['p_value']
        if p < 0.001:
            summary += "\nVERDICT: SIGNIFICANT ENVIRONMENT EFFECT DETECTED"
        elif p < 0.05:
            summary += "\nVERDICT: MARGINAL ENVIRONMENT EFFECT"
        else:
            summary += "\nVERDICT: NO SIGNIFICANT EFFECT"
    
    ax.text(0.05, 0.95, summary, transform=ax.transAxes, fontsize=10,
           verticalalignment='top', fontfamily='monospace')
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    plt.close()
    print(f"\nVisualization saved: {output_path}")


def main():
    """Main analysis pipeline."""
    print("=" * 70)
    print("TEP ENVIRONMENT TEST: SDSS 400K GALAXIES")
    print("=" * 70)
    print(f"Started: {datetime.now().isoformat()}")
    
    # Load data
    df = load_data()
    
    # Compute local density
    df = compute_local_density(df, n_neighbors=10)
    
    # Compute scaling residuals
    df = compute_scaling_residuals(df)
    
    # Test environment dependence
    env_results = test_environment_dependence(df)
    
    # Test redshift evolution
    z_results = test_redshift_evolution(df)
    
    # Test velocity anomalies
    sigma_results = test_velocity_anomaly(df)
    
    # Visualization
    fig_path = os.path.join(FIGURES_DIR, 'step_4_3_tep_environment.png')
    create_visualization(df, env_results, z_results, sigma_results, fig_path)
    
    # Save results
    results = {
        'metadata': {
            'timestamp': datetime.now().isoformat(),
            'n_galaxies': len(df),
            'z_range': [float(df['redshift'].min()), float(df['redshift'].max())],
        },
        'environment_effects': env_results,
        'redshift_evolution': z_results,
        'velocity_anomalies': sigma_results,
    }
    
    output_path = os.path.join(RESULTS_DIR, 'step_4_3_tep_environment.json')
    with open(output_path, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\nResults saved: {output_path}")
    
    # Summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    
    fj = env_results['fj_residual']
    print(f"\nFaber-Jackson Environment Effect:")
    print(f"  Void - Cluster difference: {fj['env_means']['void'] - fj['env_means']['cluster']:+.4f}")
    print(f"  p-value: {fj['p_value']:.2e}")
    print(f"  Cohen's d: {fj['cohens_d']:.3f}")
    
    if fj['p_value'] < 0.001:
        print("\n*** SIGNIFICANT ENVIRONMENT EFFECT DETECTED ***")
        if fj['cohens_d'] > 0:
            print("Void galaxies are MORE massive for their σ than cluster galaxies")
        else:
            print("Cluster galaxies are MORE massive for their σ than void galaxies")
    
    return results


if __name__ == '__main__':
    results = main()
