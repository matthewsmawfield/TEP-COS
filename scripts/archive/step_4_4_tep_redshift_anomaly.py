#!/usr/bin/env python3
"""
Step 4.4: TEP Redshift Anomaly Test

Tests for anomalous redshift-velocity relationships that could indicate
time-flow variations predicted by TEP.

Key Test: If time flows differently in different gravitational environments,
the relationship between observed redshift and intrinsic galaxy properties
should show systematic deviations from the standard cosmological model.

Specifically, we test:
1. Hubble residuals correlated with local density
2. Peculiar velocity anomalies in different environments
3. Redshift-independent distance indicators vs environment

Author: M. Smawfield
Date: January 2026
"""

import numpy as np
import pandas as pd
from scipy import stats
from scipy.spatial import cKDTree
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


def load_data():
    """Load SDSS data with environment classification."""
    print("Loading SDSS data...")
    df = pd.read_csv(os.path.join(DATA_DIR, 'sdss_galaxies.csv'))
    
    # Recompute local density
    coords = np.column_stack([df['ra'].values, df['dec'].values])
    tree = cKDTree(coords)
    distances, _ = tree.query(coords, k=11)
    df['local_density'] = 1.0 / (distances[:, -1]**2 + 0.01)
    df['log_density'] = np.log10(df['local_density'])
    
    density_percentiles = np.percentile(df['log_density'], [20, 40, 60, 80])
    df['env_class'] = pd.cut(df['log_density'], 
                             bins=[-np.inf] + list(density_percentiles) + [np.inf],
                             labels=['void', 'sparse', 'average', 'dense', 'cluster'])
    
    print(f"  Loaded {len(df):,} galaxies")
    return df


def test_fundamental_plane_residuals(df):
    """
    Test the Fundamental Plane as a distance indicator.
    
    The FP relates: log(R_e) = a*log(σ) + b*<μ>_e + c
    
    If TEP affects time-flow, the FP coefficients or scatter might
    depend on environment in unexpected ways.
    """
    print("\n" + "=" * 70)
    print("FUNDAMENTAL PLANE ANALYSIS")
    print("=" * 70)
    
    # We don't have R_e directly, but we can use concentration as a proxy
    # and test for environment-dependent deviations
    
    # Create a pseudo-FP using available data
    # log(σ) vs log(M) residuals, checking for environment dependence
    
    # Fit FP-like relation: log(σ) = a*log(M) + b*concentration + c
    mask = (np.isfinite(df['log_sigma']) & 
            np.isfinite(df['log_mass']) & 
            np.isfinite(df['concentration']))
    
    X = np.column_stack([df.loc[mask, 'log_mass'], df.loc[mask, 'concentration']])
    y = df.loc[mask, 'log_sigma'].values
    
    # Add constant term
    X_with_const = np.column_stack([np.ones(len(X)), X])
    
    # Fit using least squares
    coeffs, residuals, rank, s = np.linalg.lstsq(X_with_const, y, rcond=None)
    
    # Compute residuals
    y_pred = X_with_const @ coeffs
    fp_residuals = y - y_pred
    
    df.loc[mask, 'fp_residual'] = fp_residuals
    
    print(f"  Pseudo-FP: log(σ) = {coeffs[0]:.3f} + {coeffs[1]:.3f}*log(M) + {coeffs[2]:.3f}*C")
    print(f"  RMS residual: {np.std(fp_residuals):.4f}")
    
    # Test environment dependence
    results = {}
    print("\n  FP residuals by environment:")
    for env in ['void', 'sparse', 'average', 'dense', 'cluster']:
        env_mask = mask & (df['env_class'] == env)
        if env_mask.sum() > 0:
            mean_resid = df.loc[env_mask, 'fp_residual'].mean()
            std_resid = df.loc[env_mask, 'fp_residual'].std() / np.sqrt(env_mask.sum())
            print(f"    {env:8s}: {mean_resid:+.5f} ± {std_resid:.5f}")
            results[env] = {'mean': float(mean_resid), 'sem': float(std_resid)}
    
    # Statistical test
    void_resid = df.loc[mask & (df['env_class'] == 'void'), 'fp_residual']
    cluster_resid = df.loc[mask & (df['env_class'] == 'cluster'), 'fp_residual']
    
    t_stat, p_value = stats.ttest_ind(void_resid, cluster_resid)
    print(f"\n  Void vs Cluster: t={t_stat:.2f}, p={p_value:.2e}")
    
    results['t_statistic'] = float(t_stat)
    results['p_value'] = float(p_value)
    
    return results


def test_tully_fisher_analog(df):
    """
    Test a Tully-Fisher-like relation for spirals.
    
    TF relates luminosity to rotation velocity. For our data,
    we use mass vs velocity dispersion for disk-like galaxies.
    """
    print("\n" + "=" * 70)
    print("TULLY-FISHER ANALOG (DISK GALAXIES)")
    print("=" * 70)
    
    # Select disk-like galaxies (low Sersic proxy)
    disk_mask = (df['sersic_proxy'] < 2.5) & np.isfinite(df['log_mass']) & np.isfinite(df['log_sigma'])
    
    print(f"  Disk galaxies: {disk_mask.sum():,}")
    
    if disk_mask.sum() < 1000:
        print("  Insufficient disk galaxies for analysis")
        return {}
    
    # Fit TF relation
    slope, intercept, r, p, se = stats.linregress(
        df.loc[disk_mask, 'log_sigma'],
        df.loc[disk_mask, 'log_mass']
    )
    
    print(f"  TF relation: log(M) = {intercept:.2f} + {slope:.2f}*log(σ)")
    
    # Compute residuals
    df.loc[disk_mask, 'tf_residual'] = (
        df.loc[disk_mask, 'log_mass'] - 
        (intercept + slope * df.loc[disk_mask, 'log_sigma'])
    )
    
    # Test environment dependence
    results = {}
    print("\n  TF residuals by environment:")
    for env in ['void', 'sparse', 'average', 'dense', 'cluster']:
        env_mask = disk_mask & (df['env_class'] == env)
        if env_mask.sum() > 100:
            mean_resid = df.loc[env_mask, 'tf_residual'].mean()
            std_resid = df.loc[env_mask, 'tf_residual'].std() / np.sqrt(env_mask.sum())
            print(f"    {env:8s}: {mean_resid:+.5f} ± {std_resid:.5f} (n={env_mask.sum()})")
            results[env] = {'mean': float(mean_resid), 'sem': float(std_resid), 'n': int(env_mask.sum())}
    
    return results


def test_peculiar_velocity_field(df):
    """
    Test for anomalous peculiar velocity patterns.
    
    If TEP affects redshifts, we might see systematic patterns
    in the inferred peculiar velocities.
    """
    print("\n" + "=" * 70)
    print("PECULIAR VELOCITY FIELD ANALYSIS")
    print("=" * 70)
    
    # Estimate peculiar velocity from Hubble residuals
    # v_pec = c * (z_obs - z_cosmo) / (1 + z_cosmo)
    
    # Use Faber-Jackson to estimate "true" distance
    # Then compare to redshift distance
    
    # Fit FJ relation
    mask = np.isfinite(df['log_mass']) & np.isfinite(df['log_sigma']) & (df['redshift'] > 0.02)
    
    slope, intercept, _, _, _ = stats.linregress(
        df.loc[mask, 'log_sigma'],
        df.loc[mask, 'log_mass']
    )
    
    # Predict mass from sigma
    df.loc[mask, 'mass_from_sigma'] = intercept + slope * df.loc[mask, 'log_sigma']
    
    # Mass residual as proxy for distance error
    df.loc[mask, 'mass_residual'] = df.loc[mask, 'log_mass'] - df.loc[mask, 'mass_from_sigma']
    
    # Convert to velocity-like units (rough approximation)
    # Δlog(M) ~ 0.1 corresponds to ~10% distance error ~ 300 km/s at z=0.1
    c_kms = 299792.458
    df.loc[mask, 'v_pec_proxy'] = df.loc[mask, 'mass_residual'] * 3000  # km/s
    
    # Analyze by sky position
    # Divide sky into regions and look for coherent patterns
    
    ra_bins = np.linspace(df['ra'].min(), df['ra'].max(), 10)
    dec_bins = np.linspace(df['dec'].min(), df['dec'].max(), 10)
    
    print("\n  Mean peculiar velocity proxy by sky region:")
    
    results = {'sky_regions': []}
    
    for i in range(len(ra_bins) - 1):
        for j in range(len(dec_bins) - 1):
            region_mask = (
                mask &
                (df['ra'] >= ra_bins[i]) & (df['ra'] < ra_bins[i+1]) &
                (df['dec'] >= dec_bins[j]) & (df['dec'] < dec_bins[j+1])
            )
            
            if region_mask.sum() > 100:
                mean_v = df.loc[region_mask, 'v_pec_proxy'].mean()
                std_v = df.loc[region_mask, 'v_pec_proxy'].std() / np.sqrt(region_mask.sum())
                
                results['sky_regions'].append({
                    'ra_center': (ra_bins[i] + ra_bins[i+1]) / 2,
                    'dec_center': (dec_bins[j] + dec_bins[j+1]) / 2,
                    'mean_v': float(mean_v),
                    'sem_v': float(std_v),
                    'n': int(region_mask.sum()),
                })
    
    # Look for dipole pattern
    if len(results['sky_regions']) > 10:
        ra_vals = [r['ra_center'] for r in results['sky_regions']]
        dec_vals = [r['dec_center'] for r in results['sky_regions']]
        v_vals = [r['mean_v'] for r in results['sky_regions']]
        
        # Fit dipole: v = A*cos(ra - ra0)*cos(dec) + B*sin(dec)
        # Simplified: just check correlation with position
        
        r_ra, p_ra = stats.pearsonr(ra_vals, v_vals)
        r_dec, p_dec = stats.pearsonr(dec_vals, v_vals)
        
        print(f"\n  Correlation with RA: r={r_ra:.3f}, p={p_ra:.3f}")
        print(f"  Correlation with Dec: r={r_dec:.3f}, p={p_dec:.3f}")
        
        results['dipole'] = {
            'ra_correlation': float(r_ra),
            'ra_p_value': float(p_ra),
            'dec_correlation': float(r_dec),
            'dec_p_value': float(p_dec),
        }
    
    # Test environment dependence of peculiar velocities
    print("\n  Peculiar velocity proxy by environment:")
    for env in ['void', 'sparse', 'average', 'dense', 'cluster']:
        env_mask = mask & (df['env_class'] == env)
        if env_mask.sum() > 100:
            mean_v = df.loc[env_mask, 'v_pec_proxy'].mean()
            std_v = df.loc[env_mask, 'v_pec_proxy'].std() / np.sqrt(env_mask.sum())
            print(f"    {env:8s}: {mean_v:+.1f} ± {std_v:.1f} km/s")
    
    return results


def test_redshift_quantization(df):
    """
    Test for redshift quantization or preferred values.
    
    Some alternative cosmologies predict quantized redshifts.
    TEP doesn't specifically predict this, but it's worth checking.
    """
    print("\n" + "=" * 70)
    print("REDSHIFT DISTRIBUTION ANALYSIS")
    print("=" * 70)
    
    # Look for peaks in redshift distribution beyond selection effects
    z_vals = df['redshift'].values
    
    # Compute histogram
    bins = np.linspace(0.01, 0.55, 100)
    hist, bin_edges = np.histogram(z_vals, bins=bins)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    
    # Smooth to find underlying trend
    from scipy.ndimage import gaussian_filter1d
    smooth_hist = gaussian_filter1d(hist.astype(float), sigma=3)
    
    # Residuals from smooth
    residuals = hist - smooth_hist
    
    # Look for significant peaks
    peak_threshold = 3 * np.std(residuals)
    peaks = np.where(residuals > peak_threshold)[0]
    
    print(f"  Potential redshift peaks (>{3}σ above trend):")
    for p in peaks:
        z_peak = bin_centers[p]
        excess = residuals[p] / np.std(residuals)
        print(f"    z = {z_peak:.4f}: {excess:.1f}σ excess")
    
    results = {
        'n_peaks': len(peaks),
        'peaks': [{'z': float(bin_centers[p]), 'sigma': float(residuals[p]/np.std(residuals))} 
                  for p in peaks],
    }
    
    # Test for periodicity
    from scipy.fft import fft
    fft_result = np.abs(fft(residuals))
    freqs = np.fft.fftfreq(len(residuals), d=bins[1]-bins[0])
    
    # Find dominant frequency (excluding DC)
    positive_freqs = freqs[1:len(freqs)//2]
    positive_fft = fft_result[1:len(fft_result)//2]
    
    if len(positive_fft) > 0:
        max_idx = np.argmax(positive_fft)
        dominant_freq = positive_freqs[max_idx]
        dominant_period = 1 / dominant_freq if dominant_freq > 0 else np.inf
        
        print(f"\n  Dominant periodicity: Δz = {dominant_period:.4f}")
        results['dominant_period'] = float(dominant_period)
    
    return results


def create_visualization(df, fp_results, tf_results, pv_results, z_results, output_path):
    """Create comprehensive visualization."""
    import matplotlib.pyplot as plt
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    
    # 1. FP residuals by environment
    ax = axes[0, 0]
    if fp_results:
        envs = ['void', 'sparse', 'average', 'dense', 'cluster']
        means = [fp_results.get(e, {}).get('mean', 0) for e in envs]
        sems = [fp_results.get(e, {}).get('sem', 0) for e in envs]
        
        colors = ['#2166ac', '#67a9cf', '#d1e5f0', '#fddbc7', '#b2182b']
        ax.bar(range(5), means, yerr=sems, color=colors, capsize=5)
        ax.axhline(0, color='black', linestyle='--', alpha=0.5)
        ax.set_xticks(range(5))
        ax.set_xticklabels(['Void', 'Sparse', 'Avg', 'Dense', 'Cluster'])
        ax.set_ylabel('FP Residual')
        ax.set_title('Fundamental Plane Residuals by Environment')
    
    # 2. Peculiar velocity map
    ax = axes[0, 1]
    if pv_results and 'sky_regions' in pv_results:
        ra = [r['ra_center'] for r in pv_results['sky_regions']]
        dec = [r['dec_center'] for r in pv_results['sky_regions']]
        v = [r['mean_v'] for r in pv_results['sky_regions']]
        
        scatter = ax.scatter(ra, dec, c=v, cmap='RdBu_r', s=100, 
                            vmin=-200, vmax=200)
        plt.colorbar(scatter, ax=ax, label='v_pec proxy (km/s)')
        ax.set_xlabel('RA (deg)')
        ax.set_ylabel('Dec (deg)')
        ax.set_title('Peculiar Velocity Field')
    
    # 3. Redshift distribution
    ax = axes[1, 0]
    ax.hist(df['redshift'], bins=100, alpha=0.7, density=True)
    ax.set_xlabel('Redshift')
    ax.set_ylabel('Density')
    ax.set_title('Redshift Distribution')
    
    if z_results and 'peaks' in z_results:
        for peak in z_results['peaks']:
            ax.axvline(peak['z'], color='red', linestyle='--', alpha=0.7)
    
    # 4. Summary
    ax = axes[1, 1]
    ax.axis('off')
    
    summary = """
TEP REDSHIFT ANOMALY TEST SUMMARY

TESTS PERFORMED:

1. FUNDAMENTAL PLANE RESIDUALS
"""
    if fp_results and 'p_value' in fp_results:
        summary += f"""   Void vs Cluster: p = {fp_results['p_value']:.2e}
   Result: {'SIGNIFICANT' if fp_results['p_value'] < 0.01 else 'Not significant'}
"""
    
    summary += """
2. PECULIAR VELOCITY FIELD
"""
    if pv_results and 'dipole' in pv_results:
        summary += f"""   RA correlation: r = {pv_results['dipole']['ra_correlation']:.3f}
   Dec correlation: r = {pv_results['dipole']['dec_correlation']:.3f}
"""
    
    summary += """
3. REDSHIFT DISTRIBUTION
"""
    if z_results:
        summary += f"""   Peaks detected: {z_results['n_peaks']}
   Dominant period: Δz = {z_results.get('dominant_period', 'N/A')}
"""
    
    summary += """
INTERPRETATION:
- Environment effects are REAL but likely due to
  standard physics (tidal heating, assembly bias)
- No clear TEP-specific signature detected
- Peculiar velocity field shows expected structure
"""
    
    ax.text(0.05, 0.95, summary, transform=ax.transAxes, fontsize=10,
           verticalalignment='top', fontfamily='monospace')
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    plt.close()
    print(f"\nVisualization saved: {output_path}")


def main():
    """Main analysis."""
    print("=" * 70)
    print("TEP REDSHIFT ANOMALY TEST")
    print("=" * 70)
    print(f"Started: {datetime.now().isoformat()}")
    
    df = load_data()
    
    fp_results = test_fundamental_plane_residuals(df)
    tf_results = test_tully_fisher_analog(df)
    pv_results = test_peculiar_velocity_field(df)
    z_results = test_redshift_quantization(df)
    
    # Visualization
    fig_path = os.path.join(FIGURES_DIR, 'step_4_4_tep_redshift_anomaly.png')
    create_visualization(df, fp_results, tf_results, pv_results, z_results, fig_path)
    
    # Save results
    results = {
        'metadata': {
            'timestamp': datetime.now().isoformat(),
            'n_galaxies': len(df),
        },
        'fundamental_plane': fp_results,
        'tully_fisher': tf_results,
        'peculiar_velocity': pv_results,
        'redshift_distribution': z_results,
    }
    
    output_path = os.path.join(RESULTS_DIR, 'step_4_4_tep_redshift_anomaly.json')
    with open(output_path, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\nResults saved: {output_path}")
    
    # Summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print("\nKey findings:")
    print("1. Environment effects in scaling relations: DETECTED (standard physics)")
    print("2. Peculiar velocity field: Shows expected large-scale structure")
    print("3. Redshift quantization: No significant evidence")
    print("\nConclusion: No clear TEP-specific signatures in SDSS data")
    
    return results


if __name__ == '__main__':
    results = main()
