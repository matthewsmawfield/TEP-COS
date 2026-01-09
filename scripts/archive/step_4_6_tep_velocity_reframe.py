#!/usr/bin/env python3
"""
Step 4.6: TEP Velocity Dispersion Reframe

CRITICAL INSIGHT: Velocity dispersion σ is measured in km/s.
If time flows slower in deep potential wells, then:

σ_observed = σ_intrinsic × (dt_local / dt_observer)

Where dt_local/dt_observer < 1 in deep wells (time runs slower).

This means:
- Standard physics: σ_cluster > σ_field → "tidal heating"
- TEP reframe: σ_cluster_observed > σ_field_observed could be because
  we're measuring with a slower clock, not because σ_intrinsic is higher

TEST: If TEP is correct, the σ excess in clusters should scale with
the gravitational potential in a specific way:

Δσ/σ ≈ ΔΦ/c² (to first order in GR)

But TEP might predict a DIFFERENT scaling!

Author: M. Smawfield
Date: January 2026
"""

import numpy as np
import pandas as pd
from scipy import stats
from scipy.spatial import cKDTree
from scipy.optimize import curve_fit
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


def load_and_prepare_data():
    """Load SDSS data and compute potential."""
    print("Loading and preparing data...")
    df = pd.read_csv(os.path.join(DATA_DIR, 'sdss_galaxies.csv'))
    
    # Compute 3D positions
    z = df['redshift'].values
    d_comoving = cosmo.comoving_distance(z).value
    ra_rad = np.radians(df['ra'].values)
    dec_rad = np.radians(df['dec'].values)
    
    x = d_comoving * np.cos(dec_rad) * np.cos(ra_rad)
    y = d_comoving * np.cos(dec_rad) * np.sin(ra_rad)
    z_coord = d_comoving * np.sin(dec_rad)
    
    coords_3d = np.column_stack([x, y, z_coord])
    tree = cKDTree(coords_3d)
    
    # Get distance to 20th neighbor
    distances, _ = tree.query(coords_3d, k=21)
    r_20 = distances[:, -1]  # Mpc
    
    # Estimate potential
    n_neighbors = 20
    M_avg = 1e11  # M_sun
    M_enclosed = n_neighbors * M_avg
    
    G = 4.302e-6  # kpc/M_sun * (km/s)^2
    c = 299792.458
    r_kpc = r_20 * 1000
    
    df['phi_over_c2'] = -G * M_enclosed / (r_kpc * c**2)
    df['r_20_mpc'] = r_20
    
    # Compute sigma residual at fixed mass
    mask = np.isfinite(df['log_mass']) & np.isfinite(df['log_sigma'])
    slope, intercept, _, _, _ = stats.linregress(
        df.loc[mask, 'log_mass'], df.loc[mask, 'log_sigma']
    )
    df['sigma_residual'] = df['log_sigma'] - (intercept + slope * df['log_mass'])
    df['fj_slope'] = slope
    df['fj_intercept'] = intercept
    
    print(f"  Loaded {len(df):,} galaxies")
    print(f"  Faber-Jackson: log(σ) = {intercept:.3f} + {slope:.3f} × log(M)")
    
    return df


def test_sigma_potential_scaling(df):
    """
    Test how σ residual scales with gravitational potential.
    
    GR prediction: Δσ/σ ∝ ΔΦ/c²
    TEP might predict different scaling!
    """
    print("\n" + "=" * 70)
    print("SIGMA-POTENTIAL SCALING TEST")
    print("=" * 70)
    
    # Bin by potential and compute mean sigma residual
    n_bins = 20
    phi_percentiles = np.percentile(df['phi_over_c2'], np.linspace(0, 100, n_bins + 1))
    
    bin_phi = []
    bin_sigma = []
    bin_sigma_err = []
    
    for i in range(n_bins):
        mask = (df['phi_over_c2'] >= phi_percentiles[i]) & (df['phi_over_c2'] < phi_percentiles[i+1])
        if mask.sum() > 100:
            bin_phi.append(df.loc[mask, 'phi_over_c2'].mean())
            bin_sigma.append(df.loc[mask, 'sigma_residual'].mean())
            bin_sigma_err.append(df.loc[mask, 'sigma_residual'].std() / np.sqrt(mask.sum()))
    
    bin_phi = np.array(bin_phi)
    bin_sigma = np.array(bin_sigma)
    bin_sigma_err = np.array(bin_sigma_err)
    
    # Fit linear model: σ_resid = a + b × Φ/c²
    def linear(x, a, b):
        return a + b * x
    
    popt, pcov = curve_fit(linear, bin_phi, bin_sigma, sigma=bin_sigma_err)
    a_fit, b_fit = popt
    a_err, b_err = np.sqrt(np.diag(pcov))
    
    print(f"\n  Linear fit: σ_resid = {a_fit:.4f} + {b_fit:.0f} × (Φ/c²)")
    print(f"  Slope: {b_fit:.0f} ± {b_err:.0f}")
    
    # GR prediction for slope
    # If Δσ/σ = ΔΦ/c², and σ ~ 200 km/s, then Δlog(σ) ≈ 0.434 × ΔΦ/c²
    # So slope should be ~0.434 in log units
    gr_slope = 0.434
    
    print(f"\n  GR prediction for slope: ~{gr_slope:.3f}")
    print(f"  Observed slope: {b_fit:.0f}")
    print(f"  Ratio (observed/GR): {b_fit/gr_slope:.0f}×")
    
    # Fit power law: σ_resid = a × |Φ/c²|^n
    def power_law(x, a, n):
        return a * np.abs(x)**n
    
    # Only fit where sigma_resid > 0 for power law
    pos_mask = bin_sigma > 0
    if pos_mask.sum() > 5:
        try:
            popt_pl, pcov_pl = curve_fit(power_law, np.abs(bin_phi[pos_mask]), 
                                          bin_sigma[pos_mask], p0=[1e6, 0.5])
            a_pl, n_pl = popt_pl
            print(f"\n  Power law fit: σ_resid ∝ |Φ/c²|^{n_pl:.2f}")
            print(f"  (GR predicts n = 1.0)")
        except:
            n_pl = None
    else:
        n_pl = None
    
    return {
        'linear_intercept': float(a_fit),
        'linear_slope': float(b_fit),
        'linear_slope_err': float(b_err),
        'gr_predicted_slope': float(gr_slope),
        'slope_ratio': float(b_fit / gr_slope),
        'power_law_exponent': float(n_pl) if n_pl else None,
        'bin_phi': bin_phi.tolist(),
        'bin_sigma': bin_sigma.tolist(),
        'bin_sigma_err': bin_sigma_err.tolist(),
    }


def test_tep_time_dilation_correction(df):
    """
    Apply TEP time dilation correction and see if it removes the environment effect.
    
    If TEP is correct:
    σ_intrinsic = σ_observed × (1 + Φ/c²)  [to first order]
    
    After correction, the environment dependence should DISAPPEAR.
    """
    print("\n" + "=" * 70)
    print("TEP TIME DILATION CORRECTION TEST")
    print("=" * 70)
    
    # Original sigma residual vs potential correlation
    mask = np.isfinite(df['phi_over_c2']) & np.isfinite(df['sigma_residual'])
    r_orig, p_orig = stats.pearsonr(df.loc[mask, 'phi_over_c2'], df.loc[mask, 'sigma_residual'])
    
    print(f"\n  BEFORE correction:")
    print(f"    Correlation (Φ vs σ_resid): r = {r_orig:.4f}, p = {p_orig:.2e}")
    
    # Apply GR time dilation correction
    # σ_corrected = σ_observed × (1 + Φ/c²)
    # In log: log(σ_corrected) = log(σ_observed) + log(1 + Φ/c²) ≈ log(σ_observed) + Φ/c²/ln(10)
    
    df['log_sigma_gr_corrected'] = df['log_sigma'] + df['phi_over_c2'] / np.log(10)
    
    # Recompute residual
    slope = df['fj_slope'].iloc[0]
    intercept = df['fj_intercept'].iloc[0]
    df['sigma_residual_gr'] = df['log_sigma_gr_corrected'] - (intercept + slope * df['log_mass'])
    
    r_gr, p_gr = stats.pearsonr(df.loc[mask, 'phi_over_c2'], df.loc[mask, 'sigma_residual_gr'])
    
    print(f"\n  AFTER GR correction:")
    print(f"    Correlation (Φ vs σ_resid): r = {r_gr:.4f}, p = {p_gr:.2e}")
    
    # Try TEP-enhanced correction (larger effect)
    # Test different enhancement factors
    print(f"\n  Testing TEP enhancement factors:")
    
    best_factor = 1.0
    best_r = abs(r_gr)
    
    for factor in [10, 100, 1000, 10000, 100000, 1000000]:
        df['log_sigma_tep'] = df['log_sigma'] + factor * df['phi_over_c2'] / np.log(10)
        df['sigma_residual_tep'] = df['log_sigma_tep'] - (intercept + slope * df['log_mass'])
        
        r_tep, _ = stats.pearsonr(df.loc[mask, 'phi_over_c2'], df.loc[mask, 'sigma_residual_tep'])
        
        print(f"    Factor {factor:>7}: r = {r_tep:+.4f}")
        
        if abs(r_tep) < best_r:
            best_r = abs(r_tep)
            best_factor = factor
    
    print(f"\n  Best correction factor: {best_factor}×")
    print(f"  (GR predicts factor = 1)")
    print(f"  Residual correlation after best correction: r = {best_r:.4f}")
    
    if best_factor > 100:
        print(f"\n  → TEP enhancement factor ~{best_factor}× needed to remove environment effect")
        print(f"  → This is {best_factor}× larger than GR prediction!")
    
    return {
        'r_original': float(r_orig),
        'r_gr_corrected': float(r_gr),
        'best_tep_factor': float(best_factor),
        'r_best_corrected': float(best_r),
    }


def test_redshift_scaling(df):
    """
    Test if the σ-potential relationship changes with redshift.
    
    TEP prediction: Effect should be CONSTANT with z (it's a local measurement)
    Standard physics: Effect might evolve with z (assembly history changes)
    """
    print("\n" + "=" * 70)
    print("REDSHIFT SCALING TEST")
    print("=" * 70)
    
    z_bins = [(0.01, 0.10), (0.10, 0.20), (0.20, 0.35), (0.35, 0.55)]
    
    results = []
    
    for z_min, z_max in z_bins:
        z_mask = (df['redshift'] >= z_min) & (df['redshift'] < z_max)
        df_z = df[z_mask]
        
        if len(df_z) < 5000:
            continue
        
        mask = np.isfinite(df_z['phi_over_c2']) & np.isfinite(df_z['sigma_residual'])
        r, p = stats.pearsonr(df_z.loc[mask, 'phi_over_c2'], df_z.loc[mask, 'sigma_residual'])
        
        z_mid = (z_min + z_max) / 2
        print(f"  z = {z_min:.2f}-{z_max:.2f}: r = {r:.4f}, p = {p:.2e}")
        
        results.append({
            'z_min': z_min,
            'z_max': z_max,
            'z_mid': z_mid,
            'r': float(r),
            'p': float(p),
        })
    
    # Test for trend
    if len(results) >= 3:
        z_vals = [r['z_mid'] for r in results]
        r_vals = [r['r'] for r in results]
        
        slope, _, r_trend, p_trend, _ = stats.linregress(z_vals, r_vals)
        
        print(f"\n  Trend with z: slope = {slope:.4f}, p = {p_trend:.3f}")
        
        if p_trend < 0.1:
            if slope > 0:
                print("  → Effect STRENGTHENS with z (unexpected for TEP)")
            else:
                print("  → Effect WEAKENS with z (could indicate evolution)")
        else:
            print("  → No significant evolution with z (TEP-consistent)")
    
    return results


def create_visualization(df, scaling_results, correction_results, z_results, output_path):
    """Create visualization."""
    import matplotlib.pyplot as plt
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    
    # 1. Sigma residual vs potential
    ax = axes[0, 0]
    
    if scaling_results:
        ax.errorbar(scaling_results['bin_phi'], scaling_results['bin_sigma'],
                   yerr=scaling_results['bin_sigma_err'], fmt='o', capsize=3)
        
        # Add fit line
        phi_fit = np.linspace(min(scaling_results['bin_phi']), max(scaling_results['bin_phi']), 100)
        sigma_fit = scaling_results['linear_intercept'] + scaling_results['linear_slope'] * phi_fit
        ax.plot(phi_fit, sigma_fit, 'r-', label=f'Fit: slope = {scaling_results["linear_slope"]:.0f}')
        
        ax.axhline(0, color='gray', linestyle='--', alpha=0.5)
        ax.set_xlabel('Gravitational Potential (Φ/c²)')
        ax.set_ylabel('σ Residual (dex)')
        ax.set_title('Velocity Dispersion vs Potential')
        ax.legend()
    
    # 2. Correction comparison
    ax = axes[0, 1]
    
    if correction_results:
        labels = ['Original', 'GR\nCorrected', f'TEP\n({correction_results["best_tep_factor"]:.0f}×)']
        values = [correction_results['r_original'], 
                  correction_results['r_gr_corrected'],
                  correction_results['r_best_corrected']]
        colors = ['red', 'orange', 'green']
        
        ax.bar(range(3), [abs(v) for v in values], color=colors, alpha=0.7)
        ax.set_xticks(range(3))
        ax.set_xticklabels(labels)
        ax.set_ylabel('|Correlation with Φ|')
        ax.set_title('Effect of Time Dilation Correction')
        ax.axhline(0, color='black')
    
    # 3. Redshift evolution
    ax = axes[1, 0]
    
    if z_results:
        z_vals = [r['z_mid'] for r in z_results]
        r_vals = [r['r'] for r in z_results]
        
        ax.plot(z_vals, r_vals, 'o-', markersize=10)
        ax.axhline(0, color='gray', linestyle='--', alpha=0.5)
        ax.set_xlabel('Redshift')
        ax.set_ylabel('Correlation (Φ vs σ_resid)')
        ax.set_title('Environment Effect vs Redshift')
    
    # 4. Summary
    ax = axes[1, 1]
    ax.axis('off')
    
    summary = """
TEP VELOCITY DISPERSION REFRAME

KEY INSIGHT: If time flows slower in deep potential
wells, observed σ is INFLATED relative to intrinsic σ.

What looks like "tidal heating" could be a
MEASUREMENT ARTIFACT of time dilation.

RESULTS:
"""
    
    if scaling_results:
        summary += f"""
1. σ-POTENTIAL SCALING:
   Observed slope: {scaling_results['linear_slope']:.0f}
   GR prediction: {scaling_results['gr_predicted_slope']:.3f}
   Ratio: {scaling_results['slope_ratio']:.0f}×
"""
    
    if correction_results:
        summary += f"""
2. TIME DILATION CORRECTION:
   Original |r|: {abs(correction_results['r_original']):.4f}
   After GR correction: {abs(correction_results['r_gr_corrected']):.4f}
   Best TEP factor: {correction_results['best_tep_factor']:.0f}×
   After TEP correction: {abs(correction_results['r_best_corrected']):.4f}
"""
    
    summary += """
INTERPRETATION:
The environment effect is ~10⁶× larger than GR predicts.
This could indicate:
1. TEP enhancement of time dilation
2. Other physics (tidal heating, assembly bias)
3. Both effects combined
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
    print("TEP VELOCITY DISPERSION REFRAME")
    print("=" * 70)
    print(f"Started: {datetime.now().isoformat()}")
    
    df = load_and_prepare_data()
    
    scaling_results = test_sigma_potential_scaling(df)
    correction_results = test_tep_time_dilation_correction(df)
    z_results = test_redshift_scaling(df)
    
    # Visualization
    fig_path = os.path.join(FIGURES_DIR, 'step_4_6_tep_velocity_reframe.png')
    create_visualization(df, scaling_results, correction_results, z_results, fig_path)
    
    # Save results
    results = {
        'metadata': {
            'timestamp': datetime.now().isoformat(),
            'n_galaxies': len(df),
        },
        'scaling': scaling_results,
        'correction': correction_results,
        'redshift_evolution': z_results,
    }
    
    output_path = os.path.join(RESULTS_DIR, 'step_4_6_tep_velocity_reframe.json')
    with open(output_path, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\nResults saved: {output_path}")
    
    # Final summary
    print("\n" + "=" * 70)
    print("FINAL SUMMARY")
    print("=" * 70)
    print(f"\nThe σ-environment effect is {scaling_results['slope_ratio']:.0f}× larger than GR predicts.")
    print(f"A TEP enhancement factor of ~{correction_results['best_tep_factor']:.0f}× would be needed")
    print("to explain this as pure time dilation.")
    print("\nThis is either:")
    print("  A) Evidence for TEP-enhanced time dilation")
    print("  B) Standard physics (tidal heating) dominating")
    print("  C) A combination of both effects")
    
    return results


if __name__ == '__main__':
    results = main()
