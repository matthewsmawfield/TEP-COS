#!/usr/bin/env python3
"""
Step 4.5: TEP Isochrony Test

CRITICAL REFRAMING: Standard astrophysics assumes isochrony (uniform time flow).
If this axiom is WRONG, then:

1. "Tidal heating" in clusters could actually be TIME DILATION
   - Galaxies in deep potential wells experience slower time
   - Their internal clocks (stellar evolution, dynamics) run slower
   - They appear "younger" or "less evolved" than their redshift suggests

2. The environment effects we detected are NOT noise - they ARE the signal

3. Key TEP prediction: The effect should scale with GRAVITATIONAL POTENTIAL,
   not just local density. We need to test:
   - Does the effect correlate with estimated potential depth?
   - Does it show the correct sign (slower time = less evolved)?
   - Does it have the right magnitude for GR time dilation?

Author: M. Smawfield
Date: January 2026
"""

import numpy as np
import pandas as pd
from scipy import stats
from scipy.spatial import cKDTree
from astropy.cosmology import FlatLambdaCDM
from astropy import constants as const
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
    """Load SDSS data."""
    print("Loading SDSS data...")
    df = pd.read_csv(os.path.join(DATA_DIR, 'sdss_galaxies.csv'))
    print(f"  Loaded {len(df):,} galaxies")
    return df


def estimate_gravitational_potential(df, n_neighbors=20):
    """
    Estimate local gravitational potential for each galaxy.
    
    Under TEP, time dilation scales with potential:
    Δτ/τ ≈ ΔΦ/c²
    
    We estimate Φ from the local mass density:
    Φ ≈ -G * M_enclosed / R
    
    Where M_enclosed is estimated from neighbor count and typical mass.
    """
    print("\nEstimating gravitational potential...")
    
    # Build 3D tree using RA, Dec, and comoving distance
    z = df['redshift'].values
    d_comoving = cosmo.comoving_distance(z).value  # Mpc
    
    # Convert to Cartesian (approximate for small angles)
    ra_rad = np.radians(df['ra'].values)
    dec_rad = np.radians(df['dec'].values)
    
    x = d_comoving * np.cos(dec_rad) * np.cos(ra_rad)
    y = d_comoving * np.cos(dec_rad) * np.sin(ra_rad)
    z_coord = d_comoving * np.sin(dec_rad)
    
    coords_3d = np.column_stack([x, y, z_coord])
    
    # Build KD-tree
    tree = cKDTree(coords_3d)
    
    # Query for neighbors
    distances, indices = tree.query(coords_3d, k=n_neighbors + 1)
    
    # Distance to Nth neighbor (in Mpc)
    r_n = distances[:, -1]
    
    # Estimate enclosed mass from neighbor count
    # Assume average galaxy mass ~ 10^11 M_sun
    M_avg = 1e11  # M_sun
    M_enclosed = n_neighbors * M_avg  # M_sun
    
    # Gravitational potential: Φ = -G*M/R
    # In units where c=1, Φ/c² is dimensionless
    G = 4.302e-6  # kpc/M_sun * (km/s)^2
    c = 299792.458  # km/s
    
    # Convert r_n from Mpc to kpc
    r_kpc = r_n * 1000
    
    # Potential in units of c²
    phi_over_c2 = -G * M_enclosed / (r_kpc * c**2)
    
    df['phi_over_c2'] = phi_over_c2
    df['log_phi'] = np.log10(-phi_over_c2 + 1e-12)
    
    # Also compute local density
    df['local_density'] = n_neighbors / (4/3 * np.pi * r_n**3)
    df['log_density'] = np.log10(df['local_density'] + 1e-10)
    
    print(f"  Potential range: {phi_over_c2.min():.2e} to {phi_over_c2.max():.2e} (Φ/c²)")
    print(f"  Expected time dilation: {-phi_over_c2.max()*1e9:.1f} ns/s in deepest wells")
    
    return df


def compute_apparent_age_proxy(df):
    """
    Compute a proxy for "apparent age" of each galaxy.
    
    Under standard physics: age depends on redshift and formation time
    Under TEP: age also depends on gravitational environment
    
    Proxies for age/evolution:
    1. D4000 break strength (not available, use color proxy)
    2. Specific SFR (sSFR = SFR/M) - lower = older
    3. Concentration - higher = more evolved
    4. Velocity dispersion at fixed mass - higher = more relaxed = older
    """
    print("\nComputing apparent age proxies...")
    
    # 1. Color-based age proxy (redder = older)
    # g-r color correlates with stellar population age
    if 'petroMag_g' in df.columns and 'petroMag_r' in df.columns:
        df['g_r_color'] = df['petroMag_g'] - df['petroMag_r']
    else:
        # Use SFR as proxy (lower SFR = redder = older)
        df['g_r_color'] = -df['log_sfr'] * 0.3 + 0.7  # Rough conversion
    
    # 2. Specific SFR (log scale, more negative = older)
    df['log_ssfr'] = df['log_sfr'] - df['log_mass']
    
    # 3. Concentration already computed
    
    # 4. Sigma residual at fixed mass (Faber-Jackson residual)
    # Fit FJ relation
    mask = np.isfinite(df['log_mass']) & np.isfinite(df['log_sigma'])
    slope, intercept, _, _, _ = stats.linregress(
        df.loc[mask, 'log_mass'], df.loc[mask, 'log_sigma']
    )
    df['sigma_residual'] = df['log_sigma'] - (intercept + slope * df['log_mass'])
    
    # Combined age proxy (normalized)
    # Higher value = appears older
    df['age_proxy'] = (
        (df['g_r_color'] - df['g_r_color'].mean()) / df['g_r_color'].std() +
        (-df['log_ssfr'] - (-df['log_ssfr']).mean()) / df['log_ssfr'].std() +
        (df['concentration'] - df['concentration'].mean()) / df['concentration'].std() +
        (df['sigma_residual'] - df['sigma_residual'].mean()) / df['sigma_residual'].std()
    ) / 4
    
    print(f"  Age proxy range: {df['age_proxy'].min():.2f} to {df['age_proxy'].max():.2f}")
    
    return df


def test_tep_time_dilation(df):
    """
    Test the core TEP prediction:
    
    If time flows slower in deep gravitational wells, galaxies there
    should appear YOUNGER (less evolved) than their redshift suggests.
    
    Standard physics predicts the OPPOSITE:
    - Cluster galaxies formed earlier → should be OLDER
    - Assembly bias → cluster galaxies are more evolved
    
    TEP prediction: age_proxy should DECREASE with deeper potential (more negative Φ)
    Standard prediction: age_proxy should INCREASE with deeper potential
    """
    print("\n" + "=" * 70)
    print("TEP TIME DILATION TEST")
    print("=" * 70)
    
    # Bin by gravitational potential
    phi_percentiles = np.percentile(df['phi_over_c2'], [10, 30, 50, 70, 90])
    
    print("\nAge proxy vs gravitational potential:")
    print("(TEP predicts: deeper potential → younger appearance)")
    print("(Standard predicts: deeper potential → older appearance)")
    
    results = {'bins': []}
    
    phi_bins = [-np.inf] + list(phi_percentiles) + [0]
    labels = ['Deepest', 'Deep', 'Medium-Deep', 'Medium', 'Shallow', 'Shallowest']
    
    for i in range(len(phi_bins) - 1):
        mask = (df['phi_over_c2'] >= phi_bins[i]) & (df['phi_over_c2'] < phi_bins[i+1])
        if mask.sum() > 100:
            age = df.loc[mask, 'age_proxy'].mean()
            age_err = df.loc[mask, 'age_proxy'].std() / np.sqrt(mask.sum())
            phi_mean = df.loc[mask, 'phi_over_c2'].mean()
            
            print(f"  {labels[i]:12s} (Φ/c² ~ {phi_mean:.2e}): age = {age:+.4f} ± {age_err:.4f}")
            
            results['bins'].append({
                'label': labels[i],
                'phi_mean': float(phi_mean),
                'age_mean': float(age),
                'age_err': float(age_err),
                'n': int(mask.sum()),
            })
    
    # Correlation test
    mask = np.isfinite(df['phi_over_c2']) & np.isfinite(df['age_proxy'])
    r, p = stats.pearsonr(df.loc[mask, 'phi_over_c2'], df.loc[mask, 'age_proxy'])
    
    print(f"\n  Correlation (Φ vs age): r = {r:.4f}, p = {p:.2e}")
    
    results['correlation'] = {'r': float(r), 'p': float(p)}
    
    # Interpretation
    print("\n  INTERPRETATION:")
    if r > 0:
        print("  → Deeper potential correlates with YOUNGER appearance")
        print("  → This is CONSISTENT with TEP time dilation!")
        results['interpretation'] = 'TEP_CONSISTENT'
    else:
        print("  → Deeper potential correlates with OLDER appearance")
        print("  → This matches standard physics (assembly bias)")
        print("  → BUT: Could be masking a TEP signal underneath")
        results['interpretation'] = 'STANDARD_DOMINANT'
    
    return results


def test_redshift_dependent_dilation(df):
    """
    Test if the potential-age relationship changes with redshift.
    
    TEP prediction: The effect should be STRONGER at higher z because:
    1. More time has passed for dilation to accumulate
    2. Structures were denser in the past
    
    Standard physics: Assembly bias should be roughly constant with z
    """
    print("\n" + "=" * 70)
    print("REDSHIFT-DEPENDENT TIME DILATION TEST")
    print("=" * 70)
    
    z_bins = [(0.01, 0.08), (0.08, 0.15), (0.15, 0.25), (0.25, 0.40), (0.40, 0.60)]
    
    results = []
    
    for z_min, z_max in z_bins:
        z_mask = (df['redshift'] >= z_min) & (df['redshift'] < z_max)
        df_z = df[z_mask]
        
        if len(df_z) < 5000:
            continue
        
        # Split by potential
        phi_median = df_z['phi_over_c2'].median()
        deep_mask = df_z['phi_over_c2'] < phi_median
        shallow_mask = df_z['phi_over_c2'] >= phi_median
        
        age_deep = df_z.loc[deep_mask, 'age_proxy'].mean()
        age_shallow = df_z.loc[shallow_mask, 'age_proxy'].mean()
        
        age_diff = age_deep - age_shallow
        
        # Error on difference
        err_deep = df_z.loc[deep_mask, 'age_proxy'].std() / np.sqrt(deep_mask.sum())
        err_shallow = df_z.loc[shallow_mask, 'age_proxy'].std() / np.sqrt(shallow_mask.sum())
        age_diff_err = np.sqrt(err_deep**2 + err_shallow**2)
        
        z_mid = (z_min + z_max) / 2
        t_lookback = cosmo.lookback_time(z_mid).value
        
        print(f"  z = {z_min:.2f}-{z_max:.2f} (t = {t_lookback:.1f} Gyr): "
              f"Δage = {age_diff:+.4f} ± {age_diff_err:.4f}")
        
        results.append({
            'z_min': z_min,
            'z_max': z_max,
            'z_mid': z_mid,
            't_lookback': t_lookback,
            'age_diff': age_diff,
            'age_diff_err': age_diff_err,
        })
    
    # Test for trend with redshift
    if len(results) >= 3:
        z_vals = [r['z_mid'] for r in results]
        diff_vals = [r['age_diff'] for r in results]
        
        slope, intercept, r_value, p_value, std_err = stats.linregress(z_vals, diff_vals)
        
        print(f"\n  Trend with redshift: slope = {slope:.4f} ± {std_err:.4f}")
        print(f"  Correlation: r = {r_value:.3f}, p = {p_value:.3f}")
        
        if slope > 0 and p_value < 0.1:
            print("\n  → Effect STRENGTHENS with redshift")
            print("  → This is CONSISTENT with TEP cumulative dilation!")
        elif slope < 0 and p_value < 0.1:
            print("\n  → Effect WEAKENS with redshift")
            print("  → Inconsistent with simple TEP prediction")
        else:
            print("\n  → No significant trend with redshift")
        
        return {
            'bins': results,
            'trend_slope': float(slope),
            'trend_slope_err': float(std_err),
            'trend_r': float(r_value),
            'trend_p': float(p_value),
        }
    
    return {'bins': results}


def test_mass_independent_effect(df):
    """
    Test if the potential-age effect persists at FIXED stellar mass.
    
    This is crucial because:
    - Standard physics: More massive galaxies are older (downsizing)
    - If potential effect persists at fixed mass, it's NOT just mass-driven
    """
    print("\n" + "=" * 70)
    print("MASS-INDEPENDENT POTENTIAL EFFECT")
    print("=" * 70)
    
    mass_bins = [(9.5, 10.0), (10.0, 10.5), (10.5, 11.0), (11.0, 11.5)]
    
    results = []
    
    for m_min, m_max in mass_bins:
        mass_mask = (df['log_mass'] >= m_min) & (df['log_mass'] < m_max)
        df_m = df[mass_mask]
        
        if len(df_m) < 2000:
            continue
        
        # Correlation within mass bin
        mask = np.isfinite(df_m['phi_over_c2']) & np.isfinite(df_m['age_proxy'])
        r, p = stats.pearsonr(df_m.loc[mask, 'phi_over_c2'], df_m.loc[mask, 'age_proxy'])
        
        print(f"  M = {m_min:.1f}-{m_max:.1f}: r = {r:+.4f}, p = {p:.2e} (n={mask.sum()})")
        
        results.append({
            'mass_min': m_min,
            'mass_max': m_max,
            'r': float(r),
            'p': float(p),
            'n': int(mask.sum()),
        })
    
    # Check if effect is consistent across mass bins
    r_vals = [r['r'] for r in results]
    mean_r = np.mean(r_vals)
    std_r = np.std(r_vals)
    
    print(f"\n  Mean correlation: r = {mean_r:.4f} ± {std_r:.4f}")
    
    if np.all(np.array(r_vals) > 0):
        print("  → Effect is POSITIVE in ALL mass bins")
        print("  → Potential-age correlation is mass-independent!")
    elif np.all(np.array(r_vals) < 0):
        print("  → Effect is NEGATIVE in ALL mass bins")
    else:
        print("  → Effect changes sign across mass bins")
    
    return results


def quantify_tep_magnitude(df):
    """
    Quantify the magnitude of the potential effect and compare to GR prediction.
    
    GR time dilation: Δτ/τ = ΔΦ/c²
    
    If we see an age difference of X% between deep and shallow potential,
    and the potential difference is ΔΦ/c², then:
    
    Implied time dilation = X% / (cosmic_time_at_z)
    
    Compare to GR prediction: ΔΦ/c²
    """
    print("\n" + "=" * 70)
    print("TEP MAGNITUDE QUANTIFICATION")
    print("=" * 70)
    
    # Split into quintiles by potential
    df['phi_quintile'] = pd.qcut(df['phi_over_c2'], 5, labels=['Q1_deep', 'Q2', 'Q3', 'Q4', 'Q5_shallow'])
    
    # Age difference between deepest and shallowest
    deep = df[df['phi_quintile'] == 'Q1_deep']
    shallow = df[df['phi_quintile'] == 'Q5_shallow']
    
    age_deep = deep['age_proxy'].mean()
    age_shallow = shallow['age_proxy'].mean()
    age_diff = age_deep - age_shallow
    
    phi_deep = deep['phi_over_c2'].mean()
    phi_shallow = shallow['phi_over_c2'].mean()
    phi_diff = phi_deep - phi_shallow
    
    print(f"  Deep potential (Q1): Φ/c² = {phi_deep:.2e}, age = {age_deep:+.4f}")
    print(f"  Shallow potential (Q5): Φ/c² = {phi_shallow:.2e}, age = {age_shallow:+.4f}")
    print(f"  Difference: ΔΦ/c² = {phi_diff:.2e}, Δage = {age_diff:+.4f}")
    
    # Convert age proxy to fractional age difference
    # Age proxy is in units of standard deviations
    # Typical scatter in galaxy ages is ~2-3 Gyr out of ~10 Gyr
    age_scatter_gyr = 2.5
    cosmic_age_gyr = 10.0
    
    fractional_age_diff = age_diff * age_scatter_gyr / cosmic_age_gyr
    
    print(f"\n  Implied fractional age difference: {fractional_age_diff:.4f} ({fractional_age_diff*100:.2f}%)")
    
    # GR prediction
    gr_prediction = phi_diff  # Δτ/τ = ΔΦ/c²
    
    print(f"  GR time dilation prediction: {gr_prediction:.2e} ({gr_prediction*100:.4f}%)")
    
    # Ratio
    if abs(gr_prediction) > 1e-15:
        ratio = fractional_age_diff / gr_prediction
        print(f"\n  Observed / GR prediction: {ratio:.1f}x")
        
        if ratio > 10:
            print("  → Effect is MUCH LARGER than GR predicts")
            print("  → Could indicate TEP enhancement or other physics")
        elif ratio > 0.1:
            print("  → Effect is within order of magnitude of GR")
            print("  → Potentially consistent with TEP")
        else:
            print("  → Effect is MUCH SMALLER than GR predicts")
    
    return {
        'phi_deep': float(phi_deep),
        'phi_shallow': float(phi_shallow),
        'phi_diff': float(phi_diff),
        'age_deep': float(age_deep),
        'age_shallow': float(age_shallow),
        'age_diff': float(age_diff),
        'fractional_age_diff': float(fractional_age_diff),
        'gr_prediction': float(gr_prediction),
    }


def create_visualization(df, dilation_results, z_results, mass_results, magnitude, output_path):
    """Create comprehensive visualization."""
    import matplotlib.pyplot as plt
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    
    # 1. Age proxy vs gravitational potential
    ax = axes[0, 0]
    
    # Bin and plot
    phi_bins = np.percentile(df['phi_over_c2'], np.linspace(0, 100, 21))
    bin_centers = []
    bin_ages = []
    bin_errs = []
    
    for i in range(len(phi_bins) - 1):
        mask = (df['phi_over_c2'] >= phi_bins[i]) & (df['phi_over_c2'] < phi_bins[i+1])
        if mask.sum() > 100:
            bin_centers.append((phi_bins[i] + phi_bins[i+1]) / 2)
            bin_ages.append(df.loc[mask, 'age_proxy'].mean())
            bin_errs.append(df.loc[mask, 'age_proxy'].std() / np.sqrt(mask.sum()))
    
    ax.errorbar(bin_centers, bin_ages, yerr=bin_errs, fmt='o-', capsize=3)
    ax.set_xlabel('Gravitational Potential (Φ/c²)')
    ax.set_ylabel('Age Proxy (σ units)')
    ax.set_title('Age vs Gravitational Potential\n(TEP: deeper → younger)')
    ax.axhline(0, color='gray', linestyle='--', alpha=0.5)
    
    # Add trend line
    if dilation_results and 'correlation' in dilation_results:
        r = dilation_results['correlation']['r']
        ax.text(0.05, 0.95, f'r = {r:.4f}', transform=ax.transAxes, fontsize=12,
               verticalalignment='top', fontweight='bold')
    
    # 2. Redshift evolution
    ax = axes[0, 1]
    if z_results and 'bins' in z_results and len(z_results['bins']) > 0:
        z_vals = [r['z_mid'] for r in z_results['bins']]
        diff_vals = [r['age_diff'] for r in z_results['bins']]
        diff_errs = [r['age_diff_err'] for r in z_results['bins']]
        
        ax.errorbar(z_vals, diff_vals, yerr=diff_errs, fmt='o-', capsize=5, markersize=8)
        ax.axhline(0, color='gray', linestyle='--', alpha=0.5)
        
        if 'trend_slope' in z_results:
            z_fit = np.linspace(min(z_vals), max(z_vals), 100)
            y_fit = z_results['trend_slope'] * z_fit + (np.mean(diff_vals) - z_results['trend_slope'] * np.mean(z_vals))
            ax.plot(z_fit, y_fit, 'r--', alpha=0.7)
    
    ax.set_xlabel('Redshift')
    ax.set_ylabel('Age Difference (Deep - Shallow)')
    ax.set_title('Potential Effect vs Redshift\n(TEP: should strengthen with z)')
    
    # 3. Mass-independent test
    ax = axes[1, 0]
    if mass_results:
        masses = [(r['mass_min'] + r['mass_max'])/2 for r in mass_results]
        r_vals = [r['r'] for r in mass_results]
        
        colors = ['green' if r > 0 else 'red' for r in r_vals]
        ax.bar(range(len(masses)), r_vals, color=colors, alpha=0.7)
        ax.set_xticks(range(len(masses)))
        ax.set_xticklabels([f'{m:.1f}' for m in masses])
        ax.axhline(0, color='black', linestyle='-', alpha=0.5)
        ax.set_xlabel('log(M*/M☉)')
        ax.set_ylabel('Correlation (Φ vs Age)')
        ax.set_title('Potential-Age Correlation by Mass\n(Green = TEP-consistent)')
    
    # 4. Summary
    ax = axes[1, 1]
    ax.axis('off')
    
    summary = """
TEP ISOCHRONY TEST SUMMARY

HYPOTHESIS: If isochrony is wrong, time flows slower
in deep gravitational wells. Galaxies there should
appear YOUNGER than their redshift suggests.

KEY RESULTS:
"""
    
    if dilation_results and 'correlation' in dilation_results:
        r = dilation_results['correlation']['r']
        p = dilation_results['correlation']['p']
        interp = dilation_results.get('interpretation', 'UNKNOWN')
        
        summary += f"""
1. POTENTIAL-AGE CORRELATION:
   r = {r:.4f}, p = {p:.2e}
   Interpretation: {interp}
"""
    
    if z_results and 'trend_slope' in z_results:
        summary += f"""
2. REDSHIFT EVOLUTION:
   Slope = {z_results['trend_slope']:.4f} ± {z_results['trend_slope_err']:.4f}
   p = {z_results['trend_p']:.3f}
"""
    
    if magnitude:
        summary += f"""
3. MAGNITUDE:
   Observed age diff: {magnitude['fractional_age_diff']*100:.2f}%
   GR prediction: {magnitude['gr_prediction']*100:.4f}%
"""
    
    # Overall verdict
    if dilation_results and dilation_results.get('interpretation') == 'TEP_CONSISTENT':
        summary += "\n\nVERDICT: DATA CONSISTENT WITH TEP"
    else:
        summary += "\n\nVERDICT: STANDARD PHYSICS DOMINANT\n(But may mask underlying TEP signal)"
    
    ax.text(0.05, 0.95, summary, transform=ax.transAxes, fontsize=10,
           verticalalignment='top', fontfamily='monospace')
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    plt.close()
    print(f"\nVisualization saved: {output_path}")


def main():
    """Main analysis."""
    print("=" * 70)
    print("TEP ISOCHRONY TEST: CHALLENGING THE FUNDAMENTAL AXIOM")
    print("=" * 70)
    print(f"Started: {datetime.now().isoformat()}")
    print("\nIf standard astrophysics is WRONG about isochrony,")
    print("then time flows slower in deep gravitational wells,")
    print("and galaxies there should appear YOUNGER than expected.")
    
    df = load_data()
    df = estimate_gravitational_potential(df)
    df = compute_apparent_age_proxy(df)
    
    dilation_results = test_tep_time_dilation(df)
    z_results = test_redshift_dependent_dilation(df)
    mass_results = test_mass_independent_effect(df)
    magnitude = quantify_tep_magnitude(df)
    
    # Visualization
    fig_path = os.path.join(FIGURES_DIR, 'step_4_5_tep_isochrony.png')
    create_visualization(df, dilation_results, z_results, mass_results, magnitude, fig_path)
    
    # Save results
    results = {
        'metadata': {
            'timestamp': datetime.now().isoformat(),
            'n_galaxies': len(df),
            'hypothesis': 'If isochrony is wrong, deeper potential → younger appearance',
        },
        'time_dilation': dilation_results,
        'redshift_evolution': z_results,
        'mass_independent': mass_results,
        'magnitude': magnitude,
    }
    
    output_path = os.path.join(RESULTS_DIR, 'step_4_5_tep_isochrony.json')
    with open(output_path, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\nResults saved: {output_path}")
    
    return results


if __name__ == '__main__':
    results = main()
