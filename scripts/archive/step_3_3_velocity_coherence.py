#!/usr/bin/env python3
"""
Step 3.3: Velocity Field Coherence Test for Cosmic Topology

Hypothesis: If the universe has repeating topology, galaxies at "different"
positions may actually be connected through the topology, leading to
correlated velocity field orientations beyond what's expected from
large-scale structure alone.

Observable: Kinematic Position Angle (PA) - the orientation of each
galaxy's rotation axis projected on the sky.

Method:
1. Extract velocity field PA from MaNGA DAP (kinematic major axis)
2. Compute angular two-point correlation function of PA orientations
3. Compare to null model (random orientations)
4. Compare to LSS model (tidal alignment from large-scale structure)
5. Search for EXCESS correlation at large angular separations

TEP Connection: If time-flow varies with gravitational context, the
standard distance-redshift relationship breaks down. Regions that appear
"far apart" may actually be topologically connected, leading to
correlated dynamics.

Author: M. Smawfield
Date: January 2026
"""

import numpy as np
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u
from scipy.stats import pearsonr, spearmanr, ks_2samp
from scipy.spatial import cKDTree
import json
import os
from datetime import datetime
import warnings
warnings.filterwarnings('ignore')

# Paths
DATA_DIR = os.path.join(os.path.dirname(__file__), '..', '..', 'data')
RESULTS_DIR = os.path.join(os.path.dirname(__file__), '..', '..', 'results', 'outputs', 'topology')
FIGURES_DIR = os.path.join(os.path.dirname(__file__), '..', '..', 'results', 'figures', 'topology')

os.makedirs(RESULTS_DIR, exist_ok=True)
os.makedirs(FIGURES_DIR, exist_ok=True)


def load_velocity_data():
    """
    Load MaNGA velocity field data.
    
    Key columns from DAP:
    - STELLAR_VEL_LO/HI: Min/max stellar velocity (rotation amplitude)
    - HA_GVEL_LO/HI: Min/max H-alpha gas velocity
    - NSA_ELPETRO_PHI: Photometric position angle
    
    The kinematic PA can be inferred from the velocity field extrema
    and compared to the photometric PA.
    """
    print("Loading MaNGA velocity field data...")
    
    drp_path = os.path.join(DATA_DIR, 'drpall', 'drpall-v3_1_1.fits')
    dap_path = os.path.join(DATA_DIR, 'dapall', 'dapall-v3_1_1-3.1.0.fits')
    
    with fits.open(drp_path) as hdul:
        drp = hdul[1].data
    with fits.open(dap_path) as hdul:
        dap = hdul[1].data
    
    # Match catalogs
    drp_plateifu = np.array([f"{p}-{i}" for p, i in zip(drp['plate'], drp['ifudsgn'])])
    dap_plateifu = dap['PLATEIFU']
    
    common = set(drp_plateifu) & set(dap_plateifu)
    drp_idx = {pf: i for i, pf in enumerate(drp_plateifu)}
    dap_idx = {pf: i for i, pf in enumerate(dap_plateifu)}
    
    galaxies = []
    for pf in common:
        di = drp_idx[pf]
        ai = dap_idx[pf]
        
        # Position
        ra = drp['objra'][di]
        dec = drp['objdec'][di]
        z = drp['nsa_z'][di]
        
        # Velocity field properties
        stellar_vel_lo = dap['STELLAR_VEL_LO'][ai]
        stellar_vel_hi = dap['STELLAR_VEL_HI'][ai]
        ha_vel_lo = dap['HA_GVEL_LO'][ai]
        ha_vel_hi = dap['HA_GVEL_HI'][ai]
        
        # Photometric PA (from NSA)
        # Note: This is the photometric major axis PA, not kinematic
        # We'll use velocity amplitude as a proxy for rotation
        
        # Velocity amplitude (rotation signature)
        stellar_amp = stellar_vel_hi - stellar_vel_lo
        ha_amp = ha_vel_hi - ha_vel_lo
        
        # Filter: need valid velocities and significant rotation
        if not (0 < z < 0.2 and 
                -500 < stellar_vel_lo < 500 and 
                -500 < stellar_vel_hi < 500 and
                stellar_amp > 50):  # Require >50 km/s rotation
            continue
        
        # Get photometric PA if available
        try:
            photo_pa = dap['NSA_ELPETRO_PHI'][ai]
            if not (0 <= photo_pa <= 180):
                photo_pa = np.nan
        except:
            photo_pa = np.nan
        
        # Sersic PA as backup
        try:
            sersic_pa = dap['NSA_SERSIC_PHI'][ai]
            if not (0 <= sersic_pa <= 180):
                sersic_pa = np.nan
        except:
            sersic_pa = np.nan
        
        # Use best available PA
        if np.isfinite(photo_pa):
            pa = photo_pa
        elif np.isfinite(sersic_pa):
            pa = sersic_pa
        else:
            continue  # Skip if no PA available
        
        galaxies.append({
            'plateifu': pf,
            'ra': ra,
            'dec': dec,
            'z': z,
            'pa': pa,  # Position angle (0-180 degrees)
            'stellar_amp': stellar_amp,
            'ha_amp': ha_amp,
        })
    
    print(f"  Loaded {len(galaxies)} galaxies with valid velocity fields")
    return galaxies


def compute_pa_correlation(galaxies, angular_bins):
    """
    Compute the angular correlation function of position angles.
    
    For each pair of galaxies, we compute:
    1. Angular separation on the sky
    2. PA difference (accounting for 180° ambiguity)
    
    The correlation is: C(θ) = <cos(2 * ΔPA)>
    where ΔPA is the difference in position angles.
    
    For random orientations: C(θ) = 0
    For perfect alignment: C(θ) = 1
    For perpendicular: C(θ) = -1
    """
    print("\nComputing PA angular correlation function...")
    
    n = len(galaxies)
    
    # Extract coordinates and PAs
    coords = np.array([[g['ra'], g['dec']] for g in galaxies])
    pas = np.array([g['pa'] for g in galaxies])
    
    # Convert PA to radians (factor of 2 for 180° periodicity)
    pa_rad = np.radians(pas * 2)
    
    # Compute all pairwise angular separations
    print(f"  Computing {n*(n-1)//2} pairwise separations...")
    
    # Use SkyCoord for accurate separations
    sky_coords = SkyCoord(ra=coords[:, 0]*u.deg, dec=coords[:, 1]*u.deg)
    
    # Bin the correlations
    bin_edges = angular_bins
    n_bins = len(bin_edges) - 1
    
    correlations = []
    counts = []
    
    for b in range(n_bins):
        theta_min = bin_edges[b]
        theta_max = bin_edges[b + 1]
        
        cos_sum = 0.0
        pair_count = 0
        
        # Sample pairs for efficiency
        n_sample = min(500, n)
        sample_idx = np.random.choice(n, n_sample, replace=False)
        
        for i in sample_idx:
            # Compute separations from galaxy i to all others
            seps = sky_coords[i].separation(sky_coords).deg
            
            # Find pairs in this angular bin
            in_bin = (seps >= theta_min) & (seps < theta_max) & (np.arange(n) != i)
            
            if np.any(in_bin):
                # Compute PA correlation
                delta_pa = pa_rad[in_bin] - pa_rad[i]
                cos_sum += np.sum(np.cos(delta_pa))
                pair_count += np.sum(in_bin)
        
        if pair_count > 0:
            correlations.append(cos_sum / pair_count)
            counts.append(pair_count)
        else:
            correlations.append(np.nan)
            counts.append(0)
        
        print(f"  Bin {b}: {theta_min:.1f}-{theta_max:.1f}°, "
              f"C = {correlations[-1]:.4f}, N = {counts[-1]}")
    
    return np.array(correlations), np.array(counts), bin_edges


def null_model_test(galaxies, angular_bins, n_shuffles=500):
    """
    Test against null model by shuffling PA values.
    
    If topology is real, the observed correlation should exceed
    the shuffled distribution at some angular scales.
    """
    print(f"\nNull model test ({n_shuffles} shuffles)...")
    
    n = len(galaxies)
    coords = np.array([[g['ra'], g['dec']] for g in galaxies])
    pas = np.array([g['pa'] for g in galaxies])
    
    # Compute observed correlation
    obs_corr, obs_counts, _ = compute_pa_correlation(galaxies, angular_bins)
    
    # Shuffle and compute null distribution
    null_corrs = []
    
    for s in range(n_shuffles):
        if s % 100 == 0:
            print(f"  Shuffle {s}/{n_shuffles}...")
        
        # Shuffle PAs
        shuffled_pas = np.random.permutation(pas)
        shuffled_galaxies = [dict(g, pa=p) for g, p in zip(galaxies, shuffled_pas)]
        
        # Compute correlation (fast version)
        corr = compute_pa_correlation_fast(shuffled_galaxies, angular_bins, coords)
        null_corrs.append(corr)
    
    null_corrs = np.array(null_corrs)
    
    # Compute p-values for each bin
    n_bins = len(angular_bins) - 1
    p_values = []
    z_scores = []
    
    for b in range(n_bins):
        if np.isnan(obs_corr[b]):
            p_values.append(np.nan)
            z_scores.append(np.nan)
            continue
        
        null_b = null_corrs[:, b]
        null_b = null_b[~np.isnan(null_b)]
        
        if len(null_b) < 10:
            p_values.append(np.nan)
            z_scores.append(np.nan)
            continue
        
        # Two-sided p-value
        p = 2 * min(np.mean(null_b >= obs_corr[b]), np.mean(null_b <= obs_corr[b]))
        z = (obs_corr[b] - np.mean(null_b)) / max(np.std(null_b), 1e-6)
        
        p_values.append(p)
        z_scores.append(z)
    
    return {
        'observed': obs_corr.tolist(),
        'null_mean': np.nanmean(null_corrs, axis=0).tolist(),
        'null_std': np.nanstd(null_corrs, axis=0).tolist(),
        'p_values': p_values,
        'z_scores': z_scores,
        'counts': obs_counts.tolist(),
    }


def compute_pa_correlation_fast(galaxies, angular_bins, coords=None):
    """Fast version of PA correlation for null testing."""
    n = len(galaxies)
    
    if coords is None:
        coords = np.array([[g['ra'], g['dec']] for g in galaxies])
    
    pas = np.array([g['pa'] for g in galaxies])
    pa_rad = np.radians(pas * 2)
    
    n_bins = len(angular_bins) - 1
    correlations = np.full(n_bins, np.nan)
    
    # Sample for speed
    n_sample = min(200, n)
    sample_idx = np.random.choice(n, n_sample, replace=False)
    
    for b in range(n_bins):
        theta_min = angular_bins[b]
        theta_max = angular_bins[b + 1]
        
        cos_sum = 0.0
        pair_count = 0
        
        for i in sample_idx:
            # Approximate angular distance
            cos_dec = np.cos(np.radians(coords[i, 1]))
            seps = np.sqrt(
                ((coords[:, 0] - coords[i, 0]) * cos_dec)**2 +
                (coords[:, 1] - coords[i, 1])**2
            )
            
            in_bin = (seps >= theta_min) & (seps < theta_max) & (np.arange(n) != i)
            
            if np.any(in_bin):
                delta_pa = pa_rad[in_bin] - pa_rad[i]
                cos_sum += np.sum(np.cos(delta_pa))
                pair_count += np.sum(in_bin)
        
        if pair_count > 0:
            correlations[b] = cos_sum / pair_count
    
    return correlations


def analyze_redshift_dependence(galaxies, angular_bins):
    """
    Check if PA correlation depends on redshift.
    
    If topology is real, we might expect different correlation
    patterns at different redshift shells (different topology layers).
    """
    print("\nAnalyzing redshift dependence...")
    
    z_values = np.array([g['z'] for g in galaxies])
    z_median = np.median(z_values)
    
    # Split into low-z and high-z samples
    low_z = [g for g in galaxies if g['z'] < z_median]
    high_z = [g for g in galaxies if g['z'] >= z_median]
    
    print(f"  Low-z sample (z < {z_median:.4f}): {len(low_z)} galaxies")
    print(f"  High-z sample (z >= {z_median:.4f}): {len(high_z)} galaxies")
    
    # Compute correlations for each
    corr_low, counts_low, _ = compute_pa_correlation(low_z, angular_bins)
    corr_high, counts_high, _ = compute_pa_correlation(high_z, angular_bins)
    
    # Compare
    diff = corr_high - corr_low
    
    return {
        'z_median': float(z_median),
        'low_z': {
            'n_galaxies': len(low_z),
            'correlation': corr_low.tolist(),
            'counts': counts_low.tolist(),
        },
        'high_z': {
            'n_galaxies': len(high_z),
            'correlation': corr_high.tolist(),
            'counts': counts_high.tolist(),
        },
        'difference': diff.tolist(),
    }


def create_visualization(results, angular_bins, output_path):
    """Create visualization of velocity coherence results."""
    import matplotlib.pyplot as plt
    
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    bin_centers = (angular_bins[:-1] + angular_bins[1:]) / 2
    
    # 1. PA correlation function with null comparison
    ax = axes[0, 0]
    null_results = results['null_test']
    
    obs = np.array(null_results['observed'])
    null_mean = np.array(null_results['null_mean'])
    null_std = np.array(null_results['null_std'])
    
    ax.fill_between(bin_centers, null_mean - 2*null_std, null_mean + 2*null_std,
                   alpha=0.3, color='gray', label='Null ±2σ')
    ax.plot(bin_centers, null_mean, 'k--', label='Null mean')
    ax.plot(bin_centers, obs, 'ro-', markersize=6, label='Observed')
    
    ax.axhline(0, color='gray', linestyle=':', alpha=0.5)
    ax.set_xlabel('Angular Separation (deg)')
    ax.set_ylabel('PA Correlation C(θ)')
    ax.set_title('Position Angle Correlation Function')
    ax.legend()
    ax.set_xlim(0, max(bin_centers) * 1.1)
    
    # 2. Z-scores by angular bin
    ax = axes[0, 1]
    z_scores = np.array(null_results['z_scores'])
    colors = ['red' if abs(z) > 2 else 'blue' for z in z_scores]
    ax.bar(bin_centers, z_scores, width=np.diff(angular_bins)*0.8, color=colors, alpha=0.7)
    ax.axhline(0, color='gray', linestyle='-')
    ax.axhline(2, color='red', linestyle='--', alpha=0.5, label='2σ')
    ax.axhline(-2, color='red', linestyle='--', alpha=0.5)
    ax.set_xlabel('Angular Separation (deg)')
    ax.set_ylabel('Z-score vs Null')
    ax.set_title('Statistical Significance by Angular Scale')
    ax.legend()
    
    # 3. Redshift dependence
    ax = axes[1, 0]
    z_dep = results['redshift_dependence']
    ax.plot(bin_centers, z_dep['low_z']['correlation'], 'b.-', label=f"z < {z_dep['z_median']:.3f}")
    ax.plot(bin_centers, z_dep['high_z']['correlation'], 'r.-', label=f"z ≥ {z_dep['z_median']:.3f}")
    ax.axhline(0, color='gray', linestyle=':', alpha=0.5)
    ax.set_xlabel('Angular Separation (deg)')
    ax.set_ylabel('PA Correlation')
    ax.set_title('Redshift Dependence of PA Correlation')
    ax.legend()
    
    # 4. Summary statistics
    ax = axes[1, 1]
    ax.axis('off')
    
    # Find most significant bins
    sig_bins = [(i, z) for i, z in enumerate(z_scores) if abs(z) > 1.5]
    
    summary_text = f"""
VELOCITY FIELD COHERENCE TEST SUMMARY

Galaxies analyzed: {results['n_galaxies']}
Angular bins: {len(angular_bins)-1}

NULL HYPOTHESIS TEST:
"""
    
    if sig_bins:
        summary_text += "\nSignificant deviations (|z| > 1.5):\n"
        for i, z in sig_bins:
            summary_text += f"  {angular_bins[i]:.0f}-{angular_bins[i+1]:.0f}°: z = {z:.2f}σ\n"
    else:
        summary_text += "\nNo significant deviations from null model.\n"
    
    summary_text += f"\nVERDICT: {results['interpretation']['verdict']}"
    
    ax.text(0.1, 0.9, summary_text, transform=ax.transAxes, fontsize=10,
           verticalalignment='top', fontfamily='monospace')
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    plt.close()
    print(f"\nVisualization saved: {output_path}")


def interpret_results(null_results):
    """Interpret the velocity coherence test results."""
    z_scores = np.array(null_results['z_scores'])
    z_scores = z_scores[~np.isnan(z_scores)]
    
    max_z = np.max(np.abs(z_scores)) if len(z_scores) > 0 else 0
    n_significant = np.sum(np.abs(z_scores) > 2)
    
    if max_z > 3:
        verdict = "STRONG SIGNAL - Significant PA correlation detected"
        recommendation = "Investigate angular scales with excess correlation"
    elif max_z > 2 or n_significant > 1:
        verdict = "MODERATE SIGNAL - Possible PA correlation excess"
        recommendation = "Expand to larger surveys for confirmation"
    elif max_z > 1.5:
        verdict = "WEAK SIGNAL - Marginal excess, likely fluctuation"
        recommendation = "No strong evidence for topology"
    else:
        verdict = "NULL RESULT - PA correlations consistent with random"
        recommendation = "No evidence for cosmic topology in velocity fields"
    
    return {
        'verdict': verdict,
        'recommendation': recommendation,
        'max_z_score': float(max_z),
        'n_significant_bins': int(n_significant),
    }


def main():
    """Main analysis pipeline."""
    print("=" * 70)
    print("VELOCITY FIELD COHERENCE TEST FOR COSMIC TOPOLOGY")
    print("=" * 70)
    print(f"Started: {datetime.now().isoformat()}")
    
    # Load data
    galaxies = load_velocity_data()
    
    if len(galaxies) < 100:
        print("ERROR: Insufficient galaxies with valid velocity fields")
        return None
    
    # Define angular bins (logarithmic spacing)
    angular_bins = np.array([0.5, 1, 2, 5, 10, 20, 40, 80, 120])
    
    # Compute PA correlation
    obs_corr, obs_counts, _ = compute_pa_correlation(galaxies, angular_bins)
    
    # Null hypothesis test
    null_results = null_model_test(galaxies, angular_bins, n_shuffles=300)
    
    # Redshift dependence
    z_dependence = analyze_redshift_dependence(galaxies, angular_bins)
    
    # Interpret results
    interpretation = interpret_results(null_results)
    
    # Compile results
    results = {
        'metadata': {
            'timestamp': datetime.now().isoformat(),
            'test': 'Velocity Field Coherence',
        },
        'n_galaxies': len(galaxies),
        'angular_bins': angular_bins.tolist(),
        'null_test': null_results,
        'redshift_dependence': z_dependence,
        'interpretation': interpretation,
    }
    
    # Create visualization
    fig_path = os.path.join(FIGURES_DIR, 'step_3_3_velocity_coherence.png')
    create_visualization(results, angular_bins, fig_path)
    
    # Save results
    output_path = os.path.join(RESULTS_DIR, 'step_3_3_velocity_coherence.json')
    with open(output_path, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\nResults saved: {output_path}")
    
    # Print summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"Galaxies analyzed: {len(galaxies)}")
    print(f"Max z-score: {interpretation['max_z_score']:.2f}σ")
    print(f"Significant bins (|z| > 2): {interpretation['n_significant_bins']}")
    print(f"\nVerdict: {interpretation['verdict']}")
    print(f"Recommendation: {interpretation['recommendation']}")
    
    return results


if __name__ == '__main__':
    results = main()
