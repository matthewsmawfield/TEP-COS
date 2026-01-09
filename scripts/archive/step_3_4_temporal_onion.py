#!/usr/bin/env python3
"""
Step 3.4: Temporal Onion Hypothesis - Evolved Pattern Matching

HYPOTHESIS:
===========
The universe is layered like an onion, where each "layer" represents the same
structures at different epochs. When we look at different redshifts, we may be
seeing the SAME clusters/filaments/galaxies, but:

1. At different positions (they've moved due to cosmic expansion + peculiar velocities)
2. At different evolutionary stages (they've aged, merged, formed stars)
3. Possibly rotated or distorted

This means we might be "double-counting" structures - seeing the same galaxy
cluster at z=0.05 and z=0.10, not realizing they're the same object viewed
at different times through the temporal onion.

KEY INSIGHT:
============
Standard topology tests look for IDENTICAL patterns. This test looks for
EVOLVED patterns - the same fingerprint, but aged and transformed.

METHOD:
=======
1. For each galaxy/group at redshift z1, predict what it would look like
   at redshift z2 after cosmic evolution
2. Search for matches between observed(z2) and predicted_evolution(z1 → z2)
3. Allow for spatial displacement (the structure has moved)
4. If the universe is a temporal onion, we expect EXCESS matches

EVOLUTION MODEL:
================
- Mass: Grows via mergers and accretion (M(t) ~ M0 * (1 + growth_rate * Δt))
- Velocity dispersion: Scales with mass (σ ~ M^0.25, Faber-Jackson)
- Star formation: Declines with cosmic time (SFR ~ exp(-t/τ))
- Morphology: Evolves from disk to spheroid (Sersic n increases)

Author: M. Smawfield
Date: January 2026
"""

import numpy as np
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
from scipy.spatial import cKDTree
from scipy.stats import pearsonr, ks_2samp
from sklearn.preprocessing import StandardScaler
import json
import os
from datetime import datetime
import warnings
warnings.filterwarnings('ignore')

# Cosmology
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

# Paths
DATA_DIR = os.path.join(os.path.dirname(__file__), '..', '..', 'data')
RESULTS_DIR = os.path.join(os.path.dirname(__file__), '..', '..', 'results', 'outputs', 'topology')
FIGURES_DIR = os.path.join(os.path.dirname(__file__), '..', '..', 'results', 'figures', 'topology')

os.makedirs(RESULTS_DIR, exist_ok=True)
os.makedirs(FIGURES_DIR, exist_ok=True)


def load_galaxy_data():
    """Load MaNGA galaxy data with properties needed for evolution modeling."""
    print("Loading MaNGA galaxy data...")
    
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
        
        z = drp['nsa_z'][di]
        mass = drp['nsa_sersic_mass'][di]
        sigma = dap['STELLAR_SIGMA_1RE'][ai]
        sfr = dap['SFR_TOT'][ai]
        sersic_n = dap['NSA_SERSIC_N'][ai]
        sersic_ba = dap['NSA_SERSIC_BA'][ai]
        
        # Filter valid
        if not (0.01 < z < 0.15 and mass > 1e8 and 
                10 < sigma < 400 and -5 < sfr < 50 and 
                0.5 < sersic_n < 8 and 0.1 < sersic_ba < 1):
            continue
        
        # Lookback time
        t_lookback = cosmo.lookback_time(z).value  # Gyr
        
        galaxies.append({
            'plateifu': pf,
            'ra': drp['objra'][di],
            'dec': drp['objdec'][di],
            'z': z,
            't_lookback': t_lookback,
            'log_mass': np.log10(mass),
            'log_sigma': np.log10(sigma),
            'log_sfr': np.log10(max(sfr, 0.01)),
            'sersic_n': sersic_n,
            'sersic_ba': sersic_ba,
        })
    
    print(f"  Loaded {len(galaxies)} valid galaxies")
    return galaxies


def evolve_galaxy(galaxy, delta_t_gyr):
    """
    Predict galaxy properties after evolution by delta_t Gyr.
    
    This is a simplified evolution model based on observed scaling relations.
    Positive delta_t means forward in time (galaxy gets older).
    Negative delta_t means backward in time (galaxy was younger).
    
    Evolution rules:
    - Mass: Grows ~5% per Gyr (mergers + star formation)
    - Sigma: Scales as M^0.25 (Faber-Jackson)
    - SFR: Declines exponentially with τ ~ 3 Gyr
    - Sersic n: Increases ~0.1 per Gyr (disk → spheroid)
    - Axis ratio: Roughly constant (geometry preserved)
    """
    
    # Mass growth
    growth_rate = 0.05  # 5% per Gyr
    new_log_mass = galaxy['log_mass'] + np.log10(1 + growth_rate * delta_t_gyr)
    
    # Sigma follows mass (Faber-Jackson: σ ~ M^0.25)
    delta_log_mass = new_log_mass - galaxy['log_mass']
    new_log_sigma = galaxy['log_sigma'] + 0.25 * delta_log_mass
    
    # SFR declines exponentially
    tau_sfr = 3.0  # e-folding time in Gyr
    new_log_sfr = galaxy['log_sfr'] - delta_t_gyr / tau_sfr / np.log(10)
    new_log_sfr = max(new_log_sfr, -3)  # Floor at 0.001 Msun/yr
    
    # Sersic n increases (morphological evolution)
    sersic_rate = 0.1  # per Gyr
    new_sersic_n = galaxy['sersic_n'] + sersic_rate * delta_t_gyr
    new_sersic_n = np.clip(new_sersic_n, 0.5, 8)
    
    # Axis ratio roughly preserved
    new_sersic_ba = galaxy['sersic_ba']
    
    return {
        'log_mass': new_log_mass,
        'log_sigma': new_log_sigma,
        'log_sfr': new_log_sfr,
        'sersic_n': new_sersic_n,
        'sersic_ba': new_sersic_ba,
    }


def create_fingerprint(galaxy):
    """Create a 5D fingerprint from galaxy properties."""
    return np.array([
        galaxy['log_mass'],
        galaxy['log_sigma'],
        galaxy['log_sfr'],
        galaxy['sersic_n'],
        galaxy['sersic_ba'],
    ])


def find_evolved_matches(galaxies, z_bins, evolution_tolerance=1.5, min_angular_sep=30.0):
    """
    Search for matches between galaxies at different redshifts,
    accounting for cosmic evolution.
    
    For each galaxy at z1, we:
    1. Predict what it would look like at z2 (evolved forward/backward)
    2. Search for galaxies at z2 that match the prediction
    3. A match suggests we're seeing the same structure at different epochs
    
    evolution_tolerance: How many sigma deviation from predicted evolution
    is allowed for a match (accounts for model uncertainty)
    """
    print(f"\nSearching for evolved matches across {len(z_bins)} redshift bins...")
    
    # Organize galaxies by redshift bin
    z_values = np.array([g['z'] for g in galaxies])
    bin_indices = []
    
    for i, (z_min, z_max) in enumerate(z_bins):
        mask = (z_values >= z_min) & (z_values < z_max)
        bin_indices.append(np.where(mask)[0])
        print(f"  Bin {i}: z = {z_min:.3f} - {z_max:.3f}, N = {len(bin_indices[-1])}")
    
    # Compute fingerprints for all galaxies
    fingerprints = np.array([create_fingerprint(g) for g in galaxies])
    
    # Normalize fingerprints
    scaler = StandardScaler()
    fingerprints_norm = scaler.fit_transform(fingerprints)
    
    # Search for evolved matches across bin pairs
    matches = []
    
    for i in range(len(z_bins)):
        for j in range(i + 1, len(z_bins)):
            z_i = (z_bins[i][0] + z_bins[i][1]) / 2
            z_j = (z_bins[j][0] + z_bins[j][1]) / 2
            
            # Time difference between bins
            t_i = cosmo.lookback_time(z_i).value
            t_j = cosmo.lookback_time(z_j).value
            delta_t = t_i - t_j  # Positive if bin_i is older (higher z)
            
            print(f"\n  Comparing bin {i} (z~{z_i:.3f}) to bin {j} (z~{z_j:.3f})")
            print(f"    Time difference: {abs(delta_t):.2f} Gyr")
            
            # For each galaxy in bin_i, predict evolved state and search in bin_j
            for idx_i in bin_indices[i]:
                gal_i = galaxies[idx_i]
                
                # Predict what this galaxy would look like at z_j
                evolved = evolve_galaxy(gal_i, -delta_t)  # Negative because going to younger epoch
                evolved_fp = create_fingerprint(evolved)
                evolved_fp_norm = (evolved_fp - scaler.mean_) / scaler.scale_
                
                # Find matches in bin_j
                for idx_j in bin_indices[j]:
                    fp_j = fingerprints_norm[idx_j]
                    
                    # Distance in normalized fingerprint space
                    dist = np.linalg.norm(fp_j - evolved_fp_norm)
                    
                    if dist < evolution_tolerance:
                        gal_j = galaxies[idx_j]
                        
                        # Angular separation
                        coord_i = SkyCoord(ra=gal_i['ra']*u.deg, dec=gal_i['dec']*u.deg)
                        coord_j = SkyCoord(ra=gal_j['ra']*u.deg, dec=gal_j['dec']*u.deg)
                        ang_sep = coord_i.separation(coord_j).deg
                        
                        # CRITICAL: Require minimum angular separation
                        # to avoid matching galaxies in the same physical region
                        if ang_sep >= min_angular_sep:
                            matches.append({
                                'idx_i': int(idx_i),
                                'idx_j': int(idx_j),
                                'bin_i': i,
                                'bin_j': j,
                                'z_i': gal_i['z'],
                                'z_j': gal_j['z'],
                                'delta_t_gyr': float(delta_t),
                                'fingerprint_distance': float(dist),
                                'angular_separation_deg': float(ang_sep),
                                'ra_i': gal_i['ra'],
                                'dec_i': gal_i['dec'],
                                'ra_j': gal_j['ra'],
                                'dec_j': gal_j['dec'],
                            })
    
    print(f"\n  Total evolved matches found: {len(matches)}")
    return matches, scaler


def null_hypothesis_test(galaxies, z_bins, observed_matches, n_shuffles=300):
    """
    Test against null hypothesis by shuffling galaxy properties.
    
    If the temporal onion is real, we expect MORE evolved matches
    than random shuffling would produce.
    """
    print(f"\nNull hypothesis test ({n_shuffles} shuffles)...")
    
    observed_count = len(observed_matches)
    print(f"  Observed matches: {observed_count}")
    
    # Shuffle and count matches
    null_counts = []
    
    for s in range(n_shuffles):
        if s % 50 == 0:
            print(f"  Shuffle {s}/{n_shuffles}...")
        
        # Shuffle properties while keeping positions fixed
        # This breaks any real temporal-onion correlation
        shuffled_galaxies = []
        indices = np.random.permutation(len(galaxies))
        
        for i, g in enumerate(galaxies):
            # Keep position and redshift, shuffle properties
            source = galaxies[indices[i]]
            shuffled_galaxies.append({
                **g,
                'log_mass': source['log_mass'],
                'log_sigma': source['log_sigma'],
                'log_sfr': source['log_sfr'],
                'sersic_n': source['sersic_n'],
                'sersic_ba': source['sersic_ba'],
            })
        
        # Count matches (fast version)
        count = count_evolved_matches_fast(shuffled_galaxies, z_bins)
        null_counts.append(count)
    
    null_counts = np.array(null_counts)
    
    p_value = np.mean(null_counts >= observed_count)
    z_score = (observed_count - np.mean(null_counts)) / max(np.std(null_counts), 1)
    
    print(f"\n  Null mean: {np.mean(null_counts):.1f} ± {np.std(null_counts):.1f}")
    print(f"  p-value: {p_value:.4f}")
    print(f"  z-score: {z_score:.2f}σ")
    
    return {
        'observed': observed_count,
        'null_mean': float(np.mean(null_counts)),
        'null_std': float(np.std(null_counts)),
        'p_value': float(p_value),
        'z_score': float(z_score),
        'null_distribution': null_counts.tolist(),
    }


def count_evolved_matches_fast(galaxies, z_bins, tolerance=1.5):
    """Fast match counting for null tests."""
    z_values = np.array([g['z'] for g in galaxies])
    fingerprints = np.array([create_fingerprint(g) for g in galaxies])
    
    scaler = StandardScaler()
    fingerprints_norm = scaler.fit_transform(fingerprints)
    
    count = 0
    
    for i in range(len(z_bins)):
        for j in range(i + 1, len(z_bins)):
            mask_i = (z_values >= z_bins[i][0]) & (z_values < z_bins[i][1])
            mask_j = (z_values >= z_bins[j][0]) & (z_values < z_bins[j][1])
            
            idx_i = np.where(mask_i)[0]
            idx_j = np.where(mask_j)[0]
            
            if len(idx_i) == 0 or len(idx_j) == 0:
                continue
            
            z_i = (z_bins[i][0] + z_bins[i][1]) / 2
            z_j = (z_bins[j][0] + z_bins[j][1]) / 2
            delta_t = cosmo.lookback_time(z_i).value - cosmo.lookback_time(z_j).value
            
            # Sample for speed
            sample_i = np.random.choice(idx_i, min(100, len(idx_i)), replace=False)
            
            for ii in sample_i:
                gal = galaxies[ii]
                evolved = evolve_galaxy(gal, -delta_t)
                evolved_fp = create_fingerprint(evolved)
                evolved_fp_norm = (evolved_fp - scaler.mean_) / scaler.scale_
                
                dists = np.linalg.norm(fingerprints_norm[idx_j] - evolved_fp_norm, axis=1)
                count += np.sum(dists < tolerance)
    
    return count


def analyze_match_properties(matches, galaxies):
    """Analyze properties of the evolved matches."""
    if len(matches) == 0:
        return {}
    
    print("\nAnalyzing match properties...")
    
    ang_seps = [m['angular_separation_deg'] for m in matches]
    z_diffs = [abs(m['z_i'] - m['z_j']) for m in matches]
    fp_dists = [m['fingerprint_distance'] for m in matches]
    
    # Check if matches cluster at specific angular separations
    # (would indicate a preferred topology scale)
    
    print(f"  Angular separation: {np.median(ang_seps):.1f}° (median)")
    print(f"  Redshift difference: {np.median(z_diffs):.4f} (median)")
    print(f"  Fingerprint distance: {np.median(fp_dists):.3f} (median)")
    
    return {
        'angular_separation': {
            'median': float(np.median(ang_seps)),
            'mean': float(np.mean(ang_seps)),
            'std': float(np.std(ang_seps)),
        },
        'redshift_difference': {
            'median': float(np.median(z_diffs)),
            'mean': float(np.mean(z_diffs)),
            'std': float(np.std(z_diffs)),
        },
        'fingerprint_distance': {
            'median': float(np.median(fp_dists)),
            'mean': float(np.mean(fp_dists)),
            'std': float(np.std(fp_dists)),
        },
    }


def search_for_topology_scale(matches):
    """
    Search for a preferred angular scale in the matches.
    
    If the universe is a temporal onion with a specific "wrap" scale,
    matches should cluster at that angular separation.
    """
    if len(matches) < 50:
        return None
    
    print("\nSearching for preferred topology scale...")
    
    ang_seps = np.array([m['angular_separation_deg'] for m in matches])
    
    # Histogram of angular separations
    bins = np.linspace(0, 180, 37)  # 5-degree bins
    hist, bin_edges = np.histogram(ang_seps, bins=bins)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    
    # Expected for uniform distribution on sphere
    # (accounts for solid angle: more area at larger separations)
    expected = np.sin(np.radians(bin_centers))
    expected = expected / np.sum(expected) * len(matches)
    
    # Chi-squared test
    from scipy.stats import chisquare
    chi2, p_value = chisquare(hist, expected)
    
    # Find peaks (excess over expected)
    excess = hist - expected
    peak_idx = np.argmax(excess)
    peak_angle = bin_centers[peak_idx]
    peak_excess = excess[peak_idx]
    
    print(f"  Chi-squared: {chi2:.1f}, p = {p_value:.4f}")
    print(f"  Peak excess at: {peak_angle:.0f}° ({peak_excess:.1f} above expected)")
    
    return {
        'chi_squared': float(chi2),
        'p_value': float(p_value),
        'peak_angle_deg': float(peak_angle),
        'peak_excess': float(peak_excess),
        'histogram': hist.tolist(),
        'expected': expected.tolist(),
        'bin_centers': bin_centers.tolist(),
    }


def create_visualization(matches, galaxies, null_results, topology_scale, output_path):
    """Create visualization of temporal onion results."""
    import matplotlib.pyplot as plt
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    
    # 1. Sky distribution of matches
    ax = axes[0, 0]
    ra = [g['ra'] for g in galaxies]
    dec = [g['dec'] for g in galaxies]
    z = [g['z'] for g in galaxies]
    
    scatter = ax.scatter(ra, dec, c=z, s=2, alpha=0.3, cmap='viridis')
    plt.colorbar(scatter, ax=ax, label='Redshift')
    
    # Draw top matches
    if len(matches) > 0:
        sorted_matches = sorted(matches, key=lambda x: x['fingerprint_distance'])[:30]
        for m in sorted_matches:
            ax.plot([m['ra_i'], m['ra_j']], [m['dec_i'], m['dec_j']], 
                   'r-', alpha=0.5, linewidth=1)
    
    ax.set_xlabel('RA (deg)')
    ax.set_ylabel('Dec (deg)')
    ax.set_title('Evolved Matches Across Redshift (Top 30)')
    
    # 2. Null hypothesis test
    ax = axes[0, 1]
    if null_results:
        ax.hist(null_results['null_distribution'], bins=30, alpha=0.7, 
               color='gray', label=f'Null (μ={null_results["null_mean"]:.0f})')
        ax.axvline(null_results['observed'], color='red', linewidth=2,
                  label=f'Observed ({null_results["observed"]})')
        ax.set_xlabel('Number of Evolved Matches')
        ax.set_ylabel('Count')
        ax.set_title(f'Null Hypothesis Test\np={null_results["p_value"]:.3f}, z={null_results["z_score"]:.1f}σ')
        ax.legend()
    
    # 3. Angular separation distribution
    ax = axes[1, 0]
    if len(matches) > 0 and topology_scale:
        ax.bar(topology_scale['bin_centers'], topology_scale['histogram'], 
              width=5, alpha=0.7, label='Observed')
        ax.plot(topology_scale['bin_centers'], topology_scale['expected'], 
               'r--', linewidth=2, label='Expected (uniform)')
        ax.axvline(topology_scale['peak_angle_deg'], color='green', 
                  linestyle=':', label=f'Peak: {topology_scale["peak_angle_deg"]:.0f}°')
        ax.set_xlabel('Angular Separation (deg)')
        ax.set_ylabel('Count')
        ax.set_title('Angular Distribution of Matches')
        ax.legend()
    
    # 4. Redshift vs angular separation
    ax = axes[1, 1]
    if len(matches) > 0:
        z_diff = [abs(m['z_i'] - m['z_j']) for m in matches]
        ang_sep = [m['angular_separation_deg'] for m in matches]
        ax.scatter(z_diff, ang_sep, alpha=0.3, s=10)
        ax.set_xlabel('Redshift Difference |Δz|')
        ax.set_ylabel('Angular Separation (deg)')
        ax.set_title('Match Distribution in z-θ Space')
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    plt.close()
    print(f"\nVisualization saved: {output_path}")


def interpret_results(null_results, topology_scale):
    """Interpret the temporal onion test results."""
    if not null_results:
        return {'verdict': 'INSUFFICIENT DATA'}
    
    p = null_results['p_value']
    z = null_results['z_score']
    
    # Check for excess matches
    if p < 0.01 and z > 2.5:
        match_verdict = "STRONG SIGNAL - Excess evolved matches detected"
    elif p < 0.05 and z > 2:
        match_verdict = "MODERATE SIGNAL - Possible excess matches"
    elif p > 0.95 and z < -2:
        match_verdict = "ANTI-CORRELATION - Fewer matches than expected"
    else:
        match_verdict = "NULL - Match count consistent with random"
    
    # Check for preferred angular scale
    if topology_scale and topology_scale['p_value'] < 0.05:
        scale_verdict = f"POSSIBLE TOPOLOGY SCALE at {topology_scale['peak_angle_deg']:.0f}°"
    else:
        scale_verdict = "No preferred angular scale detected"
    
    if "STRONG" in match_verdict or "MODERATE" in match_verdict:
        overall = "TEMPORAL ONION SIGNAL DETECTED - Further investigation warranted"
        recommendation = "Investigate top matches; expand to deeper surveys"
    else:
        overall = "NULL RESULT - No evidence for temporal onion in MaNGA"
        recommendation = "Deeper surveys (SDSS, DESI) needed for definitive test"
    
    return {
        'match_verdict': match_verdict,
        'scale_verdict': scale_verdict,
        'overall': overall,
        'recommendation': recommendation,
    }


def main():
    """Main analysis pipeline."""
    print("=" * 70)
    print("TEMPORAL ONION HYPOTHESIS: EVOLVED PATTERN MATCHING")
    print("=" * 70)
    print(f"Started: {datetime.now().isoformat()}")
    print("\nHypothesis: The universe is layered like an onion, with the same")
    print("structures appearing at different redshifts in evolved forms.")
    
    # Load data
    galaxies = load_galaxy_data()
    
    # Define redshift bins
    z_bins = [
        (0.01, 0.03),
        (0.03, 0.05),
        (0.05, 0.08),
        (0.08, 0.12),
        (0.12, 0.15),
    ]
    
    # Find evolved matches
    # CRITICAL: Require 30° minimum separation to avoid matching
    # galaxies in the same physical region (which would be trivial)
    min_ang_sep = 30.0  # degrees
    matches, scaler = find_evolved_matches(galaxies, z_bins, 
                                           evolution_tolerance=1.2,
                                           min_angular_sep=min_ang_sep)
    
    # Analyze match properties
    match_properties = analyze_match_properties(matches, galaxies)
    
    # Search for topology scale
    topology_scale = search_for_topology_scale(matches)
    
    # Null hypothesis test
    null_results = null_hypothesis_test(galaxies, z_bins, matches, n_shuffles=300)
    
    # Interpret results
    interpretation = interpret_results(null_results, topology_scale)
    
    # Create visualization
    fig_path = os.path.join(FIGURES_DIR, 'step_3_4_temporal_onion.png')
    create_visualization(matches, galaxies, null_results, topology_scale, fig_path)
    
    # Compile results
    results = {
        'metadata': {
            'timestamp': datetime.now().isoformat(),
            'hypothesis': 'Temporal Onion - same structures at different epochs',
            'n_galaxies': len(galaxies),
            'n_redshift_bins': len(z_bins),
            'evolution_tolerance': 1.2,
        },
        'redshift_bins': z_bins,
        'match_statistics': {
            'total_matches': len(matches),
            'properties': match_properties,
        },
        'topology_scale': topology_scale,
        'null_hypothesis': null_results,
        'interpretation': interpretation,
        'top_matches': sorted(matches, key=lambda x: x['fingerprint_distance'])[:20],
    }
    
    # Save results
    output_path = os.path.join(RESULTS_DIR, 'step_3_4_temporal_onion.json')
    with open(output_path, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\nResults saved: {output_path}")
    
    # Print summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"Galaxies analyzed: {len(galaxies)}")
    print(f"Evolved matches found: {len(matches)}")
    print(f"Null hypothesis p-value: {null_results['p_value']:.4f}")
    print(f"Z-score vs null: {null_results['z_score']:.2f}σ")
    if topology_scale:
        print(f"Peak angular scale: {topology_scale['peak_angle_deg']:.0f}°")
    print(f"\nVerdict: {interpretation['overall']}")
    print(f"Recommendation: {interpretation['recommendation']}")
    
    return results


if __name__ == '__main__':
    results = main()
