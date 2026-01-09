#!/usr/bin/env python3
"""
Step 3.5: Strict Temporal Onion Test

The previous test found too many matches because galaxy properties follow
universal scaling relations. This test adds stricter requirements:

1. Match must be in a DIFFERENT sky region (>30° separation)
2. Properties must match the SPECIFIC predicted evolution (not just scaling relations)
3. Local neighborhood structure must also match (not just individual galaxies)

KEY INSIGHT:
============
If the universe is a temporal onion, we should see:
- The SAME cluster/group appearing at different redshifts
- With EVOLVED properties matching our predictions
- At a DIFFERENT angular position (the structure has "moved" in the onion)
- With SIMILAR local structure (neighboring galaxies also match)

This is like finding the same fingerprint, but stretched and aged.

Author: M. Smawfield
Date: January 2026
"""

import numpy as np
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
from scipy.spatial import cKDTree
from sklearn.preprocessing import StandardScaler
import json
import os
from datetime import datetime
import warnings
warnings.filterwarnings('ignore')

cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

DATA_DIR = os.path.join(os.path.dirname(__file__), '..', '..', 'data')
RESULTS_DIR = os.path.join(os.path.dirname(__file__), '..', '..', 'results', 'outputs', 'topology')
FIGURES_DIR = os.path.join(os.path.dirname(__file__), '..', '..', 'results', 'figures', 'topology')

os.makedirs(RESULTS_DIR, exist_ok=True)
os.makedirs(FIGURES_DIR, exist_ok=True)


def load_galaxy_data():
    """Load MaNGA galaxy data."""
    print("Loading MaNGA galaxy data...")
    
    drp_path = os.path.join(DATA_DIR, 'drpall', 'drpall-v3_1_1.fits')
    dap_path = os.path.join(DATA_DIR, 'dapall', 'dapall-v3_1_1-3.1.0.fits')
    
    with fits.open(drp_path) as hdul:
        drp = hdul[1].data
    with fits.open(dap_path) as hdul:
        dap = hdul[1].data
    
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
        
        if not (0.01 < z < 0.15 and mass > 1e8 and 
                10 < sigma < 400 and -5 < sfr < 50 and 
                0.5 < sersic_n < 8 and 0.1 < sersic_ba < 1):
            continue
        
        galaxies.append({
            'plateifu': pf,
            'ra': drp['objra'][di],
            'dec': drp['objdec'][di],
            'z': z,
            't_lookback': cosmo.lookback_time(z).value,
            'log_mass': np.log10(mass),
            'log_sigma': np.log10(sigma),
            'log_sfr': np.log10(max(sfr, 0.01)),
            'sersic_n': sersic_n,
            'sersic_ba': sersic_ba,
        })
    
    print(f"  Loaded {len(galaxies)} valid galaxies")
    return galaxies


def evolve_properties(props, delta_t_gyr):
    """Predict evolved properties after delta_t Gyr."""
    growth_rate = 0.05
    new_log_mass = props['log_mass'] + np.log10(1 + growth_rate * delta_t_gyr)
    delta_log_mass = new_log_mass - props['log_mass']
    
    return {
        'log_mass': new_log_mass,
        'log_sigma': props['log_sigma'] + 0.25 * delta_log_mass,
        'log_sfr': max(props['log_sfr'] - delta_t_gyr / 3.0 / np.log(10), -3),
        'sersic_n': np.clip(props['sersic_n'] + 0.1 * delta_t_gyr, 0.5, 8),
        'sersic_ba': props['sersic_ba'],
    }


def compute_local_pattern(galaxies, center_idx, n_neighbors=5):
    """
    Compute a rotation/scale-invariant descriptor of the local neighborhood.
    
    This captures the STRUCTURE around a galaxy, not just its properties.
    If the temporal onion is real, the same structure should appear elsewhere.
    """
    coords = np.array([[g['ra'], g['dec']] for g in galaxies])
    center = coords[center_idx]
    
    # Find nearest neighbors
    cos_dec = np.cos(np.radians(center[1]))
    dists = np.sqrt(((coords[:, 0] - center[0]) * cos_dec)**2 + 
                    (coords[:, 1] - center[1])**2)
    
    neighbor_idx = np.argsort(dists)[1:n_neighbors+1]  # Exclude self
    
    if len(neighbor_idx) < n_neighbors:
        return None
    
    # Compute relative positions and properties
    neighbor_dists = dists[neighbor_idx]
    neighbor_angles = np.arctan2(
        coords[neighbor_idx, 1] - center[1],
        (coords[neighbor_idx, 0] - center[0]) * cos_dec
    )
    
    # Sort by distance for consistency
    sort_idx = np.argsort(neighbor_dists)
    neighbor_dists = neighbor_dists[sort_idx]
    neighbor_angles = neighbor_angles[sort_idx]
    neighbor_idx = neighbor_idx[sort_idx]
    
    # Normalize distances (scale-invariant)
    if neighbor_dists[0] > 0:
        neighbor_dists = neighbor_dists / neighbor_dists[0]
    
    # Make rotation-invariant by using relative angles
    ref_angle = neighbor_angles[0]
    neighbor_angles = neighbor_angles - ref_angle
    neighbor_angles = np.mod(neighbor_angles + np.pi, 2*np.pi) - np.pi
    
    # Include neighbor properties (normalized)
    neighbor_props = []
    for idx in neighbor_idx:
        g = galaxies[idx]
        neighbor_props.extend([
            g['log_mass'] / 12,  # Normalize to ~1
            g['log_sigma'] / 2.5,
            g['log_sfr'] / 2,
            g['sersic_n'] / 4,
        ])
    
    # Combine into descriptor
    descriptor = np.concatenate([
        neighbor_dists[1:] / neighbor_dists[-1],  # Relative distances
        np.sort(neighbor_angles[1:]),  # Sorted relative angles
        neighbor_props,
    ])
    
    return descriptor, neighbor_idx


def find_structure_matches(galaxies, min_angular_sep=30.0, 
                          property_threshold=0.8, structure_threshold=1.5):
    """
    Find matches where BOTH individual properties AND local structure match.
    
    This is the strict test: we're looking for the same "fingerprint" 
    (galaxy + its neighborhood) appearing at a different sky position.
    """
    print(f"\nSearching for structure matches (min sep: {min_angular_sep}°)...")
    
    n = len(galaxies)
    z_values = np.array([g['z'] for g in galaxies])
    
    # Compute fingerprints and local patterns
    print("  Computing local patterns...")
    patterns = []
    pattern_indices = []
    
    for i in range(n):
        result = compute_local_pattern(galaxies, i, n_neighbors=5)
        if result is not None:
            patterns.append(result[0])
            pattern_indices.append(i)
    
    # Filter to consistent lengths
    if len(patterns) == 0:
        return []
    
    target_len = len(patterns[0])
    filtered = [(p, i) for p, i in zip(patterns, pattern_indices) if len(p) == target_len]
    patterns = np.array([f[0] for f in filtered])
    pattern_indices = [f[1] for f in filtered]
    
    print(f"  Valid patterns: {len(patterns)}")
    
    # Normalize patterns
    scaler = StandardScaler()
    patterns_norm = scaler.fit_transform(patterns)
    
    # Create property fingerprints
    props = np.array([[
        galaxies[i]['log_mass'],
        galaxies[i]['log_sigma'],
        galaxies[i]['log_sfr'],
        galaxies[i]['sersic_n'],
        galaxies[i]['sersic_ba'],
    ] for i in pattern_indices])
    
    prop_scaler = StandardScaler()
    props_norm = prop_scaler.fit_transform(props)
    
    # Search for matches across redshift
    matches = []
    
    # Define redshift bins
    z_bins = [(0.01, 0.04), (0.04, 0.08), (0.08, 0.15)]
    
    for bi in range(len(z_bins)):
        for bj in range(bi + 1, len(z_bins)):
            z_i_range = z_bins[bi]
            z_j_range = z_bins[bj]
            
            # Get indices in each bin
            idx_i = [k for k, i in enumerate(pattern_indices) 
                    if z_i_range[0] <= galaxies[i]['z'] < z_i_range[1]]
            idx_j = [k for k, i in enumerate(pattern_indices) 
                    if z_j_range[0] <= galaxies[i]['z'] < z_j_range[1]]
            
            if len(idx_i) == 0 or len(idx_j) == 0:
                continue
            
            # Time difference
            z_i_mid = (z_i_range[0] + z_i_range[1]) / 2
            z_j_mid = (z_j_range[0] + z_j_range[1]) / 2
            delta_t = cosmo.lookback_time(z_i_mid).value - cosmo.lookback_time(z_j_mid).value
            
            print(f"  Comparing z={z_i_mid:.2f} to z={z_j_mid:.2f} (Δt={abs(delta_t):.1f} Gyr)")
            
            for ki in idx_i:
                gal_i = galaxies[pattern_indices[ki]]
                
                # Predict evolved properties
                evolved = evolve_properties(gal_i, -delta_t)
                evolved_arr = np.array([
                    evolved['log_mass'],
                    evolved['log_sigma'],
                    evolved['log_sfr'],
                    evolved['sersic_n'],
                    evolved['sersic_ba'],
                ])
                evolved_norm = (evolved_arr - prop_scaler.mean_) / prop_scaler.scale_
                
                for kj in idx_j:
                    gal_j = galaxies[pattern_indices[kj]]
                    
                    # Check angular separation first (fast)
                    cos_dec = np.cos(np.radians(gal_i['dec']))
                    ang_sep = np.sqrt(
                        ((gal_j['ra'] - gal_i['ra']) * cos_dec)**2 +
                        (gal_j['dec'] - gal_i['dec'])**2
                    )
                    
                    if ang_sep < min_angular_sep:
                        continue
                    
                    # Check property match (evolved)
                    prop_dist = np.linalg.norm(props_norm[kj] - evolved_norm)
                    if prop_dist > property_threshold:
                        continue
                    
                    # Check structure match
                    struct_dist = np.linalg.norm(patterns_norm[kj] - patterns_norm[ki])
                    if struct_dist > structure_threshold:
                        continue
                    
                    # This is a candidate match!
                    matches.append({
                        'idx_i': pattern_indices[ki],
                        'idx_j': pattern_indices[kj],
                        'z_i': gal_i['z'],
                        'z_j': gal_j['z'],
                        'angular_sep_deg': float(ang_sep),
                        'property_distance': float(prop_dist),
                        'structure_distance': float(struct_dist),
                        'combined_score': float(prop_dist + struct_dist),
                        'delta_t_gyr': float(delta_t),
                        'ra_i': gal_i['ra'],
                        'dec_i': gal_i['dec'],
                        'ra_j': gal_j['ra'],
                        'dec_j': gal_j['dec'],
                    })
    
    print(f"  Total structure matches: {len(matches)}")
    return matches


def null_test(galaxies, observed_count, n_shuffles=500):
    """
    Null test by shuffling redshifts.
    
    If the temporal onion is real, the observed matches should exceed
    what we get when redshifts are randomized.
    """
    print(f"\nNull hypothesis test ({n_shuffles} shuffles)...")
    print(f"  Observed matches: {observed_count}")
    
    null_counts = []
    
    for s in range(n_shuffles):
        if s % 100 == 0:
            print(f"  Shuffle {s}/{n_shuffles}...")
        
        # Shuffle redshifts
        shuffled = [dict(g) for g in galaxies]
        z_values = [g['z'] for g in galaxies]
        np.random.shuffle(z_values)
        for g, z in zip(shuffled, z_values):
            g['z'] = z
            g['t_lookback'] = cosmo.lookback_time(z).value
        
        # Count matches (simplified for speed)
        count = count_matches_fast(shuffled)
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
    }


def count_matches_fast(galaxies, min_sep=30.0, prop_thresh=0.8, struct_thresh=1.5):
    """Fast match counting for null tests."""
    n = len(galaxies)
    z_values = np.array([g['z'] for g in galaxies])
    
    # Simple property fingerprints
    props = np.array([[
        g['log_mass'], g['log_sigma'], g['log_sfr'], g['sersic_n'], g['sersic_ba']
    ] for g in galaxies])
    
    scaler = StandardScaler()
    props_norm = scaler.fit_transform(props)
    
    count = 0
    z_bins = [(0.01, 0.04), (0.04, 0.08), (0.08, 0.15)]
    
    for bi in range(len(z_bins)):
        for bj in range(bi + 1, len(z_bins)):
            idx_i = np.where((z_values >= z_bins[bi][0]) & (z_values < z_bins[bi][1]))[0]
            idx_j = np.where((z_values >= z_bins[bj][0]) & (z_values < z_bins[bj][1]))[0]
            
            if len(idx_i) == 0 or len(idx_j) == 0:
                continue
            
            # Sample for speed
            sample_i = np.random.choice(idx_i, min(50, len(idx_i)), replace=False)
            
            for i in sample_i:
                gi = galaxies[i]
                for j in idx_j:
                    gj = galaxies[j]
                    
                    # Angular separation
                    cos_dec = np.cos(np.radians(gi['dec']))
                    ang_sep = np.sqrt(
                        ((gj['ra'] - gi['ra']) * cos_dec)**2 +
                        (gj['dec'] - gi['dec'])**2
                    )
                    
                    if ang_sep < min_sep:
                        continue
                    
                    # Property distance
                    prop_dist = np.linalg.norm(props_norm[j] - props_norm[i])
                    if prop_dist < prop_thresh:
                        count += 1
    
    return count


def create_visualization(matches, galaxies, null_results, output_path):
    """Create visualization."""
    import matplotlib.pyplot as plt
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    
    # 1. Sky map with matches
    ax = axes[0, 0]
    ra = [g['ra'] for g in galaxies]
    dec = [g['dec'] for g in galaxies]
    z = [g['z'] for g in galaxies]
    
    scatter = ax.scatter(ra, dec, c=z, s=2, alpha=0.3, cmap='viridis')
    plt.colorbar(scatter, ax=ax, label='Redshift')
    
    if len(matches) > 0:
        sorted_matches = sorted(matches, key=lambda x: x['combined_score'])[:20]
        for m in sorted_matches:
            ax.plot([m['ra_i'], m['ra_j']], [m['dec_i'], m['dec_j']], 
                   'r-', alpha=0.7, linewidth=1.5)
            ax.scatter([m['ra_i'], m['ra_j']], [m['dec_i'], m['dec_j']], 
                      c='red', s=30, zorder=5)
    
    ax.set_xlabel('RA (deg)')
    ax.set_ylabel('Dec (deg)')
    ax.set_title(f'Structure Matches (Top 20 of {len(matches)})')
    
    # 2. Match quality distribution
    ax = axes[0, 1]
    if len(matches) > 0:
        scores = [m['combined_score'] for m in matches]
        ax.hist(scores, bins=30, alpha=0.7, color='steelblue')
        ax.axvline(np.median(scores), color='red', linestyle='--', 
                  label=f'Median: {np.median(scores):.2f}')
        ax.set_xlabel('Combined Score (lower = better)')
        ax.set_ylabel('Count')
        ax.set_title('Match Quality Distribution')
        ax.legend()
    
    # 3. Angular separation distribution
    ax = axes[1, 0]
    if len(matches) > 0:
        ang_seps = [m['angular_sep_deg'] for m in matches]
        ax.hist(ang_seps, bins=20, alpha=0.7, color='steelblue')
        ax.set_xlabel('Angular Separation (deg)')
        ax.set_ylabel('Count')
        ax.set_title('Angular Separation of Matches')
    
    # 4. Summary
    ax = axes[1, 1]
    ax.axis('off')
    
    summary = f"""
STRICT TEMPORAL ONION TEST SUMMARY

Galaxies analyzed: {len(galaxies)}
Structure matches found: {len(matches)}

NULL HYPOTHESIS TEST:
Observed: {null_results['observed']}
Null mean: {null_results['null_mean']:.1f} ± {null_results['null_std']:.1f}
p-value: {null_results['p_value']:.4f}
z-score: {null_results['z_score']:.2f}σ

INTERPRETATION:
"""
    
    if null_results['z_score'] > 3:
        summary += "STRONG SIGNAL - Excess structure matches detected!\n"
        summary += "This suggests possible temporal onion topology."
    elif null_results['z_score'] > 2:
        summary += "MODERATE SIGNAL - Possible excess matches.\n"
        summary += "Further investigation warranted."
    elif null_results['z_score'] < -2:
        summary += "ANTI-CORRELATION - Fewer matches than expected.\n"
        summary += "Galaxies at different z are MORE different than random."
    else:
        summary += "NULL RESULT - Matches consistent with random.\n"
        summary += "No evidence for temporal onion in MaNGA."
    
    ax.text(0.1, 0.9, summary, transform=ax.transAxes, fontsize=10,
           verticalalignment='top', fontfamily='monospace')
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    plt.close()
    print(f"\nVisualization saved: {output_path}")


def main():
    """Main analysis."""
    print("=" * 70)
    print("STRICT TEMPORAL ONION TEST")
    print("=" * 70)
    print(f"Started: {datetime.now().isoformat()}")
    print("\nThis test requires BOTH property evolution AND structure matching.")
    
    galaxies = load_galaxy_data()
    
    # Find structure matches
    # Relax thresholds to allow for evolution uncertainty
    matches = find_structure_matches(
        galaxies, 
        min_angular_sep=30.0,
        property_threshold=1.5,  # More lenient on properties
        structure_threshold=2.5   # More lenient on structure
    )
    
    # Null test
    null_results = null_test(galaxies, len(matches), n_shuffles=500)
    
    # Visualization
    fig_path = os.path.join(FIGURES_DIR, 'step_3_5_onion_strict.png')
    create_visualization(matches, galaxies, null_results, fig_path)
    
    # Save results
    results = {
        'metadata': {
            'timestamp': datetime.now().isoformat(),
            'test': 'Strict Temporal Onion',
            'n_galaxies': len(galaxies),
        },
        'matches': {
            'count': len(matches),
            'top_20': sorted(matches, key=lambda x: x['combined_score'])[:20],
        },
        'null_hypothesis': null_results,
    }
    
    output_path = os.path.join(RESULTS_DIR, 'step_3_5_onion_strict.json')
    with open(output_path, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\nResults saved: {output_path}")
    
    # Summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"Galaxies analyzed: {len(galaxies)}")
    print(f"Structure matches: {len(matches)}")
    print(f"Null p-value: {null_results['p_value']:.4f}")
    print(f"Z-score: {null_results['z_score']:.2f}σ")
    
    if null_results['z_score'] > 2:
        print("\n*** SIGNAL DETECTED - Further investigation warranted ***")
    else:
        print("\nNo significant signal detected.")
    
    return results


if __name__ == '__main__':
    results = main()
