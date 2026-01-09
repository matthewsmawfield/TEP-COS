#!/usr/bin/env python3
"""
Step 3.1: Angular Pattern Repetition Search

Hypothesis: If the universe has closed/repeating topology, the same SPATIAL
CONFIGURATIONS of galaxies should appear at different positions on the sky.

This is fundamentally different from Step 3.0:
- Step 3.0: Same galaxy at different redshifts (requires deep surveys)
- Step 3.1: Same PATTERN of galaxies at different sky positions (works with shallow surveys)

Methodology:
1. For each galaxy, create a "local pattern" - the configuration of its N nearest neighbors
2. Encode this pattern as a rotation/scale-invariant descriptor
3. Search for statistically improbable pattern matches at DIFFERENT sky positions
4. If topology repeats, we expect excess pattern matches beyond random chance

TEP Connection: If distance-redshift relationships are non-standard, what appears
as "different regions" may actually be the same region viewed from different
light-travel paths.

Author: M. Smawfield
Date: January 2026
"""

import numpy as np
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u
from scipy.spatial import cKDTree
from scipy.stats import pearsonr, spearmanr
from sklearn.preprocessing import StandardScaler
import json
import os
from datetime import datetime
from collections import defaultdict
import warnings
warnings.filterwarnings('ignore')

# Paths
DATA_DIR = os.path.join(os.path.dirname(__file__), '..', '..', 'data')
RESULTS_DIR = os.path.join(os.path.dirname(__file__), '..', '..', 'results', 'outputs', 'topology')
FIGURES_DIR = os.path.join(os.path.dirname(__file__), '..', '..', 'results', 'figures', 'topology')

os.makedirs(RESULTS_DIR, exist_ok=True)
os.makedirs(FIGURES_DIR, exist_ok=True)


def load_manga_data():
    """Load MaNGA galaxy catalog with positions and properties."""
    print("Loading MaNGA catalogs...")
    
    drp_path = os.path.join(DATA_DIR, 'drpall', 'drpall-v3_1_1.fits')
    dap_path = os.path.join(DATA_DIR, 'dapall', 'dapall-v3_1_1-3.1.0.fits')
    
    with fits.open(drp_path) as hdul:
        drp = hdul[1].data
    with fits.open(dap_path) as hdul:
        dap = hdul[1].data
    
    # Match by PLATEIFU
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
        
        # Filter valid
        if not (0 < z < 0.2 and mass > 1e6 and 0 < sigma < 500 and 
                -10 < sfr < 100 and 0 < sersic_n < 10):
            continue
        
        gal = {
            'plateifu': pf,
            'ra': drp['objra'][di],
            'dec': drp['objdec'][di],
            'z': z,
            'mass': mass,
            'sigma': sigma,
            'sfr': sfr,
            'sersic_n': sersic_n,
        }
        galaxies.append(gal)
    
    print(f"  Loaded {len(galaxies)} valid galaxies")
    return galaxies


def compute_local_pattern(galaxies, center_idx, n_neighbors=6):
    """
    Compute a rotation/scale-invariant descriptor of the local galaxy pattern.
    
    The descriptor encodes:
    1. Relative distances between neighbors (scale-normalized)
    2. Relative angles between neighbors
    3. Property ratios between neighbors
    
    This is invariant to rotation and overall scale, allowing detection
    of the same pattern at different sky positions/orientations.
    """
    center = galaxies[center_idx]
    
    # Get angular distances to all other galaxies
    coords = np.array([[g['ra'], g['dec']] for g in galaxies])
    center_coord = np.array([center['ra'], center['dec']])
    
    # Approximate angular distance (good for small separations)
    cos_dec = np.cos(np.radians(center['dec']))
    angular_dist = np.sqrt(
        ((coords[:, 0] - center_coord[0]) * cos_dec)**2 +
        (coords[:, 1] - center_coord[1])**2
    )
    
    # Get n nearest neighbors (excluding self)
    neighbor_idx = np.argsort(angular_dist)[1:n_neighbors+1]
    
    if len(neighbor_idx) < n_neighbors:
        return None, None
    
    # Extract neighbor properties
    neighbors = [galaxies[i] for i in neighbor_idx]
    
    # Compute pattern descriptor
    # 1. Relative positions (centered, scale-normalized)
    rel_ra = np.array([n['ra'] - center['ra'] for n in neighbors]) * cos_dec
    rel_dec = np.array([n['dec'] - center['dec'] for n in neighbors])
    
    # Normalize by characteristic scale (median distance)
    scale = np.median(np.sqrt(rel_ra**2 + rel_dec**2))
    if scale < 1e-6:
        return None, None
    
    rel_ra /= scale
    rel_dec /= scale
    
    # 2. Pairwise distance ratios (rotation-invariant)
    distances = np.sqrt(rel_ra**2 + rel_dec**2)
    dist_ratios = []
    for i in range(len(distances)):
        for j in range(i+1, len(distances)):
            if distances[j] > 0:
                dist_ratios.append(distances[i] / distances[j])
    
    # 3. Pairwise angles (rotation-invariant when considering differences)
    angles = np.arctan2(rel_dec, rel_ra)
    angle_diffs = []
    for i in range(len(angles)):
        for j in range(i+1, len(angles)):
            diff = (angles[i] - angles[j]) % (2 * np.pi)
            angle_diffs.append(diff)
    
    # 4. Property ratios (mass, sigma, SFR)
    mass_ratios = []
    sigma_ratios = []
    for i, n in enumerate(neighbors):
        mass_ratios.append(np.log10(n['mass'] / center['mass']))
        sigma_ratios.append(np.log10(n['sigma'] / center['sigma']))
    
    # Combine into descriptor
    descriptor = np.concatenate([
        sorted(dist_ratios),      # Distance structure
        sorted(angle_diffs),      # Angular structure
        sorted(mass_ratios),      # Mass hierarchy
        sorted(sigma_ratios),     # Dynamical hierarchy
    ])
    
    return descriptor, {
        'center_idx': center_idx,
        'neighbor_idx': neighbor_idx.tolist(),
        'scale': float(scale),
        'center_ra': center['ra'],
        'center_dec': center['dec'],
    }


def compute_all_patterns(galaxies, n_neighbors=6):
    """Compute local pattern descriptors for all galaxies."""
    print(f"\nComputing local patterns (n_neighbors={n_neighbors})...")
    
    patterns = []
    metadata = []
    valid_indices = []
    
    for i in range(len(galaxies)):
        if i % 1000 == 0:
            print(f"  Processing galaxy {i}/{len(galaxies)}...")
        
        desc, meta = compute_local_pattern(galaxies, i, n_neighbors)
        if desc is not None:
            patterns.append(desc)
            metadata.append(meta)
            valid_indices.append(i)
    
    # Filter to consistent length descriptors
    if len(patterns) == 0:
        return np.array([]), [], []
    
    # Find most common descriptor length
    lengths = [len(p) for p in patterns]
    target_len = max(set(lengths), key=lengths.count)
    
    # Filter to matching lengths
    filtered_patterns = []
    filtered_metadata = []
    filtered_indices = []
    for p, m, i in zip(patterns, metadata, valid_indices):
        if len(p) == target_len:
            filtered_patterns.append(p)
            filtered_metadata.append(m)
            filtered_indices.append(i)
    
    patterns = np.array(filtered_patterns)
    print(f"  Valid patterns: {len(patterns)} (descriptor dim: {target_len})")
    print(f"  Descriptor dimension: {patterns.shape[1]}")
    
    return patterns, filtered_metadata, filtered_indices


def find_pattern_matches(patterns, metadata, min_angular_sep=10.0, threshold_percentile=1):
    """
    Search for matching patterns at DIFFERENT sky positions.
    
    Key constraint: Matches must be separated by at least min_angular_sep degrees
    to be considered "different regions" (not just overlapping neighborhoods).
    
    If topology repeats, we expect to find statistically improbable matches
    between distant sky regions.
    """
    print(f"\nSearching for pattern matches (min separation: {min_angular_sep}°)...")
    
    n = len(patterns)
    
    # Normalize patterns
    scaler = StandardScaler()
    patterns_norm = scaler.fit_transform(patterns)
    
    # Compute pairwise distances
    print("  Computing pairwise pattern distances...")
    
    # For efficiency, use random sampling for baseline
    n_sample = min(10000, n * (n-1) // 2)
    idx1 = np.random.randint(0, n, n_sample)
    idx2 = np.random.randint(0, n, n_sample)
    sample_distances = np.linalg.norm(patterns_norm[idx1] - patterns_norm[idx2], axis=1)
    
    threshold = np.percentile(sample_distances, threshold_percentile)
    print(f"  Distance threshold ({threshold_percentile}th percentile): {threshold:.3f}")
    
    # Find matches
    matches = []
    
    # Build spatial index for angular separation check
    coords = np.array([[m['center_ra'], m['center_dec']] for m in metadata])
    
    print("  Searching for matches...")
    for i in range(n):
        if i % 500 == 0:
            print(f"    Galaxy {i}/{n}...")
        
        # Compute distances to all other patterns
        dists = np.linalg.norm(patterns_norm - patterns_norm[i], axis=1)
        
        # Find potential matches
        potential = np.where(dists < threshold)[0]
        
        for j in potential:
            if j <= i:  # Avoid duplicates
                continue
            
            # Check angular separation
            cos_dec = np.cos(np.radians(coords[i, 1]))
            ang_sep = np.sqrt(
                ((coords[j, 0] - coords[i, 0]) * cos_dec)**2 +
                (coords[j, 1] - coords[i, 1])**2
            )
            
            if ang_sep >= min_angular_sep:
                matches.append({
                    'idx_i': int(i),
                    'idx_j': int(j),
                    'pattern_distance': float(dists[j]),
                    'angular_separation': float(ang_sep),
                    'ra_i': float(coords[i, 0]),
                    'dec_i': float(coords[i, 1]),
                    'ra_j': float(coords[j, 0]),
                    'dec_j': float(coords[j, 1]),
                })
    
    print(f"  Found {len(matches)} pattern matches")
    
    return matches, threshold


def null_hypothesis_test(patterns, metadata, min_angular_sep, threshold, n_shuffles=200):
    """
    Test against null hypothesis by shuffling pattern-position associations.
    
    If topology is real, the observed match count should exceed shuffled counts.
    """
    print(f"\nNull hypothesis test ({n_shuffles} shuffles)...")
    
    n = len(patterns)
    coords = np.array([[m['center_ra'], m['center_dec']] for m in metadata])
    
    # Normalize patterns
    scaler = StandardScaler()
    patterns_norm = scaler.fit_transform(patterns)
    
    # Count observed matches
    observed = count_matches_fast(patterns_norm, coords, min_angular_sep, threshold)
    print(f"  Observed matches: {observed}")
    
    # Shuffle and count
    null_counts = []
    for s in range(n_shuffles):
        if s % 50 == 0:
            print(f"  Shuffle {s}/{n_shuffles}...")
        
        # Shuffle pattern-position associations
        shuffle_idx = np.random.permutation(n)
        shuffled_patterns = patterns_norm[shuffle_idx]
        
        count = count_matches_fast(shuffled_patterns, coords, min_angular_sep, threshold)
        null_counts.append(count)
    
    null_counts = np.array(null_counts)
    p_value = np.mean(null_counts >= observed)
    z_score = (observed - np.mean(null_counts)) / max(np.std(null_counts), 1)
    
    print(f"\n  Null mean: {np.mean(null_counts):.1f} ± {np.std(null_counts):.1f}")
    print(f"  p-value: {p_value:.4f}")
    print(f"  z-score: {z_score:.2f}σ")
    
    return {
        'observed': int(observed),
        'null_mean': float(np.mean(null_counts)),
        'null_std': float(np.std(null_counts)),
        'p_value': float(p_value),
        'z_score': float(z_score),
    }


def count_matches_fast(patterns_norm, coords, min_angular_sep, threshold):
    """Fast match counting for null tests."""
    n = len(patterns_norm)
    count = 0
    
    # Sample-based counting for speed
    n_sample = min(2000, n)
    sample_idx = np.random.choice(n, n_sample, replace=False)
    
    for i in sample_idx:
        dists = np.linalg.norm(patterns_norm - patterns_norm[i], axis=1)
        potential = np.where(dists < threshold)[0]
        
        for j in potential:
            if j <= i:
                continue
            
            cos_dec = np.cos(np.radians(coords[i, 1]))
            ang_sep = np.sqrt(
                ((coords[j, 0] - coords[i, 0]) * cos_dec)**2 +
                (coords[j, 1] - coords[i, 1])**2
            )
            
            if ang_sep >= min_angular_sep:
                count += 1
    
    # Scale to full sample
    return count * (n / n_sample)


def analyze_best_matches(matches, galaxies, metadata, top_n=10):
    """Analyze the best pattern matches in detail."""
    if len(matches) == 0:
        return []
    
    print(f"\nAnalyzing top {top_n} matches...")
    
    # Sort by pattern distance (best matches first)
    sorted_matches = sorted(matches, key=lambda x: x['pattern_distance'])[:top_n]
    
    detailed = []
    for rank, m in enumerate(sorted_matches):
        meta_i = metadata[m['idx_i']]
        meta_j = metadata[m['idx_j']]
        
        # Get center galaxies
        gal_i = galaxies[meta_i['center_idx']]
        gal_j = galaxies[meta_j['center_idx']]
        
        # Get neighbor properties
        neighbors_i = [galaxies[idx] for idx in meta_i['neighbor_idx']]
        neighbors_j = [galaxies[idx] for idx in meta_j['neighbor_idx']]
        
        detail = {
            'rank': rank + 1,
            'pattern_distance': m['pattern_distance'],
            'angular_separation': m['angular_separation'],
            'region_i': {
                'center': {
                    'ra': gal_i['ra'],
                    'dec': gal_i['dec'],
                    'z': gal_i['z'],
                    'mass': float(gal_i['mass']),
                },
                'n_neighbors': len(neighbors_i),
                'scale_deg': meta_i['scale'],
            },
            'region_j': {
                'center': {
                    'ra': gal_j['ra'],
                    'dec': gal_j['dec'],
                    'z': gal_j['z'],
                    'mass': float(gal_j['mass']),
                },
                'n_neighbors': len(neighbors_j),
                'scale_deg': meta_j['scale'],
            },
            'z_difference': abs(gal_i['z'] - gal_j['z']),
        }
        detailed.append(detail)
        
        print(f"\n  Match {rank+1}:")
        print(f"    Pattern distance: {m['pattern_distance']:.4f}")
        print(f"    Angular separation: {m['angular_separation']:.1f}°")
        print(f"    Region I: RA={gal_i['ra']:.2f}, Dec={gal_i['dec']:.2f}, z={gal_i['z']:.4f}")
        print(f"    Region J: RA={gal_j['ra']:.2f}, Dec={gal_j['dec']:.2f}, z={gal_j['z']:.4f}")
        print(f"    Redshift difference: {detail['z_difference']:.4f}")
    
    return detailed


def create_visualization(matches, galaxies, metadata, null_results, output_path):
    """Create visualization of pattern repetition search."""
    import matplotlib.pyplot as plt
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    
    # 1. Sky distribution with match connections
    ax = axes[0, 0]
    ra = [g['ra'] for g in galaxies]
    dec = [g['dec'] for g in galaxies]
    ax.scatter(ra, dec, s=1, alpha=0.2, c='gray')
    
    # Draw top matches
    sorted_matches = sorted(matches, key=lambda x: x['pattern_distance'])[:30]
    for m in sorted_matches:
        ax.plot([m['ra_i'], m['ra_j']], [m['dec_i'], m['dec_j']], 
               'r-', alpha=0.5, linewidth=1)
        ax.scatter([m['ra_i'], m['ra_j']], [m['dec_i'], m['dec_j']], 
                  c='red', s=20, zorder=5)
    
    ax.set_xlabel('RA (deg)')
    ax.set_ylabel('Dec (deg)')
    ax.set_title('Pattern Matches Across Sky (Top 30)')
    
    # 2. Angular separation distribution
    ax = axes[0, 1]
    if len(matches) > 0:
        seps = [m['angular_separation'] for m in matches]
        ax.hist(seps, bins=30, alpha=0.7, edgecolor='black')
        ax.axvline(np.median(seps), color='r', linestyle='--', 
                  label=f'Median: {np.median(seps):.1f}°')
    ax.set_xlabel('Angular Separation (deg)')
    ax.set_ylabel('Count')
    ax.set_title('Distribution of Match Separations')
    ax.legend()
    
    # 3. Pattern distance distribution
    ax = axes[1, 0]
    if len(matches) > 0:
        dists = [m['pattern_distance'] for m in matches]
        ax.hist(dists, bins=30, alpha=0.7, edgecolor='black')
        ax.axvline(np.median(dists), color='r', linestyle='--',
                  label=f'Median: {np.median(dists):.3f}')
    ax.set_xlabel('Pattern Distance')
    ax.set_ylabel('Count')
    ax.set_title('Distribution of Pattern Similarity')
    ax.legend()
    
    # 4. Null hypothesis comparison
    ax = axes[1, 1]
    if null_results:
        # Create bar comparison
        labels = ['Observed', 'Null Mean']
        values = [null_results['observed'], null_results['null_mean']]
        errors = [0, null_results['null_std']]
        
        bars = ax.bar(labels, values, yerr=errors, capsize=5, 
                     color=['red', 'gray'], alpha=0.7)
        ax.set_ylabel('Number of Matches')
        ax.set_title(f'Null Hypothesis Test\np={null_results["p_value"]:.3f}, z={null_results["z_score"]:.1f}σ')
        
        # Add significance annotation
        if null_results['p_value'] < 0.05:
            ax.annotate('*', xy=(0, values[0]), fontsize=20, ha='center')
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    plt.close()
    print(f"\nVisualization saved: {output_path}")


def interpret_results(null_results, n_matches):
    """Interpret the pattern repetition search results."""
    if not null_results or n_matches == 0:
        return {
            'verdict': 'INSUFFICIENT DATA',
            'recommendation': 'Need more galaxies or larger survey area'
        }
    
    p = null_results['p_value']
    z = null_results['z_score']
    
    if p < 0.01 and z > 2.5:
        verdict = "STRONG SIGNAL - Excess pattern repetition detected"
        recommendation = "Investigate top matches for physical connection"
    elif p < 0.05 and z > 2:
        verdict = "MODERATE SIGNAL - Possible pattern excess"
        recommendation = "Expand analysis to larger surveys"
    elif p > 0.95 and z < -2:
        verdict = "ANTI-CORRELATION - Fewer matches than expected"
        recommendation = "Patterns are MORE unique than random - consistent with cosmic evolution"
    else:
        verdict = "NULL RESULT - Pattern matches consistent with chance"
        recommendation = "No evidence for topology in this sample"
    
    return {
        'verdict': verdict,
        'recommendation': recommendation,
        'p_value': p,
        'z_score': z,
    }


def main():
    """Main analysis pipeline."""
    print("=" * 70)
    print("COSMIC TOPOLOGY SEARCH: ANGULAR PATTERN REPETITION")
    print("=" * 70)
    print(f"Started: {datetime.now().isoformat()}")
    
    # Load data
    galaxies = load_manga_data()
    
    # Compute local patterns
    patterns, metadata, valid_indices = compute_all_patterns(galaxies, n_neighbors=6)
    
    # Find pattern matches
    matches, threshold = find_pattern_matches(
        patterns, metadata, 
        min_angular_sep=15.0,  # Require 15° separation
        threshold_percentile=0.5  # Top 0.5% most similar
    )
    
    # Null hypothesis test
    null_results = null_hypothesis_test(
        patterns, metadata, 
        min_angular_sep=15.0, 
        threshold=threshold,
        n_shuffles=200
    )
    
    # Analyze best matches
    best_matches = analyze_best_matches(matches, galaxies, metadata, top_n=10)
    
    # Create visualization
    fig_path = os.path.join(FIGURES_DIR, 'step_3_1_pattern_repetition.png')
    create_visualization(matches, galaxies, metadata, null_results, fig_path)
    
    # Compile results
    interpretation = interpret_results(null_results, len(matches))
    
    results = {
        'metadata': {
            'timestamp': datetime.now().isoformat(),
            'n_galaxies': len(galaxies),
            'n_valid_patterns': len(patterns),
            'n_neighbors': 6,
            'min_angular_separation_deg': 15.0,
            'threshold_percentile': 0.5,
        },
        'match_statistics': {
            'total_matches': len(matches),
            'pattern_threshold': float(threshold),
        },
        'null_hypothesis': null_results,
        'top_matches': best_matches,
        'interpretation': interpretation,
    }
    
    # Save results
    output_path = os.path.join(RESULTS_DIR, 'step_3_1_pattern_repetition.json')
    with open(output_path, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\nResults saved: {output_path}")
    
    # Print summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"Galaxies analyzed: {len(galaxies)}")
    print(f"Valid patterns: {len(patterns)}")
    print(f"Pattern matches found: {len(matches)}")
    print(f"Null hypothesis p-value: {null_results['p_value']:.4f}")
    print(f"Z-score vs null: {null_results['z_score']:.2f}σ")
    print(f"\nVerdict: {interpretation['verdict']}")
    print(f"Recommendation: {interpretation['recommendation']}")
    
    return results


if __name__ == '__main__':
    results = main()
