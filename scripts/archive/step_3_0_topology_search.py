#!/usr/bin/env python3
"""
Step 3.0: Cosmic Topology Search via Spectral Fingerprinting

Hypothesis: If the universe has closed/repeating topology, the same galaxies
(or their evolutionary states) may appear at different redshifts. We search
for statistically improbable matches in multi-dimensional property space.

Methodology:
1. Create N-dimensional fingerprint for each galaxy (mass, sigma, SFR, morphology)
2. Normalize to remove redshift-dependent selection effects
3. Search for matches across different redshift bins
4. For candidate matches, test if NEIGHBORS also match (topology signature)
5. Statistical validation against null hypothesis

TEP Connection: If time-flow varies with gravitational context, redshift-distance
relationships may be non-monotonic, allowing the same structure to appear at
multiple apparent distances.

Author: M. Smawfield
Date: January 2026
"""

import numpy as np
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u
from scipy.spatial import cKDTree
from scipy.stats import zscore, pearsonr, spearmanr
from sklearn.preprocessing import StandardScaler
from sklearn.neighbors import NearestNeighbors
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
    """Load and merge MaNGA DRP and DAP catalogs."""
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
    
    # Find common entries
    common = set(drp_plateifu) & set(dap_plateifu)
    print(f"  DRP entries: {len(drp)}")
    print(f"  DAP entries: {len(dap)}")
    print(f"  Common: {len(common)}")
    
    # Create lookup
    drp_idx = {pf: i for i, pf in enumerate(drp_plateifu)}
    dap_idx = {pf: i for i, pf in enumerate(dap_plateifu)}
    
    # Build merged catalog
    galaxies = []
    for pf in common:
        di = drp_idx[pf]
        ai = dap_idx[pf]
        
        gal = {
            'plateifu': pf,
            'ra': drp['objra'][di],
            'dec': drp['objdec'][di],
            'z': drp['nsa_z'][di],
            'mass': drp['nsa_sersic_mass'][di],
            'stellar_sigma': dap['STELLAR_SIGMA_1RE'][ai],
            'ha_sigma': dap['HA_GSIGMA_1RE'][ai],
            'sfr': dap['SFR_TOT'][ai],
            'sersic_n': dap['NSA_SERSIC_N'][ai],
            'sersic_ba': dap['NSA_SERSIC_BA'][ai],
        }
        galaxies.append(gal)
    
    return galaxies


def filter_valid_galaxies(galaxies):
    """Filter to galaxies with valid fingerprint properties."""
    valid = []
    for g in galaxies:
        if (0 < g['z'] < 0.2 and
            g['mass'] > 1e6 and g['mass'] < 1e14 and
            0 < g['stellar_sigma'] < 500 and
            0 < g['ha_sigma'] < 500 and
            -10 < g['sfr'] < 100 and
            0 < g['sersic_n'] < 10 and
            0 < g['sersic_ba'] <= 1):
            valid.append(g)
    
    print(f"  Valid galaxies: {len(valid)} / {len(galaxies)}")
    return valid


def create_fingerprints(galaxies):
    """
    Create multi-dimensional fingerprint for each galaxy.
    
    Fingerprint components (chosen for physical distinctiveness):
    1. log10(mass) - fundamental property
    2. log10(stellar_sigma) - dynamical state
    3. log10(ha_sigma) - gas dynamics
    4. log10(sfr + 1) - star formation activity
    5. sersic_n - morphological type
    6. sersic_ba - inclination/shape
    
    All normalized to z-scores within the sample.
    """
    print("\nCreating fingerprints...")
    
    # Extract raw features
    features = np.array([
        [np.log10(g['mass']),
         np.log10(g['stellar_sigma']),
         np.log10(max(g['ha_sigma'], 1)),
         np.log10(max(g['sfr'], 0.001) + 1),
         g['sersic_n'],
         g['sersic_ba']]
        for g in galaxies
    ])
    
    # Normalize
    scaler = StandardScaler()
    fingerprints = scaler.fit_transform(features)
    
    print(f"  Fingerprint dimensions: {fingerprints.shape[1]}")
    print(f"  Feature means: {scaler.mean_}")
    print(f"  Feature stds: {scaler.scale_}")
    
    return fingerprints, scaler


def bin_by_redshift(galaxies, n_bins=5):
    """Divide galaxies into redshift bins."""
    z_values = np.array([g['z'] for g in galaxies])
    
    # Equal-count bins
    percentiles = np.linspace(0, 100, n_bins + 1)
    bin_edges = np.percentile(z_values, percentiles)
    
    bins = []
    for i in range(n_bins):
        mask = (z_values >= bin_edges[i]) & (z_values < bin_edges[i+1])
        if i == n_bins - 1:  # Include upper edge in last bin
            mask = (z_values >= bin_edges[i]) & (z_values <= bin_edges[i+1])
        
        indices = np.where(mask)[0]
        z_range = (bin_edges[i], bin_edges[i+1])
        bins.append({
            'indices': indices,
            'z_min': z_range[0],
            'z_max': z_range[1],
            'z_mean': np.mean(z_values[mask]),
            'count': len(indices)
        })
    
    print(f"\nRedshift bins:")
    for i, b in enumerate(bins):
        print(f"  Bin {i}: z = {b['z_min']:.4f} - {b['z_max']:.4f}, N = {b['count']}")
    
    return bins


def find_cross_bin_matches(fingerprints, galaxies, z_bins, threshold_sigma=2.0):
    """
    Search for fingerprint matches ACROSS different redshift bins.
    
    A match is significant if the fingerprint distance is < threshold_sigma
    standard deviations from the mean pairwise distance.
    
    Key insight: In a non-repeating universe, matches across z-bins should
    be random. In a repeating topology, we expect EXCESS matches.
    """
    print(f"\nSearching for cross-bin matches (threshold: {threshold_sigma}σ)...")
    
    # Build KD-tree for fast neighbor search
    tree = NearestNeighbors(n_neighbors=10, metric='euclidean')
    tree.fit(fingerprints)
    
    # Calculate baseline: mean and std of random pairwise distances
    n_sample = min(5000, len(fingerprints))
    idx1 = np.random.choice(len(fingerprints), n_sample, replace=False)
    idx2 = np.random.choice(len(fingerprints), n_sample, replace=False)
    random_distances = np.linalg.norm(fingerprints[idx1] - fingerprints[idx2], axis=1)
    mean_dist = np.mean(random_distances)
    std_dist = np.std(random_distances)
    
    print(f"  Baseline distance: {mean_dist:.3f} ± {std_dist:.3f}")
    
    threshold = mean_dist - threshold_sigma * std_dist
    print(f"  Match threshold: < {threshold:.3f}")
    
    # Search for matches across bins
    matches = []
    
    for i, bin_i in enumerate(z_bins):
        for j, bin_j in enumerate(z_bins):
            if j <= i:  # Only compare different bins, avoid duplicates
                continue
            
            # For each galaxy in bin_i, find nearest neighbors in bin_j
            for idx_i in bin_i['indices']:
                fp_i = fingerprints[idx_i:idx_i+1]
                
                # Get distances to all galaxies in bin_j
                fp_j = fingerprints[bin_j['indices']]
                distances = np.linalg.norm(fp_j - fp_i, axis=1)
                
                # Find matches below threshold
                match_mask = distances < threshold
                if np.any(match_mask):
                    for k, is_match in enumerate(match_mask):
                        if is_match:
                            idx_j = bin_j['indices'][k]
                            matches.append({
                                'idx_i': int(idx_i),
                                'idx_j': int(idx_j),
                                'bin_i': i,
                                'bin_j': j,
                                'distance': float(distances[k]),
                                'sigma': float((mean_dist - distances[k]) / std_dist),
                                'z_i': galaxies[idx_i]['z'],
                                'z_j': galaxies[idx_j]['z'],
                                'ra_i': galaxies[idx_i]['ra'],
                                'dec_i': galaxies[idx_i]['dec'],
                                'ra_j': galaxies[idx_j]['ra'],
                                'dec_j': galaxies[idx_j]['dec'],
                            })
    
    print(f"  Found {len(matches)} cross-bin matches")
    
    return matches, mean_dist, std_dist


def test_neighbor_correlation(matches, galaxies, fingerprints, n_neighbors=5):
    """
    THE KEY TEST: If topology repeats, matched galaxies should have
    matched NEIGHBORS too.
    
    For each match pair (A, B), we check if A's neighbors are similar
    to B's neighbors. This is the signature of a repeating structure.
    """
    print(f"\nTesting neighbor correlation (n_neighbors={n_neighbors})...")
    
    if len(matches) == 0:
        print("  No matches to test")
        return []
    
    # Build spatial trees for neighbor finding
    coords = np.array([[g['ra'], g['dec']] for g in galaxies])
    
    # For each match, compare neighbor fingerprints
    topology_scores = []
    
    for m in matches:
        idx_i, idx_j = m['idx_i'], m['idx_j']
        
        # Find spatial neighbors of each galaxy
        coord_i = coords[idx_i]
        coord_j = coords[idx_j]
        
        # Angular distances
        dist_i = np.sqrt((coords[:, 0] - coord_i[0])**2 * np.cos(np.radians(coord_i[1]))**2 + 
                         (coords[:, 1] - coord_i[1])**2)
        dist_j = np.sqrt((coords[:, 0] - coord_j[0])**2 * np.cos(np.radians(coord_j[1]))**2 + 
                         (coords[:, 1] - coord_j[1])**2)
        
        # Get n nearest neighbors (excluding self)
        neighbors_i = np.argsort(dist_i)[1:n_neighbors+1]
        neighbors_j = np.argsort(dist_j)[1:n_neighbors+1]
        
        # Compare neighbor fingerprints
        fp_neighbors_i = fingerprints[neighbors_i]
        fp_neighbors_j = fingerprints[neighbors_j]
        
        # Compute similarity: average minimum distance between neighbor sets
        min_dists = []
        for fp_ni in fp_neighbors_i:
            dists = np.linalg.norm(fp_neighbors_j - fp_ni, axis=1)
            min_dists.append(np.min(dists))
        
        neighbor_similarity = np.mean(min_dists)
        
        topology_scores.append({
            **m,
            'neighbor_similarity': float(neighbor_similarity),
            'neighbors_i': neighbors_i.tolist(),
            'neighbors_j': neighbors_j.tolist(),
        })
    
    # Sort by neighbor similarity (lower = more similar = stronger topology signal)
    topology_scores.sort(key=lambda x: x['neighbor_similarity'])
    
    print(f"  Best neighbor similarity: {topology_scores[0]['neighbor_similarity']:.3f}")
    print(f"  Worst neighbor similarity: {topology_scores[-1]['neighbor_similarity']:.3f}")
    print(f"  Median: {topology_scores[len(topology_scores)//2]['neighbor_similarity']:.3f}")
    
    return topology_scores


def null_hypothesis_test(fingerprints, galaxies, z_bins, n_shuffles=1000):
    """
    Test against null hypothesis: shuffle redshifts and repeat analysis.
    
    If topology is real, the observed match count should exceed
    the shuffled distribution.
    """
    print(f"\nNull hypothesis test ({n_shuffles} shuffles)...")
    
    # Get observed match count
    matches, mean_dist, std_dist = find_cross_bin_matches(
        fingerprints, galaxies, z_bins, threshold_sigma=2.0
    )
    observed_count = len(matches)
    
    # Shuffle and count
    null_counts = []
    z_values = np.array([g['z'] for g in galaxies])
    
    for i in range(n_shuffles):
        if i % 100 == 0:
            print(f"  Shuffle {i}/{n_shuffles}...")
        
        # Shuffle redshifts
        shuffled_z = np.random.permutation(z_values)
        shuffled_galaxies = [dict(g, z=z) for g, z in zip(galaxies, shuffled_z)]
        
        # Re-bin
        shuffled_bins = bin_by_redshift_quiet(shuffled_galaxies, n_bins=len(z_bins))
        
        # Count matches (using pre-computed threshold)
        threshold = mean_dist - 2.0 * std_dist
        count = count_cross_bin_matches(fingerprints, shuffled_bins, threshold)
        null_counts.append(count)
    
    null_counts = np.array(null_counts)
    p_value = np.mean(null_counts >= observed_count)
    z_score = (observed_count - np.mean(null_counts)) / max(np.std(null_counts), 1)
    
    print(f"\n  Observed matches: {observed_count}")
    print(f"  Null mean: {np.mean(null_counts):.1f} ± {np.std(null_counts):.1f}")
    print(f"  p-value: {p_value:.4f}")
    print(f"  z-score: {z_score:.2f}σ")
    
    return {
        'observed': observed_count,
        'null_mean': float(np.mean(null_counts)),
        'null_std': float(np.std(null_counts)),
        'p_value': float(p_value),
        'z_score': float(z_score),
        'null_distribution': null_counts.tolist()
    }


def bin_by_redshift_quiet(galaxies, n_bins=5):
    """Silent version of bin_by_redshift for null testing."""
    z_values = np.array([g['z'] for g in galaxies])
    percentiles = np.linspace(0, 100, n_bins + 1)
    bin_edges = np.percentile(z_values, percentiles)
    
    bins = []
    for i in range(n_bins):
        mask = (z_values >= bin_edges[i]) & (z_values < bin_edges[i+1])
        if i == n_bins - 1:
            mask = (z_values >= bin_edges[i]) & (z_values <= bin_edges[i+1])
        bins.append({'indices': np.where(mask)[0]})
    
    return bins


def count_cross_bin_matches(fingerprints, z_bins, threshold):
    """Fast match counting for null tests."""
    count = 0
    for i, bin_i in enumerate(z_bins):
        for j, bin_j in enumerate(z_bins):
            if j <= i:
                continue
            for idx_i in bin_i['indices']:
                fp_i = fingerprints[idx_i:idx_i+1]
                fp_j = fingerprints[bin_j['indices']]
                distances = np.linalg.norm(fp_j - fp_i, axis=1)
                count += np.sum(distances < threshold)
    return count


def analyze_match_properties(matches, galaxies):
    """Analyze properties of matched galaxy pairs."""
    if len(matches) == 0:
        return {}
    
    print("\nAnalyzing match properties...")
    
    # Angular separations
    separations = []
    z_diffs = []
    mass_ratios = []
    
    for m in matches:
        g_i = galaxies[m['idx_i']]
        g_j = galaxies[m['idx_j']]
        
        # Angular separation
        coord_i = SkyCoord(ra=g_i['ra']*u.deg, dec=g_i['dec']*u.deg)
        coord_j = SkyCoord(ra=g_j['ra']*u.deg, dec=g_j['dec']*u.deg)
        sep = coord_i.separation(coord_j).deg
        separations.append(sep)
        
        z_diffs.append(abs(g_i['z'] - g_j['z']))
        mass_ratios.append(g_i['mass'] / g_j['mass'])
    
    separations = np.array(separations)
    z_diffs = np.array(z_diffs)
    mass_ratios = np.array(mass_ratios)
    
    print(f"  Angular separation: {np.median(separations):.1f}° (median)")
    print(f"  Redshift difference: {np.median(z_diffs):.4f} (median)")
    print(f"  Mass ratio: {np.median(mass_ratios):.2f} (median)")
    
    return {
        'angular_separation_median': float(np.median(separations)),
        'angular_separation_std': float(np.std(separations)),
        'z_diff_median': float(np.median(z_diffs)),
        'z_diff_std': float(np.std(z_diffs)),
        'mass_ratio_median': float(np.median(mass_ratios)),
        'mass_ratio_std': float(np.std(mass_ratios)),
    }


def create_visualization(matches, galaxies, topology_scores, null_results, output_path):
    """Create visualization of topology search results."""
    import matplotlib.pyplot as plt
    
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # 1. Sky distribution of matches
    ax = axes[0, 0]
    ra = [g['ra'] for g in galaxies]
    dec = [g['dec'] for g in galaxies]
    ax.scatter(ra, dec, s=1, alpha=0.3, c='gray', label='All galaxies')
    
    if len(matches) > 0:
        for m in matches[:50]:  # Top 50 matches
            ax.plot([m['ra_i'], m['ra_j']], [m['dec_i'], m['dec_j']], 
                   'r-', alpha=0.3, linewidth=0.5)
    ax.set_xlabel('RA (deg)')
    ax.set_ylabel('Dec (deg)')
    ax.set_title('Cross-Redshift Matches (Top 50)')
    
    # 2. Redshift distribution of matches
    ax = axes[0, 1]
    if len(matches) > 0:
        z_i = [m['z_i'] for m in matches]
        z_j = [m['z_j'] for m in matches]
        ax.scatter(z_i, z_j, alpha=0.5, s=10)
        ax.plot([0, 0.15], [0, 0.15], 'k--', alpha=0.3)
    ax.set_xlabel('Redshift (Galaxy 1)')
    ax.set_ylabel('Redshift (Galaxy 2)')
    ax.set_title('Redshift Pairs of Matches')
    
    # 3. Null hypothesis test
    ax = axes[1, 0]
    if null_results:
        ax.hist(null_results['null_distribution'], bins=30, alpha=0.7, 
               label=f'Null (μ={null_results["null_mean"]:.0f})')
        ax.axvline(null_results['observed'], color='r', linewidth=2,
                  label=f'Observed ({null_results["observed"]})')
        ax.set_xlabel('Number of Matches')
        ax.set_ylabel('Count')
        ax.set_title(f'Null Hypothesis Test (p={null_results["p_value"]:.3f})')
        ax.legend()
    
    # 4. Neighbor similarity distribution
    ax = axes[1, 1]
    if len(topology_scores) > 0:
        similarities = [t['neighbor_similarity'] for t in topology_scores]
        ax.hist(similarities, bins=30, alpha=0.7)
        ax.axvline(np.median(similarities), color='r', linestyle='--',
                  label=f'Median: {np.median(similarities):.2f}')
        ax.set_xlabel('Neighbor Similarity (lower = more similar)')
        ax.set_ylabel('Count')
        ax.set_title('Topology Score Distribution')
        ax.legend()
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    plt.close()
    print(f"\nVisualization saved: {output_path}")


def main():
    """Main analysis pipeline."""
    print("=" * 60)
    print("COSMIC TOPOLOGY SEARCH VIA SPECTRAL FINGERPRINTING")
    print("=" * 60)
    print(f"Started: {datetime.now().isoformat()}")
    
    # Load data
    galaxies = load_manga_data()
    galaxies = filter_valid_galaxies(galaxies)
    
    # Create fingerprints
    fingerprints, scaler = create_fingerprints(galaxies)
    
    # Bin by redshift
    z_bins = bin_by_redshift(galaxies, n_bins=5)
    
    # Find cross-bin matches
    matches, mean_dist, std_dist = find_cross_bin_matches(
        fingerprints, galaxies, z_bins, threshold_sigma=2.5
    )
    
    # Test neighbor correlation (topology signature)
    topology_scores = test_neighbor_correlation(matches, galaxies, fingerprints)
    
    # Analyze match properties
    match_properties = analyze_match_properties(matches, galaxies)
    
    # Null hypothesis test
    null_results = null_hypothesis_test(fingerprints, galaxies, z_bins, n_shuffles=500)
    
    # Create visualization
    fig_path = os.path.join(FIGURES_DIR, 'step_3_0_topology_search.png')
    create_visualization(matches, galaxies, topology_scores, null_results, fig_path)
    
    # Compile results
    results = {
        'metadata': {
            'timestamp': datetime.now().isoformat(),
            'n_galaxies': len(galaxies),
            'n_redshift_bins': len(z_bins),
            'fingerprint_dimensions': fingerprints.shape[1],
        },
        'redshift_bins': [
            {'bin': i, 'z_min': b['z_min'], 'z_max': b['z_max'], 
             'z_mean': b['z_mean'], 'count': b['count']}
            for i, b in enumerate(z_bins)
        ],
        'match_statistics': {
            'total_matches': len(matches),
            'mean_fingerprint_distance': float(mean_dist),
            'std_fingerprint_distance': float(std_dist),
            'threshold_sigma': 2.5,
        },
        'match_properties': match_properties,
        'null_hypothesis': null_results,
        'top_topology_candidates': topology_scores[:20] if topology_scores else [],
        'interpretation': interpret_results(null_results, topology_scores),
    }
    
    # Save results
    output_path = os.path.join(RESULTS_DIR, 'step_3_0_topology_search.json')
    with open(output_path, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\nResults saved: {output_path}")
    
    # Print summary
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"Galaxies analyzed: {len(galaxies)}")
    print(f"Cross-bin matches found: {len(matches)}")
    print(f"Null hypothesis p-value: {null_results['p_value']:.4f}")
    print(f"Z-score vs null: {null_results['z_score']:.2f}σ")
    print(f"\nInterpretation: {results['interpretation']['verdict']}")
    
    return results


def interpret_results(null_results, topology_scores):
    """Interpret the topology search results."""
    if not null_results:
        return {'verdict': 'INSUFFICIENT DATA'}
    
    p = null_results['p_value']
    z = null_results['z_score']
    
    if p < 0.001 and z > 3:
        verdict = "STRONG TOPOLOGY SIGNAL - Excess matches significantly exceed null expectation"
        recommendation = "Investigate top candidates with deeper spectroscopy"
    elif p < 0.05 and z > 2:
        verdict = "MODERATE TOPOLOGY SIGNAL - Possible excess matches"
        recommendation = "Expand to larger surveys (SDSS full) for confirmation"
    elif p < 0.1:
        verdict = "WEAK SIGNAL - Marginal excess, likely statistical fluctuation"
        recommendation = "No strong evidence for repeating topology in this sample"
    else:
        verdict = "NULL RESULT - Match count consistent with random expectation"
        recommendation = "No evidence for cosmic topology in MaNGA sample"
    
    return {
        'verdict': verdict,
        'recommendation': recommendation,
        'p_value': p,
        'z_score': z,
    }


if __name__ == '__main__':
    results = main()
