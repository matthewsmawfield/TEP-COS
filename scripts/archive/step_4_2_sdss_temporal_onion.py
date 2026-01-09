#!/usr/bin/env python3
"""
Step 4.2: Temporal Onion Test on SDSS Data

Tests the temporal onion hypothesis using 400,000 SDSS galaxies
spanning z = 0.01-0.75 (5.4 Gyr lookback time).

This is a much more powerful test than MaNGA because:
- 40x more galaxies
- 2.7x deeper in cosmic time
- Can see structures separated by billions of years

Author: M. Smawfield
Date: January 2026
"""

import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
from sklearn.preprocessing import StandardScaler
from scipy.spatial import cKDTree
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


def load_sdss_data():
    """Load processed SDSS galaxy catalog."""
    print("Loading SDSS galaxy catalog...")
    path = os.path.join(DATA_DIR, 'sdss_galaxies.csv')
    df = pd.read_csv(path)
    print(f"  Loaded {len(df):,} galaxies")
    print(f"  Redshift range: {df['redshift'].min():.3f} - {df['redshift'].max():.3f}")
    print(f"  Lookback time: {df['t_lookback'].min():.1f} - {df['t_lookback'].max():.1f} Gyr")
    return df


def evolve_properties(props, delta_t_gyr):
    """
    Predict evolved galaxy properties after delta_t Gyr.
    
    Evolution model:
    - Mass: +5% per Gyr (mergers + accretion)
    - Sigma: follows M^0.25 (Faber-Jackson)
    - SFR: exponential decline (τ = 3 Gyr)
    - Sersic: +0.1 per Gyr (disk → spheroid)
    """
    growth_rate = 0.05
    new_log_mass = props['log_mass'] + np.log10(1 + growth_rate * delta_t_gyr)
    delta_log_mass = new_log_mass - props['log_mass']
    
    return {
        'log_mass': new_log_mass,
        'log_sigma': props['log_sigma'] + 0.25 * delta_log_mass,
        'log_sfr': max(props['log_sfr'] - delta_t_gyr / 3.0 / np.log(10), -3),
        'sersic_proxy': np.clip(props['sersic_proxy'] + 0.1 * delta_t_gyr, 0.5, 6),
        'concentration': props['concentration'],
    }


def create_fingerprint(row):
    """Create 5D fingerprint from galaxy properties."""
    return np.array([
        row['log_mass'],
        row['log_sigma'],
        row['log_sfr'],
        row['sersic_proxy'],
        row['concentration'],
    ])


def find_evolved_matches(df, z_bins, min_angular_sep=20.0, 
                         evolution_tolerance=1.0, max_per_bin=20000):
    """
    Search for evolved matches across redshift bins.
    
    For each galaxy at z1, predict its evolved properties at z2,
    then search for matches at z2 that are at a DIFFERENT sky position.
    """
    print(f"\nSearching for evolved matches across {len(z_bins)} redshift bins...")
    print(f"  Min angular separation: {min_angular_sep}°")
    print(f"  Evolution tolerance: {evolution_tolerance}σ")
    
    # Organize by redshift bin
    bin_data = []
    for z_min, z_max in z_bins:
        mask = (df['redshift'] >= z_min) & (df['redshift'] < z_max)
        bin_df = df[mask].copy()
        if len(bin_df) > max_per_bin:
            bin_df = bin_df.sample(max_per_bin, random_state=42)
        bin_data.append(bin_df)
        t_mid = cosmo.lookback_time((z_min + z_max) / 2).value
        print(f"  Bin z={z_min:.2f}-{z_max:.2f}: {len(bin_df):,} galaxies (t~{t_mid:.1f} Gyr)")
    
    # Create fingerprints for all galaxies
    all_fps = []
    for bin_df in bin_data:
        fps = np.array([create_fingerprint(row) for _, row in bin_df.iterrows()])
        all_fps.append(fps)
    
    # Normalize fingerprints globally
    all_fps_flat = np.vstack(all_fps)
    scaler = StandardScaler()
    scaler.fit(all_fps_flat)
    
    fps_norm = [scaler.transform(fps) for fps in all_fps]
    
    # Search for matches across bin pairs
    matches = []
    
    for i in range(len(z_bins)):
        for j in range(i + 1, len(z_bins)):
            z_i = (z_bins[i][0] + z_bins[i][1]) / 2
            z_j = (z_bins[j][0] + z_bins[j][1]) / 2
            
            t_i = cosmo.lookback_time(z_i).value
            t_j = cosmo.lookback_time(z_j).value
            delta_t = t_i - t_j  # Positive if bin_i is older
            
            print(f"\n  Comparing z~{z_i:.2f} to z~{z_j:.2f} (Δt = {abs(delta_t):.1f} Gyr)...")
            
            df_i = bin_data[i]
            df_j = bin_data[j]
            fps_i = fps_norm[i]
            fps_j = fps_norm[j]
            
            # Build KD-tree for fast neighbor search in fingerprint space
            tree_j = cKDTree(fps_j)
            
            # For each galaxy in bin_i, predict evolved state
            n_matches = 0
            
            for idx_i, (_, row_i) in enumerate(df_i.iterrows()):
                # Predict evolved fingerprint
                evolved = evolve_properties(row_i.to_dict(), -delta_t)
                evolved_fp = np.array([
                    evolved['log_mass'],
                    evolved['log_sigma'],
                    evolved['log_sfr'],
                    evolved['sersic_proxy'],
                    evolved['concentration'],
                ])
                evolved_fp_norm = (evolved_fp - scaler.mean_) / scaler.scale_
                
                # Find nearby matches in fingerprint space
                nearby_idx = tree_j.query_ball_point(evolved_fp_norm, evolution_tolerance)
                
                for idx_j in nearby_idx:
                    row_j = df_j.iloc[idx_j]
                    
                    # Check angular separation
                    cos_dec = np.cos(np.radians(row_i['dec']))
                    ang_sep = np.sqrt(
                        ((row_j['ra'] - row_i['ra']) * cos_dec)**2 +
                        (row_j['dec'] - row_i['dec'])**2
                    )
                    
                    if ang_sep >= min_angular_sep:
                        fp_dist = np.linalg.norm(fps_j[idx_j] - evolved_fp_norm)
                        
                        matches.append({
                            'z_i': float(row_i['redshift']),
                            'z_j': float(row_j['redshift']),
                            't_i': float(row_i['t_lookback']),
                            't_j': float(row_j['t_lookback']),
                            'delta_t': float(delta_t),
                            'angular_sep': float(ang_sep),
                            'fp_distance': float(fp_dist),
                            'ra_i': float(row_i['ra']),
                            'dec_i': float(row_i['dec']),
                            'ra_j': float(row_j['ra']),
                            'dec_j': float(row_j['dec']),
                            'mass_i': float(row_i['log_mass']),
                            'mass_j': float(row_j['log_mass']),
                        })
                        n_matches += 1
            
            print(f"    Found {n_matches:,} matches")
    
    print(f"\n  Total matches: {len(matches):,}")
    return matches, scaler


def null_hypothesis_test(df, z_bins, observed_count, n_shuffles=200):
    """
    Test against null by shuffling redshifts.
    
    If temporal onion is real, observed matches should exceed null.
    """
    print(f"\nNull hypothesis test ({n_shuffles} shuffles)...")
    print(f"  Observed matches: {observed_count:,}")
    
    null_counts = []
    
    for s in range(n_shuffles):
        if s % 50 == 0:
            print(f"  Shuffle {s}/{n_shuffles}...")
        
        # Shuffle redshifts while keeping properties fixed
        shuffled = df.copy()
        shuffled['redshift'] = np.random.permutation(df['redshift'].values)
        shuffled['t_lookback'] = cosmo.lookback_time(shuffled['redshift'].values).value
        
        # Count matches (fast approximation)
        count = count_matches_fast(shuffled, z_bins)
        null_counts.append(count)
    
    null_counts = np.array(null_counts)
    
    p_value = np.mean(null_counts >= observed_count)
    z_score = (observed_count - np.mean(null_counts)) / max(np.std(null_counts), 1)
    
    print(f"\n  Null mean: {np.mean(null_counts):,.0f} ± {np.std(null_counts):,.0f}")
    print(f"  p-value: {p_value:.4f}")
    print(f"  z-score: {z_score:.2f}σ")
    
    return {
        'observed': observed_count,
        'null_mean': float(np.mean(null_counts)),
        'null_std': float(np.std(null_counts)),
        'p_value': float(p_value),
        'z_score': float(z_score),
    }


def count_matches_fast(df, z_bins, tolerance=1.0, min_sep=20.0, sample_size=2000):
    """Fast match counting for null tests."""
    count = 0
    
    for i in range(len(z_bins)):
        for j in range(i + 1, len(z_bins)):
            mask_i = (df['redshift'] >= z_bins[i][0]) & (df['redshift'] < z_bins[i][1])
            mask_j = (df['redshift'] >= z_bins[j][0]) & (df['redshift'] < z_bins[j][1])
            
            df_i = df[mask_i]
            df_j = df[mask_j]
            
            if len(df_i) == 0 or len(df_j) == 0:
                continue
            
            # Sample for speed
            if len(df_i) > sample_size:
                df_i = df_i.sample(sample_size, random_state=None)
            if len(df_j) > sample_size * 2:
                df_j = df_j.sample(sample_size * 2, random_state=None)
            
            # Create fingerprints
            fps_i = np.array([[r['log_mass'], r['log_sigma'], r['log_sfr'], 
                              r['sersic_proxy'], r['concentration']] 
                             for _, r in df_i.iterrows()])
            fps_j = np.array([[r['log_mass'], r['log_sigma'], r['log_sfr'],
                              r['sersic_proxy'], r['concentration']]
                             for _, r in df_j.iterrows()])
            
            # Normalize
            all_fps = np.vstack([fps_i, fps_j])
            mean = np.mean(all_fps, axis=0)
            std = np.std(all_fps, axis=0) + 1e-6
            fps_i_norm = (fps_i - mean) / std
            fps_j_norm = (fps_j - mean) / std
            
            # Count matches
            tree = cKDTree(fps_j_norm)
            
            coords_i = df_i[['ra', 'dec']].values
            coords_j = df_j[['ra', 'dec']].values
            
            for k, fp in enumerate(fps_i_norm):
                nearby = tree.query_ball_point(fp, tolerance)
                for m in nearby:
                    cos_dec = np.cos(np.radians(coords_i[k, 1]))
                    sep = np.sqrt(
                        ((coords_j[m, 0] - coords_i[k, 0]) * cos_dec)**2 +
                        (coords_j[m, 1] - coords_i[k, 1])**2
                    )
                    if sep >= min_sep:
                        count += 1
    
    return count


def analyze_matches(matches):
    """Analyze properties of matches."""
    if len(matches) == 0:
        return {}
    
    print("\nAnalyzing match properties...")
    
    ang_seps = [m['angular_sep'] for m in matches]
    delta_ts = [abs(m['delta_t']) for m in matches]
    fp_dists = [m['fp_distance'] for m in matches]
    
    print(f"  Angular separation: {np.median(ang_seps):.1f}° (median)")
    print(f"  Time difference: {np.median(delta_ts):.1f} Gyr (median)")
    print(f"  Fingerprint distance: {np.median(fp_dists):.3f} (median)")
    
    # Check for preferred angular scale
    bins = np.linspace(20, 180, 17)
    hist, _ = np.histogram(ang_seps, bins=bins)
    bin_centers = (bins[:-1] + bins[1:]) / 2
    
    # Expected for uniform on sphere
    expected = np.sin(np.radians(bin_centers))
    expected = expected / np.sum(expected) * len(matches)
    
    peak_idx = np.argmax(hist - expected)
    peak_angle = bin_centers[peak_idx]
    
    return {
        'angular_sep_median': float(np.median(ang_seps)),
        'delta_t_median': float(np.median(delta_ts)),
        'fp_distance_median': float(np.median(fp_dists)),
        'peak_angular_scale': float(peak_angle),
        'histogram': hist.tolist(),
        'expected': expected.tolist(),
        'bin_centers': bin_centers.tolist(),
    }


def create_visualization(matches, df, null_results, match_analysis, output_path):
    """Create visualization."""
    import matplotlib.pyplot as plt
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    
    # 1. Sky distribution with top matches
    ax = axes[0, 0]
    sample = df.sample(min(10000, len(df)))
    scatter = ax.scatter(sample['ra'], sample['dec'], c=sample['redshift'], 
                        s=1, alpha=0.3, cmap='viridis')
    plt.colorbar(scatter, ax=ax, label='Redshift')
    
    if len(matches) > 0:
        top_matches = sorted(matches, key=lambda x: x['fp_distance'])[:30]
        for m in top_matches:
            ax.plot([m['ra_i'], m['ra_j']], [m['dec_i'], m['dec_j']], 
                   'r-', alpha=0.5, linewidth=1)
    
    ax.set_xlabel('RA (deg)')
    ax.set_ylabel('Dec (deg)')
    ax.set_title(f'Evolved Matches (Top 30 of {len(matches):,})')
    
    # 2. Null hypothesis test
    ax = axes[0, 1]
    if null_results:
        ax.axvline(null_results['observed'], color='red', linewidth=2,
                  label=f'Observed ({null_results["observed"]:,})')
        ax.axvline(null_results['null_mean'], color='gray', linestyle='--',
                  label=f'Null mean ({null_results["null_mean"]:,.0f})')
        
        # Show ±2σ range
        ax.axvspan(null_results['null_mean'] - 2*null_results['null_std'],
                  null_results['null_mean'] + 2*null_results['null_std'],
                  alpha=0.3, color='gray', label='Null ±2σ')
        
        ax.set_xlabel('Number of Evolved Matches')
        ax.set_title(f'Null Hypothesis Test\np={null_results["p_value"]:.3f}, z={null_results["z_score"]:.1f}σ')
        ax.legend()
    
    # 3. Angular separation distribution
    ax = axes[1, 0]
    if match_analysis and 'histogram' in match_analysis:
        ax.bar(match_analysis['bin_centers'], match_analysis['histogram'],
              width=10, alpha=0.7, label='Observed')
        ax.plot(match_analysis['bin_centers'], match_analysis['expected'],
               'r--', linewidth=2, label='Expected (uniform)')
        ax.axvline(match_analysis['peak_angular_scale'], color='green',
                  linestyle=':', label=f'Peak: {match_analysis["peak_angular_scale"]:.0f}°')
        ax.set_xlabel('Angular Separation (deg)')
        ax.set_ylabel('Count')
        ax.set_title('Angular Distribution of Matches')
        ax.legend()
    
    # 4. Summary
    ax = axes[1, 1]
    ax.axis('off')
    
    summary = f"""
SDSS TEMPORAL ONION TEST SUMMARY

Galaxies analyzed: {len(df):,}
Redshift range: {df['redshift'].min():.2f} - {df['redshift'].max():.2f}
Lookback time: {df['t_lookback'].min():.1f} - {df['t_lookback'].max():.1f} Gyr

EVOLVED MATCHES:
Total found: {len(matches):,}
"""
    
    if null_results:
        summary += f"""
NULL HYPOTHESIS TEST:
Observed: {null_results['observed']:,}
Null mean: {null_results['null_mean']:,.0f} ± {null_results['null_std']:,.0f}
p-value: {null_results['p_value']:.4f}
z-score: {null_results['z_score']:.2f}σ

"""
        if null_results['z_score'] > 3:
            summary += "VERDICT: STRONG SIGNAL DETECTED"
        elif null_results['z_score'] > 2:
            summary += "VERDICT: MODERATE SIGNAL"
        elif null_results['z_score'] < -2:
            summary += "VERDICT: ANTI-CORRELATION (fewer than expected)"
        else:
            summary += "VERDICT: NULL RESULT (consistent with random)"
    
    ax.text(0.1, 0.9, summary, transform=ax.transAxes, fontsize=10,
           verticalalignment='top', fontfamily='monospace')
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    plt.close()
    print(f"\nVisualization saved: {output_path}")


def main():
    """Main analysis pipeline."""
    print("=" * 70)
    print("TEMPORAL ONION TEST: SDSS 400K GALAXIES")
    print("=" * 70)
    print(f"Started: {datetime.now().isoformat()}")
    
    # Load data
    df = load_sdss_data()
    
    # Define redshift bins spanning the full range
    z_bins = [
        (0.01, 0.05),   # t ~ 0.1-0.6 Gyr
        (0.05, 0.12),   # t ~ 0.7-1.5 Gyr
        (0.12, 0.25),   # t ~ 1.5-2.9 Gyr
        (0.25, 0.45),   # t ~ 2.9-4.5 Gyr
        (0.45, 0.75),   # t ~ 4.5-6.5 Gyr
    ]
    
    # Find evolved matches
    # Use stricter thresholds to reduce matches to manageable level
    matches, scaler = find_evolved_matches(
        df, z_bins, 
        min_angular_sep=30.0,      # Require 30° separation
        evolution_tolerance=0.5,   # Stricter tolerance (0.5σ)
        max_per_bin=8000           # Smaller sample for speed
    )
    
    # Analyze matches
    match_analysis = analyze_matches(matches)
    
    # Null hypothesis test
    null_results = null_hypothesis_test(df, z_bins, len(matches), n_shuffles=200)
    
    # Visualization
    fig_path = os.path.join(FIGURES_DIR, 'step_4_2_sdss_temporal_onion.png')
    create_visualization(matches, df, null_results, match_analysis, fig_path)
    
    # Save results
    results = {
        'metadata': {
            'timestamp': datetime.now().isoformat(),
            'n_galaxies': len(df),
            'z_range': [float(df['redshift'].min()), float(df['redshift'].max())],
            't_lookback_range': [float(df['t_lookback'].min()), float(df['t_lookback'].max())],
        },
        'z_bins': z_bins,
        'matches': {
            'count': len(matches),
            'analysis': match_analysis,
            'top_20': sorted(matches, key=lambda x: x['fp_distance'])[:20],
        },
        'null_hypothesis': null_results,
    }
    
    output_path = os.path.join(RESULTS_DIR, 'step_4_2_sdss_temporal_onion.json')
    with open(output_path, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\nResults saved: {output_path}")
    
    # Summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"Galaxies: {len(df):,}")
    print(f"Lookback time: {df['t_lookback'].max():.1f} Gyr")
    print(f"Evolved matches: {len(matches):,}")
    print(f"Null p-value: {null_results['p_value']:.4f}")
    print(f"Z-score: {null_results['z_score']:.2f}σ")
    
    if null_results['z_score'] > 2:
        print("\n*** SIGNAL DETECTED ***")
    elif null_results['z_score'] < -2:
        print("\n*** ANTI-CORRELATION: Fewer matches than random ***")
    else:
        print("\n*** NULL RESULT ***")
    
    return results


if __name__ == '__main__':
    results = main()
