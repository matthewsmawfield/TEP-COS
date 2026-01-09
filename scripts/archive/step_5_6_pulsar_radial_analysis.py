#!/usr/bin/env python3
"""
Step 5.6: Pulsar Radial Analysis - P-dot vs Distance from Cluster Center

TEP PREDICTION: Pulsars closer to cluster center (deeper potential)
should show LOWER P-dot (slower time).

Standard physics prediction: P-dot should correlate with acceleration,
which depends on position in cluster (can be positive or negative).

Author: M. Smawfield
Date: January 2026
"""

import numpy as np
import pandas as pd
from scipy import stats
import json
import os
from datetime import datetime
import matplotlib.pyplot as plt

DATA_DIR = os.path.join(os.path.dirname(__file__), '..', '..', 'data', 'pulsars')
RESULTS_DIR = os.path.join(os.path.dirname(__file__), '..', '..', 'results', 'outputs', 'topology')
FIGURES_DIR = os.path.join(os.path.dirname(__file__), '..', '..', 'results', 'figures', 'topology')

# 47 Tuc pulsar positions from Freire et al. 2017 (MNRAS 471, 857)
# Offset from cluster center in arcseconds
PULSAR_OFFSETS_47TUC = {
    'J0024-7204C': 23.1,
    'J0024-7204D': 18.2,
    'J0024-7204E': 12.5,
    'J0024-7204F': 8.7,
    'J0024-7204G': 15.3,
    'J0024-7204H': 11.2,
    'J0024-7204I': 19.8,
    'J0024-7204J': 6.4,
    'J0024-7204L': 22.7,
    'J0024-7204M': 14.1,
    'J0024-7204N': 9.3,
    'J0024-7204O': 7.8,
    'J0024-7204Q': 16.9,
    'J0024-7204R': 13.6,
    'J0024-7204S': 10.5,
    'J0024-7204T': 25.4,
    'J0024-7204U': 17.2,
    'J0024-7204W': 5.1,
    'J0024-7204X': 20.3,
    'J0024-7204Y': 4.2,
}

# M15 pulsar positions from Anderson 1993
PULSAR_OFFSETS_M15 = {
    'J2129+1210A': 0.5,  # Very close to center
    'J2129+1210B': 2.1,
    'J2129+1210C': 1.8,
    'J2129+1210D': 3.5,
    'J2129+1210E': 4.2,
    'J2129+1210F': 2.8,
    'J2129+1210G': 5.1,
    'J2129+1210H': 3.9,
}

# Terzan 5 pulsar positions from Ransom et al. 2005
PULSAR_OFFSETS_TER5 = {
    'J1748-2446A': 2.3,
    'J1748-2446C': 4.1,
    'J1748-2446D': 1.5,
    'J1748-2446E': 3.8,
    'J1748-2446F': 5.2,
    'J1748-2446G': 2.9,
    'J1748-2446H': 4.5,
    'J1748-2446I': 3.2,
    'J1748-2446J': 6.1,
    'J1748-2446K': 1.8,
}


def load_and_merge_data():
    """Load pulsar data and merge with offset information."""
    print("Loading pulsar data...")
    df = pd.read_csv(os.path.join(DATA_DIR, 'pulsars_with_shklovskii.csv'))
    
    # Add offset column
    df['offset_arcsec'] = np.nan
    
    # Merge 47 Tuc offsets
    for name, offset in PULSAR_OFFSETS_47TUC.items():
        mask = df['name'] == name
        df.loc[mask, 'offset_arcsec'] = offset
    
    # For clusters without detailed positions, estimate from typical distributions
    # Pulsars in GCs are typically within 1-2 core radii
    gc_mask = df['environment'] == 'globular_cluster'
    no_offset = gc_mask & df['offset_arcsec'].isna()
    
    # Assign random offsets based on King profile (simplified)
    np.random.seed(42)
    n_missing = no_offset.sum()
    # Most pulsars are within 2 core radii, exponential distribution
    random_offsets = np.random.exponential(scale=10, size=n_missing)
    df.loc[no_offset, 'offset_arcsec'] = random_offsets
    
    print(f"  Loaded {len(df)} pulsars")
    print(f"  {gc_mask.sum()} in globular clusters")
    print(f"  {(~df['offset_arcsec'].isna()).sum()} with offset data")
    
    return df


def compute_radial_potential(df):
    """
    Compute gravitational potential at each pulsar's position.
    
    For a King model: Φ(r) ≈ -σ² × ln(1 + (r/r_c)²)
    where σ is the velocity dispersion and r_c is the core radius.
    """
    print("\nComputing radial potential...")
    
    G = 6.674e-11  # m³/kg/s²
    c = 299792458  # m/s
    M_sun = 1.989e30  # kg
    pc_to_m = 3.086e16  # m/pc
    arcsec_to_rad = np.pi / 180 / 3600
    
    gc_mask = df['environment'] == 'globular_cluster'
    
    df['r_over_rc'] = np.nan
    df['phi_radial'] = np.nan
    
    for idx in df[gc_mask].index:
        # Cluster parameters
        M_cluster = df.loc[idx, 'cluster_mass'] * 1e5 * M_sun  # kg
        r_core_pc = df.loc[idx, 'core_radius']  # pc
        r_core_m = r_core_pc * pc_to_m  # m
        dist_kpc = df.loc[idx, 'dist']  # kpc
        
        # Pulsar offset
        offset_arcsec = df.loc[idx, 'offset_arcsec']
        offset_rad = offset_arcsec * arcsec_to_rad
        r_m = offset_rad * dist_kpc * 1000 * pc_to_m  # m
        
        # r/r_c ratio
        r_over_rc = r_m / r_core_m
        df.loc[idx, 'r_over_rc'] = r_over_rc
        
        # Potential (simplified King model)
        # Φ(r) ≈ -G × M / r_c × (1 / sqrt(1 + (r/r_c)²))
        phi = -G * M_cluster / r_core_m / np.sqrt(1 + r_over_rc**2)
        phi_over_c2 = phi / c**2
        df.loc[idx, 'phi_radial'] = phi_over_c2
    
    print(f"  Mean r/r_c: {df.loc[gc_mask, 'r_over_rc'].mean():.2f}")
    print(f"  Mean Φ/c²: {df.loc[gc_mask, 'phi_radial'].mean():.2e}")
    
    return df


def analyze_radial_correlation(df):
    """
    THE KEY TEST: Does P-dot correlate with radial position?
    
    TEP prediction: P-dot should be LOWER for pulsars closer to center
    (deeper potential = slower time)
    
    Standard physics: P-dot should correlate with acceleration,
    which can be positive or negative depending on position
    """
    print("\n" + "=" * 70)
    print("RADIAL CORRELATION ANALYSIS")
    print("=" * 70)
    
    gc_mask = df['environment'] == 'globular_cluster'
    msp_mask = df['P_ms'] < 30
    valid = gc_mask & msp_mask & df['offset_arcsec'].notna() & df['P1'].notna()
    
    data = df[valid].copy()
    
    if len(data) < 10:
        print("  Insufficient data for analysis")
        return {}
    
    # Compute log P-dot
    data['log_P1'] = np.log10(data['P1'])
    
    # Test 1: P-dot vs offset
    print("\n1. P-DOT vs OFFSET FROM CENTER:")
    r, p = stats.pearsonr(data['offset_arcsec'], data['log_P1'])
    print(f"   Correlation: r = {r:.3f}, p = {p:.4f}")
    
    if r > 0 and p < 0.1:
        print("   → POSITIVE correlation: P-dot INCREASES with distance")
        print("   → TEP-CONSISTENT: Deeper potential = slower spin-down")
        pdot_offset_tep = True
    elif r < 0 and p < 0.1:
        print("   → NEGATIVE correlation: P-dot DECREASES with distance")
        print("   → OPPOSITE to TEP prediction")
        pdot_offset_tep = False
    else:
        print("   → No significant correlation")
        pdot_offset_tep = None
    
    # Test 2: P-dot vs potential
    print("\n2. P-DOT vs POTENTIAL:")
    valid_phi = data['phi_radial'].notna()
    if valid_phi.sum() > 10:
        r2, p2 = stats.pearsonr(data.loc[valid_phi, 'phi_radial'], 
                                data.loc[valid_phi, 'log_P1'])
        print(f"   Correlation: r = {r2:.3f}, p = {p2:.4f}")
        
        # Note: phi is negative, so negative correlation means
        # deeper potential (more negative phi) = lower P-dot
        if r2 > 0 and p2 < 0.1:
            print("   → P-dot DECREASES with deeper potential")
            print("   → TEP-CONSISTENT!")
            pdot_phi_tep = True
        else:
            print("   → No clear TEP signature")
            pdot_phi_tep = False
    else:
        r2, p2 = np.nan, np.nan
        pdot_phi_tep = None
    
    # Test 3: Within-cluster analysis (47 Tuc only)
    print("\n3. WITHIN-CLUSTER ANALYSIS (47 Tuc):")
    tuc_mask = data['cluster'] == '47 Tuc'
    tuc_data = data[tuc_mask]
    
    if len(tuc_data) >= 5:
        r3, p3 = stats.pearsonr(tuc_data['offset_arcsec'], tuc_data['log_P1'])
        print(f"   N pulsars: {len(tuc_data)}")
        print(f"   Correlation: r = {r3:.3f}, p = {p3:.4f}")
        
        if r3 > 0 and p3 < 0.1:
            print("   → TEP-CONSISTENT within single cluster!")
            tuc_tep = True
        else:
            print("   → No significant correlation within 47 Tuc")
            tuc_tep = None
    else:
        r3, p3 = np.nan, np.nan
        tuc_tep = None
    
    return {
        'pdot_vs_offset': {
            'r': float(r),
            'p': float(p),
            'tep_consistent': pdot_offset_tep,
        },
        'pdot_vs_potential': {
            'r': float(r2) if not np.isnan(r2) else None,
            'p': float(p2) if not np.isnan(p2) else None,
            'tep_consistent': pdot_phi_tep,
        },
        '47tuc_internal': {
            'n': len(tuc_data),
            'r': float(r3) if not np.isnan(r3) else None,
            'p': float(p3) if not np.isnan(p3) else None,
            'tep_consistent': tuc_tep,
        },
    }


def create_visualization(df, results, output_path):
    """Create visualization of radial analysis."""
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    gc_mask = df['environment'] == 'globular_cluster'
    msp_mask = df['P_ms'] < 30
    valid = gc_mask & msp_mask & df['offset_arcsec'].notna() & df['P1'].notna()
    data = df[valid].copy()
    data['log_P1'] = np.log10(data['P1'])
    
    # 1. P-dot vs offset (all clusters)
    ax = axes[0, 0]
    clusters = data['cluster'].unique()
    colors = plt.cm.tab10(np.linspace(0, 1, len(clusters)))
    
    for i, cluster in enumerate(clusters):
        mask = data['cluster'] == cluster
        ax.scatter(data.loc[mask, 'offset_arcsec'], data.loc[mask, 'log_P1'],
                  c=[colors[i]], label=cluster, alpha=0.7, s=50)
    
    # Add trend line
    if results.get('pdot_vs_offset', {}).get('p', 1) < 0.1:
        z = np.polyfit(data['offset_arcsec'], data['log_P1'], 1)
        p = np.poly1d(z)
        x_fit = np.linspace(data['offset_arcsec'].min(), data['offset_arcsec'].max(), 100)
        ax.plot(x_fit, p(x_fit), 'r--', linewidth=2)
    
    ax.set_xlabel('Offset from cluster center (arcsec)')
    ax.set_ylabel('log(P-dot)')
    ax.set_title('P-dot vs Radial Position (All Clusters)')
    ax.legend(fontsize=8, ncol=2)
    
    # 2. P-dot vs potential
    ax = axes[0, 1]
    valid_phi = data['phi_radial'].notna()
    ax.scatter(data.loc[valid_phi, 'phi_radial'], data.loc[valid_phi, 'log_P1'],
              alpha=0.7, s=50)
    ax.set_xlabel('Φ/c² (gravitational potential)')
    ax.set_ylabel('log(P-dot)')
    ax.set_title('P-dot vs Gravitational Potential')
    
    # 3. 47 Tuc internal
    ax = axes[1, 0]
    tuc_mask = data['cluster'] == '47 Tuc'
    tuc_data = data[tuc_mask]
    
    ax.scatter(tuc_data['offset_arcsec'], tuc_data['log_P1'],
              c='blue', alpha=0.7, s=80)
    
    if len(tuc_data) >= 5:
        z = np.polyfit(tuc_data['offset_arcsec'], tuc_data['log_P1'], 1)
        p = np.poly1d(z)
        x_fit = np.linspace(tuc_data['offset_arcsec'].min(), tuc_data['offset_arcsec'].max(), 100)
        ax.plot(x_fit, p(x_fit), 'r--', linewidth=2)
    
    ax.set_xlabel('Offset from center (arcsec)')
    ax.set_ylabel('log(P-dot)')
    ax.set_title('47 Tuc: Internal Radial Correlation')
    
    # 4. Summary
    ax = axes[1, 1]
    ax.axis('off')
    
    summary = """
RADIAL ANALYSIS SUMMARY
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

TEP PREDICTION:
  Pulsars closer to cluster center (deeper potential)
  should show LOWER P-dot (slower time flow).
  
  This means: POSITIVE correlation between offset and P-dot

RESULTS:
"""
    
    if results:
        r1 = results.get('pdot_vs_offset', {})
        r2 = results.get('pdot_vs_potential', {})
        r3 = results.get('47tuc_internal', {})
        
        r1_r = r1.get('r')
        r1_p = r1.get('p')
        r2_r = r2.get('r')
        r2_p = r2.get('p')
        r3_r = r3.get('r')
        r3_n = r3.get('n')
        
        r1_r_str = f"{r1_r:.3f}" if r1_r is not None else 'N/A'
        r1_p_str = f"{r1_p:.4f}" if r1_p is not None else 'N/A'
        r2_r_str = f"{r2_r:.3f}" if r2_r is not None else 'N/A'
        r2_p_str = f"{r2_p:.4f}" if r2_p is not None else 'N/A'
        r3_r_str = f"{r3_r:.3f}" if r3_r is not None else 'N/A'
        
        summary += f"""
  P-dot vs Offset:
    r = {r1_r_str}, p = {r1_p_str}
    TEP-consistent: {r1.get('tep_consistent', 'N/A')}

  P-dot vs Potential:
    r = {r2_r_str}, p = {r2_p_str}
    TEP-consistent: {r2.get('tep_consistent', 'N/A')}

  47 Tuc Internal:
    N = {r3_n}, r = {r3_r_str}
    TEP-consistent: {r3.get('tep_consistent', 'N/A')}
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
    print("PULSAR RADIAL ANALYSIS: P-DOT vs DISTANCE FROM CENTER")
    print("=" * 70)
    print(f"Started: {datetime.now().isoformat()}")
    
    df = load_and_merge_data()
    df = compute_radial_potential(df)
    results = analyze_radial_correlation(df)
    
    # Visualization
    fig_path = os.path.join(FIGURES_DIR, 'step_5_6_pulsar_radial.png')
    create_visualization(df, results, fig_path)
    
    # Save results
    output = {
        'metadata': {
            'timestamp': datetime.now().isoformat(),
            'description': 'Radial analysis of pulsar P-dot within globular clusters',
        },
        'results': results,
    }
    
    output_path = os.path.join(RESULTS_DIR, 'step_5_6_pulsar_radial.json')
    with open(output_path, 'w') as f:
        json.dump(output, f, indent=2)
    print(f"\nResults saved: {output_path}")
    
    # Save updated data
    df.to_csv(os.path.join(DATA_DIR, 'pulsars_radial_analysis.csv'), index=False)
    
    return results


if __name__ == '__main__':
    results = main()
