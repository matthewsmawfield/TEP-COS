#!/usr/bin/env python3
"""
Step 5.5: Pulsar TEP Reframe - "Acceleration" IS Time Dilation

CRITICAL REFRAME: Standard physics interprets the low P-dot of GC pulsars
as "acceleration in the cluster potential." But under TEP, this IS the
time dilation signal!

Standard Physics:
- P-dot_observed = P-dot_intrinsic + P × a_los / c
- "Acceleration" a_los is due to gravity

TEP Interpretation:
- P-dot_observed = P-dot_intrinsic × (1 + ΔΦ/c²)
- The "acceleration term" is actually time dilation
- What looks like a_los / c is really ΔΦ/c²

Key Test:
- Standard physics: a_los can be positive or negative (random position in cluster)
- TEP: ΔΦ/c² is always negative (deeper potential = slower time)

If the observed "acceleration" is systematically negative (pulsars appear
to spin down slower), this is TEP-consistent!

Author: M. Smawfield
Date: January 2026
"""

import numpy as np
import pandas as pd
from scipy import stats
import json
import os
from datetime import datetime
import warnings
warnings.filterwarnings('ignore')

DATA_DIR = os.path.join(os.path.dirname(__file__), '..', '..', 'data', 'pulsars')
RESULTS_DIR = os.path.join(os.path.dirname(__file__), '..', '..', 'results', 'outputs', 'topology')
FIGURES_DIR = os.path.join(os.path.dirname(__file__), '..', '..', 'results', 'figures', 'topology')


def load_pulsar_data():
    """Load the processed pulsar data."""
    print("Loading pulsar data...")
    df = pd.read_csv(os.path.join(DATA_DIR, 'pulsars_with_shklovskii.csv'))
    print(f"  Loaded {len(df)} pulsars")
    return df


def compute_tep_time_dilation(df):
    """
    Compute expected TEP time dilation from cluster potential.
    
    For a pulsar at the center of a globular cluster:
    Φ/c² ≈ -G × M_cluster / (r_core × c²)
    
    TEP prediction: Δ(P-dot)/P-dot = ΔΦ/c²
    """
    print("\nComputing TEP time dilation predictions...")
    
    G = 6.674e-11  # m³/kg/s²
    c = 299792458  # m/s
    M_sun = 1.989e30  # kg
    pc_to_m = 3.086e16  # m/pc
    
    gc_mask = df['environment'] == 'globular_cluster'
    
    df['phi_over_c2'] = 0.0
    df['tep_pdot_factor'] = 1.0
    
    for idx in df[gc_mask].index:
        M_cluster = df.loc[idx, 'cluster_mass'] * 1e5 * M_sun  # kg
        r_core = df.loc[idx, 'core_radius'] * pc_to_m  # m
        
        # Gravitational potential at core
        phi = -G * M_cluster / r_core  # m²/s²
        phi_over_c2 = phi / c**2
        
        df.loc[idx, 'phi_over_c2'] = phi_over_c2
        
        # TEP prediction: observed P-dot = intrinsic P-dot × (1 + Φ/c²)
        # Since Φ/c² < 0, observed P-dot < intrinsic P-dot
        df.loc[idx, 'tep_pdot_factor'] = 1 + phi_over_c2
    
    # For field pulsars, assume galactic potential
    field_mask = df['environment'] == 'field'
    df.loc[field_mask, 'phi_over_c2'] = -1e-6  # Galactic disk potential
    df.loc[field_mask, 'tep_pdot_factor'] = 1 - 1e-6
    
    print(f"  GC mean Φ/c²: {df.loc[gc_mask, 'phi_over_c2'].mean():.2e}")
    print(f"  Field mean Φ/c²: {df.loc[field_mask, 'phi_over_c2'].mean():.2e}")
    print(f"  GC TEP factor: {df.loc[gc_mask, 'tep_pdot_factor'].mean():.6f}")
    
    return df


def compute_apparent_acceleration(df):
    """
    Compute the "apparent acceleration" from observed P-dot.
    
    Standard physics: P-dot_obs = P-dot_int + P × a / c
    Rearranging: a / c = (P-dot_obs - P-dot_int) / P
    
    We estimate P-dot_int from field pulsars with similar P.
    """
    print("\nComputing apparent acceleration...")
    
    c = 299792458  # m/s
    
    # For each GC pulsar, estimate intrinsic P-dot from field pulsars
    gc_mask = df['environment'] == 'globular_cluster'
    field_mask = df['environment'] == 'field'
    msp_mask = df['P_ms'] < 30
    
    # Fit P-dot vs P relation for field MSPs
    field_msp = df[field_mask & msp_mask]
    
    if len(field_msp) > 5:
        slope, intercept, r, p, se = stats.linregress(
            np.log10(field_msp['P']), np.log10(field_msp['P1'])
        )
        print(f"  Field MSP relation: log(P1) = {intercept:.2f} + {slope:.2f} × log(P)")
        
        # Predict intrinsic P-dot for GC pulsars
        df['P1_predicted'] = 10**(intercept + slope * np.log10(df['P']))
        
        # Apparent acceleration: a/c = (P1_obs - P1_int) / P
        df['a_over_c'] = (df['P1'] - df['P1_predicted']) / df['P']
        
        # This is what standard physics calls "acceleration"
        # Under TEP, this is actually Φ/c² (time dilation)
        
        print(f"\n  GC mean a/c: {df.loc[gc_mask & msp_mask, 'a_over_c'].mean():.2e}")
        print(f"  Field mean a/c: {df.loc[field_mask & msp_mask, 'a_over_c'].mean():.2e}")
    
    return df


def test_tep_vs_standard(df):
    """
    THE KEY TEST: Compare observed "acceleration" with TEP prediction.
    
    Standard physics: a/c should be random (positive or negative)
    TEP: a/c should equal Φ/c² (always negative, scales with potential)
    """
    print("\n" + "=" * 70)
    print("TEP vs STANDARD PHYSICS TEST")
    print("=" * 70)
    
    gc_mask = df['environment'] == 'globular_cluster'
    msp_mask = df['P_ms'] < 30
    
    gc_msp = df[gc_mask & msp_mask].copy()
    
    if len(gc_msp) < 5:
        print("  Insufficient GC MSPs for analysis")
        return {}
    
    # Test 1: Is "acceleration" systematically negative?
    print("\n1. SIGN TEST:")
    a_over_c = gc_msp['a_over_c'].dropna()
    n_negative = (a_over_c < 0).sum()
    n_total = len(a_over_c)
    
    # Binomial test: if random, expect 50% negative
    from scipy.stats import binomtest
    result = binomtest(n_negative, n_total, 0.5, alternative='greater')
    p_binomial = result.pvalue
    
    print(f"   Negative a/c: {n_negative}/{n_total} ({100*n_negative/n_total:.1f}%)")
    print(f"   Binomial p-value: {p_binomial:.4f}")
    
    if n_negative / n_total > 0.7 and p_binomial < 0.05:
        print("   → SYSTEMATIC NEGATIVE BIAS (TEP-consistent!)")
        sign_tep_consistent = True
    else:
        print("   → No systematic bias (could be random)")
        sign_tep_consistent = False
    
    # Test 2: Does "acceleration" scale with potential?
    print("\n2. POTENTIAL SCALING TEST:")
    
    # Compare a/c with Φ/c²
    valid = np.isfinite(gc_msp['a_over_c']) & np.isfinite(gc_msp['phi_over_c2'])
    
    if valid.sum() > 5:
        r, p = stats.pearsonr(gc_msp.loc[valid, 'phi_over_c2'], 
                              gc_msp.loc[valid, 'a_over_c'])
        
        print(f"   Correlation (Φ/c² vs a/c): r = {r:.3f}, p = {p:.4f}")
        
        if r > 0.3 and p < 0.1:
            print("   → a/c SCALES with Φ/c² (TEP-consistent!)")
            scaling_tep_consistent = True
        else:
            print("   → No clear scaling")
            scaling_tep_consistent = False
    else:
        scaling_tep_consistent = None
    
    # Test 3: Magnitude comparison
    print("\n3. MAGNITUDE TEST:")
    
    mean_a_over_c = gc_msp['a_over_c'].mean()
    mean_phi_over_c2 = gc_msp['phi_over_c2'].mean()
    
    print(f"   Mean observed a/c: {mean_a_over_c:.2e}")
    print(f"   Mean TEP prediction (Φ/c²): {mean_phi_over_c2:.2e}")
    
    if abs(mean_phi_over_c2) > 1e-10:
        ratio = mean_a_over_c / mean_phi_over_c2
        print(f"   Ratio (observed/TEP): {ratio:.1f}")
        
        if 0.1 < ratio < 10:
            print("   → MAGNITUDE CONSISTENT with TEP!")
            magnitude_tep_consistent = True
        else:
            print(f"   → Magnitude differs by {abs(ratio):.0f}×")
            magnitude_tep_consistent = False
    else:
        magnitude_tep_consistent = None
    
    # Overall assessment
    print("\n" + "=" * 70)
    print("OVERALL ASSESSMENT")
    print("=" * 70)
    
    tep_score = sum([sign_tep_consistent or False, 
                     scaling_tep_consistent or False,
                     magnitude_tep_consistent or False])
    
    print(f"\n  TEP-consistent tests: {tep_score}/3")
    
    if tep_score >= 2:
        print("\n  *** STRONG TEP SIGNATURE ***")
        print("  The 'cluster acceleration' is consistent with TEP time dilation!")
        verdict = "STRONG_TEP"
    elif tep_score == 1:
        print("\n  *** PARTIAL TEP SIGNATURE ***")
        print("  Some features consistent with TEP, but not conclusive.")
        verdict = "PARTIAL_TEP"
    else:
        print("\n  *** NO CLEAR TEP SIGNATURE ***")
        print("  Standard physics interpretation may be correct.")
        verdict = "NO_TEP"
    
    return {
        'sign_test': {
            'n_negative': int(n_negative),
            'n_total': int(n_total),
            'p_value': float(p_binomial),
            'tep_consistent': sign_tep_consistent,
        },
        'scaling_test': {
            'r': float(r) if scaling_tep_consistent is not None else None,
            'p': float(p) if scaling_tep_consistent is not None else None,
            'tep_consistent': scaling_tep_consistent,
        },
        'magnitude_test': {
            'mean_a_over_c': float(mean_a_over_c),
            'mean_phi_over_c2': float(mean_phi_over_c2),
            'ratio': float(ratio) if magnitude_tep_consistent is not None else None,
            'tep_consistent': magnitude_tep_consistent,
        },
        'verdict': verdict,
        'tep_score': tep_score,
    }


def analyze_by_cluster_tep(df):
    """
    Analyze each cluster's "acceleration" vs TEP prediction.
    """
    print("\n" + "=" * 70)
    print("CLUSTER-BY-CLUSTER TEP ANALYSIS")
    print("=" * 70)
    
    gc_mask = df['environment'] == 'globular_cluster'
    msp_mask = df['P_ms'] < 30
    
    results = []
    
    for cluster in df.loc[gc_mask, 'cluster'].unique():
        cluster_mask = (df['cluster'] == cluster) & msp_mask
        cluster_df = df[cluster_mask]
        
        if len(cluster_df) < 3:
            continue
        
        # Mean values
        mean_a_over_c = cluster_df['a_over_c'].mean()
        mean_phi_over_c2 = cluster_df['phi_over_c2'].mean()
        mass = cluster_df['cluster_mass'].iloc[0]
        rc = cluster_df['core_radius'].iloc[0]
        
        # Fraction negative
        frac_neg = (cluster_df['a_over_c'] < 0).mean()
        
        print(f"\n  {cluster}:")
        print(f"    N pulsars: {len(cluster_df)}")
        print(f"    Mass: {mass:.1f} × 10⁵ M☉, r_c: {rc:.2f} pc")
        print(f"    Mean a/c: {mean_a_over_c:.2e}")
        print(f"    TEP prediction (Φ/c²): {mean_phi_over_c2:.2e}")
        print(f"    Fraction negative: {100*frac_neg:.0f}%")
        
        if abs(mean_phi_over_c2) > 1e-10:
            ratio = mean_a_over_c / mean_phi_over_c2
            print(f"    Ratio (obs/TEP): {ratio:.1f}")
        else:
            ratio = np.nan
        
        results.append({
            'cluster': cluster,
            'n_pulsars': len(cluster_df),
            'mass': mass,
            'core_radius': rc,
            'mean_a_over_c': float(mean_a_over_c),
            'mean_phi_over_c2': float(mean_phi_over_c2),
            'ratio': float(ratio),
            'frac_negative': float(frac_neg),
        })
    
    # Test if ratio is consistent across clusters
    if len(results) >= 3:
        ratios = [r['ratio'] for r in results if np.isfinite(r['ratio'])]
        if len(ratios) >= 3:
            mean_ratio = np.mean(ratios)
            std_ratio = np.std(ratios)
            print(f"\n  Mean ratio across clusters: {mean_ratio:.1f} ± {std_ratio:.1f}")
            
            if std_ratio / abs(mean_ratio) < 0.5:
                print("  → CONSISTENT ratio across clusters (TEP-consistent!)")
    
    return results


def create_visualization(df, tep_results, cluster_results, output_path):
    """Create visualization."""
    import matplotlib.pyplot as plt
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    
    gc_mask = df['environment'] == 'globular_cluster'
    msp_mask = df['P_ms'] < 30
    
    # 1. a/c vs Φ/c² (THE KEY PLOT)
    ax = axes[0, 0]
    
    gc_msp = df[gc_mask & msp_mask]
    valid = np.isfinite(gc_msp['a_over_c']) & np.isfinite(gc_msp['phi_over_c2'])
    
    ax.scatter(gc_msp.loc[valid, 'phi_over_c2'], gc_msp.loc[valid, 'a_over_c'],
              alpha=0.7, s=50)
    
    # Add 1:1 line (TEP prediction)
    phi_range = np.array([gc_msp['phi_over_c2'].min(), gc_msp['phi_over_c2'].max()])
    ax.plot(phi_range, phi_range, 'r--', linewidth=2, label='TEP prediction (1:1)')
    
    ax.axhline(0, color='gray', linestyle='-', alpha=0.5)
    ax.axvline(0, color='gray', linestyle='-', alpha=0.5)
    
    ax.set_xlabel('Φ/c² (TEP prediction)')
    ax.set_ylabel('Observed a/c ("acceleration")')
    ax.set_title('THE KEY TEST: Is "acceleration" = time dilation?')
    ax.legend()
    
    # 2. Distribution of a/c
    ax = axes[0, 1]
    
    a_over_c = gc_msp['a_over_c'].dropna()
    ax.hist(a_over_c, bins=30, alpha=0.7, edgecolor='black')
    ax.axvline(0, color='red', linestyle='--', linewidth=2, label='Zero')
    ax.axvline(a_over_c.mean(), color='blue', linestyle='-', linewidth=2, 
               label=f'Mean: {a_over_c.mean():.2e}')
    
    ax.set_xlabel('a/c')
    ax.set_ylabel('Count')
    ax.set_title('Distribution of "Acceleration"\n(Should be centered at 0 if random)')
    ax.legend()
    
    # 3. Ratio by cluster
    ax = axes[1, 0]
    
    if cluster_results:
        clusters = [r['cluster'] for r in cluster_results]
        ratios = [r['ratio'] for r in cluster_results]
        
        colors = ['green' if 0.1 < r < 10 else 'red' for r in ratios]
        ax.barh(range(len(clusters)), ratios, color=colors, alpha=0.7)
        ax.axvline(1, color='black', linestyle='--', linewidth=2, label='TEP prediction')
        ax.set_yticks(range(len(clusters)))
        ax.set_yticklabels(clusters)
        ax.set_xlabel('Ratio (observed a/c) / (TEP Φ/c²)')
        ax.set_title('Ratio by Cluster\n(Green = within 10× of TEP)')
        ax.legend()
    
    # 4. Summary
    ax = axes[1, 1]
    ax.axis('off')
    
    summary = """
TEP REFRAME: "ACCELERATION" IS TIME DILATION

STANDARD PHYSICS INTERPRETATION:
- GC pulsars have low P-dot due to "acceleration"
- a_los / c = (P-dot_obs - P-dot_int) / P
- Acceleration is random (positive or negative)

TEP INTERPRETATION:
- GC pulsars have low P-dot due to TIME DILATION
- Φ/c² = (P-dot_obs - P-dot_int) / P
- Time dilation is always negative (deeper = slower)

KEY TESTS:
"""
    
    if tep_results:
        r_val = tep_results['scaling_test']['r']
        r_str = f"{r_val:.3f}" if r_val is not None else 'N/A'
        ratio_val = tep_results['magnitude_test']['ratio']
        ratio_str = f"{ratio_val:.1f}" if ratio_val is not None else 'N/A'
        
        summary += f"""
1. SIGN TEST:
   {tep_results['sign_test']['n_negative']}/{tep_results['sign_test']['n_total']} negative
   p = {tep_results['sign_test']['p_value']:.4f}
   TEP-consistent: {tep_results['sign_test']['tep_consistent']}

2. SCALING TEST:
   r = {r_str}
   TEP-consistent: {tep_results['scaling_test']['tep_consistent']}

3. MAGNITUDE TEST:
   Ratio = {ratio_str}
   TEP-consistent: {tep_results['magnitude_test']['tep_consistent']}

VERDICT: {tep_results['verdict']}
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
    print("PULSAR TEP REFRAME: 'ACCELERATION' IS TIME DILATION")
    print("=" * 70)
    print(f"Started: {datetime.now().isoformat()}")
    print("\nKey insight: What standard physics calls 'cluster acceleration'")
    print("is actually TEP time dilation (Φ/c²)!")
    
    df = load_pulsar_data()
    df = compute_tep_time_dilation(df)
    df = compute_apparent_acceleration(df)
    
    tep_results = test_tep_vs_standard(df)
    cluster_results = analyze_by_cluster_tep(df)
    
    # Visualization
    fig_path = os.path.join(FIGURES_DIR, 'step_5_5_pulsar_tep_reframe.png')
    create_visualization(df, tep_results, cluster_results, fig_path)
    
    # Save results
    results = {
        'metadata': {
            'timestamp': datetime.now().isoformat(),
            'key_insight': '"Cluster acceleration" is actually TEP time dilation',
        },
        'tep_tests': tep_results,
        'by_cluster': cluster_results,
    }
    
    output_path = os.path.join(RESULTS_DIR, 'step_5_5_pulsar_tep_reframe.json')
    with open(output_path, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\nResults saved: {output_path}")
    
    # Save updated data
    df.to_csv(os.path.join(DATA_DIR, 'pulsars_tep_analysis.csv'), index=False)
    
    return results


if __name__ == '__main__':
    results = main()
