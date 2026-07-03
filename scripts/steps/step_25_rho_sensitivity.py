#!/usr/bin/env python3
"""
Step 25: Sensitivity Analysis for Intra-Cluster Correlation (rho_intra)

The covariance-aware analysis assumes rho_intra = 0.3 (within-cluster correlation).
This script tests how sensitive the results are to this assumption by varying
rho_intra across a plausible range [0.1, 0.5].

If the GC vs Field difference remains significant across all reasonable values,
the conclusion is robust to this assumption.
"""

import numpy as np
import pandas as pd
from scipy import stats
import json
import os
from pathlib import Path

# Significance threshold (standard alpha = 0.05)
SIGNIFICANCE_THRESHOLD = 0.05
# Additional statistical thresholds
MARGINAL_THRESHOLD = 0.10  # 2x significance threshold for marginal results

def load_pulsar_data():
    """Load pulsar population control data."""
    csv_path = "results/outputs/step_02_pulsar_population_controls.csv"
    if not os.path.exists(csv_path):
        raise FileNotFoundError(f"{csv_path} not found.")
    
    df = pd.read_csv(csv_path)
    return df

def covariance_aware_ttest_varying_rho(gc_df, field_df, rho_intra):
    """
    Perform t-test accounting for within-cluster correlation with specified rho.
    """
    gc_values = gc_df['logPdot_abs'].dropna().values  # Remove NaN values (invalid pulsar measurements)
    field_values = field_df['logPdot_abs'].dropna().values  # Remove NaN values (invalid pulsar measurements)
    
    mean_gc = np.mean(gc_values)
    mean_field = np.mean(field_values)
    diff = mean_gc - mean_field
    
    # Compute effective sample size for GC
    clusters = gc_df.groupby('cluster')
    cluster_sizes = [len(c) for _, c in clusters]
    n_bar = np.mean(cluster_sizes)
    n_gc = len(gc_values)
    
    # Design effect for clustered data
    design_effect = 1 + (n_bar - 1) * rho_intra
    n_eff = n_gc / design_effect
    
    # Standard error with design effect
    var_gc = np.var(gc_values, ddof=1)
    var_field = np.var(field_values, ddof=1)
    
    se_gc = np.sqrt(var_gc * design_effect / n_gc)
    se_field = np.sqrt(var_field / len(field_values))
    se_diff = np.sqrt(se_gc**2 + se_field**2)
    
    # t-statistic with effective degrees of freedom
    t_stat = diff / se_diff
    df_eff = n_eff + len(field_values) - 2
    p_value = 2 * (1 - stats.t.cdf(abs(t_stat), df=df_eff))
    
    return {
        'mean_gc': float(mean_gc),
        'mean_field': float(mean_field),
        'difference_dex': float(diff),
        't_statistic': float(t_stat),
        'p_value': float(p_value),
        'effective_df': float(df_eff),
        'effective_sample_size_gc': float(n_eff),
        'design_effect': float(design_effect),
        'rho_intra': float(rho_intra),
        'significance_sigma': float(abs(t_stat))
    }

def main():
    """Run sensitivity analysis for within-cluster correlation (rho_intra).
    
    Tests robustness of GC vs Field difference across range of rho_intra values
    to validate result is not sensitive to correlation assumption.
    """
    print("=" * 70)
    print("STEP 5.37: Sensitivity Analysis for rho_intra")
    print("=" * 70)
    print("\nTesting robustness of GC vs Field difference to within-cluster")
    print("correlation assumption (rho_intra).")
    print()
    
    # Load data
    df = load_pulsar_data()
    gc_df = df[df['environment'] == 'globular_cluster']
    field_df = df[df['environment'] == 'field']
    
    print(f"Dataset:")
    print(f"  GC pulsars: {len(gc_df)} (in {gc_df['cluster'].nunique()} clusters)")
    print(f"  Field pulsars: {len(field_df)}")
    
    # Test range of rho_intra values
    # Expanded range including 0.0 (no correlation) to 0.7 (very high correlation)
    rho_values = [0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7]
    
    results = []
    
    print("\n" + "-" * 70)
    print("SENSITIVITY TEST RESULTS")
    print("-" * 70)
    print(f"{'rho_intra':<12} {'N_eff':<10} {'t-stat':<10} {'p-value':<12} {'sigma':<8} {'Significant?'}")
    print("-" * 70)
    
    for rho in rho_values:
        result = covariance_aware_ttest_varying_rho(gc_df, field_df, rho)
        results.append(result)
        
        significant = "YES" if result['p_value'] < SIGNIFICANCE_THRESHOLD else "NO"
        print(f"{rho:<12.2f} {result['effective_sample_size_gc']:<10.1f} "
              f"{result['t_statistic']:<10.2f} {result['p_value']:<12.2e} "
              f"{result['significance_sigma']:<8.2f} {significant}")
    
    # Summary
    min_sigma = min([r['significance_sigma'] for r in results])
    max_p = max([r['p_value'] for r in results])
    all_significant = all([r['p_value'] < SIGNIFICANCE_THRESHOLD for r in results])
    
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"Range of rho_intra tested: [{rho_values[0]:.2f}, {rho_values[-1]:.2f}]")
    print(f"Significance range: {min_sigma:.2f}σ to {results[0]['significance_sigma']:.2f}σ")
    print(f"Maximum p-value: {max_p:.2e}")
    print(f"\nRobustness Assessment:")
    if all_significant:
        print(f"  Result: The GC vs Field difference remains significant (p<{SIGNIFICANCE_THRESHOLD})")
        print(f"    across all tested values of rho_intra")
        print(f"  Even with rho_intra=0.5 (very conservative), significance = {results[-1]['significance_sigma']:.2f}σ")
        robustness = "STRONG"
    elif max_p < MARGINAL_THRESHOLD:
        print(f"  ~ MODERATE: GC vs Field difference significant at p<{SIGNIFICANCE_THRESHOLD} for most values")
        print(f"    but marginal at extreme rho_intra values")
        robustness = "MODERATE"
    else:
        print(f"  ⚠ WEAK: GC vs Field difference NOT robust to rho_intra assumption")
        robustness = "WEAK"
    
    # Save results
    output = {
        'rho_range': [float(rho_values[0]), float(rho_values[-1])],
        'rho_step': float(rho_values[1] - rho_values[0]),
        'results_by_rho': results,
        'summary': {
            'min_sigma': float(min_sigma),
            'max_p_value': float(max_p),
            'all_significant': bool(all_significant),
            'robustness_assessment': robustness,
            'baseline_rho': 0.3,
            'baseline_sigma': results[6]['significance_sigma']  # rho=0.3 is index 6
        }
    }
    
    os.makedirs('results/outputs', exist_ok=True)
    with open('results/outputs/step_25_rho_sensitivity.json', 'w') as f:
        json.dump(output, f, indent=2)
    
    print(f"\nResults saved to: results/outputs/step_25_rho_sensitivity.json")
    print("=" * 70)

if __name__ == "__main__":
    main()
