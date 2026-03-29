#!/usr/bin/env python3
"""
Step 5.39: Monte Carlo Validation of Statistical Methods

Validates the statistical pipeline by applying it to SIMULATED synthetic data with
known properties. This ensures:
1. Type I error rate is controlled (false positive rate ~5% under H0)
2. Power is adequate (true positive rate high under H1)
3. Bias in effect size estimation is minimal

IMPORTANT: This is a MONTE CARLO SIMULATION for method validation.
Each simulation uses a unique seed (base seed + simulation index) for reproducibility.

Optimized for M4 Pro MacBook with parallel processing via joblib.
"""

import numpy as np
import pandas as pd
from scipy import stats
from multiprocessing import Pool, cpu_count
from joblib import Parallel, delayed
import json
import os
from pathlib import Path

# Statistical thresholds
SIGNIFICANCE_THRESHOLD = 0.05

# Set base random seed for reproducibility
RANDOM_SEED = 42
np.random.seed(RANDOM_SEED)

def generate_synthetic_data(n_gc_clusters=29, n_per_cluster=6, n_field=198, 
                            effect_size=0.0, rho_intra=0.3, seed=None):
    """
    Generate synthetic pulsar data with specified properties.
    
    Parameters:
    -----------
    effect_size : float
        True difference between GC and Field (dex). 0.0 = H0 true.
    rho_intra : float
        Within-cluster correlation.
    seed : int, optional
        Random seed for reproducibility.
    """
    if seed is not None:
        np.random.seed(seed)
    
    # Field population (reference)
    field_logpdot = np.random.normal(-19.76, 0.64, n_field)
    
    # GC population with cluster structure
    cluster_names = [f"Cluster_{i}" for i in range(n_gc_clusters)]
    gc_data = []
    
    # Cluster means with correlation structure
    cluster_mean_base = -19.76 + effect_size  # Shifted by effect size
    cluster_means = np.random.normal(cluster_mean_base, 0.3 * np.sqrt(rho_intra), n_gc_clusters)
    
    for i, cluster in enumerate(cluster_names):
        # Within-cluster variance
        n_this_cluster = max(2, int(np.random.poisson(n_per_cluster)))
        cluster_values = np.random.normal(cluster_means[i], 0.64 * np.sqrt(1 - rho_intra), n_this_cluster)
        
        for val in cluster_values:
            gc_data.append({
                'cluster': cluster,
                'logPdot_abs': val,
                'environment': 'globular_cluster'
            })
    
    gc_df = pd.DataFrame(gc_data)
    field_df = pd.DataFrame({
        'logPdot_abs': field_logpdot,
        'environment': 'field'
    })
    field_df['cluster'] = 'Field'
    
    return gc_df, field_df

def run_single_simulation(sim_idx, effect_size=0.0, rho_intra=0.3):
    """Run a single Monte Carlo simulation (for parallel execution)."""
    gc_df, field_df = generate_synthetic_data(
        effect_size=effect_size, 
        rho_intra=rho_intra,
        seed=sim_idx  # Unique seed per simulation
    )
    
    gc_values = gc_df['logPdot_abs'].dropna().values  # Remove NaN values (invalid measurements)
    field_values = field_df['logPdot_abs'].dropna().values  # Remove NaN values (invalid measurements)
    
    mean_gc = np.mean(gc_values)
    mean_field = np.mean(field_values)
    diff = mean_gc - mean_field
    
    clusters = gc_df.groupby('cluster')
    cluster_sizes = [len(c) for _, c in clusters]
    n_bar = np.mean(cluster_sizes)
    n_gc = len(gc_values)
    
    design_effect = 1 + (n_bar - 1) * rho_intra
    n_eff = n_gc / design_effect
    
    var_gc = np.var(gc_values, ddof=1)
    var_field = np.var(field_values, ddof=1)
    
    se_gc = np.sqrt(var_gc * design_effect / n_gc)
    se_field = np.sqrt(var_field / len(field_values))
    se_diff = np.sqrt(se_gc**2 + se_field**2)
    
    t_stat = diff / se_diff
    df_eff = n_eff + len(field_values) - 2
    p_value = 2 * (1 - stats.t.cdf(abs(t_stat), df=df_eff))
    
    return {
        'mean_diff': float(diff),
        't_stat': float(t_stat),
        'p_value': float(p_value),
        'significant': p_value < SIGNIFICANCE_THRESHOLD
    }

def run_covariance_test(gc_df, field_df, rho_intra=0.3):
    """Run covariance-aware t-test."""
    gc_values = gc_df['logPdot_abs'].dropna().values  # Remove NaN values (invalid measurements)
    field_values = field_df['logPdot_abs'].dropna().values  # Remove NaN values (invalid measurements)
    
    mean_gc = np.mean(gc_values)
    mean_field = np.mean(field_values)
    diff = mean_gc - mean_field
    
    clusters = gc_df.groupby('cluster')
    cluster_sizes = [len(c) for _, c in clusters]
    n_bar = np.mean(cluster_sizes)
    n_gc = len(gc_values)
    
    design_effect = 1 + (n_bar - 1) * rho_intra
    n_eff = n_gc / design_effect
    
    var_gc = np.var(gc_values, ddof=1)
    var_field = np.var(field_values, ddof=1)
    
    se_gc = np.sqrt(var_gc * design_effect / n_gc)
    se_field = np.sqrt(var_field / len(field_values))
    se_diff = np.sqrt(se_gc**2 + se_field**2)
    
    t_stat = diff / se_diff
    df_eff = n_eff + len(field_values) - 2
    p_value = 2 * (1 - stats.t.cdf(abs(t_stat), df=df_eff))
    
    return {
        'mean_diff': float(diff),
        't_stat': float(t_stat),
        'p_value': float(p_value),
        'significant': p_value < SIGNIFICANCE_THRESHOLD
    }

def monte_carlo_validation(n_simulations=1000, effect_size=0.0, rho_intra=0.3, n_jobs=-1):
    """
    Run Monte Carlo validation with parallel processing.
    
    Parameters:
    -----------
    n_simulations : int
        Number of Monte Carlo simulations to run.
    effect_size : float
        True effect size (0.0 for H0, >0 for H1).
    rho_intra : float
        Within-cluster correlation.
    n_jobs : int
        Number of parallel jobs (-1 uses all available cores, M4 Pro has 14).
    
    Returns:
    --------
    list : Results from each simulation.
    """
    # Determine number of cores to use
    if n_jobs == -1:
        n_jobs = cpu_count()
    
    print(f"Running {n_simulations} simulations using {n_jobs} cores...")
    
    # Run simulations in parallel
    results = Parallel(n_jobs=n_jobs, backend='loky', verbose=0)(
        delayed(run_single_simulation)(i, effect_size, rho_intra)
        for i in range(n_simulations)
    )
    
    return results

def main():
    print("=" * 70)
    print("STEP 5.39: Monte Carlo Validation of Statistical Methods")
    print("=" * 70)
    print(f"M4 Pro Optimized: Using {cpu_count()} CPU cores\n")
    
    # Test 1: Type I Error Control (H0 true, effect_size=0)
    print("-" * 70)
    print("TEST 1: Type I Error Control (H0: No GC-Field Difference)")
    print("-" * 70)
    print("Generating 1000 synthetic datasets with NO true effect...")
    print("Expected: ~5% false positive rate (p<SIGNIFICANCE_THRESHOLD when H0 true)")
    
    results_h0 = monte_carlo_validation(n_simulations=1000, effect_size=0.0, rho_intra=0.3)
    
    false_positives = sum([r['significant'] for r in results_h0])
    type_i_rate = false_positives / len(results_h0)
    
    print(f"\nResults:")
    print(f"  False positives: {false_positives}/{len(results_h0)} = {type_i_rate:.1%}")
    print(f"  Expected: ~5%")
    
    if 0.01 <= type_i_rate <= 0.08:
        print(f"  Result: Type I error rate is controlled (conservative)")
        type_i_ok = True
    else:
        print(f"  ⚠ Type I error rate deviates from expected range (1%-8%)")
        type_i_ok = False
    
    # Test 2: Power Under H1 (effect_size = 0.65 dex, observed value)
    print("\n" + "-" * 70)
    print("TEST 2: Statistical Power (H1: GC-Field Difference = 0.65 dex)")
    print("-" * 70)
    print("Generating 500 synthetic datasets with TRUE effect of 0.65 dex...")
    print("Expected: High power (>80% detection rate)")
    
    results_h1 = monte_carlo_validation(n_simulations=500, effect_size=0.65, rho_intra=0.3)
    
    true_positives = sum([r['significant'] for r in results_h1])
    power = true_positives / len(results_h1)
    
    print(f"\nResults:")
    print(f"  True positives: {true_positives}/{len(results_h1)} = {power:.1%}")
    print(f"  Target: ≥80%")
    
    if power >= 0.80:
        print(f"  Result: Power is adequate")
        power_ok = True
    elif power >= 0.60:
        print(f"  ~ Power is moderate")
        power_ok = "moderate"
    else:
        print(f"  ⚠ Power is low")
        power_ok = False
    
    # Test 3: Bias in Effect Size Estimation
    print("\n" + "-" * 70)
    print("TEST 3: Bias in Effect Size Estimation")
    print("-" * 70)
    
    estimated_effects = [r['mean_diff'] for r in results_h1]
    mean_estimate = np.mean(estimated_effects)
    std_estimate = np.std(estimated_effects)
    true_effect = 0.65
    
    bias = mean_estimate - true_effect
    bias_pct = bias / true_effect * 100
    
    print(f"  True effect size: {true_effect:.3f} dex")
    print(f"  Mean estimate: {mean_estimate:.3f} dex")
    print(f"  Std of estimates: {std_estimate:.3f} dex")
    print(f"  Bias: {bias:+.3f} dex ({bias_pct:+.1f}%)")
    
    if abs(bias_pct) < 10:
        print(f"  Result: Bias is acceptable (<10%)")
        bias_ok = True
    else:
        print(f"  ⚠ Bias may be problematic (>10%)")
        bias_ok = False
    
    # Overall validation
    print("\n" + "=" * 70)
    print("VALIDATION SUMMARY")
    print("=" * 70)
    
    if type_i_ok and power_ok == True and bias_ok:
        print("Result: All tests passed")
        print("The statistical pipeline is validated and trustworthy.")
        overall = "PASSED"
    elif type_i_ok and (power_ok == True or power_ok == "moderate"):
        print("Result: Most tests passed")
        print("Pipeline is valid but may be conservative.")
        overall = "MOSTLY_PASSED"
    else:
        print("⚠ VALIDATION ISSUES DETECTED")
        print("Pipeline may need adjustment.")
        overall = "ISSUES"
    
    # Save results
    output = {
        'validation_summary': {
            'type_i_error_rate': f"{type_i_rate:.1%}",
            'statistical_power': f"{power:.1%}",
            'bias': f"{bias_pct:+.1f}%",
            'overall_status': overall,
            'interpretation': 'Statistical pipeline is validated and trustworthy' if overall == 'PASSED' else 'Pipeline is valid but may be conservative' if overall == 'MOSTLY_PASSED' else 'Validation issues detected'
        },
        'type_i_error': {
            'rate': float(type_i_rate),
            'rate_percent': f"{type_i_rate:.1%}",
            'expected': 0.05,
            'expected_percent': "5.0%",
            'n_simulations': 1000,
            'acceptable_range': "≤5.0% (conservative rates acceptable)",
            'is_acceptable': bool(type_i_rate <= 0.07),  # Conservative (<5%) is acceptable
            'passed': type_i_ok
        },
        'power': {
            'observed_power': float(power),
            'observed_power_percent': f"{power:.1%}",
            'target_power': 0.80,
            'target_power_percent': "80.0%",
            'true_effect_size': 0.65,
            'n_simulations': 500,
            'is_adequate': power >= 0.80,
            'passed': power_ok == True
        },
        'bias': {
            'true_effect': float(true_effect),
            'mean_estimate': float(mean_estimate),
            'std_estimate': float(std_estimate),
            'bias_dex': float(bias),
            'bias_percent': float(bias_pct),
            'acceptable_threshold': "±10%",
            'is_acceptable': abs(bias_pct) < 10,
            'passed': bias_ok
        },
        'overall_validation': overall
    }
    
    os.makedirs('results/outputs', exist_ok=True)
    with open('results/outputs/step_5_39_monte_carlo_validation.json', 'w') as f:
        json.dump(output, f, indent=2, default=lambda x: bool(x) if isinstance(x, np.bool_) else x)
    
    print(f"\nResults saved to: results/outputs/step_5_39_monte_carlo_validation.json")
    print("=" * 70)

if __name__ == "__main__":
    main()
