#!/usr/bin/env python3
"""
Step 5.46: Spatial Gradient Analysis
=====================================

CRITICAL N-BODY PUSHBACK PREEMPTION

Tests whether Ṗ suppression correlates with projected distance from cluster center.

N-body prediction: If mass segregation dominates, outer pulsars should show 
LESS suppression (closer to field values) than inner pulsars.

TEP prediction: Suppression should be uniform throughout cluster (field-like 
time dilation affects all pulsars regardless of position).

Methodology:
1. Load cluster pulsar data with radial positions
2. For each cluster with sufficient pulsars, test correlation between 
   projected radius and log|Ṗ|
3. Test if slope vs density relation holds across radial bins
4. Compare inner vs outer pulsar populations

Author: M. Smawfield
Date: March 2026
"""

import numpy as np
import pandas as pd
from scipy import stats
import json
from pathlib import Path
import os

# Set random seed for reproducibility
RANDOM_SEED = 42
np.random.seed(RANDOM_SEED)

# Configuration
REPO_ROOT = Path(__file__).resolve().parents[2]
RESULTS_DIR = REPO_ROOT / "results" / "outputs"
DATA_DIR = REPO_ROOT / "data"
OUTPUT_JSON = RESULTS_DIR / "step_5_46_spatial_gradient.json"
OUTPUT_MD = RESULTS_DIR / "step_5_46_spatial_gradient.md"

# Clusters with sufficient pulsars for radial analysis (min 8 for regression)
TARGET_CLUSTERS = {
    "47 Tuc": {"n_pulsars": 25, "r_half": 2.8},  # arcmin
    "Terzan 5": {"n_pulsars": 37, "r_half": 1.2},
    "M15": {"n_pulsars": 10, "r_half": 0.9},
    "M28": {"n_pulsars": 8, "r_half": 1.8},
    "M3": {"n_pulsars": 15, "r_half": 3.5},
    "Omega Cen": {"n_pulsars": 8, "r_half": 4.5},
    "M13": {"n_pulsars": 5, "r_half": 2.5},
    "M62": {"n_pulsars": 6, "r_half": 1.1},
}

# Core collapse clusters (from Harris catalog)
POST_CORE_COLLAPSE = ["M15", "M30", "M62", "NGC 6752", "NGC 6397", "Terzan 5"]


def load_pulsar_data_with_radial():
    """Load pulsar data with radial positions by parsing Freire catalog."""
    freire_path = RESULTS_DIR / "freire_GCpsr.txt"
    
    if not freire_path.exists():
        print("Error: Freire catalog not found")
        return None
    
    with open(freire_path) as f:
        lines = f.readlines()
    
    # Parse the catalog
    rows = []
    current_cluster = None
    
    for line in lines:
        line = line.rstrip('\n')
        if not line.strip() or line.startswith('#'):
            continue
        
        # Check for cluster header (no leading whitespace, not starting with J/B)
        if not line.startswith(('J', 'B')) and not line[0].isspace() and 'arcmin' not in line:
            current_cluster = line.strip()
            continue
        
        # Parse pulsar line
        if line.startswith(('J', 'B')):
            parts = line.split()
            if len(parts) >= 4:
                pulsar_name = parts[0]
                
                # Extract offset (arcmin) - 2nd column
                offset_str = parts[1] if len(parts) > 1 else '*'
                if offset_str == '*' or offset_str == 'i':
                    continue  # Skip pulsars without position data
                
                try:
                    r_arcmin = float(offset_str)
                except ValueError:
                    continue
                
                # Extract Pdot (10^-20) - 4th column
                pdot_str = parts[3] if len(parts) > 3 else '*'
                if pdot_str == '*' or pdot_str == 'i':
                    continue
                
                # Parse Pdot value (handle uncertainty in parentheses)
                import re
                pdot_match = re.match(r'^([+-]?\d+\.?\d*)\(?', pdot_str)
                if not pdot_match:
                    continue
                
                try:
                    pdot_1e20 = float(pdot_match.group(1))
                except ValueError:
                    continue
                
                # Calculate log|Pdot|
                pdot_cgs = pdot_1e20 * 1e-20
                log_pdot_abs = np.log10(abs(pdot_cgs))
                
                rows.append({
                    'cluster': current_cluster,
                    'pulsar': pulsar_name,
                    'r_arcmin': r_arcmin,
                    'pdot_1e20': pdot_1e20,
                    'logPdot_abs': log_pdot_abs
                })
    
    df = pd.DataFrame(rows)
    print(f"Parsed {len(df)} pulsars with radial positions from Freire catalog")
    return df


def load_cluster_properties():
    """Load cluster central density and core radius data."""
    # From step_5_33 or step_5_32 outputs
    cluster_data = {
        "47 Tuc": {"log_rho_c": 4.8, "r_c": 0.36, "r_half": 2.8},
        "Terzan 5": {"log_rho_c": 5.5, "r_c": 0.16, "r_half": 1.2},
        "M15": {"log_rho_c": 5.0, "r_c": 0.14, "r_half": 0.9},
        "M28": {"log_rho_c": 4.5, "r_c": 0.24, "r_half": 1.8},
        "M3": {"log_rho_c": 4.2, "r_c": 0.55, "r_half": 3.5},
        "Omega Cen": {"log_rho_c": 3.8, "r_c": 1.2, "r_half": 4.5},
        "M13": {"log_rho_c": 4.5, "r_c": 0.55, "r_half": 2.5},
        "M62": {"log_rho_c": 5.2, "r_c": 0.18, "r_half": 1.1},
        "M5": {"log_rho_c": 3.5, "r_c": 0.42, "r_half": 2.5},
        "NGC 1851": {"log_rho_c": 4.8, "r_c": 0.25, "r_half": 1.8},
        "M53": {"log_rho_c": 3.5, "r_c": 0.65, "r_half": 3.2},
        "M4": {"log_rho_c": 4.2, "r_c": 0.55, "r_half": 3.5},
    }
    return cluster_data


def analyze_cluster_spatial_gradient(cluster_name, cluster_df, field_mean_logpdot):
    """
    Analyze spatial gradient for a single cluster.
    
    Returns dict with correlation results and inner/outer comparison.
    """
    if len(cluster_df) < 5:
        return None
    
    # Check if radial data available
    if 'r_arcmin' not in cluster_df.columns and 'r_arcmin' not in cluster_df.columns:
        # Try alternative column names
        radial_col = None
        for col in cluster_df.columns:
            if 'r' in col.lower() or 'radius' in col.lower() or 'dist' in col.lower():
                radial_col = col
                break
        if radial_col is None:
            return {"error": "No radial data available", "n_pulsars": len(cluster_df)}
    else:
        radial_col = 'r_arcmin' if 'r_arcmin' in cluster_df.columns else 'r_arcmin'
    
    # Get log|Ṗ| values
    if 'logPdot_abs' in cluster_df.columns:
        logpdot_col = 'logPdot_abs'
    elif 'log_P1' in cluster_df.columns:
        logpdot_col = 'log_P1'
    else:
        return {"error": "No log|Ṗ| data available", "n_pulsars": len(cluster_df)}
    
    # Remove NaN values and convert to numeric
    valid_df = cluster_df[[radial_col, logpdot_col]].copy()
    valid_df[radial_col] = pd.to_numeric(valid_df[radial_col], errors='coerce')
    valid_df[logpdot_col] = pd.to_numeric(valid_df[logpdot_col], errors='coerce')
    valid_df = valid_df.dropna()
    
    if len(valid_df) < 5:
        return {"error": "Insufficient valid data", "n_pulsars": len(valid_df)}
    
    radii = valid_df[radial_col].values.astype(float)
    logpdots = valid_df[logpdot_col].values.astype(float)
    
    # Correlation test
    r_corr, p_corr = stats.pearsonr(radii, logpdots)
    
    # Spearman (non-parametric)
    rho_corr, p_spear = stats.spearmanr(radii, logpdots)
    
    # Split into inner/outer at half-light radius
    cluster_props = load_cluster_properties().get(cluster_name, {})
    r_half = cluster_props.get('r_half', np.median(radii))
    
    inner_mask = radii <= r_half
    outer_mask = radii > r_half
    
    inner_mean = np.mean(logpdots[inner_mask]) if np.sum(inner_mask) > 0 else None
    outer_mean = np.mean(logpdots[outer_mask]) if np.sum(outer_mask) > 0 else None
    
    inner_std = np.std(logpdots[inner_mask]) if np.sum(inner_mask) > 1 else 0
    outer_std = np.std(logpdots[outer_mask]) if np.sum(outer_mask) > 1 else 0
    
    # Compare to field mean
    inner_shift = inner_mean - field_mean_logpdot if inner_mean is not None else None
    outer_shift = outer_mean - field_mean_logpdot if outer_mean is not None else None
    
    # T-test between inner and outer
    if np.sum(inner_mask) > 1 and np.sum(outer_mask) > 1:
        t_stat, p_diff = stats.ttest_ind(
            logpdots[inner_mask], logpdots[outer_mask], 
            equal_var=False  # Welch's t-test
        )
    else:
        t_stat, p_diff = None, None
    
    return {
        "cluster": cluster_name,
        "n_pulsars": len(valid_df),
        "n_inner": int(np.sum(inner_mask)),
        "n_outer": int(np.sum(outer_mask)),
        "r_half_arcmin": float(r_half),
        "pearson_r": float(r_corr),
        "pearson_p": float(p_corr),
        "spearman_rho": float(rho_corr),
        "spearman_p": float(p_spear),
        "inner_mean_logpdot": float(inner_mean) if inner_mean else None,
        "outer_mean_logpdot": float(outer_mean) if outer_mean else None,
        "inner_shift_from_field": float(inner_shift) if inner_shift else None,
        "outer_shift_from_field": float(outer_shift) if outer_shift else None,
        "inner_std": float(inner_std),
        "outer_std": float(outer_std),
        "inner_outer_t_stat": float(t_stat) if t_stat else None,
        "inner_outer_p": float(p_diff) if p_diff else None,
    }


def test_nbody_vs_tep_predictions(results_by_cluster):
    """
    Test whether observed gradients match N-body or TEP predictions.
    
    N-body (mass segregation): Strong positive correlation expected
        - Inner pulsars: MORE acceleration (heavier, more segregated)
        - Outer pulsars: LESS acceleration (lighter, less segregated)
        - Prediction: outer_mean < inner_mean, negative r_corr
    
    TEP (field-like): No correlation expected
        - Time dilation affects all pulsars uniformly
        - Prediction: no significant correlation, similar inner/outer shifts
    """
    clusters_with_data = [r for r in results_by_cluster if 'error' not in r]
    
    if not clusters_with_data:
        return {"error": "No clusters with valid radial data"}
    
    # Count predictions
    nbody_consistent = 0
    tep_consistent = 0
    ambiguous = 0
    
    for r in clusters_with_data:
        # N-body: outer < inner (negative shift means less acceleration)
        if r['inner_shift_from_field'] is not None and r['outer_shift_from_field'] is not None:
            # Check if inner > outer (N-body prediction)
            if r['inner_shift_from_field'] > r['outer_shift_from_field'] + 0.1:
                if r['inner_outer_p'] and r['inner_outer_p'] < 0.05:
                    nbody_consistent += 1
                else:
                    ambiguous += 1
            # Check if similar (TEP prediction)
            elif abs(r['inner_shift_from_field'] - r['outer_shift_from_field']) < 0.2:
                tep_consistent += 1
            else:
                ambiguous += 1
        else:
            ambiguous += 1
    
    n_total = len(clusters_with_data)
    
    # Binomial test: is the distribution significantly different from 50/50?
    # Under null, we'd expect some mix; strong N-body signal would show
    # systematic inner > outer pattern
    
    return {
        "total_clusters_analyzed": n_total,
        "nbody_consistent_clusters": nbody_consistent,
        "tep_consistent_clusters": tep_consistent,
        "ambiguous_clusters": ambiguous,
        "nbody_fraction": nbody_consistent / n_total if n_total > 0 else None,
        "tep_fraction": tep_consistent / n_total if n_total > 0 else None,
        "interpretation": (
            "N-body mass segregation favored" if nbody_consistent > tep_consistent + 2 
            else "TEP uniform suppression favored" if tep_consistent > nbody_consistent + 2
            else "Inconclusive / mixed evidence"
        )
    }


def main_analysis():
    """Main spatial gradient analysis."""
    print("=" * 70)
    print("STEP 5.46: SPATIAL GRADIENT ANALYSIS")
    print("=" * 70)
    print("\nPurpose: Test if Ṗ suppression varies with cluster radius")
    print("N-body prediction: Mass segregation → inner pulsars show MORE suppression")
    print("TEP prediction: Uniform field-like time dilation throughout cluster")
    print()
    
    # Load data
    df = load_pulsar_data_with_radial()
    if df is None:
        print("Error: Could not load pulsar data")
        return None
    
    print(f"Loaded {len(df)} pulsar records")
    
    # Load field mean from step_5_10 - REQUIRED
    field_json = RESULTS_DIR / "step_5_10_pulsar_population_controls.json"
    if not field_json.exists():
        print(f"ERROR: Required input file not found: {field_json}")
        print(f"Spatial gradient analysis requires population control results from step_5_10.")
        raise RuntimeError("Missing required input: step_5_10_pulsar_population_controls.json")
    
    with open(field_json) as f:
        field_data = json.load(f)
    field_mean_logpdot = field_data['base_log10_abs_pdot']['field_mean']
    print(f"Field mean log|Ṗ|: {field_mean_logpdot:.3f}")
    
    # Analyze each cluster
    cluster_results = []
    
    for cluster_name in TARGET_CLUSTERS.keys():
        # Match cluster name variations
        cluster_df = df[df['cluster'].str.contains(
            cluster_name.replace(' ', '\\s*'), 
            case=False, regex=True, na=False
        )]
        
        if len(cluster_df) == 0:
            # Try exact match
            cluster_df = df[df['cluster'] == cluster_name]
        
        if len(cluster_df) > 0:
            print(f"\n  Analyzing {cluster_name}: {len(cluster_df)} pulsars")
            result = analyze_cluster_spatial_gradient(
                cluster_name, cluster_df, field_mean_logpdot
            )
            if result:
                cluster_results.append(result)
                if 'error' not in result:
                    print(f"    Pearson r = {result['pearson_r']:.3f} (p = {result['pearson_p']:.3f})")
                    inner_str = f"{result['inner_shift_from_field']:.3f}" if result['inner_shift_from_field'] is not None else "N/A"
                    outer_str = f"{result['outer_shift_from_field']:.3f}" if result['outer_shift_from_field'] is not None else "N/A"
                    print(f"    Inner shift: {inner_str} dex")
                    print(f"    Outer shift: {outer_str} dex")
                else:
                    print(f"    {result['error']}")
        else:
            print(f"\n  {cluster_name}: No data found")
    
    # Test predictions
    prediction_test = test_nbody_vs_tep_predictions(cluster_results)
    
    print(f"\n{'='*70}")
    print("PREDICTION TEST RESULTS")
    print(f"{'='*70}")
    print(f"Clusters analyzed: {prediction_test.get('total_clusters_analyzed', 0)}")
    print(f"N-body consistent: {prediction_test.get('nbody_consistent_clusters', 0)}")
    print(f"TEP consistent: {prediction_test.get('tep_consistent_clusters', 0)}")
    print(f"Ambiguous: {prediction_test.get('ambiguous_clusters', 0)}")
    print(f"\nInterpretation: {prediction_test.get('interpretation', 'N/A')}")
    
    # Statistical summary
    valid_clusters = [r for r in cluster_results if 'error' not in r]
    mean_corr = None
    t_stat = None
    p_mean = None
    positive_corr = 0
    negative_corr = 0
    
    if valid_clusters:
        all_correlations = [r['pearson_r'] for r in valid_clusters]
        mean_corr = np.mean(all_correlations)
        
        # One-sample t-test: is mean correlation significantly different from 0?
        if len(all_correlations) > 1:
            t_stat, p_mean = stats.ttest_1samp(all_correlations, 0)
        
        # Check for systematic trend
        positive_corr = sum(1 for r in all_correlations if r > 0)
        negative_corr = sum(1 for r in all_correlations if r < 0)
        
        print(f"\nMean correlation across clusters: {mean_corr:.3f}" if mean_corr else "\nNo valid correlation data")
        if t_stat is not None:
            print(f"Significance vs zero: t = {t_stat:.2f}, p = {p_mean:.3f}")
        print(f"Positive correlations: {positive_corr}, Negative: {negative_corr}")
    
    # Save results
    output = {
        "timestamp": pd.Timestamp.now().isoformat(),
        "method": "Spatial gradient analysis: testing correlation between projected radius and log|Ṗ|",
        "field_mean_logpdot": float(field_mean_logpdot),
        "cluster_results": cluster_results,
        "prediction_test": prediction_test,
        "statistical_summary": {
            "mean_correlation": float(mean_corr) if valid_clusters else None,
            "t_statistic_vs_zero": float(t_stat) if t_stat else None,
            "p_value_vs_zero": float(p_mean) if p_mean else None,
            "n_positive_corr": positive_corr if valid_clusters else None,
            "n_negative_corr": negative_corr if valid_clusters else None,
        },
        "nbody_prediction": "Mass segregation: inner pulsars show MORE acceleration (higher Ṗ)",
        "tep_prediction": "Uniform suppression: no radial gradient in Ṗ",
        "conclusion": prediction_test.get('interpretation', 'N/A')
    }
    
    os.makedirs(RESULTS_DIR, exist_ok=True)
    with open(OUTPUT_JSON, 'w') as f:
        json.dump(output, f, indent=2)
    
    # Generate markdown report
    mean_corr_str = f"{mean_corr:.3f}" if mean_corr is not None else "N/A"
    
    md_content = f"""# Spatial Gradient Analysis Report

## Purpose
Test whether Ṗ suppression correlates with projected distance from cluster center.

- **N-body prediction**: Mass segregation → outer pulsars show LESS suppression
- **TEP prediction**: Uniform field-like time dilation → no radial gradient

## Results Summary

| Metric | Value |
|--------|-------|
| Clusters analyzed | {prediction_test.get('total_clusters_analyzed', 0)} |
| N-body consistent | {prediction_test.get('nbody_consistent_clusters', 0)} |
| TEP consistent | {prediction_test.get('tep_consistent_clusters', 0)} |
| Ambiguous | {prediction_test.get('ambiguous_clusters', 0)} |
| Mean correlation | {mean_corr_str} |

## Conclusion

{prediction_test.get('interpretation', 'N/A')}

## Individual Cluster Results

| Cluster | N | Pearson r | p-value | Inner Shift | Outer Shift | Status |
|---------|---|-----------|---------|-------------|-------------|--------|
"""
    
    for r in cluster_results:
        if 'error' not in r:
            inner_val = r.get('inner_shift_from_field')
            outer_val = r.get('outer_shift_from_field')
            inner_str = f"{inner_val:.3f}" if inner_val is not None else "N/A"
            outer_str = f"{outer_val:.3f}" if outer_val is not None else "N/A"
            status = "N/A"
            if inner_val is not None and outer_val is not None:
                if inner_val > outer_val + 0.1:
                    status = "N-body"
                elif abs(inner_val - outer_val) < 0.2:
                    status = "TEP"
                else:
                    status = "Ambiguous"
            else:
                status = "Ambiguous"
            md_content += f"| {r['cluster']} | {r['n_pulsars']} | {r['pearson_r']:.3f} | {r['pearson_p']:.3f} | {inner_str} | {outer_str} | {status} |\n"
        else:
            md_content += f"| {r.get('cluster', 'Unknown')} | {r.get('n_pulsars', 0)} | N/A | N/A | N/A | N/A | {r['error']} |\n"
    
    md_content += f"""

## Implications for N-Body Pushback

This analysis directly addresses the "messy dynamics" critique by testing whether
observed suppression patterns match N-body predictions. If mass segregation dominated:

1. Inner pulsars (massive, segregated) should show HIGHER Ṗ
2. Outer pulsars (lighter, unsegregated) should show LOWER Ṗ
3. Strong negative correlation between radius and Ṗ expected

The observed pattern provides a quantitative test of this alternative explanation.

---

*Report generated by step_5_46_spatial_gradient.py*
"""
    
    with open(OUTPUT_MD, 'w') as f:
        f.write(md_content)
    
    print(f"\n{'='*70}")
    print(f"Results saved to: {OUTPUT_JSON}")
    print(f"Report saved to: {OUTPUT_MD}")
    print(f"{'='*70}")
    
    return output


if __name__ == "__main__":
    main_analysis()
