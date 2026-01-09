#!/usr/bin/env python3
"""
Step 5.31: Per-Cluster Controlled Residuals

For each cluster, compute the CONTROLLED residual by matching
GC pulsars to field pulsars with similar period and B-proxy.

This tests whether the 0.13 dex residual is truly constant across clusters.
"""

import numpy as np
import pandas as pd
from scipy import stats
import json
from pathlib import Path
from datetime import datetime, timezone

REPO_ROOT = Path(__file__).resolve().parents[2]
RESULTS_DIR = REPO_ROOT / "results" / "outputs"
OUT_JSON = RESULTS_DIR / "step_5_31_per_cluster_controlled_residuals.json"
OUT_MD = RESULTS_DIR / "step_5_31_per_cluster_controlled_residuals.md"


def match_to_field(gc_pulsars, field_pulsars, n_matches=5):
    """
    For each GC pulsar, find the n_matches closest field pulsars
    by period and B-proxy, and compute the residual.
    """
    residuals = []
    
    for _, gc in gc_pulsars.iterrows():
        # Distance in (logP, log_b_proxy) space
        dP = (field_pulsars['logP'].values - gc['logP'])**2
        dB = (field_pulsars['log_b_proxy'].values - gc['log_b_proxy'])**2
        dist = np.sqrt(dP + dB)
        
        # Find closest matches
        closest_idx = np.argsort(dist)[:n_matches]
        closest = field_pulsars.iloc[closest_idx]
        
        # Compute residual: GC - mean(matched field)
        field_mean = closest['logPdot_abs'].mean()
        residual = gc['logPdot_abs'] - field_mean
        residuals.append(residual)
    
    return np.array(residuals)


def compute_cluster_controlled_residual(cluster_name, gc_cluster, field_pulsars, n_matches=5):
    """
    Compute the controlled residual for a single cluster.
    """
    if len(gc_cluster) < 2:
        return None
    
    residuals = match_to_field(gc_cluster, field_pulsars, n_matches)
    
    return {
        'cluster': cluster_name,
        'n_pulsars': len(gc_cluster),
        'controlled_residual': float(np.mean(residuals)),
        'residual_std': float(np.std(residuals)),
        'residual_sem': float(np.std(residuals) / np.sqrt(len(residuals))),
    }


def main():
    print("="*70)
    print("PER-CLUSTER CONTROLLED RESIDUALS")
    print("="*70)
    
    # Load data
    csv_path = REPO_ROOT / "results" / "outputs" / "step_5_10_pulsar_population_controls.csv"
    df = pd.read_csv(csv_path)
    
    field = df[df['environment'] == 'field'].copy()
    gc = df[df['environment'] == 'globular_cluster'].copy()
    
    print(f"\nField MSPs: {len(field)}")
    print(f"GC MSPs: {len(gc)}")
    
    # Compute distance column for matching
    field['dist'] = 0  # Placeholder
    
    # Get per-cluster controlled residuals
    results = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "method": "Period + B-proxy matching (5 nearest neighbors)",
        "clusters": {},
    }
    
    cluster_names = gc['cluster'].value_counts().index
    
    print(f"\n{'Cluster':<30} {'N':>4} {'Controlled Residual':>20} {'SEM':>10}")
    print("-"*70)
    
    for cluster_name in cluster_names:
        gc_cluster = gc[gc['cluster'] == cluster_name]
        
        result = compute_cluster_controlled_residual(cluster_name, gc_cluster, field)
        
        if result:
            results["clusters"][cluster_name] = result
            print(f"{cluster_name:<30} {result['n_pulsars']:>4} {result['controlled_residual']:>+18.3f} dex {result['residual_sem']:>8.3f}")
    
    # Summary statistics
    controlled_residuals = [r['controlled_residual'] for r in results['clusters'].values()]
    n_pulsars = [r['n_pulsars'] for r in results['clusters'].values()]
    
    # Weighted mean
    weighted_mean = np.average(controlled_residuals, weights=n_pulsars)
    
    # Is it constant? Check std
    residual_std = np.std(controlled_residuals)
    residual_mean = np.mean(controlled_residuals)
    
    print(f"\n{'='*70}")
    print("SUMMARY")
    print(f"{'='*70}")
    print(f"\nNumber of clusters: {len(controlled_residuals)}")
    print(f"Mean controlled residual: {residual_mean:+.3f} dex")
    print(f"Weighted mean: {weighted_mean:+.3f} dex")
    print(f"Std across clusters: {residual_std:.3f} dex")
    print(f"Range: {min(controlled_residuals):+.3f} to {max(controlled_residuals):+.3f} dex")
    
    # Test: are controlled residuals constant (no density correlation)?
    density_map = {
        'Terzan 5': 5.5, '47 Tuc (NGC 104)': 4.88, 'NGC 6517': 5.8,
        'M28 (NGC 6626)': 4.52, 'M62 (NGC 6266)': 5.16, 'M13 (NGC 6205)': 3.79,
        'M15 (NGC 7078)': 5.05, 'M5 (NGC 5904)': 3.53, 'NGC 6752': 4.30,
        'M2 (NGC 7089)': 4.15, 'Terzan 1': 5.0, 'Omega Centauri (NGC 5139)': 3.12,
        'M53 (NGC 5024)': 2.96, 'M71 (NGC 6838)': 2.29, 'M3 (NGC 5272)': 3.68,
        'NGC 6397': 5.68, 'NGC 1851': 5.09, 'NGC 6522': 5.50, 'NGC 6544': 5.20,
        'NGC 6624': 5.60, 'NGC 6760': 3.80, 'M22 (NGC 6656)': 2.97,
        'M80 (NGC 6093)': 4.79, 'M92 (NGC 6341)': 4.30, 'NGC 6712': 3.70,
        'NGC 6652': 4.50, 'M14 (NGC 6402)': 3.44, 'NGC 6539': 3.30, 'M4 (NGC 6121)': 2.85,
    }
    
    # Get density and controlled residual for each cluster
    densities = []
    ctrl_residuals = []
    for name, data in results['clusters'].items():
        if name in density_map:
            densities.append(density_map[name])
            ctrl_residuals.append(data['controlled_residual'])
    
    if len(densities) >= 5:
        r_ctrl, p_ctrl = stats.pearsonr(densities, ctrl_residuals)
        
        print(f"\n--- DENSITY CORRELATION TEST ---")
        print(f"Controlled residual vs log(ρc):")
        print(f"  Pearson r = {r_ctrl:.3f}")
        print(f"  p-value   = {p_ctrl:.4f}")
        
        if abs(r_ctrl) < 0.3 and p_ctrl > 0.05:
            print(f"\n✓ CONTROLLED residuals do NOT correlate with density!")
            print(f"  This confirms the Universality Constraint.")
        else:
            print(f"\n⚠️ CONTROLLED residuals still show density correlation")
            print(f"  The Universality Constraint may not hold.")
        
        results["summary"] = {
            "n_clusters": len(controlled_residuals),
            "mean_controlled_residual": float(residual_mean),
            "weighted_mean": float(weighted_mean),
            "std_across_clusters": float(residual_std),
            "density_correlation_r": float(r_ctrl),
            "density_correlation_p": float(p_ctrl),
        }
    
    # Save
    with open(OUT_JSON, 'w') as f:
        json.dump(results, f, indent=2)
    
    # Markdown
    md = f"""# Per-Cluster Controlled Residuals

**Generated:** {results['timestamp_utc']}

## Method
For each GC pulsar, find the 5 closest field MSPs by period and B-proxy,
then compute the residual (GC - mean of matched field).

## Results

| Cluster | N | Controlled Residual |
|---------|---|---------------------|
"""
    
    for name, data in sorted(results['clusters'].items(), key=lambda x: -x[1]['n_pulsars']):
        md += f"| {name} | {data['n_pulsars']} | {data['controlled_residual']:+.3f} |\n"
    
    md += f"""
## Summary

- Mean controlled residual: {residual_mean:+.3f} dex
- Std across clusters: {residual_std:.3f} dex
- Density correlation: r = {results.get('summary', {}).get('density_correlation_r', 'N/A'):.3f}
"""
    
    with open(OUT_MD, 'w') as f:
        f.write(md)
    
    print(f"\nResults saved to:")
    print(f"  {OUT_JSON}")
    print(f"  {OUT_MD}")


if __name__ == "__main__":
    main()
