#!/usr/bin/env python3
"""
Step 55: Cluster Bootstrap

Bootstrap by cluster, not by pulsar, to ensure significance is not inflated by
within-cluster dependence.  Pulsars within the same cluster share environment,
so they are not independent draws.

Author: M. Smawfield
Date: 2026-06-23
"""

import csv
import json
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from scipy import stats

REPO_ROOT = Path(__file__).resolve().parents[2]
RESULTS_DIR = REPO_ROOT / "results" / "outputs"
STEP_5_10_CSV = RESULTS_DIR / "step_02_pulsar_population_controls.csv"

OUT_JSON = RESULTS_DIR / "step_55_cluster_bootstrap.json"
OUT_MD = RESULTS_DIR / "step_55_cluster_bootstrap.md"

CLUSTER_PARAMS = {
    "Terzan 5": {"rho_c": 5.50}, "47 Tuc (NGC 104)": {"rho_c": 4.88},
    "NGC 6517": {"rho_c": 5.80}, "M28 (NGC 6626)": {"rho_c": 4.52},
    "M62 (NGC 6266)": {"rho_c": 5.16}, "M13 (NGC 6205)": {"rho_c": 3.79},
    "M15 (NGC 7078)": {"rho_c": 5.05}, "M5 (NGC 5904)": {"rho_c": 3.53},
    "Terzan 1": {"rho_c": 5.00}, "NGC 6752": {"rho_c": 4.30},
    "M2 (NGC 7089)": {"rho_c": 4.15}, "Omega Centauri (NGC 5139)": {"rho_c": 3.12},
    "M53 (NGC 5024)": {"rho_c": 2.96}, "M3 (NGC 5272)": {"rho_c": 3.68},
    "M71 (NGC 6838)": {"rho_c": 2.29}, "NGC 6397": {"rho_c": 5.68},
    "NGC 1851": {"rho_c": 5.09}, "NGC 6522": {"rho_c": 5.50},
    "NGC 6544": {"rho_c": 5.20}, "NGC 6624": {"rho_c": 5.60},
    "NGC 6760": {"rho_c": 3.80}, "M22 (NGC 6656)": {"rho_c": 2.97},
    "M80 (NGC 6093)": {"rho_c": 4.79}, "M92 (NGC 6341)": {"rho_c": 4.30},
    "NGC 6712": {"rho_c": 3.70}, "NGC 6652": {"rho_c": 4.50},
    "M14 (NGC 6402)": {"rho_c": 3.44}, "NGC 6539": {"rho_c": 3.30},
    "M4 (NGC 6121)": {"rho_c": 2.85}, "NGC 6440": {"rho_c": 5.10},
    "NGC 6441": {"rho_c": 5.00}, "NGC 6316": {"rho_c": 4.80},
    "M30 (NGC 7099)": {"rho_c": 4.20},
}


def load_data():
    gc_rows = []
    field_rows = []
    with open(STEP_5_10_CSV, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                r = {
                    'environment': row['environment'],
                    'cluster': row['cluster'],
                    'logPdot_abs': float(row['logPdot_abs']),
                    'logP': float(row['logP']),
                    'log_b_proxy': float(row['log_b_proxy']),
                }
                if r['environment'] == 'globular_cluster':
                    gc_rows.append(r)
                elif r['environment'] == 'field':
                    field_rows.append(r)
            except (ValueError, KeyError):
                continue
    return gc_rows, field_rows


def cluster_bootstrap_raw(gc_rows, field_rows, n_bootstrap=5000, seed=42):
    """Bootstrap clusters for raw mean difference."""
    rng = np.random.default_rng(seed)
    clusters = defaultdict(list)
    for r in gc_rows:
        clusters[r['cluster']].append(r)
    cluster_names = list(clusters.keys())
    n_clusters = len(cluster_names)

    field_vals = np.array([r['logPdot_abs'] for r in field_rows])
    field_mean = np.mean(field_vals)

    diffs = []
    for _ in range(n_bootstrap):
        sampled_names = rng.choice(cluster_names, size=n_clusters, replace=True)
        sampled_gc = []
        for name in sampled_names:
            sampled_gc.extend(clusters[name])
        gc_mean = np.mean([r['logPdot_abs'] for r in sampled_gc])
        diffs.append(gc_mean - field_mean)

    diffs = np.array(diffs)
    return {
        'mean_diff': float(np.mean(diffs)),
        'std_diff': float(np.std(diffs)),
        'ci_lower': float(np.percentile(diffs, 2.5)),
        'ci_upper': float(np.percentile(diffs, 97.5)),
        'p_value': float(2 * min(np.mean(diffs <= 0), np.mean(diffs >= 0))) if np.mean(diffs) != 0 else 1.0,
    }


def cluster_bootstrap_matched(gc_rows, field_rows, match_keys, n_matches=5, n_bootstrap=5000, seed=43):
    """Bootstrap clusters for matched residual."""
    rng = np.random.default_rng(seed)
    clusters = defaultdict(list)
    for r in gc_rows:
        clusters[r['cluster']].append(r)
    cluster_names = list(clusters.keys())
    n_clusters = len(cluster_names)

    field_vals = np.array([r['logPdot_abs'] for r in field_rows])
    field_features = np.array([[r[k] for k in match_keys] for r in field_rows])
    combined = np.vstack([field_features, field_features])  # dummy for std
    means = np.mean(field_features, axis=0)
    stds = np.std(field_features, axis=0)
    stds[stds == 0] = 1.0
    field_std = (field_features - means) / stds

    residuals = []
    for _ in range(n_bootstrap):
        sampled_names = rng.choice(cluster_names, size=n_clusters, replace=True)
        sampled_gc = []
        for name in sampled_names:
            sampled_gc.extend(clusters[name])

        gc_vals = np.array([r['logPdot_abs'] for r in sampled_gc])
        gc_features = np.array([[r[k] for k in match_keys] for r in sampled_gc])
        gc_std = (gc_features - means) / stds

        resid = []
        for i, gc_f in enumerate(gc_std):
            dist = np.sqrt(np.sum((field_std - gc_f) ** 2, axis=1))
            closest = np.argsort(dist)[:n_matches]
            field_mean = field_vals[closest].mean()
            resid.append(gc_vals[i] - field_mean)
        residuals.append(np.mean(resid))

    residuals = np.array(residuals)
    return {
        'mean_residual': float(np.mean(residuals)),
        'std_residual': float(np.std(residuals)),
        'ci_lower': float(np.percentile(residuals, 2.5)),
        'ci_upper': float(np.percentile(residuals, 97.5)),
        'p_value': float(2 * min(np.mean(residuals <= 0), np.mean(residuals >= 0))) if np.mean(residuals) != 0 else 1.0,
    }


def cluster_bootstrap_gamma(gc_rows, field_rows, match_keys=None, n_matches=5, n_bootstrap=5000, seed=44):
    """Bootstrap clusters for density scaling slope."""
    rng = np.random.default_rng(seed)
    clusters = defaultdict(list)
    for r in gc_rows:
        clusters[r['cluster']].append(r)
    cluster_names = list(clusters.keys())
    n_clusters = len(cluster_names)

    field_vals = np.array([r['logPdot_abs'] for r in field_rows])
    field_features = None
    if match_keys:
        field_features = np.array([[r[k] for k in match_keys] for r in field_rows])
        means = np.mean(field_features, axis=0)
        stds = np.std(field_features, axis=0)
        stds[stds == 0] = 1.0
        field_std = (field_features - means) / stds

    slopes = []
    for _ in range(n_bootstrap):
        sampled_names = rng.choice(cluster_names, size=n_clusters, replace=True)

        cluster_means = []
        densities = []
        for name in sampled_names:
            if name not in CLUSTER_PARAMS:
                continue
            gc_list = clusters[name]
            if match_keys:
                gc_vals = np.array([r['logPdot_abs'] for r in gc_list])
                gc_features = np.array([[r[k] for k in match_keys] for r in gc_list])
                gc_std = (gc_features - means) / stds
                resid = []
                for i, gc_f in enumerate(gc_std):
                    dist = np.sqrt(np.sum((field_std - gc_f) ** 2, axis=1))
                    closest = np.argsort(dist)[:n_matches]
                    field_mean = field_vals[closest].mean()
                    resid.append(gc_vals[i] - field_mean)
                cmean = np.mean(resid)
            else:
                cmean = np.mean([r['logPdot_abs'] for r in gc_list])
            cluster_means.append(cmean)
            densities.append(CLUSTER_PARAMS[name]['rho_c'])

        if len(densities) >= 3:
            slope, _, _, _, _ = stats.linregress(densities, cluster_means)
            slopes.append(slope)

    slopes = np.array(slopes)
    return {
        'mean_gamma': float(np.mean(slopes)),
        'std_gamma': float(np.std(slopes)),
        'ci_lower': float(np.percentile(slopes, 2.5)),
        'ci_upper': float(np.percentile(slopes, 97.5)),
        'p_value': float(2 * min(np.mean(slopes <= 0), np.mean(slopes >= 0))) if np.mean(slopes) != 0 else 1.0,
    }


def main():
    print("=" * 78)
    print("STEP 5.75: CLUSTER BOOTSTRAP")
    print("=" * 78)
    print()

    gc_rows, field_rows = load_data()
    print(f"Loaded {len(gc_rows)} GC pulsars, {len(field_rows)} field pulsars")
    print()

    # 1. Raw mean difference
    print("--- Raw mean difference (cluster bootstrap) ---")
    raw_boot = cluster_bootstrap_raw(gc_rows, field_rows, n_bootstrap=5000)
    print(f"  Mean diff: {raw_boot['mean_diff']:+.3f} dex")
    print(f"  Std:       {raw_boot['std_diff']:.3f} dex")
    print(f"  95% CI:    [{raw_boot['ci_lower']:+.3f}, {raw_boot['ci_upper']:+.3f}]")
    print(f"  p-value:   {raw_boot['p_value']:.4f}")
    print()

    # 2. Period-only matched residual
    print("--- Period-only matched residual (cluster bootstrap) ---")
    po_boot = cluster_bootstrap_matched(gc_rows, field_rows, ['logP'], n_matches=5, n_bootstrap=5000)
    print(f"  Mean resid: {po_boot['mean_residual']:+.3f} dex")
    print(f"  Std:        {po_boot['std_residual']:.3f} dex")
    print(f"  95% CI:     [{po_boot['ci_lower']:+.3f}, {po_boot['ci_upper']:+.3f}]")
    print(f"  p-value:    {po_boot['p_value']:.4f}")
    print()

    # 3. Period+B-field matched residual
    print("--- Period+B-field matched residual (cluster bootstrap) ---")
    pb_boot = cluster_bootstrap_matched(gc_rows, field_rows, ['logP', 'log_b_proxy'], n_matches=5, n_bootstrap=5000)
    print(f"  Mean resid: {pb_boot['mean_residual']:+.3f} dex")
    print(f"  Std:        {pb_boot['std_residual']:.3f} dex")
    print(f"  95% CI:     [{pb_boot['ci_lower']:+.3f}, {pb_boot['ci_upper']:+.3f}]")
    print(f"  p-value:    {pb_boot['p_value']:.4f}")
    print()

    # 4. Gamma raw
    print("--- Gamma raw (cluster bootstrap) ---")
    gamma_raw_boot = cluster_bootstrap_gamma(gc_rows, field_rows, match_keys=None, n_bootstrap=5000)
    print(f"  Mean Gamma: {gamma_raw_boot['mean_gamma']:.3f}")
    print(f"  Std:        {gamma_raw_boot['std_gamma']:.3f}")
    print(f"  95% CI:     [{gamma_raw_boot['ci_lower']:.3f}, {gamma_raw_boot['ci_upper']:.3f}]")
    print(f"  p-value:    {gamma_raw_boot['p_value']:.4f}")
    print()

    # 5. Gamma period-only
    print("--- Gamma period-only (cluster bootstrap) ---")
    gamma_po_boot = cluster_bootstrap_gamma(gc_rows, field_rows, match_keys=['logP'], n_matches=5, n_bootstrap=5000)
    print(f"  Mean Gamma: {gamma_po_boot['mean_gamma']:.3f}")
    print(f"  Std:        {gamma_po_boot['std_gamma']:.3f}")
    print(f"  95% CI:     [{gamma_po_boot['ci_lower']:.3f}, {gamma_po_boot['ci_upper']:.3f}]")
    print(f"  p-value:    {gamma_po_boot['p_value']:.4f}")
    print()

    results = {
        'timestamp_utc': datetime.now(timezone.utc).isoformat(),
        'raw': raw_boot,
        'period_only': po_boot,
        'period_bfield': pb_boot,
        'gamma_raw': gamma_raw_boot,
        'gamma_period_only': gamma_po_boot,
    }

    with open(OUT_JSON, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"Saved JSON: {OUT_JSON}")

    md = f"""# Step 55: Cluster Bootstrap

**Generated:** {results['timestamp_utc']}

Cluster-level bootstrap (resampling clusters with replacement) to account for
within-cluster dependence.

## Results

| Test | Mean | Std | 95% CI Lower | 95% CI Upper | p-value |
|------|------|-----|--------------|--------------|---------|
| Raw diff (dex) | {raw_boot['mean_diff']:+.3f} | {raw_boot['std_diff']:.3f} | {raw_boot['ci_lower']:+.3f} | {raw_boot['ci_upper']:+.3f} | {raw_boot['p_value']:.4f} |
| Period-only resid (dex) | {po_boot['mean_residual']:+.3f} | {po_boot['std_residual']:.3f} | {po_boot['ci_lower']:+.3f} | {po_boot['ci_upper']:+.3f} | {po_boot['p_value']:.4f} |
| Period+B resid (dex) | {pb_boot['mean_residual']:+.3f} | {pb_boot['std_residual']:.3f} | {pb_boot['ci_lower']:+.3f} | {pb_boot['ci_upper']:+.3f} | {pb_boot['p_value']:.4f} |
| Gamma raw | {gamma_raw_boot['mean_gamma']:.3f} | {gamma_raw_boot['std_gamma']:.3f} | {gamma_raw_boot['ci_lower']:.3f} | {gamma_raw_boot['ci_upper']:.3f} | {gamma_raw_boot['p_value']:.4f} |
| Gamma period-only | {gamma_po_boot['mean_gamma']:.3f} | {gamma_po_boot['std_gamma']:.3f} | {gamma_po_boot['ci_lower']:.3f} | {gamma_po_boot['ci_upper']:.3f} | {gamma_po_boot['p_value']:.4f} |

## Interpretation

- If cluster-bootstrap p-values are **larger** than pulsar-bootstrap p-values,
  the original analysis was inflating significance by treating within-cluster
  pulsars as independent.
- The 95% CI should be used for robust inference.
"""

    with open(OUT_MD, 'w') as f:
        f.write(md)
    print(f"Saved Markdown: {OUT_MD}")
    print("Done.")


if __name__ == '__main__':
    main()
