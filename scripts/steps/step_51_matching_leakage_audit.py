#!/usr/bin/env python3
"""
Step 51: Matching Leakage Audit

Compare raw, period-only, period+non-outcome, period+B-field, hybrid expanded,
and cluster-level residualized outputs to quantify how much each control
procedure attenuates the true signal.

This follows the injection-recovery finding that period+B-field matching
suppresses ~50% of the injected effect because B_proxy ~ sqrt(P*Pdot) partly
conditions on the outcome variable.

Author: M. Smawfield
Date: 2026-06-23
"""

import csv
import json
import math
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from scipy import stats

REPO_ROOT = Path(__file__).resolve().parents[2]
RESULTS_DIR = REPO_ROOT / "results" / "outputs"
STEP_5_10_CSV = RESULTS_DIR / "step_02_pulsar_population_controls.csv"

OUT_JSON = RESULTS_DIR / "step_51_matching_leakage_audit.json"
OUT_MD = RESULTS_DIR / "step_51_matching_leakage_audit.md"

CLUSTER_PARAMS = {
    "Terzan 5":         {"rho_c": 5.50}, "47 Tuc (NGC 104)": {"rho_c": 4.88},
    "NGC 6517":         {"rho_c": 5.80}, "M28 (NGC 6626)":   {"rho_c": 4.52},
    "M62 (NGC 6266)":   {"rho_c": 5.16}, "M13 (NGC 6205)":   {"rho_c": 3.79},
    "M15 (NGC 7078)":   {"rho_c": 5.05}, "M5 (NGC 5904)":    {"rho_c": 3.53},
    "Terzan 1":         {"rho_c": 5.00}, "NGC 6752":         {"rho_c": 4.30},
    "M2 (NGC 7089)":    {"rho_c": 4.15}, "Omega Centauri (NGC 5139)": {"rho_c": 3.12},
    "M53 (NGC 5024)":   {"rho_c": 2.96}, "M3 (NGC 5272)":    {"rho_c": 3.68},
    "M71 (NGC 6838)":   {"rho_c": 2.29}, "NGC 6397":         {"rho_c": 5.68},
    "NGC 1851":         {"rho_c": 5.09}, "NGC 6522":         {"rho_c": 5.50},
    "NGC 6544":         {"rho_c": 5.20}, "NGC 6624":         {"rho_c": 5.60},
    "NGC 6760":         {"rho_c": 3.80}, "M22 (NGC 6656)":   {"rho_c": 2.97},
    "M80 (NGC 6093)":   {"rho_c": 4.79}, "M92 (NGC 6341)":   {"rho_c": 4.30},
    "NGC 6712":         {"rho_c": 3.70}, "NGC 6652":         {"rho_c": 4.50},
    "M14 (NGC 6402)":   {"rho_c": 3.44}, "NGC 6539":         {"rho_c": 3.30},
    "M4 (NGC 6121)":    {"rho_c": 2.85}, "NGC 6440":         {"rho_c": 5.10},
    "NGC 6441":         {"rho_c": 5.00}, "NGC 6316":         {"rho_c": 4.80},
    "M30 (NGC 7099)":   {"rho_c": 4.20},
}


def load_data():
    gc_rows = []
    field_rows = []
    with open(STEP_5_10_CSV, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                r = {
                    'source': row['source'],
                    'environment': row['environment'],
                    'cluster': row['cluster'],
                    'name': row['name'],
                    'P0_s': float(row['P0_s']),
                    'P1_sps': float(row['P1_sps']),
                    'logP': float(row['logP']),
                    'logPdot_abs': float(row['logPdot_abs']),
                    'log_b_proxy': float(row['log_b_proxy']),
                    'log_tau_c': float(row['log_tau_c']),
                }
                if r['environment'] == 'globular_cluster':
                    gc_rows.append(r)
                elif r['environment'] == 'field':
                    field_rows.append(r)
            except (ValueError, KeyError):
                continue
    return gc_rows, field_rows


def raw_difference(gc_rows, field_rows):
    gc_vals = np.array([r['logPdot_abs'] for r in gc_rows])
    field_vals = np.array([r['logPdot_abs'] for r in field_rows])
    return float(np.mean(gc_vals) - np.mean(field_vals))


def matched_residual(gc_rows, field_rows, match_keys, n_matches=5):
    if not gc_rows or not field_rows:
        return 0.0
    gc_vals = np.array([r['logPdot_abs'] for r in gc_rows])
    field_vals = np.array([r['logPdot_abs'] for r in field_rows])

    gc_features = np.array([[r[k] for k in match_keys] for r in gc_rows])
    field_features = np.array([[r[k] for k in match_keys] for r in field_rows])

    combined = np.vstack([gc_features, field_features])
    means = np.mean(combined, axis=0)
    stds = np.std(combined, axis=0)
    stds[stds == 0] = 1.0

    gc_std = (gc_features - means) / stds
    field_std = (field_features - means) / stds

    residuals = []
    for i, gc_f in enumerate(gc_std):
        dist = np.sqrt(np.sum((field_std - gc_f) ** 2, axis=1))
        closest = np.argsort(dist)[:n_matches]
        field_mean = field_vals[closest].mean()
        residuals.append(gc_vals[i] - field_mean)

    return float(np.mean(residuals))


def cluster_residualized(gc_rows, field_rows, match_keys, n_matches=5):
    """Per-cluster controlled residual, then average across clusters."""
    clusters = defaultdict(list)
    for r in gc_rows:
        clusters[r['cluster']].append(r)

    cluster_residuals = []
    for cluster_name, gc_list in clusters.items():
        if len(gc_list) < 2:
            continue
        resid = matched_residual(gc_list, field_rows, match_keys, n_matches)
        cluster_residuals.append(resid)

    if not cluster_residuals:
        return 0.0
    return float(np.mean(cluster_residuals))


def density_scaling_slope(gc_rows, field_rows, residual_mode='raw', match_keys=None, n_matches=5):
    clusters = defaultdict(list)
    for r in gc_rows:
        clusters[r['cluster']].append(r)

    cluster_means = []
    densities = []
    for cluster_name, gc_list in clusters.items():
        if cluster_name not in CLUSTER_PARAMS:
            continue
        rho_c = CLUSTER_PARAMS[cluster_name]['rho_c']
        if residual_mode == 'raw':
            vals = np.array([r['logPdot_abs'] for r in gc_list])
            cmean = np.mean(vals)
        elif residual_mode == 'matched':
            cmean = matched_residual(gc_list, field_rows, match_keys, n_matches)
        elif residual_mode == 'cluster_residualized':
            cmean = cluster_residualized(gc_list, field_rows, match_keys, n_matches)
        else:
            raise ValueError(f"Unknown residual_mode: {residual_mode}")
        cluster_means.append(cmean)
        densities.append(rho_c)

    if len(densities) < 3:
        return None, None, None

    slope, intercept, r_value, p_value, std_err = stats.linregress(densities, cluster_means)
    return float(slope), float(std_err), float(r_value ** 2)


def run_variant(name, gc_rows, field_rows, mode, match_keys=None, n_matches=5):
    if mode == 'raw':
        resid = raw_difference(gc_rows, field_rows)
        gamma, gamma_err, gamma_r2 = density_scaling_slope(gc_rows, field_rows, 'raw')
    elif mode == 'matched':
        resid = matched_residual(gc_rows, field_rows, match_keys, n_matches)
        gamma, gamma_err, gamma_r2 = density_scaling_slope(gc_rows, field_rows, 'matched', match_keys, n_matches)
    elif mode == 'cluster_residualized':
        resid = cluster_residualized(gc_rows, field_rows, match_keys, n_matches)
        gamma, gamma_err, gamma_r2 = density_scaling_slope(gc_rows, field_rows, 'cluster_residualized', match_keys, n_matches)
    else:
        raise ValueError(f"Unknown mode: {mode}")

    return {
        'name': name,
        'mode': mode,
        'match_keys': match_keys,
        'n_matches': n_matches,
        'residual_dex': resid,
        'gamma': gamma,
        'gamma_err': gamma_err,
        'gamma_r2': gamma_r2,
    }


def main():
    print("=" * 78)
    print("STEP 5.71: MATCHING LEAKAGE AUDIT")
    print("=" * 78)
    print()

    gc_rows, field_rows = load_data()
    print(f"Loaded {len(gc_rows)} GC pulsars, {len(field_rows)} field pulsars")
    print()

    variants = [
        ('Raw (no matching)', 'raw', None, 5),
        ('Period-only matching (5 NN)', 'matched', ['logP'], 5),
        ('Period + log(tau_c) matching (5 NN)', 'matched', ['logP', 'log_tau_c'], 5),
        ('Period + B-proxy matching (5 NN)', 'matched', ['logP', 'log_b_proxy'], 5),
        ('Period + B-proxy matching (15 NN)', 'matched', ['logP', 'log_b_proxy'], 15),
        ('Period + B-proxy cluster-residualized (5 NN)', 'cluster_residualized', ['logP', 'log_b_proxy'], 5),
    ]

    results = {
        'timestamp_utc': datetime.now(timezone.utc).isoformat(),
        'variants': [],
    }

    print(f"{'Variant':<45} {'Residual':>12} {'Gamma':>10} {'Gamma err':>10}")
    print("-" * 78)

    for name, mode, keys, n in variants:
        variant = run_variant(name, gc_rows, field_rows, mode, keys, n)
        results['variants'].append(variant)
        gamma_str = f"{variant['gamma']:.3f}" if variant['gamma'] is not None else "N/A"
        gamma_err_str = f"+/- {variant['gamma_err']:.3f}" if variant['gamma_err'] is not None else "N/A"
        print(f"{name:<45} {variant['residual_dex']:>+11.3f} dex {gamma_str:>10} {gamma_err_str:>10}")

    # Compute attenuation relative to raw
    raw_resid = results['variants'][0]['residual_dex']
    raw_gamma = results['variants'][0]['gamma']

    print()
    print("Attenuation relative to raw:")
    print(f"{'Variant':<45} {'Resid atten':>12} {'Gamma atten':>12}")
    print("-" * 78)
    for v in results['variants']:
        resid_atten = v['residual_dex'] / raw_resid if raw_resid != 0 else None
        gamma_atten = (v['gamma'] / raw_gamma) if (raw_gamma and v['gamma'] is not None) else None
        print(f"{v['name']:<45} {resid_atten:>11.2%} {gamma_atten:>11.2%}" if resid_atten is not None else f"{v['name']:<45} {'N/A':>12} {gamma_atten:>11.2%}")

    # Save
    with open(OUT_JSON, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\nSaved JSON: {OUT_JSON}")

    md = f"""# Step 51: Matching Leakage Audit

**Generated:** {results['timestamp_utc']}

## Results

| Variant | Residual (dex) | Gamma | Gamma err |
|---------|----------------|-------|-----------|
"""
    for v in results['variants']:
        gamma_str = f"{v['gamma']:.3f}" if v['gamma'] is not None else "N/A"
        gamma_err_str = f"+/- {v['gamma_err']:.3f}" if v['gamma_err'] is not None else "N/A"
        md += f"| {v['name']} | {v['residual_dex']:+.3f} | {gamma_str} | {gamma_err_str} |\n"

    md += """
## Interpretation

- **Raw** gives the unattenuated signal.
- **Period-only** should recover nearly the full amplitude if matching is clean.
- **Period + B-proxy** will attenuate if B-proxy leaks outcome information.
- **Expanded NN (15)** may over-smooth and further suppress amplitude.
- **Cluster-residualized** subtracts the field mean per cluster; if the field mean
carries density-correlated structure, this will attenuate the slope.
"""

    with open(OUT_MD, 'w') as f:
        f.write(md)
    print(f"Saved Markdown: {OUT_MD}")

    print("\nDone.")


if __name__ == '__main__':
    main()
