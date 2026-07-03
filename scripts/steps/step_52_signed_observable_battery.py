#!/usr/bin/env python3
"""
Step 52: Signed Observable Battery

The analysis uses log|Ṗ| as the main observable, which can compress or distort
a signed acceleration signal.  This script runs the full diagnostic battery on:

  1. log|Ṗ|                (current headline)
  2. log|Ṗ/P|              (direct acceleration observable)
  3. signed Ṗ/P            (preserves acceleration direction)
  4. positive-only Ṗ       (intrinsic spin-down branch)
  5. negative-only Ṗ       (acceleration-dominated branch)
  6. distribution width / MAD (broadening instead of mean shift)

If TEP is a gradient/timing-response effect, the strongest signal may live in
signed Ṗ/P and distributional broadening rather than mean log|Ṗ|.

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

OUT_JSON = RESULTS_DIR / "step_52_signed_observable_battery.json"
OUT_MD = RESULTS_DIR / "step_52_signed_observable_battery.md"

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

C = 299792458.0


def load_data():
    gc_rows = []
    field_rows = []
    with open(STEP_5_10_CSV, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                p0 = float(row['P0_s'])
                p1 = float(row['P1_sps'])
                r = {
                    'source': row['source'],
                    'environment': row['environment'],
                    'cluster': row['cluster'],
                    'name': row['name'],
                    'P0_s': p0,
                    'P1_sps': p1,
                    'logP': float(row['logP']),
                    'logPdot_abs': float(row['logPdot_abs']),
                    'log_b_proxy': float(row['log_b_proxy']),
                    'log_tau_c': float(row['log_tau_c']),
                    'pdot_over_p': p1 / p0 if p0 > 0 else 0.0,
                    'log_abs_pdot_over_p': math.log10(abs(p1 / p0)) if p0 > 0 and p1 != 0 else -99.0,
                    'pdot_sign': 1 if p1 > 0 else (-1 if p1 < 0 else 0),
                }
                if r['environment'] == 'globular_cluster':
                    gc_rows.append(r)
                elif r['environment'] == 'field':
                    field_rows.append(r)
            except (ValueError, KeyError):
                continue
    return gc_rows, field_rows


def compute_observable(rows, mode):
    """
    Return a list of numeric values for the given observable mode.
    """
    out = []
    for r in rows:
        if mode in ('log_abs_pdot', 'logPdot_abs'):
            out.append(r['logPdot_abs'])
        elif mode == 'log_abs_pdot_over_p':
            out.append(r['log_abs_pdot_over_p'])
        elif mode == 'signed_pdot_over_p':
            out.append(r['pdot_over_p'])
        elif mode == 'positive_only_pdot':
            if r['P1_sps'] > 0:
                out.append(r['logPdot_abs'])
        elif mode == 'negative_only_pdot':
            if r['P1_sps'] < 0:
                out.append(r['logPdot_abs'])
        elif mode == 'width_mad':
            out.append(r['logPdot_abs'])
    return np.array(out) if out else np.array([])


def raw_summary(gc_rows, field_rows, mode):
    gc_vals = compute_observable(gc_rows, mode)
    field_vals = compute_observable(field_rows, mode)

    if len(gc_vals) == 0 or len(field_vals) == 0:
        return {}

    if mode == 'width_mad':
        # For width, compare MAD (median absolute deviation) and std
        gc_mad = float(np.median(np.abs(gc_vals - np.median(gc_vals))))
        field_mad = float(np.median(np.abs(field_vals - np.median(field_vals))))
        gc_std = float(np.std(gc_vals, ddof=1))
        field_std = float(np.std(field_vals, ddof=1))
        return {
            'gc_mad': gc_mad,
            'field_mad': field_mad,
            'delta_mad': gc_mad - field_mad,
            'gc_std': gc_std,
            'field_std': field_std,
            'delta_std': gc_std - field_std,
            'n_gc': len(gc_vals),
            'n_field': len(field_vals),
        }

    if mode == 'signed_pdot_over_p':
        # For signed, use mean and median of signed values
        gc_mean = float(np.mean(gc_vals))
        field_mean = float(np.mean(field_vals))
        gc_median = float(np.median(gc_vals))
        field_median = float(np.median(field_vals))
        # Also report negative fraction
        gc_neg = float(np.sum(gc_vals < 0) / len(gc_vals))
        field_neg = float(np.sum(field_vals < 0) / len(field_vals))
        return {
            'gc_mean': gc_mean,
            'field_mean': field_mean,
            'delta_mean': gc_mean - field_mean,
            'gc_median': gc_median,
            'field_median': field_median,
            'delta_median': gc_median - field_median,
            'gc_negative_fraction': gc_neg,
            'field_negative_fraction': field_neg,
            'delta_negative_fraction': gc_neg - field_neg,
            'n_gc': len(gc_vals),
            'n_field': len(field_vals),
        }

    # Default: mean difference for log observables
    gc_mean = float(np.mean(gc_vals))
    field_mean = float(np.mean(field_vals))
    return {
        'gc_mean': gc_mean,
        'field_mean': field_mean,
        'delta_mean': gc_mean - field_mean,
        'n_gc': len(gc_vals),
        'n_field': len(field_vals),
    }


def matched_residual(gc_rows, field_rows, match_keys, n_matches=5, observable='logPdot_abs'):
    if not gc_rows or not field_rows:
        return 0.0

    def get_val(r):
        if observable == 'logPdot_abs':
            return r['logPdot_abs']
        elif observable == 'log_abs_pdot_over_p':
            return r['log_abs_pdot_over_p']
        elif observable == 'signed_pdot_over_p':
            return r['pdot_over_p']
        else:
            return r['logPdot_abs']

    gc_vals = np.array([get_val(r) for r in gc_rows])
    field_vals = np.array([get_val(r) for r in field_rows])

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


def density_scaling_slope(gc_rows, field_rows, mode, observable='logPdot_abs', match_keys=None):
    clusters = defaultdict(list)
    for r in gc_rows:
        clusters[r['cluster']].append(r)

    def get_val(r):
        if observable == 'logPdot_abs':
            return r['logPdot_abs']
        elif observable == 'log_abs_pdot_over_p':
            return r['log_abs_pdot_over_p']
        elif observable == 'signed_pdot_over_p':
            return r['pdot_over_p']
        else:
            return r['logPdot_abs']

    cluster_means = []
    densities = []
    for cluster_name, gc_list in clusters.items():
        if cluster_name not in CLUSTER_PARAMS:
            continue
        rho_c = CLUSTER_PARAMS[cluster_name]['rho_c']

        if mode == 'raw':
            vals = np.array([get_val(r) for r in gc_list])
            cmean = np.mean(vals)
        elif mode == 'matched':
            cmean = matched_residual(gc_list, field_rows, match_keys, n_matches=5, observable=observable)
        else:
            raise ValueError(f"Unknown mode: {mode}")

        cluster_means.append(cmean)
        densities.append(rho_c)

    if len(densities) < 3:
        return None, None, None

    slope, intercept, r_value, p_value, std_err = stats.linregress(densities, cluster_means)
    return float(slope), float(std_err), float(r_value ** 2)


def run_battery(name, gc_rows, field_rows, mode, match_keys=None, observable='logPdot_abs'):
    result = {'name': name, 'mode': mode, 'observable': observable, 'match_keys': match_keys}

    if observable in ('log_abs_pdot_over_p', 'logPdot_abs'):
        raw = raw_summary(gc_rows, field_rows, observable.replace('log_abs_pdot_over_p', 'log_abs_pdot_over_p').replace('logPdot_abs', 'log_abs_pdot'))
    elif observable == 'signed_pdot_over_p':
        raw = raw_summary(gc_rows, field_rows, 'signed_pdot_over_p')
    else:
        raw = {}

    result['raw'] = raw

    if mode == 'matched' and match_keys:
        resid = matched_residual(gc_rows, field_rows, match_keys, observable=observable)
        gamma, gamma_err, gamma_r2 = density_scaling_slope(gc_rows, field_rows, 'matched', observable, match_keys)
    else:
        resid = raw.get('delta_mean', 0.0) if raw else 0.0
        gamma, gamma_err, gamma_r2 = density_scaling_slope(gc_rows, field_rows, 'raw', observable)

    result['residual_dex'] = resid
    result['gamma'] = gamma
    result['gamma_err'] = gamma_err
    result['gamma_r2'] = gamma_r2

    return result


def main():
    print("=" * 78)
    print("STEP 5.72: SIGNED OBSERVABLE BATTERY")
    print("=" * 78)
    print()

    gc_rows, field_rows = load_data()
    print(f"Loaded {len(gc_rows)} GC pulsars, {len(field_rows)} field pulsars")
    print()

    # Raw comparisons for all observables
    observables = [
        ('log|Ṗ|', 'logPdot_abs'),
        ('log|Ṗ/P|', 'log_abs_pdot_over_p'),
        ('signed Ṗ/P', 'signed_pdot_over_p'),
        ('positive-only Ṗ', 'positive_only_pdot'),
        ('negative-only Ṗ', 'negative_only_pdot'),
        ('width / MAD', 'width_mad'),
    ]

    print("Raw comparisons:")
    print(f"{'Observable':<25} {'GC mean':>12} {'Field mean':>12} {'Delta':>12} {'N_GC':>6} {'N_Field':>8}")
    print("-" * 78)

    raw_results = []
    for label, mode in observables:
        raw = raw_summary(gc_rows, field_rows, mode)
        raw_results.append((label, mode, raw))
        if not raw:
            print(f"{label:<25} {'N/A':>12} {'N/A':>12} {'N/A':>12} {'N/A':>6} {'N/A':>8}")
        elif mode == 'width_mad':
            print(f"{label:<25} {raw.get('gc_mad', 0):>12.3f} {raw.get('field_mad', 0):>12.3f} {raw.get('delta_mad', 0):>+11.3f} {raw.get('n_gc', 0):>6} {raw.get('n_field', 0):>8}")
        elif mode == 'signed_pdot_over_p':
            print(f"{label:<25} {raw.get('gc_mean', 0):>12.3e} {raw.get('field_mean', 0):>12.3e} {raw.get('delta_mean', 0):>+11.3e} {raw.get('n_gc', 0):>6} {raw.get('n_field', 0):>8}")
        else:
            print(f"{label:<25} {raw.get('gc_mean', 0):>12.3f} {raw.get('field_mean', 0):>12.3f} {raw.get('delta_mean', 0):>+11.3f} {raw.get('n_gc', 0):>6} {raw.get('n_field', 0):>8}")

    print()
    print("Matched residuals by observable:")
    print(f"{'Observable / Matching':<45} {'Residual':>12} {'Gamma':>10} {'Gamma err':>10}")
    print("-" * 78)

    # Matched residuals for key observables
    battery_results = []
    for label, observable in [
        ('log|Ṗ|  raw', 'logPdot_abs'),
        ('log|Ṗ|  period-only', 'logPdot_abs'),
        ('log|Ṗ|  period+B-proxy', 'logPdot_abs'),
        ('log|Ṗ/P| raw', 'log_abs_pdot_over_p'),
        ('log|Ṗ/P| period-only', 'log_abs_pdot_over_p'),
        ('log|Ṗ/P| period+B-proxy', 'log_abs_pdot_over_p'),
        ('signed Ṗ/P raw', 'signed_pdot_over_p'),
        ('signed Ṗ/P period-only', 'signed_pdot_over_p'),
        ('signed Ṗ/P period+B-proxy', 'signed_pdot_over_p'),
    ]:
        if 'raw' in label:
            mode = 'raw'
            keys = None
        elif 'period-only' in label:
            mode = 'matched'
            keys = ['logP']
        else:
            mode = 'matched'
            keys = ['logP', 'log_b_proxy']

        result = run_battery(label, gc_rows, field_rows, mode, keys, observable)
        battery_results.append(result)
        gamma_str = f"{result['gamma']:.3f}" if result['gamma'] is not None else "N/A"
        gamma_err_str = f"+/- {result['gamma_err']:.3f}" if result['gamma_err'] is not None else "N/A"
        print(f"{label:<45} {result['residual_dex']:>+11.3f} {gamma_str:>10} {gamma_err_str:>10}")

    # Save
    results = {
        'timestamp_utc': datetime.now(timezone.utc).isoformat(),
        'raw_comparisons': {label: raw for label, mode, raw in raw_results},
        'battery': battery_results,
    }

    with open(OUT_JSON, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\nSaved JSON: {OUT_JSON}")

    md = f"""# Step 52: Signed Observable Battery

**Generated:** {results['timestamp_utc']}

## Raw Comparisons

| Observable | GC | Field | Delta | N_GC | N_Field |
|------------|----|-------|-------|------|---------|
"""
    for label, mode, raw in raw_results:
        if mode == 'width_mad':
            md += f"| {label} | {raw.get('gc_mad', 0):.3f} | {raw.get('field_mad', 0):.3f} | {raw.get('delta_mad', 0):+.3f} | {raw.get('n_gc', 0)} | {raw.get('n_field', 0)} |\n"
        elif mode == 'signed_pdot_over_p':
            md += f"| {label} | {raw.get('gc_mean', 0):.3e} | {raw.get('field_mean', 0):.3e} | {raw.get('delta_mean', 0):+.3e} | {raw.get('n_gc', 0)} | {raw.get('n_field', 0)} |\n"
        else:
            md += f"| {label} | {raw.get('gc_mean', 0):.3f} | {raw.get('field_mean', 0):.3f} | {raw.get('delta_mean', 0):+.3f} | {raw.get('n_gc', 0)} | {raw.get('n_field', 0)} |\n"

    md += """
## Matched Battery

| Variant | Residual | Gamma | Gamma err |
|---------|----------|-------|-----------|
"""
    for b in battery_results:
        gamma_str = f"{b['gamma']:.3f}" if b['gamma'] is not None else "N/A"
        gamma_err_str = f"+/- {b['gamma_err']:.3f}" if b['gamma_err'] is not None else "N/A"
        md += f"| {b['name']} | {b['residual_dex']:+.3f} | {gamma_str} | {gamma_err_str} |\n"

    md += """
## Interpretation

- **log|Ṗ/P|** is the direct acceleration observable; if TEP acts on Ṗ/P, this
  should show the cleanest signal.
- **signed Ṗ/P** preserves direction; an excess of negative values in GCs
  indicates acceleration-dominated line-of-sight components.
- **positive-only / negative-only** isolate the intrinsic and acceleration
  branches respectively.
- **width / MAD** tests whether the effect manifests as broadening rather than
  a mean shift (expected if TEP is a stochastic or multi-component process).
"""

    with open(OUT_MD, 'w') as f:
        f.write(md)
    print(f"Saved Markdown: {OUT_MD}")

    print("\nDone.")


if __name__ == '__main__':
    main()
