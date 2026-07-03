#!/usr/bin/env python3
"""
Step 50: Injection-Recovery Test

Inject known TEP-like signals into the observed catalog, then run the full
analysis pipeline unchanged.  This is the single most decisive bug test: if the
pipeline cannot recover a known injected effect at full strength, the current
"weak" results are almost certainly methodological attenuation.

Injected signals:
  1. +0.60 dex uniform GC excess
  2. Gamma = 0.39 density-flattened signal
  3. Gamma = 0.75 Newtonian scaling
  4. Binary shielding -0.33 dex
  5. Signed acceleration mixture (known positive/negative fractions)

Recovery metrics:
  - Raw mean amplitude (dex)
  - Period-only matched residual (dex)
  - Period+B-field matched residual (dex)
  - Density-scaling slope Gamma
  - Binary sign fraction (GC negative %)

Author: M. Smawfield
Date: 2026-06-23
"""

import csv
import json
import math
import re
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from scipy import stats

REPO_ROOT = Path(__file__).resolve().parents[2]
RESULTS_DIR = REPO_ROOT / "results" / "outputs"
DATA_DIR = REPO_ROOT / "data"
FREIRE_TXT = DATA_DIR / "freire_GCpsr.txt"

OUT_JSON = RESULTS_DIR / "step_50_injection_recovery.json"
OUT_MD = RESULTS_DIR / "step_50_injection_recovery.md"

STEP_5_10_CSV = RESULTS_DIR / "step_02_pulsar_population_controls.csv"

# ---------------------------------------------------------------------------
# Cluster parameters (from step_11 and step_07)
# ---------------------------------------------------------------------------
CLUSTER_PARAMS = {
    "Terzan 5":         {"M": 2.0e6, "Rc": 0.16, "Rt": 5.0,  "rho_c": 5.50},
    "47 Tuc (NGC 104)": {"M": 1.0e6, "Rc": 0.36, "Rt": 42.0, "rho_c": 4.88},
    "NGC 6517":         {"M": 2.0e5, "Rc": 0.06, "Rt": 5.0,  "rho_c": 5.80},
    "M28 (NGC 6626)":   {"M": 5.0e5, "Rc": 0.24, "Rt": 12.0, "rho_c": 4.52},
    "M62 (NGC 6266)":   {"M": 1.0e6, "Rc": 0.18, "Rt": 8.0,  "rho_c": 5.16},
    "M13 (NGC 6205)":   {"M": 6.0e5, "Rc": 0.62, "Rt": 25.0, "rho_c": 3.79},
    "M15 (NGC 7078)":   {"M": 5.0e5, "Rc": 0.14, "Rt": 21.0, "rho_c": 5.05},
    "M5 (NGC 5904)":    {"M": 5.0e5, "Rc": 0.42, "Rt": 28.0, "rho_c": 3.53},
    "Terzan 1":         {"M": 1.5e5, "Rc": 0.10, "Rt": 4.0,  "rho_c": 5.00},
    "NGC 6752":         {"M": 3.0e5, "Rc": 0.17, "Rt": 25.0, "rho_c": 4.30},
    "M2 (NGC 7089)":    {"M": 6.0e5, "Rc": 0.32, "Rt": 21.0, "rho_c": 4.15},
    "Omega Centauri (NGC 5139)": {"M": 4.0e6, "Rc": 2.37, "Rt": 57.0, "rho_c": 3.12},
    "M53 (NGC 5024)":   {"M": 3.0e5, "Rc": 0.65, "Rt": 22.0, "rho_c": 2.96},
    "M3 (NGC 5272)":    {"M": 5.0e5, "Rc": 0.37, "Rt": 38.0, "rho_c": 3.68},
    "M71 (NGC 6838)":   {"M": 2.0e4, "Rc": 0.63, "Rt": 8.0,  "rho_c": 2.29},
    "NGC 6397":         {"M": 1.0e5, "Rc": 0.05, "Rt": 15.0, "rho_c": 5.68},
    "NGC 1851":         {"M": 3.0e5, "Rc": 0.09, "Rt": 11.0, "rho_c": 5.09},
    "NGC 6522":         {"M": 2.0e5, "Rc": 0.05, "Rt": 5.0,  "rho_c": 5.50},
    "NGC 6544":         {"M": 5.0e4, "Rc": 0.05, "Rt": 3.0,  "rho_c": 5.20},
    "NGC 6624":         {"M": 2.0e5, "Rc": 0.06, "Rt": 6.0,  "rho_c": 5.60},
    "NGC 6760":         {"M": 2.0e5, "Rc": 0.34, "Rt": 8.0,  "rho_c": 3.80},
    "M22 (NGC 6656)":   {"M": 5.0e5, "Rc": 1.33, "Rt": 32.0, "rho_c": 2.97},
    "M80 (NGC 6093)":   {"M": 4.0e5, "Rc": 0.15, "Rt": 13.0, "rho_c": 4.79},
    "M92 (NGC 6341)":   {"M": 3.0e5, "Rc": 0.26, "Rt": 15.0, "rho_c": 4.30},
    "NGC 6712":         {"M": 1.5e5, "Rc": 0.33, "Rt": 7.0,  "rho_c": 3.70},
    "NGC 6652":         {"M": 1.0e5, "Rc": 0.10, "Rt": 5.0,  "rho_c": 4.50},
    "M14 (NGC 6402)":   {"M": 1.0e6, "Rc": 0.78, "Rt": 18.0, "rho_c": 3.44},
    "NGC 6539":         {"M": 3.0e5, "Rc": 0.60, "Rt": 10.0, "rho_c": 3.30},
    "M4 (NGC 6121)":    {"M": 1.0e5, "Rc": 0.83, "Rt": 33.0, "rho_c": 2.85},
    "NGC 6440":         {"M": 1.5e5, "Rc": 0.12, "Rt": 5.5,  "rho_c": 5.10},
    "NGC 6441":         {"M": 8.0e5, "Rc": 0.20, "Rt": 12.0, "rho_c": 5.00},
    "NGC 6316":         {"M": 1.2e5, "Rc": 0.15, "Rt": 5.0,  "rho_c": 4.80},
    "M30 (NGC 7099)":   {"M": 2.5e5, "Rc": 0.25, "Rt": 18.0, "rho_c": 4.20},
}

_NUM_RE = re.compile(r'^[+-]?\d+\.?\d*(?:[eE][+-]?\d+)?$')


def _parse_numeric(val: str):
    if not val or val in ('*', 'i'):
        return None
    val = re.sub(r'\([^)]*\)', '', val).strip().lstrip('<>')
    if _NUM_RE.match(val):
        return float(val)
    return None


def parse_freire_catalog(text: str) -> list[dict]:
    rows = []
    current_cluster = None
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        if not line.startswith('J') and not line.startswith('B') and not line.startswith('*'):
            if '(' in line or 'NGC' in line or 'Terzan' in line or 'Omega' in line:
                current_cluster = line.strip()
                continue
            parts = line.split()
            if len(parts) <= 3 and not any(c.isdigit() and '.' in line for c in line):
                current_cluster = line.strip()
                continue
        if current_cluster is None:
            continue
        parts = line.split()
        if len(parts) < 5:
            continue
        name = parts[0]
        if not (name.startswith('J') or name.startswith('B')):
            continue
        try:
            idx = 1
            _parse_numeric(parts[idx]) if idx < len(parts) else None
            idx += 1
            _parse_numeric(parts[idx]) if idx < len(parts) else None
            idx += 1
            _parse_numeric(parts[idx]) if idx < len(parts) else None
            idx += 1
            _parse_numeric(parts[idx]) if idx < len(parts) else None
            idx += 1
            Pb_days = _parse_numeric(parts[idx]) if idx < len(parts) else None
            idx += 1
            _parse_numeric(parts[idx]) if idx < len(parts) else None
            idx += 1
            _parse_numeric(parts[idx]) if idx < len(parts) else None
            idx += 1
            _parse_numeric(parts[idx]) if idx < len(parts) else None
            is_binary = Pb_days is not None and Pb_days > 0
            rows.append({
                'cluster': current_cluster,
                'name': name,
                'is_binary': is_binary,
            })
        except (IndexError, ValueError):
            continue
    return rows


def load_step_02():
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
                }
                if r['environment'] == 'globular_cluster':
                    gc_rows.append(r)
                elif r['environment'] == 'field':
                    field_rows.append(r)
            except (ValueError, KeyError):
                continue
    return gc_rows, field_rows


def get_binary_flags(gc_rows):
    if not FREIRE_TXT.exists():
        return {}
    with open(FREIRE_TXT, 'r', errors='ignore') as f:
        freire_rows = parse_freire_catalog(f.read())
    binary_map = {}
    for r in freire_rows:
        binary_map[r['name']] = r['is_binary']
    out = {}
    for row in gc_rows:
        name = row['name']
        if name in binary_map:
            out[name] = binary_map[name]
            continue
        base = re.sub(r'[A-Za-z]+$', '', name)
        for fname, flag in binary_map.items():
            if fname.startswith(base) and len(fname) == len(name):
                out[name] = flag
                break
    return out


# ---------------------------------------------------------------------------
# Core analysis functions (replicated from existing steps)
# ---------------------------------------------------------------------------

def raw_mean_difference(gc_rows, field_rows):
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


def density_scaling_slope(gc_rows, field_rows, residual_mode='raw'):
    clusters = defaultdict(list)
    for r in gc_rows:
        clusters[r['cluster']].append(r)

    field_arr = field_rows

    cluster_means = []
    densities = []
    for cluster_name, gc_list in clusters.items():
        if cluster_name not in CLUSTER_PARAMS:
            continue
        rho_c = CLUSTER_PARAMS[cluster_name]['rho_c']
        if residual_mode == 'raw':
            vals = np.array([r['logPdot_abs'] for r in gc_list])
            cmean = np.mean(vals)
        elif residual_mode == 'period_only':
            cmean = matched_residual(gc_list, field_arr, ['logP'], n_matches=5)
        elif residual_mode == 'period_bfield':
            cmean = matched_residual(gc_list, field_arr, ['logP', 'log_b_proxy'], n_matches=5)
        else:
            raise ValueError(f"Unknown residual_mode: {residual_mode}")
        cluster_means.append(cmean)
        densities.append(rho_c)

    if len(densities) < 3:
        return None, None, None

    slope, intercept, r_value, p_value, std_err = stats.linregress(densities, cluster_means)
    return float(slope), float(std_err), float(r_value ** 2)


def binary_sign_fraction(gc_rows, binary_flags):
    signed = [r for r in gc_rows if r['P1_sps'] != 0]
    if not signed:
        return 0.0, 0, 0
    neg = sum(1 for r in signed if r['P1_sps'] < 0)
    n = len(signed)
    return float(neg / n), neg, n


def binary_fraction_by_type(gc_rows, binary_flags):
    signed = [r for r in gc_rows if r['P1_sps'] != 0]
    binary_neg = 0
    binary_n = 0
    iso_neg = 0
    iso_n = 0
    for r in signed:
        is_bin = binary_flags.get(r['name'], False)
        if is_bin:
            binary_n += 1
            if r['P1_sps'] < 0:
                binary_neg += 1
        else:
            iso_n += 1
            if r['P1_sps'] < 0:
                iso_neg += 1
    return {
        'binary_negative_fraction': float(binary_neg / binary_n) if binary_n else 0.0,
        'binary_n': binary_n,
        'isolated_negative_fraction': float(iso_neg / iso_n) if iso_n else 0.0,
        'isolated_n': iso_n,
    }


# ---------------------------------------------------------------------------
# Injection generators
# ---------------------------------------------------------------------------

def inject_null(gc_rows, field_rows):
    return [dict(r) for r in gc_rows], [dict(r) for r in field_rows]


def inject_uniform_excess(gc_rows, field_rows, dex=0.60):
    gc_out = []
    for r in gc_rows:
        r2 = dict(r)
        sign = 1.0 if r2['P1_sps'] >= 0 else -1.0
        abs_p1 = abs(r2['P1_sps'])
        new_abs = abs_p1 * (10 ** dex)
        r2['P1_sps'] = sign * new_abs
        r2['logPdot_abs'] = math.log10(new_abs)
        r2['log_b_proxy'] = math.log10(math.sqrt(r2['P0_s'] * new_abs))
        gc_out.append(r2)
    return gc_out, [dict(r) for r in field_rows]


def inject_density_scaling(gc_rows, field_rows, gamma_target=0.39, baseline_dex=0.0):
    gc_out = []
    for r in gc_rows:
        r2 = dict(r)
        cluster = r2['cluster']
        rho_c = CLUSTER_PARAMS.get(cluster, {}).get('rho_c', 4.0)
        extra_dex = baseline_dex + gamma_target * rho_c
        sign = 1.0 if r2['P1_sps'] >= 0 else -1.0
        abs_p1 = abs(r2['P1_sps'])
        new_abs = abs_p1 * (10 ** extra_dex)
        r2['P1_sps'] = sign * new_abs
        r2['logPdot_abs'] = math.log10(new_abs)
        r2['log_b_proxy'] = math.log10(math.sqrt(r2['P0_s'] * new_abs))
        gc_out.append(r2)
    return gc_out, [dict(r) for r in field_rows]


def inject_binary_shielding(gc_rows, field_rows, binary_flags, dex=-0.33):
    gc_out = []
    for r in gc_rows:
        r2 = dict(r)
        is_bin = binary_flags.get(r2['name'], False)
        if is_bin:
            sign = 1.0 if r2['P1_sps'] >= 0 else -1.0
            abs_p1 = abs(r2['P1_sps'])
            new_abs = abs_p1 * (10 ** dex)
            r2['P1_sps'] = sign * new_abs
            r2['logPdot_abs'] = math.log10(new_abs)
            r2['log_b_proxy'] = math.log10(math.sqrt(r2['P0_s'] * new_abs))
        gc_out.append(r2)
    return gc_out, [dict(r) for r in field_rows]


def inject_signed_mixture(gc_rows, field_rows, target_negative_fraction=0.40, amplification_dex=1.0):
    rng = np.random.default_rng(42)
    n_gc = len(gc_rows)
    n_neg = int(round(n_gc * target_negative_fraction))
    neg_indices = set(rng.choice(n_gc, size=n_neg, replace=False))

    gc_out = []
    for i, r in enumerate(gc_rows):
        r2 = dict(r)
        abs_p1 = abs(r2['P1_sps'])
        new_abs = abs_p1 * (10 ** amplification_dex)
        if i in neg_indices:
            r2['P1_sps'] = -new_abs
        else:
            r2['P1_sps'] = +new_abs
        r2['logPdot_abs'] = math.log10(new_abs)
        r2['log_b_proxy'] = math.log10(math.sqrt(r2['P0_s'] * new_abs))
        gc_out.append(r2)
    return gc_out, [dict(r) for r in field_rows]


# ---------------------------------------------------------------------------
# Recovery metrics
# ---------------------------------------------------------------------------

def compute_recovery(injected, baseline, injected_params):
    recovery = {}
    delta_injected = injected['raw_amplitude_dex'] - baseline['raw_amplitude_dex']
    target_raw = injected_params.get('raw_dex', 0.0)
    if target_raw != 0:
        recovery['raw_amplitude_recovery_fraction'] = float(delta_injected / target_raw)
    recovery['raw_amplitude_delta'] = float(delta_injected)

    delta_po = injected['period_only_residual_dex'] - baseline['period_only_residual_dex']
    target_po = injected_params.get('period_only_dex', target_raw)
    if target_po != 0:
        recovery['period_only_recovery_fraction'] = float(delta_po / target_po)
    recovery['period_only_delta'] = float(delta_po)

    delta_pb = injected['period_bfield_residual_dex'] - baseline['period_bfield_residual_dex']
    target_pb = injected_params.get('period_bfield_dex', target_raw)
    if target_pb != 0:
        recovery['period_bfield_recovery_fraction'] = float(delta_pb / target_pb)
    recovery['period_bfield_delta'] = float(delta_pb)

    target_gamma = injected_params.get('gamma', None)
    if target_gamma is not None and baseline['gamma_raw'] is not None and injected['gamma_raw'] is not None:
        delta_gamma = injected['gamma_raw'] - baseline['gamma_raw']
        recovery['gamma_raw_delta'] = float(delta_gamma)
        recovery['gamma_raw_injected'] = float(target_gamma)
        recovery['gamma_raw_recovered'] = float(injected['gamma_raw'])
        if abs(target_gamma) > 0.01:
            recovery['gamma_raw_recovery_fraction'] = float(delta_gamma / target_gamma)

    if target_gamma is not None and baseline['gamma_period_only'] is not None and injected['gamma_period_only'] is not None:
        delta_gamma_po = injected['gamma_period_only'] - baseline['gamma_period_only']
        recovery['gamma_period_only_delta'] = float(delta_gamma_po)
        if abs(target_gamma) > 0.01:
            recovery['gamma_period_only_recovery_fraction'] = float(delta_gamma_po / target_gamma)

    target_neg = injected_params.get('negative_fraction', None)
    if target_neg is not None:
        delta_neg = injected['negative_fraction'] - baseline['negative_fraction']
        recovery['negative_fraction_delta'] = float(delta_neg)
        recovery['negative_fraction_target'] = float(target_neg)
        recovery['negative_fraction_recovered'] = float(injected['negative_fraction'])
        recovery['negative_fraction_match'] = float(injected['negative_fraction'])

    return recovery


def run_analysis_battery(gc_rows, field_rows, binary_flags, label=""):
    raw_diff = raw_mean_difference(gc_rows, field_rows)
    period_only_resid = matched_residual(gc_rows, field_rows, ['logP'], n_matches=5)
    period_bfield_resid = matched_residual(gc_rows, field_rows, ['logP', 'log_b_proxy'], n_matches=5)

    gamma_raw, gamma_raw_err, gamma_raw_r2 = density_scaling_slope(gc_rows, field_rows, 'raw')
    gamma_po, gamma_po_err, gamma_po_r2 = density_scaling_slope(gc_rows, field_rows, 'period_only')
    gamma_pb, gamma_pb_err, gamma_pb_r2 = density_scaling_slope(gc_rows, field_rows, 'period_bfield')

    neg_frac, neg_n, total_n = binary_sign_fraction(gc_rows, binary_flags)
    bin_type = binary_fraction_by_type(gc_rows, binary_flags)

    return {
        'label': label,
        'n_gc': len(gc_rows),
        'n_field': len(field_rows),
        'raw_amplitude_dex': raw_diff,
        'period_only_residual_dex': period_only_resid,
        'period_bfield_residual_dex': period_bfield_resid,
        'gamma_raw': gamma_raw,
        'gamma_raw_err': gamma_raw_err,
        'gamma_raw_r2': gamma_raw_r2,
        'gamma_period_only': gamma_po,
        'gamma_period_only_err': gamma_po_err,
        'gamma_period_only_r2': gamma_po_r2,
        'gamma_period_bfield': gamma_pb,
        'gamma_period_bfield_err': gamma_pb_err,
        'gamma_period_bfield_r2': gamma_pb_r2,
        'negative_fraction': neg_frac,
        'negative_n': neg_n,
        'total_n': total_n,
        'binary_negative_fraction': bin_type['binary_negative_fraction'],
        'binary_n': bin_type['binary_n'],
        'isolated_negative_fraction': bin_type['isolated_negative_fraction'],
        'isolated_n': bin_type['isolated_n'],
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print("=" * 78)
    print("STEP 5.70: INJECTION-RECOVERY TEST")
    print("=" * 78)
    print()

    gc_rows_orig, field_rows_orig = load_step_02()
    binary_flags = get_binary_flags(gc_rows_orig)

    print(f"Loaded {len(gc_rows_orig)} GC pulsars, {len(field_rows_orig)} field pulsars")
    print(f"Binary flags resolved for {len(binary_flags)} GC pulsars")
    print()

    # Baseline (no injection)
    print("--- BASELINE (no injection) ---")
    baseline = run_analysis_battery(gc_rows_orig, field_rows_orig, binary_flags, label="baseline")
    print(f"  Raw amplitude:         {baseline['raw_amplitude_dex']:+.3f} dex")
    print(f"  Period-only residual:  {baseline['period_only_residual_dex']:+.3f} dex")
    print(f"  Period+B residual:     {baseline['period_bfield_residual_dex']:+.3f} dex")
    if baseline['gamma_raw'] is not None:
        print(f"  Gamma (raw):           {baseline['gamma_raw']:.3f} +/- {baseline['gamma_raw_err']:.3f}")
    print(f"  Negative fraction:     {baseline['negative_fraction']:.1%} ({baseline['negative_n']}/{baseline['total_n']})")
    print()

    scenarios = [
        {
            'name': 'Uniform +0.60 dex GC excess',
            'injector': lambda gc, field: inject_uniform_excess(gc, field, dex=0.60),
            'params': {'raw_dex': 0.60, 'period_only_dex': 0.60, 'period_bfield_dex': 0.60},
        },
        {
            'name': 'Gamma = 0.39 density-flattened',
            'injector': lambda gc, field: inject_density_scaling(gc, field, gamma_target=0.39),
            'params': {'gamma': 0.39},
        },
        {
            'name': 'Gamma = 0.75 Newtonian scaling',
            'injector': lambda gc, field: inject_density_scaling(gc, field, gamma_target=0.75),
            'params': {'gamma': 0.75},
        },
        {
            'name': 'Binary shielding -0.33 dex',
            'injector': lambda gc, field: inject_binary_shielding(gc, field, binary_flags, dex=-0.33),
            'params': {'raw_dex': -0.33 * (baseline['binary_n'] / baseline['total_n']) if baseline['total_n'] else 0.0},
        },
        {
            'name': 'Signed acceleration mixture (40% negative)',
            'injector': lambda gc, field: inject_signed_mixture(gc, field, target_negative_fraction=0.40, amplification_dex=1.0),
            'params': {'negative_fraction': 0.40, 'raw_dex': 1.0},
        },
    ]

    results = {
        'timestamp_utc': datetime.now(timezone.utc).isoformat(),
        'baseline': baseline,
        'scenarios': [],
    }

    for scenario in scenarios:
        name = scenario['name']
        print(f"--- INJECTION: {name} ---")
        gc_inj, field_inj = scenario['injector'](gc_rows_orig, field_rows_orig)
        injected = run_analysis_battery(gc_inj, field_inj, binary_flags, label=name)

        print(f"  Raw amplitude:         {injected['raw_amplitude_dex']:+.3f} dex")
        print(f"  Period-only residual:  {injected['period_only_residual_dex']:+.3f} dex")
        print(f"  Period+B residual:     {injected['period_bfield_residual_dex']:+.3f} dex")
        if injected['gamma_raw'] is not None:
            print(f"  Gamma (raw):           {injected['gamma_raw']:.3f} +/- {injected['gamma_raw_err']:.3f}")
        print(f"  Negative fraction:     {injected['negative_fraction']:.1%}")

        recovery = compute_recovery(injected, baseline, scenario['params'])
        print()
        print("  RECOVERY:")
        for k, v in recovery.items():
            if 'fraction' in k and isinstance(v, float):
                print(f"    {k}: {v:.2%}")
            else:
                print(f"    {k}: {v}")
        print()

        results['scenarios'].append({
            'name': name,
            'injected_params': scenario['params'],
            'injected_results': injected,
            'recovery': recovery,
        })

    # Save JSON
    with open(OUT_JSON, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"Saved JSON: {OUT_JSON}")

    # Save Markdown
    md = f"""# Step 50: Injection-Recovery Test

**Generated:** {results['timestamp_utc']}

## Baseline (no injection)

| Metric | Value |
|--------|-------|
| Raw amplitude | {baseline['raw_amplitude_dex']:+.3f} dex |
| Period-only residual | {baseline['period_only_residual_dex']:+.3f} dex |
| Period+B-field residual | {baseline['period_bfield_residual_dex']:+.3f} dex |
| Gamma (raw) | {baseline['gamma_raw']:.3f} +/- {baseline['gamma_raw_err']:.3f} |
| Negative fraction | {baseline['negative_fraction']:.1%} ({baseline['negative_n']}/{baseline['total_n']}) |

## Injection Scenarios

"""
    for sc in results['scenarios']:
        md += f"### {sc['name']}\n\n"
        md += "| Metric | Injected | Recovered | Recovery |\n"
        md += "|--------|----------|-----------|----------|\n"
        md += f"| Raw amplitude (dex) | {sc['injected_params'].get('raw_dex', 'N/A')} | {sc['injected_results']['raw_amplitude_dex']:+.3f} | {sc['recovery'].get('raw_amplitude_recovery_fraction', 'N/A')} |\n"
        md += f"| Period-only residual (dex) | {sc['injected_params'].get('period_only_dex', 'N/A')} | {sc['injected_results']['period_only_residual_dex']:+.3f} | {sc['recovery'].get('period_only_recovery_fraction', 'N/A')} |\n"
        md += f"| Period+B residual (dex) | {sc['injected_params'].get('period_bfield_dex', 'N/A')} | {sc['injected_results']['period_bfield_residual_dex']:+.3f} | {sc['recovery'].get('period_bfield_recovery_fraction', 'N/A')} |\n"
        md += f"| Gamma (raw) | {sc['injected_params'].get('gamma', 'N/A')} | {sc['injected_results']['gamma_raw']:.3f} | {sc['recovery'].get('gamma_raw_recovery_fraction', 'N/A')} |\n"
        md += f"| Negative fraction | {sc['injected_params'].get('negative_fraction', 'N/A')} | {sc['injected_results']['negative_fraction']:.1%} | {sc['recovery'].get('negative_fraction_delta', 'N/A')} |\n"
        md += "\n"

    with open(OUT_MD, 'w') as f:
        f.write(md)
    print(f"Saved Markdown: {OUT_MD}")

    print()
    print("=" * 78)
    print("SUMMARY")
    print("=" * 78)
    print()
    print("Recovery fractions (1.0 = perfect recovery, <0.7 = attenuation warning):")
    print()
    for sc in results['scenarios']:
        r = sc['recovery']
        print(f"  {sc['name']}:")
        if 'raw_amplitude_recovery_fraction' in r:
            print(f"    Raw amplitude:       {r['raw_amplitude_recovery_fraction']:.2%}")
        if 'period_only_recovery_fraction' in r:
            print(f"    Period-only resid:   {r['period_only_recovery_fraction']:.2%}")
        if 'period_bfield_recovery_fraction' in r:
            print(f"    Period+B resid:      {r['period_bfield_recovery_fraction']:.2%}")
        if 'gamma_raw_recovery_fraction' in r:
            print(f"    Gamma (raw):         {r['gamma_raw_recovery_fraction']:.2%}")
        if 'negative_fraction_delta' in r:
            print(f"    Neg fraction delta:  {r['negative_fraction_delta']:+.3f}")
        print()

    print("=" * 78)
    print("Done.")
    print("=" * 78)


if __name__ == '__main__':
    main()
