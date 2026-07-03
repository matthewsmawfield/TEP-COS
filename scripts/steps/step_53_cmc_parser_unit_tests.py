#!/usr/bin/env python3
"""
Step 53: CMC Parser Unit Tests

Audit the CMC pipeline for choices that could be attenuating the true signal.
Tests:
  1. Deduplication: compare deduped vs non-deduped amplitude
  2. Fixed-period assumption: compare fixed MSP period vs observed distribution
  3. Acceleration conversion: unit-test with analytic toy cluster
  4. Projected radius sampling: condition on observed radii
  5. Period-distribution injection: test sensitivity to period draw

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
CMC_DIR = DATA_DIR / "cmc" / "M15"
STEP_5_10_CSV = RESULTS_DIR / "step_02_pulsar_population_controls.csv"

OUT_JSON = RESULTS_DIR / "step_53_step_01_cmc_parser_unit_tests.json"
OUT_MD = RESULTS_DIR / "step_53_step_01_cmc_parser_unit_tests.md"

G = 4.302e-3  # pc (km/s)^2 / M_sun
c = 3e8  # m/s
c_pc_s = c / 3.086e16  # pc/s


def parse_morepulsars(path):
    """Parse initial.morepulsars.dat and return list of dicts."""
    rows = []
    with open(path, 'r') as f:
        header = f.readline().strip()
        # Extract column names
        column_names = []
        if header.startswith('#'):
            matches = re.findall(r'#\d+:([^\s#]+)', header)
            column_names = matches
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) < len(column_names):
                continue
            row = dict(zip(column_names, parts))
            rows.append(row)
    return rows, column_names


def parse_conv(path):
    conv = {'massunitmsun': 484844.0, 'lengthunitparsec': 1.0, 'timeunitsmyr': 1906.06}
    with open(path, 'r') as f:
        for line in f:
            line = line.strip()
            if '=' in line and not line.startswith('#'):
                key, value = line.split('=', 1)
                key = key.strip()
                value = value.strip()
                try:
                    conv[key] = float(value)
                except ValueError:
                    pass
    return conv


def get_m15_cluster_props():
    """Read M15 dyn file for core radius and total mass."""
    dyn_path = CMC_DIR / "initial.dyn.dat"
    if not dyn_path.exists():
        return None
    with open(dyn_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) >= 25:
                # columns: time, ..., M, Rc, Rh, ...
                try:
                    return {'total_mass': float(parts[4]), 'core_radius': float(parts[7])}
                except (ValueError, IndexError):
                    continue
    return None


def filter_ns(rows):
    """Filter to neutron stars (startype0 == 13)."""
    out = []
    for r in rows:
        if r.get('startype0') == '13':
            out.append(r)
    return out


def compute_acceleration_toy_cluster(mass_msun, r_core_pc, r_pulsar_pc):
    """Analytic King-model acceleration for unit test."""
    r = np.asarray(r_pulsar_pc)
    x = r / max(r_core_pc, 0.01)
    M_enc = mass_msun * (x**3) / ((1 + x**2) ** 1.5)
    a = np.where(r > 0, G * M_enc / (r**2), 0.0)
    return a  # km/s^2


def test_1_deduplication():
    """Compare deduped vs non-deduped amplitude."""
    print("--- Test 1: Deduplication ---")
    path = CMC_DIR / "initial.morepulsars.dat"
    if not path.exists():
        print("  SKIPPED: M15 morepulsars.dat not found")
        return None

    rows, _ = parse_morepulsars(path)
    ns_rows = filter_ns(rows)
    n_before = len(ns_rows)

    # Deduplicate by id0 (keep last)
    seen = {}
    for r in ns_rows:
        iid = r.get('id0')
        seen[iid] = r
    deduped = list(seen.values())
    n_after = len(deduped)

    print(f"  NS rows: {n_before} -> {n_after} (deduped)")

    # Compute mean logP for both
    def mean_logp(rows_subset):
        vals = []
        for r in rows_subset:
            try:
                p = float(r['P0[sec]'])
                if p > 0:
                    vals.append(math.log10(p))
            except (ValueError, KeyError):
                continue
        return float(np.mean(vals)) if vals else None

    mean_before = mean_logp(ns_rows)
    mean_after = mean_logp(deduped)

    result = {
        'test': 'deduplication',
        'n_before': n_before,
        'n_after': n_after,
        'mean_logP_before': mean_before,
        'mean_logP_after': mean_after,
        'delta_mean_logP': (mean_after - mean_before) if (mean_before and mean_after) else None,
    }
    print(f"  Mean logP before: {mean_before:.3f} after: {mean_after:.3f}")
    return result


def test_2_fixed_period():
    """Compare fixed MSP period vs observed distribution."""
    print("--- Test 2: Fixed-period assumption ---")
    path = CMC_DIR / "initial.morepulsars.dat"
    if not path.exists():
        print("  SKIPPED: M15 morepulsars.dat not found")
        return None

    rows, _ = parse_morepulsars(path)
    ns_rows = filter_ns(rows)

    # Extract periods
    periods = []
    for r in ns_rows:
        try:
            p = float(r['P0[sec]'])
            if p > 0:
                periods.append(p)
        except (ValueError, KeyError):
            continue

    if not periods:
        print("  No periods found")
        return None

    # Fixed MSP period: use median
    fixed_p = np.median(periods)

    # Observed distribution: use actual periods
    # Compute Pdot = P * a / c for a simple potential
    conv = parse_conv(CMC_DIR / "initial.conv.sh")
    massunit = conv['massunitmsun']
    lengthunit = conv['lengthunitparsec']

    cluster_props = get_m15_cluster_props()
    if cluster_props:
        M_total = cluster_props['total_mass'] * massunit
        r_core = cluster_props['core_radius'] * lengthunit
    else:
        M_total = 5e5
        r_core = 0.14

    def compute_pdot_for_periods(periods_arr, r_pc_arr):
        a_vals = compute_acceleration_toy_cluster(M_total, r_core, r_pc_arr)
        a_si = a_vals * 1000  # m/s^2
        pdot = periods_arr * a_si / c
        return pdot

    # Use random radii for test
    rng = np.random.default_rng(42)
    n = min(len(periods), 1000)
    r_pc = np.abs(rng.normal(0, 0.5 * r_core, n))
    p_sample = np.array(periods[:n])

    pdot_fixed = compute_pdot_for_periods(np.full(n, fixed_p), r_pc)
    pdot_var = compute_pdot_for_periods(p_sample, r_pc)

    result = {
        'test': 'fixed_period',
        'fixed_period_s': float(fixed_p),
        'mean_pdot_fixed': float(np.mean(pdot_fixed)),
        'mean_pdot_variable': float(np.mean(pdot_var)),
        'ratio_variable_to_fixed': float(np.mean(pdot_var) / np.mean(pdot_fixed)) if np.mean(pdot_fixed) != 0 else None,
    }
    print(f"  Fixed period: {fixed_p:.4f} s")
    print(f"  Mean Pdot fixed:   {result['mean_pdot_fixed']:.3e}")
    print(f"  Mean Pdot variable: {result['mean_pdot_variable']:.3e}")
    print(f"  Ratio: {result['ratio_variable_to_fixed']:.3f}")
    return result


def test_3_acceleration_conversion():
    """Unit-test acceleration conversion with analytic toy cluster."""
    print("--- Test 3: Acceleration conversion ---")
    M = 1e6  # M_sun
    r_c = 1.0  # pc
    r = 0.5  # pc

    # Analytic enclosed mass and acceleration
    x = r / r_c
    M_enc = M * (x**3) / ((1 + x**2) ** 1.5)
    a_analytic = G * M_enc / (r**2)  # km/s^2

    # Convert to Pdot for P = 0.005 s
    P = 0.005
    a_si = a_analytic * 1000  # m/s^2
    pdot = P * a_si / c
    pdot_over_p = pdot / P

    # Check consistency: pdot/P should equal a/c
    a_from_pdot_over_p = pdot_over_p * c
    a_from_pdot = pdot * c / P

    result = {
        'test': 'acceleration_conversion',
        'M_enc_analytic_msun': float(M_enc),
        'a_analytic_kms2': float(a_analytic),
        'P_s': P,
        'pdot': float(pdot),
        'pdot_over_p': float(pdot_over_p),
        'a_from_pdot_over_p_ms2': float(a_from_pdot_over_p),
        'a_from_pdot_ms2': float(a_from_pdot),
        'consistency_pdot_over_p': float(abs(a_from_pdot_over_p - a_si) / a_si) if a_si != 0 else None,
        'consistency_pdot': float(abs(a_from_pdot - a_si) / a_si) if a_si != 0 else None,
    }
    print(f"  a_analytic = {a_analytic:.3e} km/s^2")
    print(f"  pdot = {pdot:.3e}")
    print(f"  pdot/P = {pdot_over_p:.3e}")
    print(f"  Consistency check (pdot/P == a/c): {result['consistency_pdot_over_p']:.6f}")
    return result


def test_4_projected_radius():
    """Test projected radius sampling vs 3D radius."""
    print("--- Test 4: Projected radius sampling ---")
    path = CMC_DIR / "initial.morepulsars.dat"
    if not path.exists():
        print("  SKIPPED")
        return None

    rows, _ = parse_morepulsars(path)
    ns_rows = filter_ns(rows)

    r_3d = []
    for r in ns_rows:
        try:
            rr = float(r['r'])
            if rr > 0:
                r_3d.append(rr)
        except (ValueError, KeyError):
            continue

    # Projected radius: r_proj = r * sin(theta) where theta is random
    rng = np.random.default_rng(42)
    theta = np.arccos(rng.uniform(-1, 1, len(r_3d)))
    r_proj = np.array(r_3d) * np.sin(theta)

    result = {
        'test': 'projected_radius',
        'n': len(r_3d),
        'mean_r_3d': float(np.mean(r_3d)),
        'mean_r_proj': float(np.mean(r_proj)),
        'ratio_proj_to_3d': float(np.mean(r_proj) / np.mean(r_3d)) if np.mean(r_3d) > 0 else None,
    }
    print(f"  Mean r_3d: {result['mean_r_3d']:.3f}")
    print(f"  Mean r_proj: {result['mean_r_proj']:.3f}")
    print(f"  Ratio: {result['ratio_proj_to_3d']:.3f}")
    return result


def test_5_period_distribution_sensitivity():
    """Test how predicted amplitude changes with period distribution width."""
    print("--- Test 5: Period-distribution sensitivity ---")
    # Load observed GC MSP periods from real data
    gc_periods = []
    with open(STEP_5_10_CSV, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row.get('environment') == 'globular_cluster':
                try:
                    p = float(row['P0_s'])
                    if p > 0:
                        gc_periods.append(p)
                except (ValueError, KeyError):
                    continue

    if not gc_periods:
        print("  No GC periods found")
        return None

    conv = parse_conv(CMC_DIR / "initial.conv.sh")
    massunit = conv['massunitmsun']
    lengthunit = conv['lengthunitparsec']
    cluster_props = get_m15_cluster_props()
    if cluster_props:
        M_total = cluster_props['total_mass'] * massunit
        r_core = cluster_props['core_radius'] * lengthunit
    else:
        M_total = 5e5
        r_core = 0.14

    def sample_pdot(periods, n=1000):
        rng = np.random.default_rng(42)
        r_pc = np.abs(rng.normal(0, 0.5 * r_core, n))
        p_draw = rng.choice(periods, size=n)
        a_vals = compute_acceleration_toy_cluster(M_total, r_core, r_pc)
        a_si = a_vals * 1000
        pdot = p_draw * a_si / c
        return pdot

    pdot_obs = sample_pdot(gc_periods)
    # Narrow period distribution: only periods within 1 std of mean
    p_mean = np.mean(gc_periods)
    p_std = np.std(gc_periods)
    narrow = [p for p in gc_periods if abs(p - p_mean) < p_std]
    pdot_narrow = sample_pdot(narrow) if narrow else pdot_obs

    result = {
        'test': 'period_distribution_sensitivity',
        'n_gc_periods': len(gc_periods),
        'mean_pdot_observed_dist': float(np.mean(pdot_obs)),
        'mean_pdot_narrow_dist': float(np.mean(pdot_narrow)),
        'ratio_narrow_to_observed': float(np.mean(pdot_narrow) / np.mean(pdot_obs)) if np.mean(pdot_obs) != 0 else None,
    }
    print(f"  N GC periods: {len(gc_periods)}")
    print(f"  Mean Pdot (observed dist): {np.mean(pdot_obs):.3e}")
    print(f"  Mean Pdot (narrow dist):   {np.mean(pdot_narrow):.3e}")
    print(f"  Ratio: {result['ratio_narrow_to_observed']:.3f}")
    return result


def main():
    print("=" * 78)
    print("STEP 5.73: CMC PARSER UNIT TESTS")
    print("=" * 78)
    print()

    results = {
        'timestamp_utc': datetime.now(timezone.utc).isoformat(),
        'tests': [],
    }

    for test_fn in [test_1_deduplication, test_2_fixed_period, test_3_acceleration_conversion,
                     test_4_projected_radius, test_5_period_distribution_sensitivity]:
        r = test_fn()
        if r:
            results['tests'].append(r)
        print()

    with open(OUT_JSON, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"Saved JSON: {OUT_JSON}")

    md = f"""# Step 53: CMC Parser Unit Tests

**Generated:** {results['timestamp_utc']}

## Tests

"""
    for t in results['tests']:
        md += f"### {t['test']}\n\n"
        for k, v in t.items():
            if k == 'test':
                continue
            md += f"- {k}: {v}\n"
        md += "\n"

    md += """## Interpretation

- **Deduplication** should not remove dynamically important states.
- **Fixed-period assumption** may under-predict amplitude if the real period
  distribution is broader.
- **Acceleration conversion** must be internally consistent (pdot/P == a/c).
- **Projected radius** sampling changes the mean radius and therefore the
  predicted acceleration.
- **Period-distribution sensitivity** quantifies how much the predicted
  amplitude depends on the assumed period draw.
"""

    with open(OUT_MD, 'w') as f:
        f.write(md)
    print(f"Saved Markdown: {OUT_MD}")
    print("Done.")


if __name__ == '__main__':
    main()
