#!/usr/bin/env python3
"""
Step 09: Pdot/P Parallel Observable Analysis

Addresses the review concern that Pdot/P is the more direct physical observable
because acceleration enters naturally through:

    Pdot/P ~ a_parallel / c

This script runs the full primary test battery on log|Pdot/P| in parallel with
the existing log|Pdot| results, and reports side-by-side comparisons.

Tests performed:
  1. Base GC vs field comparison (t-test, Mann-Whitney)
  2. Period-matched bootstrap control (without replacement)
  3. Period+B-proxy matched bootstrap control (without replacement, standardized)
  4. Hybrid maximum analysis (using expanded field sample from step_06)
"""

import csv
import json
import math
from pathlib import Path
from datetime import datetime, timezone

import numpy as np
from scipy import stats

REPO_ROOT = Path(__file__).resolve().parents[2]
RESULTS_DIR = REPO_ROOT / "results" / "outputs"

OUT_JSON = RESULTS_DIR / "step_09_pdot_over_p_analysis.json"
OUT_CSV = RESULTS_DIR / "step_09_pdot_over_p_analysis.csv"
OUT_MD = RESULTS_DIR / "step_09_pdot_over_p_analysis.md"

STEP_5_10_CSV = RESULTS_DIR / "step_02_pulsar_population_controls.csv"
STEP_5_27_CSV = RESULTS_DIR / "step_06_hybrid_pulsar_sample.csv"


def load_step_02_data():
    """Load robust step_02 parsed data."""
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


def load_step_06_data():
    """Load hybrid maximum sample. Note: step_06 CSV uses P_ms, not P0_s."""
    gc_rows = []
    field_rows = []
    with open(STEP_5_27_CSV, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                # step_06 uses P_ms (period in ms), convert to P0_s (seconds)
                p_ms = float(row['P_ms'])
                p0_s = p_ms / 1000.0
                r = {
                    'source': row['source'],
                    'environment': row['environment'],
                    'cluster': row.get('cluster', ''),
                    'name': row['name'],
                    'P0_s': p0_s,
                    'P_ms': p_ms,
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


def compute_pdot_over_p(rows):
    """Compute log|Pdot/P| for each pulsar."""
    out = []
    for r in rows:
        p0 = r['P0_s']
        p1 = r['P1_sps']
        if p0 <= 0:
            continue
        pdot_over_p = p1 / p0
        if pdot_over_p == 0:
            continue
        log_abs = math.log10(abs(pdot_over_p))
        r2 = dict(r)
        r2['pdot_over_p'] = float(pdot_over_p)
        r2['log_abs_pdot_over_p'] = float(log_abs)
        out.append(r2)
    return out


def _base_comparison(gc_rows, field_rows):
    """Standard GC vs field comparison."""
    gc_y = np.array([r['log_abs_pdot_over_p'] for r in gc_rows])
    field_y = np.array([r['log_abs_pdot_over_p'] for r in field_rows])
    t_stat, t_p = stats.ttest_ind(gc_y, field_y, equal_var=False)
    mw_u, mw_p = stats.mannwhitneyu(gc_y, field_y, alternative='two-sided')
    return {
        'gc_mean': float(np.mean(gc_y)),
        'gc_std': float(np.std(gc_y, ddof=1)),
        'field_mean': float(np.mean(field_y)),
        'field_std': float(np.std(field_y, ddof=1)),
        'diff_dex': float(np.mean(gc_y) - np.mean(field_y)),
        't_stat': float(t_stat),
        't_p': float(t_p),
        'mw_u': float(mw_u),
        'mw_p': float(mw_p),
        'gc_n': int(len(gc_y)),
        'field_n': int(len(field_y)),
    }


def _period_matched_bootstrap(gc_rows, field_rows, n_boot=2000, seed=42):
    """Bootstrap period-matched comparison on log|Pdot/P| WITHOUT replacement."""
    rng = np.random.default_rng(seed)
    gc_logp = np.array([r['logP'] for r in gc_rows])
    gc_y = np.array([r['log_abs_pdot_over_p'] for r in gc_rows])
    field_logp = np.array([r['logP'] for r in field_rows])
    field_y = np.array([r['log_abs_pdot_over_p'] for r in field_rows])

    n_gc = len(gc_rows)
    n_field = len(field_rows)
    distances = np.zeros((n_gc, n_field))
    for i in range(n_gc):
        distances[i, :] = np.abs(field_logp - gc_logp[i])

    diffs = []
    for _ in range(n_boot):
        idx_gc = rng.integers(0, n_gc, size=n_gc)
        used_field = set()
        f_sel = []
        matched_gc_idx = []
        order = rng.permutation(len(idx_gc))
        for idx in order:
            i = idx_gc[idx]
            sorted_indices = np.argsort(distances[i, :])
            for j in sorted_indices:
                if j not in used_field:
                    used_field.add(j)
                    f_sel.append(field_y[j])
                    matched_gc_idx.append(idx)
                    break
        if len(f_sel) == 0:
            continue
        f_sel = np.array(f_sel)
        gc_matched_y = gc_y[idx_gc[matched_gc_idx]]
        diffs.append(float(np.mean(gc_matched_y) - np.mean(f_sel)))

    diffs = np.array(diffs)
    return {
        'n_boot': int(n_boot),
        'diff_mean': float(np.mean(diffs)),
        'diff_ci16': float(np.quantile(diffs, 0.16)),
        'diff_ci84': float(np.quantile(diffs, 0.84)),
        'p_two_sided': float(2 * min(np.mean(diffs <= 0), np.mean(diffs >= 0))),
    }


def _two_dim_match_bootstrap(gc_rows, field_rows, n_boot=2000, seed=42):
    """Bootstrap matching in standardized (logP, log_b_proxy) on log|Pdot/P|."""
    rng = np.random.default_rng(seed)
    gc_x = np.array([[r['logP'], r['log_b_proxy']] for r in gc_rows])
    gc_y = np.array([r['log_abs_pdot_over_p'] for r in gc_rows])
    field_x = np.array([[r['logP'], r['log_b_proxy']] for r in field_rows])
    field_y = np.array([r['log_abs_pdot_over_p'] for r in field_rows])

    combined_x = np.vstack([gc_x, field_x])
    means = np.mean(combined_x, axis=0)
    stds = np.std(combined_x, axis=0)
    gc_x_std = (gc_x - means) / stds
    field_x_std = (field_x - means) / stds

    n_gc = len(gc_rows)
    n_field = len(field_rows)
    distances = np.zeros((n_gc, n_field))
    for i in range(n_gc):
        dx = field_x_std[:, 0] - gc_x_std[i, 0]
        dy = field_x_std[:, 1] - gc_x_std[i, 1]
        distances[i, :] = np.sqrt(dx * dx + dy * dy)

    diffs = []
    for _ in range(n_boot):
        idx_gc = rng.integers(0, n_gc, size=n_gc)
        used_field = set()
        f_sel = []
        matched_gc_idx = []
        order = rng.permutation(len(idx_gc))
        for idx in order:
            i = idx_gc[idx]
            sorted_indices = np.argsort(distances[i, :])
            for j in sorted_indices:
                if j not in used_field:
                    used_field.add(j)
                    f_sel.append(field_y[j])
                    matched_gc_idx.append(idx)
                    break
        if len(f_sel) == 0:
            continue
        f_sel = np.array(f_sel)
        gc_matched_y = gc_y[idx_gc[matched_gc_idx]]
        diffs.append(float(np.mean(gc_matched_y) - np.mean(f_sel)))

    diffs = np.array(diffs)
    return {
        'n_boot': int(n_boot),
        'diff_mean': float(np.mean(diffs)),
        'diff_ci16': float(np.quantile(diffs, 0.16)),
        'diff_ci84': float(np.quantile(diffs, 0.84)),
        'p_two_sided': float(2 * min(np.mean(diffs <= 0), np.mean(diffs >= 0))),
    }


def run_analysis():
    print("Step 09: Pdot/P Parallel Observable Analysis")
    print("=" * 60)

    # --- Base sample (step_02) ---
    gc_base, field_base = load_step_02_data()
    gc_base = compute_pdot_over_p(gc_base)
    field_base = compute_pdot_over_p(field_base)
    print(f"  Base sample: {len(gc_base)} GC, {len(field_base)} field")

    base_comparison = _base_comparison(gc_base, field_base)
    print(f"  Base comparison: diff = {base_comparison['diff_dex']:.3f} dex, t-p = {base_comparison['t_p']:.2e}")

    period_matched = _period_matched_bootstrap(gc_base, field_base, n_boot=2000, seed=42)
    print(f"  Period-matched: diff = {period_matched['diff_mean']:.3f} dex, p = {period_matched['p_two_sided']:.4f}")

    bproxy_matched = _two_dim_match_bootstrap(gc_base, field_base, n_boot=2000, seed=42)
    print(f"  B-proxy matched: diff = {bproxy_matched['diff_mean']:.3f} dex, p = {bproxy_matched['p_two_sided']:.4f}")

    # --- Hybrid sample (step_06) ---
    gc_hybrid, field_hybrid = load_step_06_data()
    gc_hybrid = compute_pdot_over_p(gc_hybrid)
    field_hybrid = compute_pdot_over_p(field_hybrid)
    print(f"  Hybrid sample: {len(gc_hybrid)} GC, {len(field_hybrid)} field")

    hybrid_comparison = _base_comparison(gc_hybrid, field_hybrid)
    print(f"  Hybrid comparison: diff = {hybrid_comparison['diff_dex']:.3f} dex, t-p = {hybrid_comparison['t_p']:.2e}")

    hybrid_matched = _period_matched_bootstrap(gc_hybrid, field_hybrid, n_boot=5000, seed=42)
    print(f"  Hybrid period-matched: diff = {hybrid_matched['diff_mean']:.3f} dex, p = {hybrid_matched['p_two_sided']:.4f}")

    result = {
        'meta': {
            'timestamp_utc': datetime.now(timezone.utc).isoformat(),
            'note': 'Pdot/P parallel observable. All tests rerun on log|Pdot/P| alongside existing log|Pdot|.',
        },
        'base_sample': {
            'counts': {'gc_n': len(gc_base), 'field_n': len(field_base)},
            'base_comparison': base_comparison,
            'controls': {
                'period_matched': period_matched,
                'period_and_bproxy_matched': bproxy_matched,
            }
        },
        'hybrid_sample': {
            'counts': {'gc_n': len(gc_hybrid), 'field_n': len(field_hybrid)},
            'base_comparison': hybrid_comparison,
            'controls': {
                'period_matched': hybrid_matched,
            }
        },
        'manuscript_summary': {
            'base_raw_difference_dex': round(base_comparison['diff_dex'], 2),
            'base_period_matched_residual_dex': round(period_matched['diff_mean'], 2),
            'base_bproxy_matched_residual_dex': round(bproxy_matched['diff_mean'], 2),
            'hybrid_raw_difference_dex': round(hybrid_comparison['diff_dex'], 2),
            'hybrid_period_matched_residual_dex': round(hybrid_matched['diff_mean'], 2),
            'hybrid_p_value': hybrid_matched['p_two_sided'],
        }
    }

    with open(OUT_JSON, 'w') as f:
        json.dump(result, f, indent=2)
    print(f"  Saved JSON: {OUT_JSON}")

    # Write CSV
    with open(OUT_CSV, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['name', 'environment', 'cluster', 'P0_s', 'P1_sps',
                         'pdot_over_p', 'log_abs_pdot_over_p', 'logP', 'log_b_proxy'])
        for r in gc_base + field_base:
            writer.writerow([
                r['name'], r['environment'], r.get('cluster', ''),
                r['P0_s'], r['P1_sps'], r['pdot_over_p'],
                r['log_abs_pdot_over_p'], r['logP'], r['log_b_proxy']
            ])
    print(f"  Saved CSV: {OUT_CSV}")

    # Markdown summary
    md = f"""# Step 09: Pdot/P Parallel Observable Analysis

**Timestamp:** {result['meta']['timestamp_utc']}

## Base Sample (step_02)

Counts: {result['base_sample']['counts']['gc_n']} GC MSPs, {result['base_sample']['counts']['field_n']} field MSPs

### Base Comparison (log|Pdot/P|)

| Statistic | GC | Field | Difference |
|-----------|-----|-------|------------|
| Mean (dex) | {base_comparison['gc_mean']:.3f} | {base_comparison['field_mean']:.3f} | {base_comparison['diff_dex']:.3f} |
| Std (dex) | {base_comparison['gc_std']:.3f} | {base_comparison['field_std']:.3f} | — |
| t-test | — | — | t = {base_comparison['t_stat']:.2f}, p = {base_comparison['t_p']:.2e} |
| Mann-Whitney | — | — | p = {base_comparison['mw_p']:.2e} |

### Controls on log|Pdot/P|

| Control | Diff (dex) | 16th–84th CI | p (two-sided) |
|---------|------------|--------------|---------------|
| Period-matched | {period_matched['diff_mean']:.3f} | [{period_matched['diff_ci16']:.3f}, {period_matched['diff_ci84']:.3f}] | {period_matched['p_two_sided']:.4f} |
| Period+B-proxy | {bproxy_matched['diff_mean']:.3f} | [{bproxy_matched['diff_ci16']:.3f}, {bproxy_matched['diff_ci84']:.3f}] | {bproxy_matched['p_two_sided']:.4f} |

## Hybrid Sample (step_06)

Counts: {result['hybrid_sample']['counts']['gc_n']} GC MSPs, {result['hybrid_sample']['counts']['field_n']} field MSPs

### Base Comparison (log|Pdot/P|)

| Statistic | GC | Field | Difference |
|-----------|-----|-------|------------|
| Mean (dex) | {hybrid_comparison['gc_mean']:.3f} | {hybrid_comparison['field_mean']:.3f} | {hybrid_comparison['diff_dex']:.3f} |
| Std (dex) | {hybrid_comparison['gc_std']:.3f} | {hybrid_comparison['field_std']:.3f} | — |
| t-test | — | — | t = {hybrid_comparison['t_stat']:.2f}, p = {hybrid_comparison['t_p']:.2e} |
| Mann-Whitney | — | — | p = {hybrid_comparison['mw_p']:.2e} |

### Hybrid Period-Matched Control

- Mean difference: {hybrid_matched['diff_mean']:.3f} dex
- 16th–84th percentile CI: [{hybrid_matched['diff_ci16']:.3f}, {hybrid_matched['diff_ci84']:.3f}]
- Two-sided p: {hybrid_matched['p_two_sided']:.4f}

## Manuscript Summary

| Observable | Raw diff (dex) | Controlled residual (dex) | p-value |
|------------|----------------|---------------------------|---------|
| log|Pdot| (existing) | 0.59 | 0.606 (period-matched) | < 10⁻¹³ |
| log|Pdot/P| (this work) | {result['manuscript_summary']['base_raw_difference_dex']:.2f} | {result['manuscript_summary']['base_period_matched_residual_dex']:.2f} (period-matched) | {period_matched['p_two_sided']:.4f} |
| log|Pdot/P| hybrid | {result['manuscript_summary']['hybrid_raw_difference_dex']:.2f} | {result['manuscript_summary']['hybrid_period_matched_residual_dex']:.2f} (period-matched) | {result['manuscript_summary']['hybrid_p_value']:.4f} |
"""
    with open(OUT_MD, 'w') as f:
        f.write(md)
    print(f"  Saved MD: {OUT_MD}")
    print("Done.")


if __name__ == "__main__":
    run_analysis()
