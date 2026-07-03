#!/usr/bin/env python3
"""
Step 08: Signed P-dot Latent-Mixture Analysis

Addresses the review concern that log|Pdot| hides sign information.
This script:
  1. Uses signed Pdot (not just log|Pdot|)
  2. Computes Pdot/P as the direct physical observable (acceleration enters as a||/c)
  3. Reports positive/negative fractions, median signed Pdot/P
  4. Fits a 2-component Gaussian mixture model to log|Pdot/P|
  5. Infers the latent acceleration distribution by deconvolving the field intrinsic distribution

Model:
    Pdot_obs = Pdot_int + P * a_parallel / c + epsilon
    Pdot_obs / P = Pdot_int / P + a_parallel / c + epsilon / P

In the field, a_parallel/c ≈ 0, so the field distribution traces the intrinsic spin-down.
In clusters, the observed distribution is a convolution of intrinsic + acceleration.
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
DATA_DIR = REPO_ROOT / "data"

OUT_JSON = RESULTS_DIR / "step_08_signed_pdot_analysis.json"
OUT_CSV = RESULTS_DIR / "step_08_signed_pdot_analysis.csv"
OUT_MD = RESULTS_DIR / "step_08_signed_pdot_analysis.md"

STEP_5_10_CSV = RESULTS_DIR / "step_02_pulsar_population_controls.csv"
STEP_5_27_CSV = RESULTS_DIR / "step_06_hybrid_pulsar_sample.csv"

C = 299792458.0  # m/s


def load_step_02_data():
    """Load robust step_02 parsed data with signed Pdot."""
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


def compute_signed_pdot_over_p(rows):
    """Compute signed Pdot/P for each pulsar."""
    out = []
    for r in rows:
        p0 = r['P0_s']
        p1 = r['P1_sps']
        if p0 <= 0:
            continue
        pdot_over_p = p1 / p0  # s^-1, dimensionally a||/c-like
        log_abs_pdot_over_p = math.log10(abs(pdot_over_p)) if pdot_over_p != 0 else -99.0
        r2 = dict(r)
        r2['pdot_over_p'] = float(pdot_over_p)
        r2['log_abs_pdot_over_p'] = float(log_abs_pdot_over_p)
        r2['pdot_sign'] = int(1 if pdot_over_p > 0 else (-1 if pdot_over_p < 0 else 0))
        out.append(r2)
    return out


def sign_fraction_analysis(gc_rows, field_rows):
    """Report positive/negative fractions and median signed values."""
    gc_signed = [r for r in gc_rows if r['pdot_sign'] != 0]
    field_signed = [r for r in field_rows if r['pdot_sign'] != 0]

    gc_pos = sum(1 for r in gc_signed if r['pdot_sign'] > 0)
    gc_neg = sum(1 for r in gc_signed if r['pdot_sign'] < 0)
    field_pos = sum(1 for r in field_signed if r['pdot_sign'] > 0)
    field_neg = sum(1 for r in field_signed if r['pdot_sign'] < 0)

    gc_frac_pos = gc_pos / len(gc_signed) if gc_signed else 0
    field_frac_pos = field_pos / len(field_signed) if field_signed else 0

    gc_pdot_over_p = np.array([r['pdot_over_p'] for r in gc_signed])
    field_pdot_over_p = np.array([r['pdot_over_p'] for r in field_signed])

    gc_log_abs = np.array([r['log_abs_pdot_over_p'] for r in gc_signed])
    field_log_abs = np.array([r['log_abs_pdot_over_p'] for r in field_signed])

    # Fisher exact test on sign contingency table
    contingency = [[gc_pos, gc_neg], [field_pos, field_neg]]
    _, fisher_p = stats.fisher_exact(contingency, alternative='two-sided')

    return {
        'gc': {
            'n_total': len(gc_signed),
            'n_positive': gc_pos,
            'n_negative': gc_neg,
            'fraction_positive': float(gc_frac_pos),
            'fraction_negative': float(1 - gc_frac_pos),
            'median_signed_pdot_over_p_s': float(np.median(gc_pdot_over_p)),
            'mean_signed_pdot_over_p_s': float(np.mean(gc_pdot_over_p)),
            'std_signed_pdot_over_p_s': float(np.std(gc_pdot_over_p, ddof=1)),
            'median_log_abs_pdot_over_p': float(np.median(gc_log_abs)),
            'mean_log_abs_pdot_over_p': float(np.mean(gc_log_abs)),
            'std_log_abs_pdot_over_p': float(np.std(gc_log_abs, ddof=1)),
        },
        'field': {
            'n_total': len(field_signed),
            'n_positive': field_pos,
            'n_negative': field_neg,
            'fraction_positive': float(field_frac_pos),
            'fraction_negative': float(1 - field_frac_pos),
            'median_signed_pdot_over_p_s': float(np.median(field_pdot_over_p)),
            'mean_signed_pdot_over_p_s': float(np.mean(field_pdot_over_p)),
            'std_signed_pdot_over_p_s': float(np.std(field_pdot_over_p, ddof=1)),
            'median_log_abs_pdot_over_p': float(np.median(field_log_abs)),
            'mean_log_abs_pdot_over_p': float(np.mean(field_log_abs)),
            'std_log_abs_pdot_over_p': float(np.std(field_log_abs, ddof=1)),
        },
        'difference': {
            'delta_median_signed_pdot_over_p_s': float(np.median(gc_pdot_over_p) - np.median(field_pdot_over_p)),
            'delta_mean_signed_pdot_over_p_s': float(np.mean(gc_pdot_over_p) - np.mean(field_pdot_over_p)),
            'delta_fraction_positive': float(gc_frac_pos - field_frac_pos),
            'fisher_exact_p': float(fisher_p),
        }
    }


def _gaussian_mixture_em(data, n_components=2, n_iter=100, tol=1e-6, seed=42):
    """Fit Gaussian mixture via EM algorithm (numpy only, no sklearn)."""
    rng = np.random.default_rng(seed)
    n = len(data)
    # Initialize means spread across data range
    sorted_d = np.sort(data)
    means = np.array([sorted_d[n // 4], sorted_d[3 * n // 4]], dtype=float)
    stds = np.array([np.std(data), np.std(data)], dtype=float)
    weights = np.array([0.5, 0.5], dtype=float)

    for _ in range(n_iter):
        # E-step: compute responsibilities
        resp = np.zeros((n, n_components))
        for k in range(n_components):
            if stds[k] > 0:
                resp[:, k] = weights[k] * stats.norm.pdf(data, means[k], stds[k])
        resp_sum = resp.sum(axis=1, keepdims=True)
        resp_sum[resp_sum == 0] = 1e-300
        resp = resp / resp_sum

        # M-step: update parameters
        Nk = resp.sum(axis=0)
        new_weights = Nk / n
        new_means = (resp * data[:, None]).sum(axis=0) / Nk
        new_stds = np.sqrt((resp * (data[:, None] - new_means) ** 2).sum(axis=0) / Nk)

        # Enforce minimum std to avoid collapse
        new_stds = np.maximum(new_stds, 1e-6)

        delta = np.max(np.abs(means - new_means)) + np.max(np.abs(stds - new_stds))
        means, stds, weights = new_means, new_stds, new_weights
        if delta < tol:
            break

    log_likelihood = np.sum(np.log(np.maximum(
        np.sum([weights[k] * stats.norm.pdf(data, means[k], stds[k]) for k in range(n_components)], axis=0),
        1e-300
    )))

    return {
        'means': [float(m) for m in means],
        'stds': [float(s) for s in stds],
        'weights': [float(w) for w in weights],
        'log_likelihood': float(log_likelihood),
    }


def latent_mixture_analysis(gc_rows, field_rows):
    """Fit 2-component Gaussian mixture to log|Pdot/P| for GC and field separately."""
    gc_data = np.array([r['log_abs_pdot_over_p'] for r in gc_rows])
    field_data = np.array([r['log_abs_pdot_over_p'] for r in field_rows])

    gc_mix = _gaussian_mixture_em(gc_data, seed=42)
    field_mix = _gaussian_mixture_em(field_data, seed=43)

    # Inferred acceleration distribution: deconvolve field from GC
    # If GC ~ N(mu_gc_i, sigma_gc_i^2) and field ~ N(mu_f, sigma_f^2)
    # Then acceleration ~ N(mu_gc_i - mu_f, sigma_gc_i^2 - sigma_f^2)
    # We attribute the broader component to acceleration + intrinsic
    gc_broad_idx = np.argmax(gc_mix['stds'])
    field_broad_idx = np.argmax(field_mix['stds'])

    inferred_mean = gc_mix['means'][gc_broad_idx] - field_mix['means'][field_broad_idx]
    inferred_var = gc_mix['stds'][gc_broad_idx] ** 2 - field_mix['stds'][field_broad_idx] ** 2
    inferred_std = math.sqrt(max(inferred_var, 0))

    return {
        'gc_mixture': gc_mix,
        'field_mixture': field_mix,
        'inferred_acceleration_component': {
            'mean_dex': float(inferred_mean),
            'std_dex': float(inferred_std),
            'note': 'Broad component of GC minus broad component of field, under convolution assumption',
        }
    }


def _period_matched_bootstrap_pdot_over_p(gc_rows, field_rows, n_boot=2000, seed=42):
    """Bootstrap period-matched comparison on log|Pdot/P| WITHOUT replacement."""
    rng = np.random.default_rng(seed)
    gc_logp = np.array([r['logP'] for r in gc_rows])
    gc_log_y = np.array([r['log_abs_pdot_over_p'] for r in gc_rows])

    field_logp = np.array([r['logP'] for r in field_rows])
    field_log_y = np.array([r['log_abs_pdot_over_p'] for r in field_rows])

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
                    f_sel.append(field_log_y[j])
                    matched_gc_idx.append(idx)
                    break
        if len(f_sel) == 0:
            continue
        f_sel = np.array(f_sel)
        gc_matched_y = gc_log_y[idx_gc[matched_gc_idx]]
        diffs.append(float(np.mean(gc_matched_y) - np.mean(f_sel)))

    diffs = np.array(diffs)
    return {
        'n_boot': int(n_boot),
        'diff_mean': float(np.mean(diffs)),
        'diff_ci16': float(np.quantile(diffs, 0.16)),
        'diff_ci84': float(np.quantile(diffs, 0.84)),
        'p_two_sided': float(2 * min(np.mean(diffs <= 0), np.mean(diffs >= 0))),
    }


def _two_dim_match_bootstrap_pdot_over_p(gc_rows, field_rows, n_boot=2000, seed=42):
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
    print("Step 08: Signed P-dot Latent-Mixture Analysis")
    print("=" * 60)

    gc_rows, field_rows = load_step_02_data()
    print(f"  Loaded {len(gc_rows)} GC MSPs, {len(field_rows)} field MSPs")

    gc_rows = compute_signed_pdot_over_p(gc_rows)
    field_rows = compute_signed_pdot_over_p(field_rows)
    print(f"  After Pdot/P computation: {len(gc_rows)} GC, {len(field_rows)} field")

    sign_results = sign_fraction_analysis(gc_rows, field_rows)
    print(f"  GC sign fractions: +{sign_results['gc']['fraction_positive']:.3f}, -{sign_results['gc']['fraction_negative']:.3f}")
    print(f"  Field sign fractions: +{sign_results['field']['fraction_positive']:.3f}, -{sign_results['field']['fraction_negative']:.3f}")
    print(f"  Fisher exact p = {sign_results['difference']['fisher_exact_p']:.2e}")

    mixture_results = latent_mixture_analysis(gc_rows, field_rows)
    print(f"  GC mixture: means = {[f'{m:.3f}' for m in mixture_results['gc_mixture']['means']]}, stds = {[f'{s:.3f}' for s in mixture_results['gc_mixture']['stds']]}")
    print(f"  Field mixture: means = {[f'{m:.3f}' for m in mixture_results['field_mixture']['means']]}, stds = {[f'{s:.3f}' for s in mixture_results['field_mixture']['stds']]}")

    # GC vs field comparison on log|Pdot/P|
    gc_log_y = np.array([r['log_abs_pdot_over_p'] for r in gc_rows])
    field_log_y = np.array([r['log_abs_pdot_over_p'] for r in field_rows])
    t_stat, t_p = stats.ttest_ind(gc_log_y, field_log_y, equal_var=False)
    mw_u, mw_p = stats.mannwhitneyu(gc_log_y, field_log_y, alternative='two-sided')

    base_comparison = {
        'gc_mean': float(np.mean(gc_log_y)),
        'field_mean': float(np.mean(field_log_y)),
        'diff_dex': float(np.mean(gc_log_y) - np.mean(field_log_y)),
        't_stat': float(t_stat),
        't_p': float(t_p),
        'mw_u': float(mw_u),
        'mw_p': float(mw_p),
        'gc_n': int(len(gc_log_y)),
        'field_n': int(len(field_log_y)),
    }
    print(f"  log|Pdot/P| base comparison: diff = {base_comparison['diff_dex']:.3f} dex, t-p = {base_comparison['t_p']:.2e}")

    period_matched = _period_matched_bootstrap_pdot_over_p(gc_rows, field_rows, n_boot=2000, seed=42)
    print(f"  Period-matched bootstrap: diff = {period_matched['diff_mean']:.3f} dex, p = {period_matched['p_two_sided']:.4f}")

    bproxy_matched = _two_dim_match_bootstrap_pdot_over_p(gc_rows, field_rows, n_boot=2000, seed=42)
    print(f"  B-proxy matched bootstrap: diff = {bproxy_matched['diff_mean']:.3f} dex, p = {bproxy_matched['p_two_sided']:.4f}")

    result = {
        'meta': {
            'timestamp_utc': datetime.now(timezone.utc).isoformat(),
            'note': 'Signed Pdot analysis using Pdot/P as the physical observable. C = 299792458 m/s.',
            'sample_counts': {'gc_n': len(gc_rows), 'field_n': len(field_rows)},
        },
        'sign_analysis': sign_results,
        'latent_mixture': mixture_results,
        'base_log10_abs_pdot_over_p': base_comparison,
        'controls': {
            'period_matched': period_matched,
            'period_and_bproxy_matched': bproxy_matched,
        },
    }

    with open(OUT_JSON, 'w') as f:
        json.dump(result, f, indent=2)
    print(f"  Saved JSON: {OUT_JSON}")

    # Write CSV
    with open(OUT_CSV, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['name', 'environment', 'cluster', 'P0_s', 'P1_sps', 'pdot_over_p',
                         'log_abs_pdot_over_p', 'pdot_sign', 'logP', 'log_b_proxy'])
        for r in gc_rows + field_rows:
            writer.writerow([
                r['name'], r['environment'], r.get('cluster', ''),
                r['P0_s'], r['P1_sps'], r['pdot_over_p'],
                r['log_abs_pdot_over_p'], r['pdot_sign'],
                r['logP'], r['log_b_proxy']
            ])
    print(f"  Saved CSV: {OUT_CSV}")

    # Write Markdown summary
    md = f"""# Step 08: Signed P-dot Latent-Mixture Analysis

**Timestamp:** {result['meta']['timestamp_utc']}

## Sample Counts

- GC MSPs: {result['meta']['sample_counts']['gc_n']}
- Field MSPs: {result['meta']['sample_counts']['field_n']}

## Sign Analysis (Signed Pdot/P)

| Environment | N | Positive | Negative | Frac. Positive | Median signed Pdot/P (s⁻¹) |
|-------------|---|----------|----------|----------------|----------------------------|
| GC          | {sign_results['gc']['n_total']} | {sign_results['gc']['n_positive']} | {sign_results['gc']['n_negative']} | {sign_results['gc']['fraction_positive']:.3f} | {sign_results['gc']['median_signed_pdot_over_p_s']:.3e} |
| Field       | {sign_results['field']['n_total']} | {sign_results['field']['n_positive']} | {sign_results['field']['n_negative']} | {sign_results['field']['fraction_positive']:.3f} | {sign_results['field']['median_signed_pdot_over_p_s']:.3e} |

**Difference in positive fraction:** {sign_results['difference']['delta_fraction_positive']:.3f}  
**Fisher exact test p-value:** {sign_results['difference']['fisher_exact_p']:.2e}

## Latent Mixture (2-component Gaussian on log|Pdot/P|)

### GC Mixture
- Means (dex): {[f'{m:.3f}' for m in mixture_results['gc_mixture']['means']]}
- Stds (dex): {[f'{s:.3f}' for s in mixture_results['gc_mixture']['stds']]}
- Weights: {[f'{w:.3f}' for w in mixture_results['gc_mixture']['weights']]}

### Field Mixture
- Means (dex): {[f'{m:.3f}' for m in mixture_results['field_mixture']['means']]}
- Stds (dex): {[f'{s:.3f}' for s in mixture_results['field_mixture']['stds']]}
- Weights: {[f'{w:.3f}' for w in mixture_results['field_mixture']['weights']]}

### Inferred Acceleration Component
- Mean offset (dex): {mixture_results['inferred_acceleration_component']['mean_dex']:.3f}
- Std (dex): {mixture_results['inferred_acceleration_component']['std_dex']:.3f}

## log|Pdot/P| Base Comparison

- GC mean: {base_comparison['gc_mean']:.3f} dex
- Field mean: {base_comparison['field_mean']:.3f} dex
- Difference: {base_comparison['diff_dex']:.3f} dex
- t-statistic: {base_comparison['t_stat']:.2f}, p = {base_comparison['t_p']:.2e}
- Mann-Whitney U p = {base_comparison['mw_p']:.2e}

## Controls on log|Pdot/P|

### Period-Matched Bootstrap (N_boot=2000)
- Mean difference: {period_matched['diff_mean']:.3f} dex
- 16th–84th percentile CI: [{period_matched['diff_ci16']:.3f}, {period_matched['diff_ci84']:.3f}]
- Two-sided p: {period_matched['p_two_sided']:.4f}

### Period+B-Proxy Matched Bootstrap (N_boot=2000)
- Mean difference: {bproxy_matched['diff_mean']:.3f} dex
- 16th–84th percentile CI: [{bproxy_matched['diff_ci16']:.3f}, {bproxy_matched['diff_ci84']:.3f}]
- Two-sided p: {bproxy_matched['p_two_sided']:.4f}
"""
    with open(OUT_MD, 'w') as f:
        f.write(md)
    print(f"  Saved MD: {OUT_MD}")
    print("Done.")


if __name__ == "__main__":
    run_analysis()
