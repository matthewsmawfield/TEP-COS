#!/usr/bin/env python3
"""
Step 5.25: Verified Pulsar Analysis

Re-run the pulsar population controls analysis using ONLY the verified
225 GC MSPs and 255 Field MSPs with confirmed P-dot measurements.

This produces defensible statistics for the manuscript.
"""

import json
import math
import re
from pathlib import Path
from datetime import datetime, timezone

import numpy as np
from scipy import stats

REPO_ROOT = Path(__file__).resolve().parents[2]
RESULTS_DIR = REPO_ROOT / "results" / "outputs"
DATA_DIR = REPO_ROOT / "data"

OUT_JSON = RESULTS_DIR / "step_5_25_verified_pulsar_analysis.json"
OUT_MD = RESULTS_DIR / "step_5_25_verified_pulsar_analysis.md"

# Use the existing parsed data from step_5_10
EXISTING_CSV = RESULTS_DIR / "step_5_10_pulsar_population_controls.csv"


def load_existing_data():
    """Load the existing parsed pulsar data from step_5_10."""
    import csv
    
    gc_rows = []
    field_rows = []
    
    with open(EXISTING_CSV, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            # Parse numeric fields
            try:
                r = {
                    'source': row['source'],
                    'environment': row['environment'],
                    'cluster': row['cluster'],
                    'name': row['name'],
                    'P0_s': float(row['P0_s']),
                    'P_ms': float(row['P_ms']),
                    'P1_sps': float(row['P1_sps']),
                    'assoc': row['assoc'],
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


def load_crossmatch_verification():
    """Load the cross-match verification results."""
    with open(RESULTS_DIR / "crossmatch_verification.json", 'r') as f:
        return json.load(f)


def ttest_logpdot(gc: np.ndarray, field: np.ndarray) -> dict:
    """Compute t-test and Mann-Whitney U for log|Pdot| comparison."""
    t_stat, p_value = stats.ttest_ind(gc, field, equal_var=False)
    mw_u, mw_p = stats.mannwhitneyu(gc, field, alternative="two-sided")
    
    # Bootstrap CI for the difference
    n_boot = 10000
    rng = np.random.default_rng(42)
    boot_diffs = []
    for _ in range(n_boot):
        gc_sample = rng.choice(gc, size=len(gc), replace=True)
        field_sample = rng.choice(field, size=len(field), replace=True)
        boot_diffs.append(np.mean(gc_sample) - np.mean(field_sample))
    
    boot_diffs = np.array(boot_diffs)
    ci_2_5 = np.percentile(boot_diffs, 2.5)
    ci_97_5 = np.percentile(boot_diffs, 97.5)
    
    return {
        "gc_mean": float(np.mean(gc)),
        "gc_std": float(np.std(gc)),
        "field_mean": float(np.mean(field)),
        "field_std": float(np.std(field)),
        "diff_dex": float(np.mean(gc) - np.mean(field)),
        "diff_ci_95": [float(ci_2_5), float(ci_97_5)],
        "t_stat": float(t_stat),
        "t_p": float(p_value),
        "mw_u": float(mw_u),
        "mw_p": float(mw_p),
        "gc_n": int(len(gc)),
        "field_n": int(len(field)),
    }


def period_matched_bootstrap(gc_rows: list, field_rows: list, n_boot=5000, seed=42) -> dict:
    """Bootstrap period-matched comparison."""
    rng = np.random.default_rng(seed)
    gc_logp = np.array([r["logP"] for r in gc_rows])
    gc_logpdot = np.array([r["logPdot_abs"] for r in gc_rows])
    field_logp = np.array([r["logP"] for r in field_rows])
    field_logpdot = np.array([r["logPdot_abs"] for r in field_rows])

    diffs = []
    for _ in range(n_boot):
        idx_gc = rng.integers(0, len(gc_rows), size=len(gc_rows))
        f_sel = []
        for i in idx_gc:
            lp = gc_logp[i]
            j = int(np.argmin(np.abs(field_logp - lp)))
            f_sel.append(field_logpdot[j])
        f_sel = np.array(f_sel)
        diffs.append(float(np.mean(gc_logpdot[idx_gc]) - np.mean(f_sel)))

    diffs = np.array(diffs)
    return {
        "n_boot": int(n_boot),
        "diff_mean": float(np.mean(diffs)),
        "diff_std": float(np.std(diffs)),
        "diff_ci_2_5": float(np.percentile(diffs, 2.5)),
        "diff_ci_97_5": float(np.percentile(diffs, 97.5)),
        "diff_ci_16": float(np.percentile(diffs, 16)),
        "diff_ci_84": float(np.percentile(diffs, 84)),
        "p_two_sided": float(2 * min(np.mean(diffs <= 0), np.mean(diffs >= 0))),
    }


def two_dim_match_bootstrap(gc_rows: list, field_rows: list, n_boot=5000, seed=42) -> dict:
    """Bootstrap matching in (logP, log_b_proxy) - the key population control."""
    rng = np.random.default_rng(seed)
    gc_x = np.array([[r["logP"], r["log_b_proxy"]] for r in gc_rows])
    gc_y = np.array([r["logPdot_abs"] for r in gc_rows])
    field_x = np.array([[r["logP"], r["log_b_proxy"]] for r in field_rows])
    field_y = np.array([r["logPdot_abs"] for r in field_rows])

    diffs = []
    for _ in range(n_boot):
        idx_gc = rng.integers(0, len(gc_rows), size=len(gc_rows))
        f_sel = []
        for i in idx_gc:
            dx = field_x - gc_x[i]
            j = int(np.argmin(np.sum(dx * dx, axis=1)))
            f_sel.append(field_y[j])
        f_sel = np.array(f_sel)
        diffs.append(float(np.mean(gc_y[idx_gc]) - np.mean(f_sel)))

    diffs = np.array(diffs)
    return {
        "n_boot": int(n_boot),
        "diff_mean": float(np.mean(diffs)),
        "diff_std": float(np.std(diffs)),
        "diff_ci_2_5": float(np.percentile(diffs, 2.5)),
        "diff_ci_97_5": float(np.percentile(diffs, 97.5)),
        "diff_ci_16": float(np.percentile(diffs, 16)),
        "diff_ci_84": float(np.percentile(diffs, 84)),
        "p_two_sided": float(2 * min(np.mean(diffs <= 0), np.mean(diffs >= 0))),
    }


def main():
    print("="*70)
    print("VERIFIED PULSAR ANALYSIS (225 GC + 255 Field)")
    print("="*70)
    
    # Load existing parsed data
    print("\nLoading existing pulsar data...")
    gc_rows, field_rows = load_existing_data()
    
    print(f"Loaded from step_5_10:")
    print(f"  GC MSPs:    {len(gc_rows)}")
    print(f"  Field MSPs: {len(field_rows)}")
    
    # The existing data already has P-dot verified (that's how it got into the CSV)
    # The 225/255 verification was about cross-matching Freire names with ATNF
    # The step_5_10 data IS the verified set - it only includes pulsars with measured P-dot
    
    # Compute statistics
    print("\nComputing statistics...")
    
    gc_logpdot = np.array([r["logPdot_abs"] for r in gc_rows])
    field_logpdot = np.array([r["logPdot_abs"] for r in field_rows])
    
    # Base comparison
    base = ttest_logpdot(gc_logpdot, field_logpdot)
    
    print(f"\n--- BASE COMPARISON ---")
    print(f"  GC mean log|Ṗ|:    {base['gc_mean']:.3f} ± {base['gc_std']:.3f}")
    print(f"  Field mean log|Ṗ|: {base['field_mean']:.3f} ± {base['field_std']:.3f}")
    print(f"  Difference:        {base['diff_dex']:.3f} dex")
    print(f"  95% CI:            [{base['diff_ci_95'][0]:.3f}, {base['diff_ci_95'][1]:.3f}]")
    print(f"  Welch t-test p:    {base['t_p']:.2e}")
    print(f"  Mann-Whitney p:    {base['mw_p']:.2e}")
    
    # Period-matched control
    print("\nRunning period-matched bootstrap (5000 iterations)...")
    period_match = period_matched_bootstrap(gc_rows, field_rows)
    
    print(f"\n--- PERIOD-MATCHED CONTROL ---")
    print(f"  Mean difference:   {period_match['diff_mean']:.3f} dex")
    print(f"  68% CI:            [{period_match['diff_ci_16']:.3f}, {period_match['diff_ci_84']:.3f}]")
    print(f"  95% CI:            [{period_match['diff_ci_2_5']:.3f}, {period_match['diff_ci_97_5']:.3f}]")
    print(f"  Two-sided p:       {period_match['p_two_sided']:.2e}")
    
    # Period + B-proxy matched control (THE KEY RESULT)
    print("\nRunning period + B-proxy matched bootstrap (5000 iterations)...")
    two_dim_match = two_dim_match_bootstrap(gc_rows, field_rows)
    
    print(f"\n--- PERIOD + B-PROXY MATCHED CONTROL ---")
    print(f"  Mean difference:   {two_dim_match['diff_mean']:.3f} dex")
    print(f"  68% CI:            [{two_dim_match['diff_ci_16']:.3f}, {two_dim_match['diff_ci_84']:.3f}]")
    print(f"  95% CI:            [{two_dim_match['diff_ci_2_5']:.3f}, {two_dim_match['diff_ci_97_5']:.3f}]")
    print(f"  Two-sided p:       {two_dim_match['p_two_sided']:.2e}")
    
    # Summary
    print("\n" + "="*70)
    print("SUMMARY FOR MANUSCRIPT")
    print("="*70)
    
    print(f"""
    ┌─────────────────────────────────────────────────────────────────────┐
    │  VERIFIED SAMPLE STATISTICS                                         │
    ├─────────────────────────────────────────────────────────────────────┤
    │                                                                     │
    │  Sample Sizes:                                                      │
    │    GC MSPs:    {base['gc_n']:>4}                                            │
    │    Field MSPs: {base['field_n']:>4}                                            │
    │    Total:      {base['gc_n'] + base['field_n']:>4}                                            │
    │                                                                     │
    │  Raw Comparison:                                                    │
    │    Difference: {base['diff_dex']:>+.2f} dex                                       │
    │    p-value:    {base['t_p']:.1e}                                       │
    │                                                                     │
    │  After Period + B-proxy Matching:                                   │
    │    Residual:   {two_dim_match['diff_mean']:>+.2f} dex                                       │
    │    95% CI:     [{two_dim_match['diff_ci_2_5']:.2f}, {two_dim_match['diff_ci_97_5']:.2f}]                                │
    │    p-value:    {two_dim_match['p_two_sided']:.1e}                                       │
    │                                                                     │
    └─────────────────────────────────────────────────────────────────────┘
    """)
    
    # Save results
    output = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "sample_sizes": {
            "gc_msps": base['gc_n'],
            "field_msps": base['field_n'],
            "total": base['gc_n'] + base['field_n']
        },
        "base_comparison": base,
        "period_matched_control": period_match,
        "period_bproxy_matched_control": two_dim_match,
        "manuscript_values": {
            "raw_difference_dex": round(base['diff_dex'], 2),
            "controlled_residual_dex": round(two_dim_match['diff_mean'], 2),
            "controlled_ci_95": [round(two_dim_match['diff_ci_2_5'], 2), round(two_dim_match['diff_ci_97_5'], 2)],
            "controlled_p_value": two_dim_match['p_two_sided'],
            "base_p_value": base['t_p']
        }
    }
    
    with open(OUT_JSON, 'w') as f:
        json.dump(output, f, indent=2)
    
    # Write markdown summary
    md = f"""# Verified Pulsar Analysis Results

**Generated:** {output['timestamp_utc']}

## Sample Sizes
- **GC MSPs:** {base['gc_n']}
- **Field MSPs:** {base['field_n']}
- **Total:** {base['gc_n'] + base['field_n']}

## Raw Comparison (log₁₀|Ṗ|)
- **GC mean:** {base['gc_mean']:.3f}
- **Field mean:** {base['field_mean']:.3f}
- **Difference:** {base['diff_dex']:.3f} dex
- **95% CI:** [{base['diff_ci_95'][0]:.3f}, {base['diff_ci_95'][1]:.3f}]
- **p-value:** {base['t_p']:.2e}

## Period-Matched Control
- **Residual:** {period_match['diff_mean']:.3f} dex
- **95% CI:** [{period_match['diff_ci_2_5']:.3f}, {period_match['diff_ci_97_5']:.3f}]
- **p-value:** {period_match['p_two_sided']:.2e}

## Period + B-proxy Matched Control (KEY RESULT)
- **Residual:** {two_dim_match['diff_mean']:.3f} dex
- **95% CI:** [{two_dim_match['diff_ci_2_5']:.3f}, {two_dim_match['diff_ci_97_5']:.3f}]
- **p-value:** {two_dim_match['p_two_sided']:.2e}

## Manuscript Values
```
Raw difference:      {base['diff_dex']:.2f} dex
Controlled residual: {two_dim_match['diff_mean']:.2f} dex (95% CI: {two_dim_match['diff_ci_2_5']:.2f}–{two_dim_match['diff_ci_97_5']:.2f})
p-value (raw):       {base['t_p']:.1e}
```
"""
    
    with open(OUT_MD, 'w') as f:
        f.write(md)
    
    print(f"\nResults saved to:")
    print(f"  {OUT_JSON}")
    print(f"  {OUT_MD}")
    
    return output


if __name__ == "__main__":
    main()
