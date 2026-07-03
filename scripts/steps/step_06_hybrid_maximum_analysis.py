#!/usr/bin/env python3
"""
Step 06: Hybrid Maximum Analysis

Start with robust step_02 data, then expand by:
1. Cross-matching Freire pulsars without P-dot against ATNF
2. Adding field MSPs from ATNF that weren't in the original sample

This gives us MAXIMUM pulsars with RELIABLE P-dot values.
"""

import csv
import json
import math
import re
from pathlib import Path
from datetime import datetime, timezone
from collections import defaultdict

import numpy as np
from scipy import stats

REPO_ROOT = Path(__file__).resolve().parents[2]
RESULTS_DIR = REPO_ROOT / "results" / "outputs"
DATA_DIR = REPO_ROOT / "data"

OUT_JSON = RESULTS_DIR / "step_06_hybrid_maximum_analysis.json"
OUT_CSV = RESULTS_DIR / "step_06_hybrid_pulsar_sample.csv"
OUT_MD = RESULTS_DIR / "step_06_hybrid_maximum_analysis.md"

# Existing data
STEP_5_10_CSV = RESULTS_DIR / "step_02_pulsar_population_controls.csv"
FREIRE_TXT = DATA_DIR / "freire_GCpsr.txt"


def load_step_02_data():
    """Load the robust step_02 parsed data."""
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


def parse_freire_names_without_pdot():
    """Get Freire MSP names that DON'T have P-dot in Freire (marked with *)."""
    
    missing_pdot = []
    current_cluster = None
    
    with open(FREIRE_TXT, 'r', errors='ignore') as f:
        lines = f.readlines()
    
    for line in lines:
        line_stripped = line.strip()
        if not line_stripped:
            continue
        
        # Cluster header
        if not line_stripped.startswith('J') and not line_stripped.startswith('B'):
            if any(x in line_stripped for x in ['NGC', 'Terzan', 'M ', 'M1', 'M2', 'M3', 'M4', 'M5', 
                                        'M6', 'M7', 'M8', 'M9', 'Omega', 'Tuc', 'Pal', 'Liller',
                                        '47 ', 'IC ', 'Djorg']):
                current_cluster = line_stripped.split('(')[0].strip() if '(' in line_stripped else line_stripped.strip()
            continue
        
        # Parse pulsar line - tab-separated
        parts = [p for p in line.split('\t') if p.strip()]
        if len(parts) < 4:
            parts = line.split()
        if len(parts) < 4:
            continue
        
        name = parts[0].strip()
        
        # Check period (column 2 or 3)
        period_ms = None
        pdot_col = None
        
        for i, part in enumerate(parts[1:], 1):
            part = part.strip()
            if part in ('*', 'i', ''):
                continue
            
            # Try to parse as number
            clean = re.sub(r'\([^)]*\)', '', part)
            try:
                val = float(clean)
                if period_ms is None and 0.5 < val < 1000:
                    period_ms = val
                    # P-dot is typically the next column
                    if i + 1 < len(parts):
                        pdot_col = parts[i + 1].strip() if i + 1 < len(parts) else None
                    break
            except (ValueError, TypeError):
                pass
        
        # Check if P-dot is missing (marked with * or i)
        is_msp = period_ms is not None and period_ms < 30
        pdot_missing = pdot_col in ('*', 'i', None, '')
        
        if is_msp and pdot_missing:
            missing_pdot.append({
                'name': name,
                'cluster': current_cluster,
                'period_ms': period_ms
            })
    
    return missing_pdot


def parse_atnf_db():
    """Parse ATNF psrcat.db for P-dot values."""
    
    db_path = DATA_DIR / "atnf_psrcat.db"
    if not db_path.exists():
        print(f"Warning: ATNF db not found at {db_path}")
        return {}
    
    atnf = {}
    current = {}
    
    with open(db_path, 'r', errors='ignore') as f:
        content = f.read()
    
    for line in content.splitlines():
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        
        if line.startswith('@'):
            if current:
                name = current.get('PSRJ') or current.get('PSRB')
                p0 = current.get('P0')
                p1 = current.get('P1')
                f0 = current.get('F0')
                f1 = current.get('F1')
                
                if p0 is None and f0:
                    p0 = 1.0 / f0
                if p1 is None and f0 and f1:
                    p1 = -f1 / (f0 * f0)
                
                if name and p0 and p1 is not None:
                    atnf[name] = {
                        'period_s': p0,
                        'period_ms': p0 * 1000,
                        'pdot': p1
                    }
            current = {}
            continue
        
        parts = line.split()
        if len(parts) >= 2:
            key = parts[0]
            val = parts[1]
            
            if key in ('PSRJ', 'PSRB'):
                current[key] = val
            elif key in ('P0', 'P1', 'F0', 'F1'):
                try:
                    current[key] = float(val)
                except (ValueError, TypeError):
                    pass
    
    return atnf


def crossmatch_missing(missing_pdot, atnf, existing_gc_names):
    """Find P-dot from ATNF for Freire pulsars missing P-dot."""
    
    found = []
    
    for pulsar in missing_pdot:
        name = pulsar['name']
        
        # Skip if already in existing sample
        if name in existing_gc_names:
            continue
        
        # Try exact match
        atnf_data = atnf.get(name)
        
        # Try variations
        if not atnf_data:
            # Try base name without letter suffix
            base = re.sub(r'[A-Za-z]+$', '', name)
            for atnf_name, data in atnf.items():
                if atnf_name.startswith(base) and atnf_name.endswith(name[-1]):
                    atnf_data = data
                    break
        
        if atnf_data and atnf_data['pdot'] is not None and atnf_data['pdot'] != 0:
            p0 = atnf_data['period_s']
            p1 = atnf_data['pdot']
            p1_abs = abs(p1)
            if p1_abs <= 0 or p0 <= 0:
                continue
            
            found.append({
                'source': 'atnf_crossmatch',
                'environment': 'globular_cluster',
                'cluster': pulsar['cluster'],
                'name': name,
                'P0_s': p0,
                'P_ms': p0 * 1000,
                'P1_sps': p1,
                'assoc': pulsar['cluster'],
                'logP': math.log10(p0),
                'logPdot_abs': math.log10(p1_abs),
                'log_b_proxy': math.log10(math.sqrt(p0 * p1_abs)),
                'log_tau_c': math.log10(p0 / (2 * p1_abs)),
            })
    
    return found


def get_additional_field_msps(atnf, existing_field_names, gc_names):
    """Get additional field MSPs from ATNF not in original sample."""
    
    additional = []
    
    for name, data in atnf.items():
        # Skip if already in sample
        if name in existing_field_names:
            continue
        
        # Skip if GC pulsar
        if name in gc_names:
            continue
        
        # Skip if looks like GC (letter suffix + GC-like name)
        if re.search(r'[A-Z]$', name):
            continue
        
        # MSP cut
        if data['period_ms'] >= 30:
            continue
        
        # Must have P-dot
        if data['pdot'] is None:
            continue
        
        p0 = data['period_s']
        p1 = data['pdot']
        p1_abs = abs(p1)
        
        additional.append({
            'source': 'atnf_additional',
            'environment': 'field',
            'cluster': None,
            'name': name,
            'P0_s': p0,
            'P_ms': p0 * 1000,
            'P1_sps': p1,
            'assoc': '',
            'logP': math.log10(p0),
            'logPdot_abs': math.log10(p1_abs),
            'log_b_proxy': math.log10(math.sqrt(p0 * p1_abs)),
            'log_tau_c': math.log10(p0 / (2 * p1_abs)),
        })
    
    return additional


def ttest_logpdot(gc: np.ndarray, field: np.ndarray) -> dict:
    """Compute statistics."""
    t_stat, p_value = stats.ttest_ind(gc, field, equal_var=False)
    mw_u, mw_p = stats.mannwhitneyu(gc, field, alternative="two-sided")
    
    n_boot = 10000
    rng = np.random.default_rng(42)
    boot_diffs = []
    for _ in range(n_boot):
        gc_sample = rng.choice(gc, size=len(gc), replace=True)
        field_sample = rng.choice(field, size=len(field), replace=True)
        boot_diffs.append(np.mean(gc_sample) - np.mean(field_sample))
    
    boot_diffs = np.array(boot_diffs)
    
    return {
        "gc_mean": float(np.mean(gc)),
        "gc_std": float(np.std(gc)),
        "field_mean": float(np.mean(field)),
        "field_std": float(np.std(field)),
        "diff_dex": float(np.mean(gc) - np.mean(field)),
        "diff_ci_95": [float(np.percentile(boot_diffs, 2.5)), float(np.percentile(boot_diffs, 97.5))],
        "t_p": float(p_value),
        "mw_p": float(mw_p),
        "gc_n": int(len(gc)),
        "field_n": int(len(field)),
    }


def two_dim_match_bootstrap(gc_rows: list, field_rows: list, n_boot=5000, seed=42) -> dict:
    """Bootstrap matching in (logP, log_b_proxy) WITHOUT replacement.
    
    Each GC pulsar is matched to a unique field pulsar to avoid bias from overmatching.
    """
    rng = np.random.default_rng(seed)
    gc_x = np.array([[r["logP"], r["log_b_proxy"]] for r in gc_rows])
    gc_y = np.array([r["logPdot_abs"] for r in gc_rows])
    field_x = np.array([[r["logP"], r["log_b_proxy"]] for r in field_rows])
    field_y = np.array([r["logPdot_abs"] for r in field_rows])

    # Pre-compute all pairwise distances
    n_gc = len(gc_rows)
    n_field = len(field_rows)
    distances = np.zeros((n_gc, n_field))
    for i in range(n_gc):
        dx = field_x[:, 0] - gc_x[i, 0]
        dy = field_x[:, 1] - gc_x[i, 1]
        distances[i, :] = np.sqrt(dx*dx + dy*dy)

    diffs = []
    for _ in range(n_boot):
        # Resample GC pulsars with replacement (bootstrap)
        idx_gc = rng.integers(0, n_gc, size=n_gc)
        
        # Match WITHOUT replacement
        used_field = set()
        f_sel = []
        
        # Randomize order to avoid systematic bias
        order = rng.permutation(len(idx_gc))
        
        for idx in order:
            i = idx_gc[idx]
            # Find nearest unused field pulsar
            sorted_indices = np.argsort(distances[i, :])
            for j in sorted_indices:
                if j not in used_field:
                    used_field.add(j)
                    f_sel.append(field_y[j])
                    break
        
        f_sel = np.array(f_sel)
        diffs.append(float(np.mean(gc_y[idx_gc]) - np.mean(f_sel)))

    diffs = np.array(diffs)
    # Compute p-value with floor at 1/n_boot (minimum resolvable for bootstrap)
    p_left = np.mean(diffs <= 0)
    p_right = np.mean(diffs >= 0)
    p_two_sided_raw = 2 * min(p_left, p_right)
    p_floor = 1.0 / n_boot  # Minimum resolvable p-value
    p_two_sided = max(p_two_sided_raw, p_floor)  # Ensure we don't report p=0.0
    
    return {
        "n_boot": int(n_boot),
        "diff_mean": float(np.mean(diffs)),
        "diff_ci_2_5": float(np.percentile(diffs, 2.5)),
        "diff_ci_97_5": float(np.percentile(diffs, 97.5)),
        "p_two_sided": float(p_two_sided),
        "p_two_sided_note": f"min_resolvable={p_floor:.1e}" if p_two_sided_raw < p_floor else None,
    }


def main():
    """Run hybrid maximum pulsar analysis combining base and cross-matched samples.
    
    Loads robust base data, cross-matches with ATNF for additional pulsars,
    and computes both raw and population-controlled statistics.
    """
    print("="*70)
    print("HYBRID MAXIMUM PULSAR ANALYSIS")
    print("="*70)
    
    # Load robust base data
    print("\nLoading step_02 robust data...")
    gc_base, field_base = load_step_02_data()
    print(f"  Base GC MSPs:    {len(gc_base)}")
    print(f"  Base Field MSPs: {len(field_base)}")
    
    existing_gc_names = set(r['name'] for r in gc_base)
    existing_field_names = set(r['name'] for r in field_base)
    
    # Find Freire pulsars without P-dot
    print("\nFinding Freire MSPs without P-dot in Freire...")
    missing_pdot = parse_freire_names_without_pdot()
    print(f"  Freire MSPs missing P-dot: {len(missing_pdot)}")
    
    # Parse ATNF
    print("\nParsing ATNF catalog...")
    atnf = parse_atnf_db()
    print(f"  ATNF entries with P-dot: {len(atnf)}")
    
    # Cross-match to find P-dot
    print("\nCross-matching to find P-dot from ATNF...")
    gc_additional = crossmatch_missing(missing_pdot, atnf, existing_gc_names)
    print(f"  Additional GC MSPs from cross-match: {len(gc_additional)}")
    
    # Get additional field MSPs
    print("\nFinding additional field MSPs...")
    all_gc_names = existing_gc_names | set(r['name'] for r in gc_additional)
    field_additional = get_additional_field_msps(atnf, existing_field_names, all_gc_names)
    print(f"  Additional Field MSPs: {len(field_additional)}")
    
    # Combine samples
    gc_all = gc_base + gc_additional
    field_all = field_base + field_additional
    
    print("\n" + "="*70)
    print("MAXIMUM SAMPLE")
    print("="*70)
    print(f"  GC MSPs:    {len(gc_all)} (base {len(gc_base)} + additional {len(gc_additional)})")
    print(f"  Field MSPs: {len(field_all)} (base {len(field_base)} + additional {len(field_additional)})")
    print(f"  TOTAL:      {len(gc_all) + len(field_all)}")
    
    # Compute statistics
    print("\nComputing statistics...")
    
    gc_logpdot = np.array([r["logPdot_abs"] for r in gc_all])
    field_logpdot = np.array([r["logPdot_abs"] for r in field_all])
    
    base = ttest_logpdot(gc_logpdot, field_logpdot)
    
    print(f"\n--- RAW COMPARISON ---")
    print(f"  GC mean log|Ṗ|:    {base['gc_mean']:.3f}")
    print(f"  Field mean log|Ṗ|: {base['field_mean']:.3f}")
    print(f"  Difference:        {base['diff_dex']:.3f} dex")
    print(f"  95% CI:            [{base['diff_ci_95'][0]:.3f}, {base['diff_ci_95'][1]:.3f}]")
    print(f"  p-value:           {base['t_p']:.2e}")
    
    print("\nRunning population-controlled bootstrap...")
    controlled = two_dim_match_bootstrap(gc_all, field_all)
    
    print(f"\n--- CONTROLLED RESIDUAL ---")
    print(f"  Residual:          {controlled['diff_mean']:.3f} dex")
    print(f"  95% CI:            [{controlled['diff_ci_2_5']:.3f}, {controlled['diff_ci_97_5']:.3f}]")
    print(f"  p-value:           {controlled['p_two_sided']:.2e}")
    
    # Summary
    total = len(gc_all) + len(field_all)
    print("\n" + "="*70)
    print("FINAL MAXIMUM SAMPLE")
    print("="*70)
    print(f"""
    ┌─────────────────────────────────────────────────────────────────────┐
    │  MAXIMUM PULSAR SAMPLE (Hybrid Approach)                            │
    ├─────────────────────────────────────────────────────────────────────┤
    │                                                                     │
    │  GC MSPs:             {len(gc_all):>4}  ({len(gc_base)} base + {len(gc_additional)} ATNF cross-match)     │
    │  Field MSPs:          {len(field_all):>4}  ({len(field_base)} base + {len(field_additional)} additional)        │
    │  ─────────────────────────────                                      │
    │  TOTAL:               {total:>4}                                            │
    │                                                                     │
    │  Raw difference:      {base['diff_dex']:>+.2f} dex                                      │
    │  Controlled residual: {controlled['diff_mean']:>+.2f} dex (95% CI: {controlled['diff_ci_2_5']:.2f}–{controlled['diff_ci_97_5']:.2f})       │
    │  p-value:             {controlled['p_two_sided']:.1e}                                      │
    │                                                                     │
    └─────────────────────────────────────────────────────────────────────┘
    """)
    
    # Save results
    output = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "sample_sizes": {
            "gc_base": len(gc_base),
            "gc_additional": len(gc_additional),
            "gc_total": len(gc_all),
            "field_base": len(field_base),
            "field_additional": len(field_additional),
            "field_total": len(field_all),
            "total": total
        },
        "base_comparison": base,
        "controlled_residual": controlled,
        "manuscript_values": {
            "gc_count": len(gc_all),
            "field_count": len(field_all),
            "total_count": total,
            "raw_difference_dex": round(base['diff_dex'], 2),
            "controlled_residual_dex": round(controlled['diff_mean'], 2),
            "controlled_ci_95": [round(controlled['diff_ci_2_5'], 2), round(controlled['diff_ci_97_5'], 2)],
            "p_value": controlled['p_two_sided']
        }
    }
    
    with open(OUT_JSON, 'w') as f:
        json.dump(output, f, indent=2)
    
    # Save CSV
    with open(OUT_CSV, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['source', 'environment', 'cluster', 'name', 'P_ms', 'P1_sps', 
                    'logP', 'logPdot_abs', 'log_b_proxy'])
        for r in gc_all + field_all:
            w.writerow([r['source'], r['environment'], r.get('cluster', ''), r['name'],
                       f"{r['P_ms']:.4f}", f"{r['P1_sps']:.4e}",
                       f"{r['logP']:.4f}", f"{r['logPdot_abs']:.4f}", f"{r['log_b_proxy']:.4f}"])
    
    print(f"\nResults saved to:")
    print(f"  {OUT_JSON}")
    print(f"  {OUT_CSV}")
    
    return output


if __name__ == "__main__":
    main()
