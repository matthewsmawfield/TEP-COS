#!/usr/bin/env python3
"""
Step 5.26: MAXIMUM Pulsar Analysis

Cross-match Freire GC pulsars with ATNF to get P-dot values for ALL possible
MSPs, then compute statistics on the expanded sample.

Goal: Maximize the sample while maintaining rigor.
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

OUT_JSON = RESULTS_DIR / "step_5_26_maximum_pulsar_analysis.json"
OUT_CSV = RESULTS_DIR / "step_5_26_maximum_pulsar_sample.csv"
OUT_MD = RESULTS_DIR / "step_5_26_maximum_pulsar_analysis.md"

# Existing data files
FREIRE_TXT = DATA_DIR / "freire_gcpsr_2025.txt"
ATNF_HTML = DATA_DIR / "atnf_full_catalog.html"


def parse_freire_all():
    """Parse ALL Freire pulsars, including those without P-dot in Freire."""
    
    pulsars = []
    current_cluster = None
    
    with open(FREIRE_TXT, 'r', errors='ignore') as f:
        lines = f.readlines()
    
    for line in lines:
        line = line.strip()
        if not line:
            continue
        
        # Cluster header detection
        if not line.startswith('J') and not line.startswith('B'):
            if any(x in line for x in ['NGC', 'Terzan', 'M ', 'M1', 'M2', 'M3', 'M4', 'M5', 
                                        'M6', 'M7', 'M8', 'M9', 'Omega', 'Tuc', 'Pal', 'Liller',
                                        '47 ', 'IC ', 'Djorg']):
                current_cluster = line.split('(')[0].strip() if '(' in line else line.strip()
            continue
        
        # Parse pulsar line
        parts = line.split()
        if len(parts) < 3:
            continue
        
        name = parts[0]
        
        # Find period (ms range)
        period_ms = None
        pdot_freire = None
        
        for i, part in enumerate(parts[1:], 1):
            if part == '*' or part == 'i':
                continue
            
            # Try to parse as number
            clean = re.sub(r'\([^)]*\)', '', part)
            try:
                val = float(clean)
                if period_ms is None and 0.5 < val < 1000:
                    period_ms = val
                elif period_ms is not None and pdot_freire is None:
                    # This could be P-dot (in 10^-20 units)
                    if abs(val) < 1000:  # Reasonable P-dot range
                        pdot_freire = val * 1e-20
                    break
            except:
                pass
        
        if period_ms is not None:
            pulsars.append({
                'name': name,
                'cluster': current_cluster,
                'period_ms': period_ms,
                'period_s': period_ms / 1000.0,
                'pdot_freire': pdot_freire,
                'has_pdot_freire': pdot_freire is not None
            })
    
    return pulsars


def parse_atnf_for_crossmatch():
    """Parse ATNF catalog and build lookup by name."""
    
    if not ATNF_HTML.exists():
        # Try the psrcat.db file
        db_path = RESULTS_DIR / "atnf_psrcat.db"
        if db_path.exists():
            return parse_atnf_db(db_path)
        return {}
    
    with open(ATNF_HTML, 'r', errors='ignore') as f:
        content = f.read()
    
    # Remove HTML tags
    text = re.sub(r'<[^>]+>', ' ', content)
    
    # Pattern for pulsar entries
    pattern = r'([JB]\d{4}[+-]\d{2,4}[A-Za-z]*)\s+(\d+\.\d+(?:e[+-]?\d+)?)\s+([+-]?\d+\.\d+(?:e[+-]?\d+)?|\*)'
    
    atnf = {}
    for match in re.finditer(pattern, text):
        name, p0_str, p1_str = match.groups()
        try:
            p0 = float(p0_str)
            has_pdot = p1_str != '*'
            pdot = float(p1_str) if has_pdot else None
            
            atnf[name] = {
                'period_s': p0,
                'period_ms': p0 * 1000,
                'pdot': pdot,
                'has_pdot': has_pdot
            }
        except:
            pass
    
    return atnf


def parse_atnf_db(db_path):
    """Parse ATNF psrcat.db format."""
    
    with open(db_path, 'r', errors='ignore') as f:
        content = f.read()
    
    atnf = {}
    current = {}
    
    for line in content.splitlines():
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        
        if line.startswith('@'):
            # Flush current record
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
                
                if name and p0:
                    atnf[name] = {
                        'period_s': p0,
                        'period_ms': p0 * 1000,
                        'pdot': p1,
                        'has_pdot': p1 is not None
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
                except:
                    pass
    
    return atnf


def normalize_name(name):
    """Normalize pulsar name for matching."""
    name = name.upper().strip()
    # Handle J vs B prefixes and coordinate variations
    return name


def crossmatch_and_expand(freire_pulsars, atnf_lookup):
    """Cross-match Freire pulsars with ATNF to get P-dot where available."""
    
    expanded = []
    
    for fp in freire_pulsars:
        name = fp['name']
        
        # Try exact match
        atnf_data = atnf_lookup.get(name)
        
        # Try without letter suffix (e.g., J0024-7204A -> J0024-7204)
        if not atnf_data:
            base_name = re.sub(r'[A-Za-z]+$', '', name)
            for atnf_name in atnf_lookup:
                if atnf_name.startswith(base_name):
                    atnf_data = atnf_lookup[atnf_name]
                    break
        
        # Determine P-dot source
        pdot = fp['pdot_freire']
        pdot_source = 'freire' if pdot else None
        
        if not pdot and atnf_data and atnf_data['has_pdot']:
            pdot = atnf_data['pdot']
            pdot_source = 'atnf'
        
        if pdot and fp['period_ms'] < 30:  # MSP cut
            p0 = fp['period_s']
            p1_abs = abs(pdot)
            
            expanded.append({
                'name': name,
                'cluster': fp['cluster'],
                'period_s': p0,
                'period_ms': fp['period_ms'],
                'pdot': pdot,
                'pdot_source': pdot_source,
                'logP': math.log10(p0),
                'logPdot_abs': math.log10(p1_abs),
                'b_proxy': math.sqrt(p0 * p1_abs),
                'log_b_proxy': math.log10(math.sqrt(p0 * p1_abs)),
                'environment': 'globular_cluster'
            })
    
    return expanded


def get_field_msps(atnf_lookup, gc_names):
    """Get field MSPs from ATNF excluding GC pulsars."""
    
    field = []
    
    # Common GC name patterns
    gc_patterns = ['NGC', 'M ', 'M1', 'M2', 'M3', 'M4', 'M5', 'M6', 'M7', 'M8', 'M9',
                   'Terzan', 'Omega', 'Tuc', '47 Tuc', 'Pal', 'Liller', 'Djorg']
    
    for name, data in atnf_lookup.items():
        # Skip if in GC list
        if name in gc_names:
            continue
        
        # Skip if looks like GC pulsar (letter suffix common in GCs)
        if re.search(r'[A-Z]$', name) and any(p in name for p in gc_patterns):
            continue
        
        # MSP cut and P-dot required
        if data['period_ms'] < 30 and data['has_pdot']:
            p0 = data['period_s']
            p1_abs = abs(data['pdot'])
            
            field.append({
                'name': name,
                'cluster': None,
                'period_s': p0,
                'period_ms': data['period_ms'],
                'pdot': data['pdot'],
                'pdot_source': 'atnf',
                'logP': math.log10(p0),
                'logPdot_abs': math.log10(p1_abs),
                'b_proxy': math.sqrt(p0 * p1_abs),
                'log_b_proxy': math.log10(math.sqrt(p0 * p1_abs)),
                'environment': 'field'
            })
    
    return field


def ttest_logpdot(gc: np.ndarray, field: np.ndarray) -> dict:
    """Compute statistics for log|Pdot| comparison."""
    t_stat, p_value = stats.ttest_ind(gc, field, equal_var=False)
    mw_u, mw_p = stats.mannwhitneyu(gc, field, alternative="two-sided")
    
    # Bootstrap CI
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
        "t_stat": float(t_stat),
        "t_p": float(p_value),
        "mw_u": float(mw_u),
        "mw_p": float(mw_p),
        "gc_n": int(len(gc)),
        "field_n": int(len(field)),
    }


def two_dim_match_bootstrap(gc_rows: list, field_rows: list, n_boot=5000, seed=42) -> dict:
    """Bootstrap matching in (logP, log_b_proxy)."""
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
    print("MAXIMUM PULSAR ANALYSIS")
    print("="*70)
    
    # Parse Freire - ALL pulsars
    print("\nParsing Freire catalog (all pulsars)...")
    freire_pulsars = parse_freire_all()
    freire_msps = [p for p in freire_pulsars if p['period_ms'] < 30]
    print(f"  Total Freire pulsars: {len(freire_pulsars)}")
    print(f"  Freire MSPs (P<30ms): {len(freire_msps)}")
    print(f"  With P-dot in Freire: {sum(1 for p in freire_msps if p['has_pdot_freire'])}")
    
    # Parse ATNF for cross-match
    print("\nParsing ATNF catalog...")
    atnf_lookup = parse_atnf_for_crossmatch()
    print(f"  ATNF entries: {len(atnf_lookup)}")
    
    # Cross-match to expand GC sample
    print("\nCross-matching to expand GC sample...")
    gc_expanded = crossmatch_and_expand(freire_msps, atnf_lookup)
    print(f"  GC MSPs with P-dot (expanded): {len(gc_expanded)}")
    print(f"    From Freire: {sum(1 for p in gc_expanded if p['pdot_source'] == 'freire')}")
    print(f"    From ATNF:   {sum(1 for p in gc_expanded if p['pdot_source'] == 'atnf')}")
    
    # Get field MSPs
    print("\nGetting field MSPs...")
    gc_names = set(p['name'] for p in gc_expanded)
    field_msps = get_field_msps(atnf_lookup, gc_names)
    print(f"  Field MSPs with P-dot: {len(field_msps)}")
    
    # Compute statistics
    print("\n" + "="*70)
    print("COMPUTING STATISTICS")
    print("="*70)
    
    gc_logpdot = np.array([r["logPdot_abs"] for r in gc_expanded])
    field_logpdot = np.array([r["logPdot_abs"] for r in field_msps])
    
    base = ttest_logpdot(gc_logpdot, field_logpdot)
    
    print(f"\n--- RAW COMPARISON ---")
    print(f"  GC mean log|Ṗ|:    {base['gc_mean']:.3f} ± {base['gc_std']:.3f}")
    print(f"  Field mean log|Ṗ|: {base['field_mean']:.3f} ± {base['field_std']:.3f}")
    print(f"  Difference:        {base['diff_dex']:.3f} dex")
    print(f"  95% CI:            [{base['diff_ci_95'][0]:.3f}, {base['diff_ci_95'][1]:.3f}]")
    print(f"  p-value:           {base['t_p']:.2e}")
    
    # Population-controlled analysis
    print("\nRunning period + B-proxy matched bootstrap (5000 iterations)...")
    two_dim_match = two_dim_match_bootstrap(gc_expanded, field_msps)
    
    print(f"\n--- CONTROLLED RESIDUAL (Period + B-proxy matched) ---")
    print(f"  Residual:          {two_dim_match['diff_mean']:.3f} dex")
    print(f"  95% CI:            [{two_dim_match['diff_ci_2_5']:.3f}, {two_dim_match['diff_ci_97_5']:.3f}]")
    print(f"  p-value:           {two_dim_match['p_two_sided']:.2e}")
    
    # Final summary
    print("\n" + "="*70)
    print("MAXIMUM SAMPLE SUMMARY")
    print("="*70)
    
    total = len(gc_expanded) + len(field_msps)
    print(f"""
    ┌─────────────────────────────────────────────────────────────────────┐
    │  MAXIMUM PULSAR SAMPLE                                              │
    ├─────────────────────────────────────────────────────────────────────┤
    │                                                                     │
    │  GC MSPs with Ṗ:      {len(gc_expanded):>4}  (Freire + ATNF cross-match)        │
    │  Field MSPs with Ṗ:   {len(field_msps):>4}  (ATNF non-GC)                       │
    │  ─────────────────────────────                                      │
    │  TOTAL:               {total:>4}                                            │
    │                                                                     │
    │  Raw difference:      {base['diff_dex']:>+.2f} dex                                      │
    │  Controlled residual: {two_dim_match['diff_mean']:>+.2f} dex (95% CI: {two_dim_match['diff_ci_2_5']:.2f}–{two_dim_match['diff_ci_97_5']:.2f})       │
    │  p-value:             {two_dim_match['p_two_sided']:.1e}                                      │
    │                                                                     │
    └─────────────────────────────────────────────────────────────────────┘
    """)
    
    # Save results
    output = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "sample_sizes": {
            "gc_msps": len(gc_expanded),
            "gc_from_freire": sum(1 for p in gc_expanded if p['pdot_source'] == 'freire'),
            "gc_from_atnf": sum(1 for p in gc_expanded if p['pdot_source'] == 'atnf'),
            "field_msps": len(field_msps),
            "total": total
        },
        "base_comparison": base,
        "controlled_residual": two_dim_match,
        "manuscript_values": {
            "gc_count": len(gc_expanded),
            "field_count": len(field_msps),
            "total_count": total,
            "raw_difference_dex": round(base['diff_dex'], 2),
            "controlled_residual_dex": round(two_dim_match['diff_mean'], 2),
            "controlled_ci_95": [round(two_dim_match['diff_ci_2_5'], 2), round(two_dim_match['diff_ci_97_5'], 2)],
            "controlled_p_value": two_dim_match['p_two_sided'],
            "base_p_value": base['t_p']
        }
    }
    
    with open(OUT_JSON, 'w') as f:
        json.dump(output, f, indent=2)
    
    # Save CSV
    with open(OUT_CSV, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['name', 'cluster', 'environment', 'period_ms', 'pdot', 'pdot_source', 
                    'logP', 'logPdot_abs', 'log_b_proxy'])
        for r in gc_expanded + field_msps:
            w.writerow([r['name'], r.get('cluster', ''), r['environment'], 
                       f"{r['period_ms']:.4f}", f"{r['pdot']:.4e}", r['pdot_source'],
                       f"{r['logP']:.4f}", f"{r['logPdot_abs']:.4f}", f"{r['log_b_proxy']:.4f}"])
    
    # Save markdown
    md = f"""# Maximum Pulsar Analysis Results

**Generated:** {output['timestamp_utc']}

## Sample Sizes (MAXIMUM)
- **GC MSPs:** {len(gc_expanded)} (Freire + ATNF cross-match)
  - From Freire: {output['sample_sizes']['gc_from_freire']}
  - From ATNF: {output['sample_sizes']['gc_from_atnf']}
- **Field MSPs:** {len(field_msps)} (ATNF non-GC)
- **TOTAL:** {total}

## Raw Comparison
- **Difference:** {base['diff_dex']:.2f} dex (GC higher)
- **95% CI:** [{base['diff_ci_95'][0]:.2f}, {base['diff_ci_95'][1]:.2f}]
- **p-value:** {base['t_p']:.2e}

## Controlled Residual (Period + B-proxy matched)
- **Residual:** {two_dim_match['diff_mean']:.2f} dex
- **95% CI:** [{two_dim_match['diff_ci_2_5']:.2f}, {two_dim_match['diff_ci_97_5']:.2f}]
- **p-value:** {two_dim_match['p_two_sided']:.2e}
"""
    
    with open(OUT_MD, 'w') as f:
        f.write(md)
    
    print(f"\nResults saved to:")
    print(f"  {OUT_JSON}")
    print(f"  {OUT_CSV}")
    print(f"  {OUT_MD}")
    
    return output


if __name__ == "__main__":
    main()
