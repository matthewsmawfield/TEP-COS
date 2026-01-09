#!/usr/bin/env python3
"""
Step 5.11: Binary Pulsar Analysis for TEP

This script performs two key analyses:
1. Binary vs Isolated MSP comparison within GCs
2. Orbital parameter analysis for GC binaries

The goal is to identify cleaner probes of TEP that are less contaminated
by line-of-sight acceleration than raw Ṗ measurements.

Author: M. Smawfield
Date: 2026-01-03
"""

import json
import hashlib
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional
import numpy as np
from scipy import stats

# Paths
SCRIPT_DIR = Path(__file__).parent
PROJECT_ROOT = SCRIPT_DIR.parent.parent
RESULTS_DIR = PROJECT_ROOT / "results" / "outputs"
FREIRE_CACHE = RESULTS_DIR / "freire_GCpsr.txt"

# Constants
MSP_PERIOD_CUT_MS = 30.0  # P < 30 ms defines MSP

# Regex for parsing numeric values with uncertainties
_NUM_RE = re.compile(r'^[+-]?\d+\.?\d*(?:[eE][+-]?\d+)?$')


def _parse_numeric(val: str) -> Optional[float]:
    """Parse a numeric value, handling uncertainties in parentheses."""
    if not val or val == '*' or val == 'i':
        return None
    # Remove uncertainty in parentheses: e.g., "24.599(2)" -> "24.599"
    val = re.sub(r'\([^)]*\)', '', val)
    # Remove leading/trailing whitespace
    val = val.strip()
    # Handle < or > prefixes
    val = val.lstrip('<>')
    if _NUM_RE.match(val):
        return float(val)
    return None


def parse_freire_catalog(text: str) -> list[dict]:
    """
    Parse the Freire GCpsr catalog, extracting both spin and binary parameters.
    
    Returns a list of dicts with keys:
    - cluster, name, offset_arcmin, P_ms, Pdot_1e20, DM
    - Pb_days, x_s, e, m2_msun (binary parameters, None if isolated)
    - is_binary: bool
    """
    rows = []
    current_cluster = None
    
    for line in text.splitlines():
        line = line.strip()
        
        # Skip empty lines and comments
        if not line or line.startswith('#'):
            continue
        
        # Check if this is a cluster header (no tabs, contains NGC or name)
        # Cluster headers are lines that don't start with 'J' or 'B' and contain cluster names
        if not line.startswith('J') and not line.startswith('B') and not line.startswith('*'):
            # This might be a cluster header
            if '(' in line or 'NGC' in line or 'M' in line.split()[0] if line.split() else False:
                current_cluster = line.strip()
                continue
            # Check if it's a simple cluster name
            parts = line.split()
            if len(parts) <= 3 and not any(c.isdigit() and '.' in line for c in line):
                current_cluster = line.strip()
                continue
        
        if current_cluster is None:
            continue
        
        # Parse pulsar line
        # Format: Name Offset Period Pdot DM Pb x e m2
        # Fields may be separated by tabs or multiple spaces
        parts = line.split()
        if len(parts) < 5:
            continue
        
        # First field is name (starts with J or B, or * for unmeasured offset)
        name = parts[0]
        if name == '*':
            # Offset is unmeasured, name is in next field
            if len(parts) < 6:
                continue
            offset = None
            name = parts[0]  # Actually this is tricky - need to handle differently
            # Skip lines where we can't parse properly
            continue
        
        if not (name.startswith('J') or name.startswith('B')):
            continue
        
        # Parse remaining fields
        try:
            idx = 1
            
            # Offset (arcmin) - may be * or numeric
            offset_str = parts[idx] if idx < len(parts) else '*'
            offset = _parse_numeric(offset_str)
            idx += 1
            
            # Period (ms)
            P_ms = _parse_numeric(parts[idx]) if idx < len(parts) else None
            idx += 1
            
            # Pdot (10^-20)
            Pdot_1e20 = _parse_numeric(parts[idx]) if idx < len(parts) else None
            idx += 1
            
            # DM
            DM = _parse_numeric(parts[idx]) if idx < len(parts) else None
            idx += 1
            
            # Binary parameters (optional)
            Pb_days = _parse_numeric(parts[idx]) if idx < len(parts) else None
            idx += 1
            x_s = _parse_numeric(parts[idx]) if idx < len(parts) else None
            idx += 1
            e = _parse_numeric(parts[idx]) if idx < len(parts) else None
            idx += 1
            m2_msun = _parse_numeric(parts[idx]) if idx < len(parts) else None
            
            # Determine if binary
            # 'i' means isolated, '*' means unmeasured
            is_binary = Pb_days is not None and Pb_days > 0
            
            rows.append({
                'cluster': current_cluster,
                'name': name,
                'offset_arcmin': offset,
                'P_ms': P_ms,
                'Pdot_1e20': Pdot_1e20,
                'DM': DM,
                'Pb_days': Pb_days,
                'x_s': x_s,
                'e': e,
                'm2_msun': m2_msun,
                'is_binary': is_binary,
            })
            
        except (IndexError, ValueError):
            continue
    
    return rows


def compute_derived_quantities(rows: list[dict]) -> list[dict]:
    """Add derived quantities to each pulsar row."""
    for r in rows:
        P_ms = r.get('P_ms')
        Pdot_1e20 = r.get('Pdot_1e20')
        
        if P_ms is not None and P_ms > 0:
            r['P_s'] = P_ms / 1000.0
            r['logP'] = np.log10(r['P_s'])
        else:
            r['P_s'] = None
            r['logP'] = None
        
        if Pdot_1e20 is not None:
            r['Pdot_sps'] = Pdot_1e20 * 1e-20
            if Pdot_1e20 != 0:
                r['logPdot_abs'] = np.log10(abs(Pdot_1e20 * 1e-20))
            else:
                r['logPdot_abs'] = None
            r['Pdot_sign'] = 'positive' if Pdot_1e20 > 0 else 'negative' if Pdot_1e20 < 0 else 'zero'
        else:
            r['Pdot_sps'] = None
            r['logPdot_abs'] = None
            r['Pdot_sign'] = None
        
        # Characteristic age proxy (if Pdot > 0)
        if r['P_s'] is not None and r['Pdot_sps'] is not None and r['Pdot_sps'] > 0:
            r['tau_c_yr'] = r['P_s'] / (2 * r['Pdot_sps']) / (365.25 * 24 * 3600)
            r['log_tau_c'] = np.log10(r['tau_c_yr'])
        else:
            r['tau_c_yr'] = None
            r['log_tau_c'] = None
        
        # Surface B-field proxy (if Pdot > 0)
        if r['P_s'] is not None and r['Pdot_sps'] is not None and r['Pdot_sps'] > 0:
            r['B_proxy'] = 3.2e19 * np.sqrt(r['P_s'] * r['Pdot_sps'])
            r['log_B_proxy'] = np.log10(r['B_proxy'])
        else:
            r['B_proxy'] = None
            r['log_B_proxy'] = None
    
    return rows


def is_msp(r: dict) -> bool:
    """Check if a pulsar is an MSP (P < 30 ms)."""
    P_ms = r.get('P_ms')
    return P_ms is not None and P_ms < MSP_PERIOD_CUT_MS


def binary_vs_isolated_comparison(rows: list[dict]) -> dict:
    """
    Compare binary vs isolated MSPs within GCs.
    
    If the low |Ṗ| effect is due to cluster acceleration, both should be affected equally.
    If there's a systematic difference, it could indicate a population or selection effect.
    """
    # Filter to MSPs with measured Pdot
    msps = [r for r in rows if is_msp(r) and r.get('logPdot_abs') is not None]
    
    binary_msps = [r for r in msps if r['is_binary']]
    isolated_msps = [r for r in msps if not r['is_binary']]
    
    binary_logPdot = np.array([r['logPdot_abs'] for r in binary_msps])
    isolated_logPdot = np.array([r['logPdot_abs'] for r in isolated_msps])
    
    result = {
        'n_binary': len(binary_msps),
        'n_isolated': len(isolated_msps),
        'binary_mean_logPdot': float(np.mean(binary_logPdot)) if len(binary_logPdot) > 0 else None,
        'isolated_mean_logPdot': float(np.mean(isolated_logPdot)) if len(isolated_logPdot) > 0 else None,
        'binary_std_logPdot': float(np.std(binary_logPdot)) if len(binary_logPdot) > 0 else None,
        'isolated_std_logPdot': float(np.std(isolated_logPdot)) if len(isolated_logPdot) > 0 else None,
    }
    
    if len(binary_logPdot) >= 3 and len(isolated_logPdot) >= 3:
        # Welch t-test
        t_stat, t_p = stats.ttest_ind(binary_logPdot, isolated_logPdot, equal_var=False)
        # Mann-Whitney U test
        mw_u, mw_p = stats.mannwhitneyu(binary_logPdot, isolated_logPdot, alternative='two-sided')
        
        result['diff_dex'] = result['binary_mean_logPdot'] - result['isolated_mean_logPdot']
        result['t_stat'] = float(t_stat)
        result['t_p'] = float(t_p)
        result['mw_u'] = float(mw_u)
        result['mw_p'] = float(mw_p)
    
    # Sign distribution comparison
    binary_signs = [r['Pdot_sign'] for r in binary_msps if r['Pdot_sign']]
    isolated_signs = [r['Pdot_sign'] for r in isolated_msps if r['Pdot_sign']]
    
    result['binary_sign_dist'] = {
        'positive': sum(1 for s in binary_signs if s == 'positive'),
        'negative': sum(1 for s in binary_signs if s == 'negative'),
    }
    result['isolated_sign_dist'] = {
        'positive': sum(1 for s in isolated_signs if s == 'positive'),
        'negative': sum(1 for s in isolated_signs if s == 'negative'),
    }
    
    # Fraction negative
    if len(binary_signs) > 0:
        result['binary_frac_negative'] = result['binary_sign_dist']['negative'] / len(binary_signs)
    if len(isolated_signs) > 0:
        result['isolated_frac_negative'] = result['isolated_sign_dist']['negative'] / len(isolated_signs)
    
    # ---------------------------------------------------------
    # SPATIAL STRATIFICATION (Reviewer Response)
    # ---------------------------------------------------------
    # Compare Binary vs Isolated within radial bins to control for mass segregation
    
    # Filter to those with valid offsets
    msps_spatial = [r for r in msps if r.get('logPdot_abs') is not None and r['offset_arcmin'] is not None]
    
    if len(msps_spatial) >= 10:
        offsets = [r['offset_arcmin'] for r in msps_spatial]
        median_offset = float(np.median(offsets))
        
        inner_msps = [r for r in msps_spatial if r['offset_arcmin'] <= median_offset]
        outer_msps = [r for r in msps_spatial if r['offset_arcmin'] > median_offset]
        
        result['spatial_control'] = {
            'median_offset_arcmin': median_offset,
            'inner': _compare_subsets(inner_msps),
            'outer': _compare_subsets(outer_msps)
        }
    
    return result


def _compare_subsets(subset: list[dict]) -> dict:
    """Helper to compare binary vs isolated in a subset."""
    binaries = [r for r in subset if r['is_binary']]
    isolated = [r for r in subset if not r['is_binary']]
    
    b_vals = [r['logPdot_abs'] for r in binaries]
    i_vals = [r['logPdot_abs'] for r in isolated]
    
    res = {
        'n_binary': len(binaries),
        'n_isolated': len(isolated),
        'mean_binary': float(np.mean(b_vals)) if b_vals else None,
        'mean_isolated': float(np.mean(i_vals)) if i_vals else None,
    }
    
    if len(b_vals) >= 2 and len(i_vals) >= 2:
        res['diff'] = res['mean_binary'] - res['mean_isolated']
        _, p = stats.ttest_ind(b_vals, i_vals, equal_var=False)
        res['p_value'] = float(p)
    else:
        res['diff'] = None
        res['p_value'] = None
        
    return res


def orbital_parameter_analysis(rows: list[dict]) -> dict:
    """
    Analyze orbital parameters of GC binaries.
    
    Look for correlations between orbital parameters and spin-down.
    """
    # Filter to binary MSPs with measured Pdot and orbital parameters
    binaries = [r for r in rows if is_msp(r) and r['is_binary'] and r.get('logPdot_abs') is not None]
    
    result = {
        'n_binaries_with_pdot': len(binaries),
    }
    
    if len(binaries) < 5:
        result['note'] = 'Insufficient binaries with measured Pdot for orbital analysis'
        return result
    
    # Extract arrays
    logPdot = np.array([r['logPdot_abs'] for r in binaries])
    Pb = np.array([r['Pb_days'] for r in binaries if r['Pb_days'] is not None])
    e = np.array([r['e'] for r in binaries if r['e'] is not None])
    m2 = np.array([r['m2_msun'] for r in binaries if r['m2_msun'] is not None])
    
    # Correlations
    if len(Pb) >= 5:
        logPb = np.log10(Pb)
        # Match indices
        logPdot_Pb = np.array([r['logPdot_abs'] for r in binaries if r['Pb_days'] is not None])
        if len(logPdot_Pb) == len(logPb):
            r_Pb, p_Pb = stats.pearsonr(logPb, logPdot_Pb)
            result['Pb_correlation'] = {
                'r': float(r_Pb),
                'p': float(p_Pb),
                'n': len(logPb),
            }
    
    if len(e) >= 5:
        loge = np.log10(np.maximum(e, 1e-10))  # Avoid log(0)
        logPdot_e = np.array([r['logPdot_abs'] for r in binaries if r['e'] is not None])
        if len(logPdot_e) == len(loge):
            r_e, p_e = stats.pearsonr(loge, logPdot_e)
            result['e_correlation'] = {
                'r': float(r_e),
                'p': float(p_e),
                'n': len(loge),
            }
    
    if len(m2) >= 5:
        logm2 = np.log10(np.maximum(m2, 0.001))
        logPdot_m2 = np.array([r['logPdot_abs'] for r in binaries if r['m2_msun'] is not None])
        if len(logPdot_m2) == len(logm2):
            r_m2, p_m2 = stats.pearsonr(logm2, logPdot_m2)
            result['m2_correlation'] = {
                'r': float(r_m2),
                'p': float(p_m2),
                'n': len(logm2),
            }
    
    # Summary statistics
    result['Pb_stats'] = {
        'median_days': float(np.median(Pb)) if len(Pb) > 0 else None,
        'min_days': float(np.min(Pb)) if len(Pb) > 0 else None,
        'max_days': float(np.max(Pb)) if len(Pb) > 0 else None,
    }
    
    return result


def cluster_summary(rows: list[dict]) -> dict:
    """Generate per-cluster summary of binary vs isolated MSPs."""
    clusters = {}
    
    for r in rows:
        if not is_msp(r):
            continue
        
        cluster = r['cluster']
        if cluster not in clusters:
            clusters[cluster] = {
                'n_binary': 0,
                'n_isolated': 0,
                'n_with_pdot': 0,
                'binary_with_pdot': 0,
                'isolated_with_pdot': 0,
            }
        
        if r['is_binary']:
            clusters[cluster]['n_binary'] += 1
            if r.get('logPdot_abs') is not None:
                clusters[cluster]['binary_with_pdot'] += 1
        else:
            clusters[cluster]['n_isolated'] += 1
            if r.get('logPdot_abs') is not None:
                clusters[cluster]['isolated_with_pdot'] += 1
        
        if r.get('logPdot_abs') is not None:
            clusters[cluster]['n_with_pdot'] += 1
    
    return clusters


def main():
    # Load Freire catalog
    if not FREIRE_CACHE.exists():
        raise FileNotFoundError(f"Freire catalog not found at {FREIRE_CACHE}")
    
    freire_text = FREIRE_CACHE.read_text()
    freire_sha256 = hashlib.sha256(freire_text.encode()).hexdigest()
    
    # Parse catalog
    rows = parse_freire_catalog(freire_text)
    rows = compute_derived_quantities(rows)
    
    print(f"Parsed {len(rows)} pulsars from Freire catalog")
    
    # Filter to MSPs
    msps = [r for r in rows if is_msp(r)]
    print(f"MSPs (P < 30 ms): {len(msps)}")
    
    # Binary vs isolated comparison
    binary_isolated = binary_vs_isolated_comparison(rows)
    print(f"\nBinary vs Isolated MSPs in GCs:")
    print(f"  Binary MSPs: {binary_isolated['n_binary']}")
    print(f"  Isolated MSPs: {binary_isolated['n_isolated']}")
    if binary_isolated.get('diff_dex') is not None:
        print(f"  Difference (binary - isolated): {binary_isolated['diff_dex']:.3f} dex")
        print(f"  t-test p-value: {binary_isolated['t_p']:.4g}")
    
    # Orbital parameter analysis
    orbital = orbital_parameter_analysis(rows)
    print(f"\nOrbital Parameter Analysis:")
    print(f"  Binaries with measured Pdot: {orbital['n_binaries_with_pdot']}")
    
    # Cluster summary
    cluster_stats = cluster_summary(rows)
    
    # Build output
    output = {
        'meta': {
            'timestamp_utc': datetime.now(timezone.utc).isoformat(),
            'source': {
                'freire_gcpsr': {
                    'path': str(FREIRE_CACHE),
                    'sha256': freire_sha256,
                }
            },
            'selection': {
                'msp_cut': f'P < {MSP_PERIOD_CUT_MS} ms',
            },
            'counts': {
                'total_pulsars': len(rows),
                'total_msps': len(msps),
                'binary_msps': binary_isolated['n_binary'],
                'isolated_msps': binary_isolated['n_isolated'],
            }
        },
        'binary_vs_isolated': binary_isolated,
        'orbital_analysis': orbital,
        'cluster_summary': cluster_stats,
    }
    
    # Write outputs
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    
    json_path = RESULTS_DIR / "step_5_11_binary_pulsar_analysis.json"
    with open(json_path, 'w') as f:
        json.dump(output, f, indent=2)
    print(f"\nWrote: {json_path}")
    
    # Write markdown summary
    md_path = RESULTS_DIR / "step_5_11_binary_pulsar_analysis.md"
    md_lines = [
        "# Binary Pulsar Analysis (Freire GCpsr)",
        f"**Source:** {FREIRE_CACHE}",
        f"**SHA256:** `{freire_sha256[:16]}...`",
        "",
        "## Sample Sizes",
        f"- **Total GC pulsars parsed:** {len(rows)}",
        f"- **GC MSPs (P < 30 ms):** {len(msps)}",
        f"- **Binary MSPs:** {binary_isolated['n_binary']}",
        f"- **Isolated MSPs:** {binary_isolated['n_isolated']}",
        "",
        "## Binary vs Isolated Comparison",
        "",
        "| Metric | Binary MSPs | Isolated MSPs |",
        "| --- | --- | --- |",
        f"| N | {binary_isolated['n_binary']} | {binary_isolated['n_isolated']} |",
    ]
    
    if binary_isolated.get('binary_mean_logPdot') is not None:
        md_lines.append(f"| Mean log|Ṗ| | {binary_isolated['binary_mean_logPdot']:.3f} | {binary_isolated['isolated_mean_logPdot']:.3f} |")
        md_lines.append(f"| Std log|Ṗ| | {binary_isolated['binary_std_logPdot']:.3f} | {binary_isolated['isolated_std_logPdot']:.3f} |")
    
    if binary_isolated.get('binary_frac_negative') is not None:
        md_lines.append(f"| Fraction negative Ṗ | {binary_isolated['binary_frac_negative']:.1%} | {binary_isolated.get('isolated_frac_negative', 0):.1%} |")
    
    md_lines.append("")
    
    if binary_isolated.get('diff_dex') is not None:
        md_lines.extend([
            "### Statistical Tests",
            f"- **Difference (binary - isolated):** {binary_isolated['diff_dex']:.3f} dex",
            f"- **Welch t-test p:** {binary_isolated['t_p']:.4g}",
            f"- **Mann-Whitney p:** {binary_isolated['mw_p']:.4g}",
            "",
        ])
    
    md_lines.extend([
        "## Interpretation",
        "",
        "If the low |Ṗ| effect in GC pulsars were purely due to cluster acceleration, we would expect:",
        "- Binary and isolated MSPs to be affected equally (same line-of-sight acceleration)",
        "- No significant difference in log|Ṗ| between the two populations",
        "",
        "Any significant difference would suggest:",
        "- Population/selection effects (binary MSPs may have different intrinsic properties)",
        "- Or a TEP-like effect that couples differently to binary vs isolated systems",
        "",
    ])
    
    if binary_isolated.get('diff_dex') is not None:
        if abs(binary_isolated['t_p']) < 0.05:
            md_lines.append(f"**Result:** Significant difference detected (p = {binary_isolated['t_p']:.4g}). This warrants further investigation.")
        else:
            md_lines.append(f"**Result:** No significant difference (p = {binary_isolated['t_p']:.2f}). Binary and isolated MSPs show similar log|Ṗ| distributions, consistent with cluster acceleration affecting both equally.")
    
    with open(md_path, 'w') as f:
        f.write('\n'.join(md_lines))
    print(f"Wrote: {md_path}")


if __name__ == "__main__":
    main()
