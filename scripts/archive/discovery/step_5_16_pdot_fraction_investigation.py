#!/usr/bin/env python3
"""
Step 5.16: Investigate P-dot Fraction Change

The manuscript reports 45% negative P-dot in GC pulsars, but the updated
catalog shows only 22%. This script investigates the discrepancy.

Possible explanations:
1. Different filtering criteria (P < 30ms vs different cutoff)
2. New pulsars discovered in outer regions (less acceleration)
3. Parsing differences in the catalog format
4. Updated measurements changing sign for some pulsars
"""

import re
import numpy as np
from pathlib import Path
from collections import defaultdict

DATA_DIR = Path("/Users/matthewsmawfield/www/TEP-COS/data")

def parse_freire_detailed(filepath):
    """Parse Freire catalog with detailed P-dot extraction."""
    pulsars = []
    current_cluster = None
    
    with open(filepath, 'r') as f:
        lines = f.readlines()
    
    for line in lines:
        original_line = line
        line = line.strip()
        
        if not line or line.startswith('#'):
            continue
        
        # Cluster header detection
        if not line.startswith('J') and not line.startswith('B'):
            if any(x in line for x in ['NGC', 'Terzan', 'M', 'Omega', 'Tuc']):
                current_cluster = line.split('(')[0].strip() if '(' in line else line.strip()
                continue
        
        # Parse pulsar line
        if line.startswith('J') or line.startswith('B'):
            parts = line.split()
            if len(parts) < 3:
                continue
            
            name = parts[0]
            
            # Look for P-dot value - it's typically the 4th numeric column
            # Format: Name, offset, Period, dP/dt, DM, ...
            # P-dot is in units of 10^-20 and can have parenthetical errors
            
            pdot = None
            pdot_raw = None
            period = None
            
            for i, part in enumerate(parts[1:], 1):
                if part == '*':
                    continue
                    
                # Clean the part
                clean = re.sub(r'\([^)]*\)', '', part)  # Remove parenthetical errors
                
                try:
                    val = float(clean)
                    
                    # First numeric is offset (arcmin), second is period, third is P-dot
                    if period is None and 0.5 < val < 10000:  # Period in ms
                        period = val
                    elif period is not None and pdot is None:
                        # This should be P-dot
                        pdot = val
                        pdot_raw = part
                        break
                except:
                    continue
            
            if period is not None:
                pulsar = {
                    'name': name,
                    'cluster': current_cluster,
                    'period_ms': period,
                    'pdot': pdot,
                    'pdot_raw': pdot_raw,
                    'is_msp': period < 30,
                    'has_measured_pdot': pdot is not None
                }
                pulsars.append(pulsar)
    
    return pulsars

def analyze_pdot_distribution(pulsars):
    """Analyze P-dot distribution in detail."""
    
    # Filter to MSPs with measured P-dot
    msps = [p for p in pulsars if p['is_msp'] and p['has_measured_pdot'] and p['pdot'] is not None]
    
    print(f"\n{'='*60}")
    print("P-DOT DISTRIBUTION ANALYSIS")
    print(f"{'='*60}")
    
    print(f"\nTotal MSPs (P < 30ms) with measured P-dot: {len(msps)}")
    
    # Count positive/negative
    positive = [p for p in msps if p['pdot'] > 0]
    negative = [p for p in msps if p['pdot'] < 0]
    zero = [p for p in msps if p['pdot'] == 0]
    
    print(f"\nP-dot sign distribution:")
    print(f"  Positive (spin-down): {len(positive)} ({100*len(positive)/len(msps):.1f}%)")
    print(f"  Negative (spin-up):   {len(negative)} ({100*len(negative)/len(msps):.1f}%)")
    print(f"  Zero:                 {len(zero)}")
    
    # By cluster
    print(f"\n{'='*60}")
    print("P-DOT SIGN BY CLUSTER (MSPs only)")
    print(f"{'='*60}")
    
    cluster_stats = defaultdict(lambda: {'pos': 0, 'neg': 0, 'total': 0})
    
    for p in msps:
        cluster = p['cluster'] or 'Unknown'
        cluster_stats[cluster]['total'] += 1
        if p['pdot'] > 0:
            cluster_stats[cluster]['pos'] += 1
        elif p['pdot'] < 0:
            cluster_stats[cluster]['neg'] += 1
    
    # Sort by total count
    sorted_clusters = sorted(cluster_stats.items(), key=lambda x: x[1]['total'], reverse=True)
    
    print(f"\n{'Cluster':<20} {'Total':>6} {'Pos':>6} {'Neg':>6} {'%Neg':>8}")
    print("-"*50)
    
    for cluster, stats in sorted_clusters[:20]:
        if stats['total'] > 0:
            frac_neg = 100 * stats['neg'] / stats['total']
            print(f"{cluster:<20} {stats['total']:>6} {stats['pos']:>6} {stats['neg']:>6} {frac_neg:>7.1f}%")
    
    # High negative-fraction clusters (like in manuscript)
    high_neg_clusters = [(c, s) for c, s in sorted_clusters if s['total'] >= 5 and s['neg']/s['total'] > 0.3]
    
    print(f"\n{'='*60}")
    print("CLUSTERS WITH >30% NEGATIVE P-DOT (N≥5)")
    print(f"{'='*60}")
    
    if high_neg_clusters:
        for cluster, stats in high_neg_clusters:
            frac = 100 * stats['neg'] / stats['total']
            print(f"  {cluster}: {stats['neg']}/{stats['total']} = {frac:.1f}% negative")
    else:
        print("  None found!")
    
    # Sample some negative P-dot pulsars
    print(f"\n{'='*60}")
    print("SAMPLE NEGATIVE P-DOT MSPs")
    print(f"{'='*60}")
    
    for p in negative[:10]:
        print(f"  {p['name']:<16} P={p['period_ms']:.2f}ms  Pdot={p['pdot']:.4f}  [{p['cluster']}]")
    
    # Check if manuscript used different criteria
    print(f"\n{'='*60}")
    print("COMPARISON WITH MANUSCRIPT VALUES")
    print(f"{'='*60}")
    
    print(f"\nManuscript claims:")
    print(f"  - Field MSPs: 2% negative P-dot")
    print(f"  - GC MSPs: 45% negative P-dot")
    print(f"  - Sample size: 181 GC MSPs")
    
    print(f"\nCurrent catalog:")
    print(f"  - GC MSPs with P-dot: {len(msps)}")
    print(f"  - Negative fraction: {100*len(negative)/len(msps):.1f}%")
    
    # The discrepancy suggests different data sources or filtering
    print(f"\n{'='*60}")
    print("DIAGNOSIS")
    print(f"{'='*60}")
    
    print("""
The 45% negative P-dot fraction in the manuscript likely came from:
1. The ATNF catalog filtering, not Freire catalog directly
2. Different sample selection (only pulsars with well-measured P-dot)
3. The manuscript sample included pulsars with large acceleration terms

The Freire catalog may include newer pulsars discovered in outer regions
with smaller acceleration contamination, reducing the negative fraction.

RECOMMENDATION: 
- Cross-check with ATNF catalog for consistency
- The key TEP result (0.13 dex residual) is based on population controls,
  not raw negative fraction
- Update manuscript to reflect new sample size but note methodology
""")
    
    return {
        'total_msps': len(msps),
        'n_positive': len(positive),
        'n_negative': len(negative),
        'frac_negative': len(negative) / len(msps) if msps else 0,
        'cluster_stats': dict(cluster_stats)
    }

def main():
    catalog_path = DATA_DIR / "freire_gcpsr_2025.txt"
    
    if not catalog_path.exists():
        print(f"Error: Catalog not found at {catalog_path}")
        return
    
    pulsars = parse_freire_detailed(catalog_path)
    print(f"Parsed {len(pulsars)} pulsars")
    
    results = analyze_pdot_distribution(pulsars)
    
    return results

if __name__ == "__main__":
    main()
