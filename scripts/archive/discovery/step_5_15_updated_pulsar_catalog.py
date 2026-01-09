#!/usr/bin/env python3
"""
Step 5.15: Updated Pulsar Catalog Analysis

Parses the updated Freire GCpsr catalog (2025-01-16) and compares to 
the previous analysis to identify new pulsars with measured P-dot values.
"""

import re
import json
import numpy as np
from pathlib import Path
from collections import defaultdict

DATA_DIR = Path("/Users/matthewsmawfield/www/TEP-COS/data")
OUTPUT_DIR = Path("/Users/matthewsmawfield/www/TEP-COS/results/outputs")

def parse_freire_catalog(filepath):
    """Parse the Freire GCpsr catalog format."""
    pulsars = []
    current_cluster = None
    
    with open(filepath, 'r') as f:
        lines = f.readlines()
    
    for line in lines:
        line = line.strip()
        
        # Skip empty lines and comments
        if not line or line.startswith('#'):
            continue
        
        # Check for cluster header (lines without tabs at start that aren't pulsar names)
        if not line.startswith('J') and not line.startswith('B') and not line.startswith('\t'):
            # This might be a cluster name
            if '(' in line or 'NGC' in line or 'Terzan' in line or 'Omega' in line or 'M' in line:
                current_cluster = line.split('(')[0].strip() if '(' in line else line.strip()
                continue
        
        # Parse pulsar line
        if line.startswith('J') or line.startswith('B'):
            parts = line.split()
            if len(parts) < 2:
                continue
            
            name = parts[0]
            
            # Try to extract period and P-dot
            # Format varies but typically: Name, offset, Period, dP/dt, DM, ...
            period = None
            pdot = None
            offset = None
            
            for i, part in enumerate(parts[1:], 1):
                # Skip asterisks (unmeasured values)
                if part == '*':
                    continue
                
                # Try to identify period (typically 1-100 ms for MSPs)
                try:
                    val = float(part.replace('(', '').split(')')[0])
                    if period is None and 0.5 < val < 1000:  # Likely period in ms
                        # Check if previous was offset
                        if i > 1 and offset is None:
                            try:
                                offset = float(parts[i-1])
                            except:
                                pass
                        period = val
                    elif period is not None and pdot is None:
                        # This might be P-dot (can be positive or negative)
                        if abs(val) < 1000 and val != 0:  # P-dot in units of 10^-20
                            pdot = val
                            break
                except:
                    continue
            
            if period is not None:
                pulsar = {
                    'name': name,
                    'cluster': current_cluster,
                    'period_ms': period,
                    'pdot': pdot,
                    'offset_arcmin': offset,
                    'has_pdot': pdot is not None and str(pdot) != '*'
                }
                pulsars.append(pulsar)
    
    return pulsars

def analyze_catalog(pulsars):
    """Analyze the parsed catalog."""
    results = {
        'total_pulsars': len(pulsars),
        'with_pdot': sum(1 for p in pulsars if p['has_pdot']),
        'msps': sum(1 for p in pulsars if p['period_ms'] and p['period_ms'] < 30),
        'msps_with_pdot': sum(1 for p in pulsars if p['period_ms'] and p['period_ms'] < 30 and p['has_pdot']),
        'clusters': defaultdict(lambda: {'total': 0, 'with_pdot': 0, 'msps': 0})
    }
    
    for p in pulsars:
        cluster = p['cluster'] or 'Unknown'
        results['clusters'][cluster]['total'] += 1
        if p['has_pdot']:
            results['clusters'][cluster]['with_pdot'] += 1
        if p['period_ms'] and p['period_ms'] < 30:
            results['clusters'][cluster]['msps'] += 1
    
    # Convert defaultdict to regular dict
    results['clusters'] = dict(results['clusters'])
    
    return results

def compute_pdot_statistics(pulsars):
    """Compute P-dot statistics for MSPs."""
    msps_with_pdot = [p for p in pulsars if p['period_ms'] and p['period_ms'] < 30 and p['has_pdot'] and p['pdot'] is not None]
    
    if not msps_with_pdot:
        return None
    
    pdots = [p['pdot'] for p in msps_with_pdot]
    abs_pdots = [abs(p) for p in pdots]
    log_abs_pdots = [np.log10(abs(p)) for p in pdots if abs(p) > 0]
    
    # Count positive vs negative
    n_positive = sum(1 for p in pdots if p > 0)
    n_negative = sum(1 for p in pdots if p < 0)
    
    stats = {
        'n_msps_with_pdot': len(msps_with_pdot),
        'n_positive_pdot': n_positive,
        'n_negative_pdot': n_negative,
        'frac_negative': n_negative / len(pdots) if pdots else 0,
        'mean_log_abs_pdot': np.mean(log_abs_pdots) if log_abs_pdots else None,
        'std_log_abs_pdot': np.std(log_abs_pdots) if log_abs_pdots else None,
        'median_abs_pdot': np.median(abs_pdots) if abs_pdots else None
    }
    
    return stats

def main():
    print("="*60)
    print("UPDATED FREIRE GCPSR CATALOG ANALYSIS (2025-01-16)")
    print("="*60)
    
    catalog_path = DATA_DIR / "freire_gcpsr_2025.txt"
    
    if not catalog_path.exists():
        print(f"Error: Catalog not found at {catalog_path}")
        return
    
    # Parse catalog
    pulsars = parse_freire_catalog(catalog_path)
    print(f"\nParsed {len(pulsars)} pulsars from catalog")
    
    # Analyze
    results = analyze_catalog(pulsars)
    
    print(f"\n--- CATALOG SUMMARY ---")
    print(f"Total pulsars: {results['total_pulsars']}")
    print(f"With measured P-dot: {results['with_pdot']}")
    print(f"MSPs (P < 30 ms): {results['msps']}")
    print(f"MSPs with P-dot: {results['msps_with_pdot']}")
    
    print(f"\n--- TOP CLUSTERS BY MSP COUNT ---")
    cluster_counts = [(c, d['msps']) for c, d in results['clusters'].items()]
    cluster_counts.sort(key=lambda x: x[1], reverse=True)
    for cluster, count in cluster_counts[:15]:
        pdot_count = results['clusters'][cluster]['with_pdot']
        print(f"  {cluster}: {count} MSPs ({pdot_count} with P-dot)")
    
    # P-dot statistics
    stats = compute_pdot_statistics(pulsars)
    if stats:
        print(f"\n--- P-DOT STATISTICS (MSPs) ---")
        print(f"MSPs with P-dot: {stats['n_msps_with_pdot']}")
        print(f"Positive P-dot: {stats['n_positive_pdot']} ({100*stats['n_positive_pdot']/stats['n_msps_with_pdot']:.1f}%)")
        print(f"Negative P-dot: {stats['n_negative_pdot']} ({100*stats['frac_negative']:.1f}%)")
        if stats['mean_log_abs_pdot']:
            print(f"Mean log|P-dot|: {stats['mean_log_abs_pdot']:.2f}")
            print(f"Std log|P-dot|: {stats['std_log_abs_pdot']:.2f}")
    
    # Compare to previous (manuscript values)
    print(f"\n--- COMPARISON TO MANUSCRIPT ---")
    print(f"Previous GC MSPs (manuscript): 181")
    print(f"Updated GC MSPs (this catalog): {results['msps_with_pdot']}")
    print(f"Change: +{results['msps_with_pdot'] - 181} MSPs with P-dot")
    
    # Save results
    output = {
        'catalog_date': '2025-01-16',
        'source': 'https://www3.mpifr-bonn.mpg.de/staff/pfreire/GCpsr.txt',
        'summary': results,
        'pdot_statistics': stats,
        'comparison': {
            'previous_gc_msps': 181,
            'updated_gc_msps': results['msps_with_pdot'],
            'new_pulsars': results['msps_with_pdot'] - 181
        },
        'pulsars': pulsars
    }
    
    output_path = OUTPUT_DIR / "updated_pulsar_catalog_analysis.json"
    
    # Convert numpy types for JSON
    def convert(obj):
        if isinstance(obj, (np.integer, np.floating)):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        return obj
    
    with open(output_path, 'w') as f:
        json.dump(output, f, indent=2, default=convert)
    
    print(f"\nResults saved to: {output_path}")
    
    return output

if __name__ == "__main__":
    main()
