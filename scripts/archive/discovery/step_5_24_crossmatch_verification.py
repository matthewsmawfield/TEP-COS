#!/usr/bin/env python3
"""
Step 5.24: Cross-Match Verification

Rigorous cross-match of Freire GC pulsars against ATNF to determine
how many ACTUALLY have measured P-dot values.
"""

import re
from pathlib import Path
from collections import defaultdict

DATA_DIR = Path("/Users/matthewsmawfield/www/TEP-COS/data")
OUTPUT_DIR = Path("/Users/matthewsmawfield/www/TEP-COS/results/outputs")

def parse_freire_catalog():
    """Parse Freire catalog and extract pulsar names and P-dot status."""
    
    filepath = DATA_DIR / "freire_gcpsr_2025.txt"
    
    pulsars = []
    current_cluster = None
    
    with open(filepath, 'r', errors='ignore') as f:
        lines = f.readlines()
    
    for line in lines:
        line = line.strip()
        if not line:
            continue
        
        # Check for cluster header
        if not line.startswith('J') and not line.startswith('B'):
            # Could be cluster name
            if any(x in line for x in ['NGC', 'Terzan', 'M ', 'M1', 'M2', 'M3', 'M4', 'M5', 
                                        'M6', 'M7', 'M8', 'M9', 'Omega', 'Tuc', 'Pal', 'Liller',
                                        '47 ', 'IC ', 'Djorg']):
                current_cluster = line.split('(')[0].strip() if '(' in line else line.strip()
                continue
        
        # Parse pulsar line
        if line.startswith('J') or line.startswith('B'):
            parts = line.split()
            if len(parts) < 2:
                continue
            
            name = parts[0]
            
            # Check if P-dot is present (look for scientific notation or explicit value)
            # P-dot is typically in column 3 or 4 (after offset and period)
            has_pdot = False
            pdot_value = None
            
            # Look for P-dot pattern: +/- number with parentheses or scientific notation
            for i, part in enumerate(parts[1:], 1):
                # Skip if it's clearly period (ms range) or DM
                if part == '*':
                    continue
                    
                # P-dot indicators: starts with +/-, has parentheses, very small number
                if re.match(r'^[+-]?\d+\.?\d*\([^)]+\)', part):
                    # This looks like a measured value with uncertainty
                    has_pdot = True
                    pdot_value = part
                    break
                elif re.match(r'^[+-]?\d+\.\d+e[+-]?\d+', part, re.IGNORECASE):
                    has_pdot = True
                    pdot_value = part
                    break
            
            # Determine if MSP (period < 30ms)
            period_ms = None
            for part in parts[1:]:
                if part == '*':
                    continue
                try:
                    val = float(re.sub(r'\([^)]*\)', '', part))
                    if 0.5 < val < 1000:  # Likely period in ms
                        period_ms = val
                        break
                except:
                    continue
            
            is_msp = period_ms is not None and period_ms < 30
            
            pulsars.append({
                'name': name,
                'cluster': current_cluster,
                'period_ms': period_ms,
                'is_msp': is_msp,
                'has_pdot_freire': has_pdot,
                'pdot_freire': pdot_value
            })
    
    return pulsars

def parse_atnf_html():
    """Parse ATNF HTML and extract pulsars with P-dot."""
    
    filepath = DATA_DIR / "atnf_full_catalog.html"
    
    if not filepath.exists():
        print(f"ATNF file not found: {filepath}")
        return {}
    
    with open(filepath, 'r', errors='ignore') as f:
        content = f.read()
    
    # Remove HTML tags
    text = re.sub(r'<[^>]+>', ' ', content)
    
    # Find pulsar entries with P-dot
    # Pattern: Jxxxx+xxxx or Bxxxx+xx followed by period and P-dot
    pattern = r'([JB]\d{4}[+-]\d{2,4}[A-Za-z]*)\s+(\d+\.\d+(?:e[+-]?\d+)?)\s+([+-]?\d+\.\d+(?:e[+-]?\d+)?|\*)'
    
    atnf_pulsars = {}
    matches = re.findall(pattern, text)
    
    for name, p0_str, p1_str in matches:
        try:
            p0 = float(p0_str)
            has_pdot = p1_str != '*'
            pdot = float(p1_str) if has_pdot else None
            
            if p0 < 0.030:  # MSP
                atnf_pulsars[name] = {
                    'period_s': p0,
                    'period_ms': p0 * 1000,
                    'has_pdot': has_pdot,
                    'pdot': pdot
                }
        except:
            pass
    
    return atnf_pulsars

def normalize_name(name):
    """Normalize pulsar name for matching."""
    # Remove common variations
    name = name.upper()
    name = re.sub(r'[+-]', '', name)
    return name

def crossmatch(freire_pulsars, atnf_pulsars):
    """Cross-match Freire pulsars against ATNF."""
    
    # Build normalized ATNF lookup
    atnf_normalized = {}
    for name, data in atnf_pulsars.items():
        norm = normalize_name(name)
        atnf_normalized[norm] = (name, data)
    
    results = []
    
    for fp in freire_pulsars:
        norm_name = normalize_name(fp['name'])
        
        # Try exact match first
        atnf_match = None
        if norm_name in atnf_normalized:
            atnf_match = atnf_normalized[norm_name]
        else:
            # Try partial match (handle A/B/C suffixes)
            base = norm_name.rstrip('ABCDEFGHIJKLMNOPQRSTUVWXYZ')
            for an, data in atnf_normalized.items():
                if an.startswith(base) or base.startswith(an[:10]):
                    atnf_match = (an, data[1])
                    break
        
        has_pdot_atnf = False
        pdot_atnf = None
        if atnf_match:
            has_pdot_atnf = atnf_match[1]['has_pdot']
            pdot_atnf = atnf_match[1]['pdot']
        
        # Determine final P-dot status
        has_pdot_anywhere = fp['has_pdot_freire'] or has_pdot_atnf
        
        results.append({
            'name': fp['name'],
            'cluster': fp['cluster'],
            'is_msp': fp['is_msp'],
            'period_ms': fp['period_ms'],
            'has_pdot_freire': fp['has_pdot_freire'],
            'has_pdot_atnf': has_pdot_atnf,
            'has_pdot_anywhere': has_pdot_anywhere,
            'in_atnf': atnf_match is not None
        })
    
    return results

def main():
    print("="*70)
    print("CROSS-MATCH VERIFICATION: Freire vs ATNF")
    print("="*70)
    
    # Parse both catalogs
    print("\nParsing Freire catalog...")
    freire_pulsars = parse_freire_catalog()
    print(f"Found {len(freire_pulsars)} pulsars in Freire")
    
    print("\nParsing ATNF catalog...")
    atnf_pulsars = parse_atnf_html()
    print(f"Found {len(atnf_pulsars)} MSPs with data in ATNF")
    
    # Cross-match
    print("\nCross-matching...")
    results = crossmatch(freire_pulsars, atnf_pulsars)
    
    # Analyze results
    print("\n" + "="*70)
    print("RESULTS")
    print("="*70)
    
    # Filter to MSPs only
    msps = [r for r in results if r['is_msp']]
    
    n_msps = len(msps)
    n_pdot_freire = sum(1 for r in msps if r['has_pdot_freire'])
    n_pdot_atnf = sum(1 for r in msps if r['has_pdot_atnf'])
    n_pdot_anywhere = sum(1 for r in msps if r['has_pdot_anywhere'])
    n_in_atnf = sum(1 for r in msps if r['in_atnf'])
    
    print(f"""
    FREIRE GC MSPs (P < 30ms):
    ──────────────────────────
    Total MSPs in Freire:              {n_msps}
    With P-dot in Freire:              {n_pdot_freire}
    With P-dot in ATNF:                {n_pdot_atnf}
    With P-dot in EITHER:              {n_pdot_anywhere}
    Found in ATNF:                     {n_in_atnf}
    """)
    
    # Breakdown by cluster
    print("\nTOP CLUSTERS (MSPs with P-dot):")
    print("─"*50)
    
    cluster_counts = defaultdict(lambda: {'total': 0, 'with_pdot': 0})
    for r in msps:
        if r['cluster']:
            cluster_counts[r['cluster']]['total'] += 1
            if r['has_pdot_anywhere']:
                cluster_counts[r['cluster']]['with_pdot'] += 1
    
    sorted_clusters = sorted(cluster_counts.items(), key=lambda x: x[1]['with_pdot'], reverse=True)[:15]
    for cluster, counts in sorted_clusters:
        print(f"  {cluster:25} {counts['with_pdot']:3}/{counts['total']:3} have P-dot")
    
    # FINAL VERDICT
    print("\n" + "="*70)
    print("FINAL VERDICT")
    print("="*70)
    
    print(f"""
    ┌─────────────────────────────────────────────────────────────────────┐
    │                     VERIFIED GC MSP COUNT                           │
    ├─────────────────────────────────────────────────────────────────────┤
    │                                                                     │
    │   Total GC MSPs in Freire:           {n_msps:4}                         │
    │   With measured P-dot (verified):    {n_pdot_anywhere:4}                         │
    │                                                                     │
    │   This is the DEFENSIBLE number for the manuscript.                 │
    │                                                                     │
    └─────────────────────────────────────────────────────────────────────┘
    
    COMPARISON:
    - Original manuscript:        181 GC MSPs
    - Claimed (inflated):         398 GC MSPs  
    - Verified (defensible):      {n_pdot_anywhere} GC MSPs
    
    True improvement: {n_pdot_anywhere}/181 = {n_pdot_anywhere/181:.1%} ({n_pdot_anywhere - 181:+d})
    """)
    
    # Save results
    import json
    output = {
        'freire_total_msps': n_msps,
        'pdot_in_freire': n_pdot_freire,
        'pdot_in_atnf': n_pdot_atnf,
        'pdot_anywhere': n_pdot_anywhere,
        'in_atnf': n_in_atnf,
        'defensible_gc_count': n_pdot_anywhere,
        'original_count': 181,
        'improvement_factor': n_pdot_anywhere / 181
    }
    
    with open(OUTPUT_DIR / "crossmatch_verification.json", 'w') as f:
        json.dump(output, f, indent=2)
    
    print(f"\nResults saved to: {OUTPUT_DIR / 'crossmatch_verification.json'}")
    
    return output

if __name__ == "__main__":
    main()
