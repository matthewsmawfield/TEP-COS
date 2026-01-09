#!/usr/bin/env python3
"""
Step 5.18: Parse Full ATNF Catalog

Extract all MSPs with measured P-dot from the downloaded ATNF HTML.
"""

import re
from pathlib import Path
from collections import defaultdict
import json

DATA_DIR = Path("/Users/matthewsmawfield/www/TEP-COS/data")
OUTPUT_DIR = Path("/Users/matthewsmawfield/www/TEP-COS/results/outputs")

def parse_atnf_html(filepath):
    """Parse ATNF catalog HTML to extract pulsar data."""
    
    with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
        content = f.read()
    
    pulsars = []
    
    # Find the data table - look for patterns like "J0023+0923" followed by numbers
    # ATNF format: NAME  P0(s)  P1  ASSOC
    
    # Split by table rows
    rows = re.findall(r'<tr[^>]*>(.*?)</tr>', content, re.DOTALL | re.IGNORECASE)
    
    for row in rows:
        # Extract cells
        cells = re.findall(r'<td[^>]*>(.*?)</td>', row, re.DOTALL | re.IGNORECASE)
        
        if len(cells) < 3:
            continue
        
        # Clean cell contents
        cells = [re.sub(r'<[^>]+>', '', c).strip() for c in cells]
        
        # Look for pulsar name in first cell
        name_match = re.search(r'([JB]\d{4}[+-]\d{2,4}[A-Za-z]*)', cells[0])
        if not name_match:
            continue
        
        name = name_match.group(1)
        
        # Try to find P0 and P1
        p0 = None
        p1 = None
        assoc = None
        
        for i, cell in enumerate(cells[1:], 1):
            cell_clean = cell.replace('*', '').strip()
            
            # Check if it's a number
            try:
                val = float(cell_clean)
                if p0 is None and 0.0001 < val < 100:  # Period in seconds
                    p0 = val
                elif p0 is not None and p1 is None:
                    p1 = val
            except:
                # Check for association
                if 'GC' in cell or 'NGC' in cell or 'Ter' in cell or 'M ' in cell:
                    assoc = cell
                elif any(x in cell.lower() for x in ['gc', 'globular']):
                    assoc = 'GC'
        
        if p0 is not None:
            pulsars.append({
                'name': name,
                'P0_s': p0,
                'P0_ms': p0 * 1000,
                'P1': p1,
                'assoc': assoc,
                'is_msp': p0 < 0.030,
                'is_gc': assoc is not None and ('GC' in str(assoc) or 'NGC' in str(assoc) or 'Ter' in str(assoc))
            })
    
    return pulsars

def alternative_parse(filepath):
    """Alternative regex-based parsing."""
    
    with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
        content = f.read()
    
    pulsars = []
    
    # Look for lines with pulsar names and numeric data
    # Pattern: pulsar name, followed by period (small number), then P-dot (scientific notation)
    
    # Find all pulsar names first
    names = re.findall(r'([JB]\d{4}[+-]\d{2,4}[A-Za-z]*)', content)
    unique_names = list(dict.fromkeys(names))  # Remove duplicates, preserve order
    
    print(f"Found {len(unique_names)} unique pulsar names")
    
    # For each name, try to find associated data
    for name in unique_names:
        # Find context around this name
        pattern = re.escape(name) + r'[^<]*?(\d+\.?\d*(?:e[+-]?\d+)?)[^<]*?([+-]?\d+\.?\d*(?:e[+-]?\d+)?)?'
        match = re.search(pattern, content, re.IGNORECASE)
        
        if match:
            try:
                p0 = float(match.group(1))
                p1 = float(match.group(2)) if match.group(2) else None
                
                if 0.0001 < p0 < 100:  # Valid period range
                    pulsars.append({
                        'name': name,
                        'P0_s': p0,
                        'P0_ms': p0 * 1000,
                        'P1': p1,
                        'is_msp': p0 < 0.030
                    })
            except:
                pass
    
    return pulsars

def simple_count(filepath):
    """Simple line-by-line extraction."""
    
    with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
        content = f.read()
    
    # Remove HTML tags
    text = re.sub(r'<[^>]+>', ' ', content)
    
    # Find pulsar entries
    # Pattern: Jxxxx+xxxx or Bxxxx+xx followed by a period value
    
    pattern = r'([JB]\d{4}[+-]\d{2,4}[A-Za-z]*)\s+(\d+\.\d+(?:e[+-]?\d+)?)\s+([+-]?\d+\.\d+(?:e[+-]?\d+)?|\*)'
    
    matches = re.findall(pattern, text)
    
    pulsars = []
    for name, p0_str, p1_str in matches:
        try:
            p0 = float(p0_str)
            p1 = float(p1_str) if p1_str != '*' else None
            
            pulsars.append({
                'name': name,
                'P0_s': p0,
                'P0_ms': p0 * 1000,
                'P1': p1,
                'is_msp': p0 < 0.030,
                'has_pdot': p1 is not None
            })
        except:
            pass
    
    return pulsars

def main():
    print("="*70)
    print("PARSING FULL ATNF CATALOG")
    print("="*70)
    
    filepath = DATA_DIR / "atnf_full_catalog.html"
    
    if not filepath.exists():
        print(f"ERROR: {filepath} not found")
        return
    
    # Try simple extraction
    pulsars = simple_count(filepath)
    
    print(f"\nExtracted {len(pulsars)} pulsars with period data")
    
    # Filter to MSPs
    msps = [p for p in pulsars if p['is_msp']]
    msps_with_pdot = [p for p in msps if p.get('has_pdot', False)]
    
    print(f"MSPs (P < 30ms): {len(msps)}")
    print(f"MSPs with measured P-dot: {len(msps_with_pdot)}")
    
    # Check for GC associations in the raw content
    with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
        content = f.read()
    
    gc_keywords = ['GlobClust', 'NGC', 'Terzan', '47Tuc', 'M28', 'M15', 'M62']
    gc_count = 0
    for kw in gc_keywords:
        gc_count += content.count(kw)
    
    print(f"\nGC-related keywords found: {gc_count}")
    
    # Summary
    print("\n" + "="*70)
    print("ATNF CATALOG SUMMARY")
    print("="*70)
    
    print(f"""
   Total pulsars extracted: {len(pulsars)}
   MSPs (P < 30ms):         {len(msps)}
   MSPs with P-dot:         {len(msps_with_pdot)}
   
   Combined with Freire:
      Freire GC MSPs:       333 (with 202 having P-dot)
      ATNF total MSPs:      {len(msps_with_pdot)}
      
   After deduplication (estimate):
      GC MSPs:              ~333 (from Freire, authoritative for GCs)
      Field MSPs:           ~{len(msps_with_pdot) - 200} (ATNF minus GC overlap)
      
   TOTAL AVAILABLE:         ~{333 + max(0, len(msps_with_pdot) - 200)} MSPs with P-dot
""")
    
    # Sample of field MSPs (non-GC)
    print("\nSample field MSPs (first 10):")
    field_msps = [p for p in msps_with_pdot][:10]
    for p in field_msps:
        print(f"   {p['name']}: P = {p['P0_ms']:.3f} ms, P-dot = {p['P1']:.2e}")
    
    # Save results
    output = {
        'total_extracted': len(pulsars),
        'total_msps': len(msps),
        'msps_with_pdot': len(msps_with_pdot),
        'source': 'ATNF v2.5.1',
        'sample_msps': [{'name': p['name'], 'P_ms': p['P0_ms'], 'Pdot': p['P1']} 
                        for p in msps_with_pdot[:50]]
    }
    
    output_path = OUTPUT_DIR / "atnf_msp_extraction.json"
    with open(output_path, 'w') as f:
        json.dump(output, f, indent=2)
    
    print(f"\nResults saved to: {output_path}")
    
    return msps_with_pdot

if __name__ == "__main__":
    main()
