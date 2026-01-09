#!/usr/bin/env python3
"""
Binary Spatial Distribution Figure

Plots the cumulative radial distribution of Binary vs Isolated MSPs in GCs
to test for mass segregation effects.

Author: M. Smawfield
Date: 2026-01-09
"""

import re
import json
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from scipy import stats

# Configuration
STYLE_CONFIG = {
    "figsize": (9, 6),
    "dpi": 300,
    "font_family": "serif",
    "font_size": 12,
    "label_size": 14,
    "colors": {
        "binary": "#1f77b4",  # Blue
        "isolated": "#7f7f7f" # Gray
    }
}

# Paths
SCRIPT_DIR = Path(__file__).parent
PROJECT_ROOT = SCRIPT_DIR.parent.parent
DATA_PATH = PROJECT_ROOT / "data" / "freire_gcpsr_2025.txt"
OUTPUT_DIR = PROJECT_ROOT / "results" / "figures" / "manuscript"
STATS_DIR = PROJECT_ROOT / "results" / "outputs"

# Constants
MSP_PERIOD_CUT_MS = 30.0
_NUM_RE = re.compile(r'^[+-]?\d+\.?\d*(?:[eE][+-]?\d+)?$')

def _parse_numeric(val: str):
    if not val or val == '*' or val == 'i':
        return None
    val = re.sub(r'\([^)]*\)', '', val).strip().lstrip('<>')
    if _NUM_RE.match(val):
        return float(val)
    return None

def parse_catalog(text: str):
    rows = []
    current_cluster = None
    
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith('#'):
            continue
            
        # Cluster detection (simplified from step_5_11)
        if not line.startswith('J') and not line.startswith('B') and not line.startswith('*'):
            # Heuristic for cluster headers
            if any(x in line for x in ['NGC', 'Ter', 'M', '47 Tuc']):
                current_cluster = line.strip()
                continue
                
        if current_cluster is None:
            continue
            
        parts = line.split()
        if len(parts) < 5:
            continue
            
        name = parts[0]
        if not (name.startswith('J') or name.startswith('B')):
            continue
            
        try:
            # Flexible parsing due to variable whitespace
            # Name Offset P Pdot ...
            offset_str = parts[1]
            offset = _parse_numeric(offset_str)
            
            p_ms_str = parts[2]
            p_ms = _parse_numeric(p_ms_str)
            
            # Binary info is further down.
            # Standard columns: Name, Offset, P, Pdot, DM, Pb, x, e, m2
            # But sometimes fields are missing or merged.
            # Using the logic that Pb is the 6th column (index 5) usually
            # But let's look at line structure:
            # J0024-7205E 0.65 3.53 +9.85 24.2 2.25 ...
            
            is_binary = False
            pb = None
            
            # Attempt to find Pb (Period binary)
            # It usually appears after DM.
            # Fields: Name(0) Offset(1) P(2) Pdot(3) DM(4) Pb(5)...
            if len(parts) > 5:
                pb_str = parts[5]
                pb = _parse_numeric(pb_str)
                if pb is not None and pb > 0:
                    is_binary = True
            
            if p_ms is not None:
                rows.append({
                    'cluster': current_cluster,
                    'name': name,
                    'offset': offset,
                    'p_ms': p_ms,
                    'is_binary': is_binary
                })
                
        except Exception:
            continue
            
    return rows

def main():
    # Setup
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    STATS_DIR.mkdir(parents=True, exist_ok=True)
    
    # Load Data
    if not DATA_PATH.exists():
        raise FileNotFoundError(f"Data not found at {DATA_PATH}")
        
    rows = parse_catalog(DATA_PATH.read_text())
    
    # Filter MSPs with valid offsets
    msps = [r for r in rows if r['p_ms'] < MSP_PERIOD_CUT_MS and r['offset'] is not None]
    
    # Separate populations
    binaries = [r['offset'] for r in msps if r['is_binary']]
    isolated = [r['offset'] for r in msps if not r['is_binary']]
    
    print(f"Binary MSPs: {len(binaries)}")
    print(f"Isolated MSPs: {len(isolated)}")
    
    # Normalize offsets? 
    # The prompt asks for "spatial distribution", usually normalized by core radius 
    # would be better, but "cumulative radial fraction" of raw offsets is a good first step 
    # if we assume the mix of clusters is similar.
    # However, different clusters have different sizes. 
    # A raw offset comparison might be biased if binaries prefer certain clusters.
    # But usually r/r_c is used. I don't have r_c handy in this file.
    # I will plot raw offset distribution first as requested ("radial fraction" usually implies normalized, 
    # but "cumulative radial fraction of the samples" might just mean CDF of offsets).
    # Given the prompt "Plot the spatial distribution... of the samples in the catalog", 
    # and the advice "If binaries are actually more centrally concentrated...".
    
    # Statistics
    ks_stat, ks_p = stats.ks_2samp(binaries, isolated)
    
    stats_out = {
        "n_binary": len(binaries),
        "n_isolated": len(isolated),
        "median_offset_binary": float(np.median(binaries)),
        "median_offset_isolated": float(np.median(isolated)),
        "mean_offset_binary": float(np.mean(binaries)),
        "mean_offset_isolated": float(np.mean(isolated)),
        "ks_stat": float(ks_stat),
        "ks_p": float(ks_p),
        "interpretation": "Binaries are more centrally concentrated" if np.median(binaries) < np.median(isolated) else "No concentration difference or Isolated more concentrated"
    }
    
    with open(STATS_DIR / "binary_spatial_stats.json", "w") as f:
        json.dump(stats_out, f, indent=2)
        
    # Plotting
    plt.style.use('default')
    plt.rcParams['font.family'] = STYLE_CONFIG['font_family']
    plt.rcParams['font.size'] = STYLE_CONFIG['font_size']
    
    fig, ax = plt.subplots(figsize=STYLE_CONFIG['figsize'])
    
    # Sorted data for CDF
    bin_sorted = np.sort(binaries)
    iso_sorted = np.sort(isolated)
    
    y_bin = np.arange(1, len(bin_sorted) + 1) / len(bin_sorted)
    y_iso = np.arange(1, len(iso_sorted) + 1) / len(iso_sorted)
    
    ax.plot(bin_sorted, y_bin, label=f'Binary MSPs (N={len(binaries)})', 
            color=STYLE_CONFIG['colors']['binary'], lw=2)
    ax.plot(iso_sorted, y_iso, label=f'Isolated MSPs (N={len(isolated)})', 
            color=STYLE_CONFIG['colors']['isolated'], lw=2, linestyle='--')
    
    # Styling
    ax.set_xlabel('Projected Offset (arcmin)', fontsize=STYLE_CONFIG['label_size'])
    ax.set_ylabel('Cumulative Fraction', fontsize=STYLE_CONFIG['label_size'])
    ax.set_title('Spatial Distribution: Binary vs Isolated MSPs', fontsize=STYLE_CONFIG['label_size'], pad=15)
    ax.legend(frameon=False)
    ax.grid(True, alpha=0.3)
    ax.set_xscale('log') # Log scale helps see the core concentration
    ax.set_xlim(0.01, 20)
    
    # Annotate stats
    stats_text = (f"KS Test p-value: {ks_p:.1e}\n"
                  f"Median Offset (Bin): {np.median(binaries):.2f}'\n"
                  f"Median Offset (Iso): {np.median(isolated):.2f}'")
    
    ax.text(0.05, 0.4, stats_text, transform=ax.transAxes, 
            bbox=dict(facecolor='white', alpha=0.8, edgecolor='none'))
            
    plt.tight_layout()
    
    out_png = OUTPUT_DIR / "binary_spatial_distribution.png"
    plt.savefig(out_png, dpi=STYLE_CONFIG['dpi'])
    print(f"Saved figure to {out_png}")

if __name__ == "__main__":
    main()
