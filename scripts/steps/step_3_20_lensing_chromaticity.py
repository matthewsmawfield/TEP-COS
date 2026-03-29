#!/usr/bin/env python3
"""
Step 3.20: Lensing Chromaticity Test (FIXED VERSION)
====================================================

CRITICAL FIX: Removed simulation-based analysis. 
Now performs actual data availability check and reports honest status.

Purpose: Distinguish TEP (achromatic) from microlensing (chromatic)
TEP predicts: ΔΓ = 0 across bands (achromatic)
Microlensing predicts: ΔΓ ≠ 0 (chromatic, wavelength-dependent)

STATUS: Real multi-band analysis requires archival data processing.
Previous versions used np.random.normal() simulations which were misleading.

Author: M. Smawfield
Date: March 2026 (Fixed)
"""

import numpy as np
import pandas as pd
from scipy import stats
import json
from pathlib import Path
import os

# Configuration
REPO_ROOT = Path(__file__).resolve().parents[2]
RESULTS_DIR = REPO_ROOT / "results" / "outputs"
DATA_DIR = REPO_ROOT / "data" / "cosmograil"
OUTPUT_JSON = RESULTS_DIR / "step_3_20_lensing_chromaticity.json"

# Random seed for reproducibility
RANDOM_SEED = 42
np.random.seed(RANDOM_SEED)

# Multi-band systems from COSMOGRAIL
MULTIBAND_SYSTEMS = {
    "DESJ0408": {
        "bands": ["R"],
        "z_lens": 0.60,
        "z_source": 2.21,
        "has_multiband": False,
        "note": "Requires archival multi-band data"
    },
    "RXJ1131": {
        "bands": ["R", "V"],
        "z_lens": 0.30,
        "z_source": 0.66,
        "has_multiband": True,
        "note": "RXJ1131-1231 has R and V band"
    },
    "Q2237": {
        "bands": ["R", "V", "I"],
        "z_lens": 0.04,
        "z_source": 1.69,
        "has_multiband": True,
        "note": "Einstein Cross - multi-epoch multiband"
    },
    "HE0435": {
        "bands": ["R", "V"],
        "z_lens": 0.45,
        "z_source": 1.69,
        "has_multiband": True,
        "note": "H0LiCOW system"
    },
    "PG1115": {
        "bands": ["R", "I"],
        "z_lens": 0.31,
        "z_source": 1.72,
        "has_multiband": True,
        "note": "Multi-band available"
    }
}


def simulate_chromaticity_test(system_name, bands, n_bootstrap=100):
    """
    OBSOLETE: Simulation-based analysis removed.
    
    Previous versions used np.random.normal() to simulate gamma measurements.
    This was misleading as it generated fake results that looked like real data.
    
    Use analyze_chromaticity_from_real_data() instead.
    """
    raise NotImplementedError(
        "Simulation-based analysis removed. "
        "Use check_data_availability() and analyze from real data."
    )


def check_data_availability():
    """Check what real multi-band data is actually available."""
    available_data = {}
    if not DATA_DIR.exists():
        return {}
    for f in DATA_DIR.glob("**/*"):
        fname = f.name.lower()
        for system in ["desj0408", "rxj1131", "q2237", "he0435", "pg1115"]:
            if system in fname:
                if system not in available_data:
                    available_data[system] = {"files": [], "bands": set()}
                available_data[system]["files"].append(str(f.name))
                for band in ["g", "r", "i", "z", "R", "V", "I", "B"]:
                    if f"_{band.lower()}." in fname or f"_{band.upper()}." in fname:
                        available_data[system]["bands"].add(band.upper())
    return available_data


def analyze_chromaticity_from_real_data(system_name, available_data):
    """Check data availability for chromaticity analysis."""
    system_key = system_name.lower()
    if system_key not in available_data:
        return {"system": system_name, "status": "NO_DATA", "n_bands": 0,
                "message": "No data files found"}
    data = available_data[system_key]
    bands = list(data["bands"])
    n_bands = len(bands)
    if n_bands < 2:
        return {"system": system_name, "status": "INSUFFICIENT_BANDS", 
                "n_bands": n_bands, "bands": bands,
                "message": f"Only {n_bands} band(s). Need ≥2 for chromaticity test."}
    return {"system": system_name, "status": "DATA_AVAILABLE_ANALYSIS_PENDING",
            "n_bands": n_bands, "bands": bands, "n_files": len(data["files"]),
            "message": f"Multi-band data available ({n_bands} bands). Analysis pending."}


def analyze_all_systems():
    """Main analysis: Check real data availability (NO SIMULATION)."""
    print("=" * 70)
    print("STEP 3.20: LENSING CHROMATICITY TEST (FIXED)")
    print("=" * 70)
    print("\nCRITICAL: Simulation removed. Reporting actual data status.")
    print()
    
    available_data = check_data_availability()
    print("Data availability check:")
    print("-" * 50)
    
    results = {}
    data_available = 0
    insufficient = 0
    no_data = 0
    
    for system_name in MULTIBAND_SYSTEMS:
        result = analyze_chromaticity_from_real_data(system_name, available_data)
        results[system_name] = result
        print(f"\n{system_name}: {result['status']}")
        if result['status'] == "DATA_AVAILABLE_ANALYSIS_PENDING":
            data_available += 1
            print(f"  Bands: {result['bands']}, Files: {result['n_files']}")
        elif result['status'] == "INSUFFICIENT_BANDS":
            insufficient += 1
        else:
            no_data += 1
    
    print(f"\n{'='*70}")
    print("HONEST ASSESSMENT (NO SIMULATION)")
    print(f"{'='*70}")
    print(f"Multi-band data available: {data_available}/{len(MULTIBAND_SYSTEMS)}")
    print(f"Insufficient bands: {insufficient}/{len(MULTIBAND_SYSTEMS)}")
    print(f"No data: {no_data}/{len(MULTIBAND_SYSTEMS)}")
    
    if data_available == 0:
        overall_verdict = "DATA_REQUIRED"
        print("\nVerdict: DATA_REQUIRED")
        print("Real multi-band analysis not performed.")
    else:
        overall_verdict = "ANALYSIS_PENDING"
        print(f"\nVerdict: ANALYSIS_PENDING")
        print(f"{data_available} systems have multi-band data.")
    
    print("\nWARNING: Previous simulation results are INVALID.")
    
    output = {
        "timestamp": pd.Timestamp.now().isoformat(),
        "method": "Real data availability check (NO SIMULATION)",
        "note": "Simulation removed. Reports actual data status only.",
        "systems": results,
        "summary": {
            "total_systems": len(MULTIBAND_SYSTEMS),
            "data_available": data_available,
            "insufficient_bands": insufficient,
            "no_data": no_data,
            "verdict": overall_verdict
        },
        "interpretation": {
            "tep_prediction": "Achromatic (ΔΓ = 0)",
            "microlensing_prediction": "Chromatic (ΔΓ ≠ 0)",
            "current_status": "Real analysis pending - simulation removed"
        },
        "recommendations": [
            "URGENT: Perform real multi-band temporal shear analysis",
            "Q2237: 13 files with B,G,I,R,V - analyze each band",
            "HE0435: 4 files with R,V - analyze each band",
            "Remove all claims based on previous simulation"
        ]
    }
    
    os.makedirs(RESULTS_DIR, exist_ok=True)
    with open(OUTPUT_JSON, 'w') as f:
        json.dump(output, f, indent=2)
    
    print(f"\n{'='*70}")
    print(f"Results saved to: {OUTPUT_JSON}")
    print(f"{'='*70}")
    return output


if __name__ == "__main__":
    analyze_all_systems()
