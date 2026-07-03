#!/usr/bin/env python3
"""
Step 20: Field Binary Pulsar Analysis

Tests the "Field Binary Study" prediction:
"Compare binary vs. isolated pulsars in the field (where cluster acceleration is absent).
If the difference observed in clusters disappears in the field, it supports the TEP environmental interpretation."

This script:
1. Loads the ATNF pulsar catalog (cached from Step 20).
2. Filters for Field MSPs (P < 30 ms, not in GCs).
3. Separates them into Binary (P1 present, Binary type) vs Isolated.
4. Compares their Pdot distributions.

Author: TEP Collaboration
"""

import csv
import io
import json
import math
import re
import tarfile
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
from scipy import stats

# Paths
SCRIPT_DIR = Path(__file__).parent
PROJECT_ROOT = SCRIPT_DIR.parent.parent
RESULTS_DIR = PROJECT_ROOT / "results" / "outputs"
DATA_DIR = PROJECT_ROOT / "data"

# Inputs (from Step 20)
ATNF_DB_PATH = DATA_DIR / "atnf_psrcat.db"
ATNF_TGZ_PATH = DATA_DIR / "atnf_psrcat_pkg.tar.gz"

# Outputs
OUT_JSON = RESULTS_DIR / "step_20_field_binary_analysis.json"
OUT_MD = RESULTS_DIR / "step_20_field_binary_analysis.md"

# Constants
MSP_PERIOD_MAX_S = 0.030  # 30 ms

_NUM_RE = re.compile(r"^[+-]?(?:\d+\.?\d*|\d*\.?\d+)(?:[eE][+-]?\d+)?$")

def _load_atnf_db_content() -> str:
    """Load ATNF DB content from .db file or .tar.gz."""
    # Try .db file first (extracted)
    if ATNF_DB_PATH.exists():
        return ATNF_DB_PATH.read_text(errors="replace")
    
    # Try .tar.gz
    if ATNF_TGZ_PATH.exists():
        with tarfile.open(ATNF_TGZ_PATH, mode="r:gz") as tf:
            # Try common names
            for member in tf.getmembers():
                base = Path(member.name).name
                if base.lower() in {"psrcat.db", "psrcat.db.txt"}:
                    f = tf.extractfile(member)
                    if f is None:
                        continue
                    return f.read().decode("utf-8", errors="replace")
            
            # Fallback: first file with .db extension
            for member in tf.getmembers():
                if member.name.lower().endswith(".db"):
                    f = tf.extractfile(member)
                    if f is None:
                        continue
                    return f.read().decode("utf-8", errors="replace")

    raise FileNotFoundError(
        f"ATNF database not found at {ATNF_DB_PATH} or {ATNF_TGZ_PATH}. "
        "Run Step 20 first."
    )

def _parse_atnf_psrcat_db(db_text: str) -> List[Dict]:
    """Parse ATNF psrcat.db (ASCII block format)."""
    rows = []
    current = {}

    def flush():
        nonlocal current
        if not current:
            return

        name = current.get("PSRJ") or current.get("PSRB")
        
        # Binary status
        binary_type = current.get("BINARY")
        is_binary = binary_type is not None

        # P0, P1
        p0 = current.get("P0")
        p1 = current.get("P1")
        f0 = current.get("F0")
        f1 = current.get("F1")
        
        if p0 is None and f0 is not None and f0 != 0:
            p0 = 1.0 / f0
        if p1 is None and f0 is not None and f1 is not None and f0 != 0:
            p1 = -f1 / (f0 * f0)

        if name is None or p0 is None:
            current = {}
            return

        # Pdot can be None (if not measured)
        p1_val = p1 if p1 is not None else None

        rows.append({
            "name": name,
            "P0_s": p0,
            "P1_sps": p1_val,
            "assoc": current.get("ASSOC", ""),
            "type": current.get("TYPE", ""),
            "is_binary": is_binary,
            "binary_type": binary_type
        })
        current = {}

    for raw_line in db_text.splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith("#"):
            continue
        if line.startswith("@"):  # record separator
            flush()
            continue

        parts = line.split()
        if len(parts) < 2:
            continue

        key = parts[0].strip()
        val0 = parts[1].strip()
        val_rest = " ".join(parts[1:]).strip()

        if key in {"P0", "P1", "F0", "F1"}:
            if _NUM_RE.match(val0):
                current[key] = float(val0)
            continue

        if key in {"PSRJ", "PSRB", "BINARY"}:
            current[key] = val0
        elif key in {"ASSOC", "TYPE"}:
            current[key] = val_rest

    flush()
    return rows

def analyze_field_binaries():
    print("Loading ATNF database...")
    db_text = _load_atnf_db_content()
    rows = _parse_atnf_psrcat_db(db_text)
    print(f"Parsed {len(rows)} pulsars.")

    # Filter for Field MSPs
    field_msps = []
    for r in rows:
        # MSP Cut
        if r["P0_s"] >= MSP_PERIOD_MAX_S:
            continue
            
        # GC Exclusion
        assoc = r["assoc"].lower()
        if "gc" in assoc or "globular" in assoc:
            continue
        if "gc" in str(r["type"]).lower():
            continue
            
        # Must have Pdot measured
        if r["P1_sps"] is None or r["P1_sps"] <= 0:
             continue
             
        field_msps.append(r)

    print(f"Found {len(field_msps)} Field MSPs with measured Pdot.")

    # Split Binary vs Isolated
    binaries = [r for r in field_msps if r["is_binary"]]
    isolated = [r for r in field_msps if not r["is_binary"]]

    print(f"  Binaries: {len(binaries)}")
    print(f"  Isolated: {len(isolated)}")

    # Extract log10(|Pdot|)
    def get_log_pdot(plist):
        return np.array([math.log10(abs(r["P1_sps"] / 1e-20)) for r in plist])

    bin_logpdot = get_log_pdot(binaries)
    iso_logpdot = get_log_pdot(isolated)

    # Statistics
    stats_res = {}
    
    mean_bin = np.mean(bin_logpdot)
    std_bin = np.std(bin_logpdot)
    mean_iso = np.mean(iso_logpdot)
    std_iso = np.std(iso_logpdot)
    
    diff = mean_bin - mean_iso
    
    # T-test
    t_stat, p_val = stats.ttest_ind(bin_logpdot, iso_logpdot, equal_var=False)
    
    # Mann-Whitney
    u_stat, u_p = stats.mannwhitneyu(bin_logpdot, iso_logpdot, alternative='two-sided')

    stats_res = {
        "binary_n": len(binaries),
        "isolated_n": len(isolated),
        "binary_mean_logpdot": float(mean_bin),
        "binary_std_logpdot": float(std_bin),
        "isolated_mean_logpdot": float(mean_iso),
        "isolated_std_logpdot": float(std_iso),
        "diff_dex": float(diff),
        "t_stat": float(t_stat),
        "t_p_value": float(p_val),
        "mw_u_stat": float(u_stat),
        "mw_p_value": float(u_p)
    }

    # Interpretation
    # If p < 0.05, they are different. 
    # TEP prediction: In the field, there should be NO difference (unlike in clusters).
    # So we WANT p > 0.05 (failure to reject null) to support TEP environmental interpretation for the GC signal.
    # However, standard binary evolution theory suggests binaries might be slightly recycled differently, 
    # but "acceleration" shouldn't be a factor in the field.
    
    if p_val > 0.05:
        interpretation = "CONSISTENT: No significant difference between Binary and Isolated Field MSPs. This supports the hypothesis that the difference observed in GCs is environmental (e.g. acceleration or TEP cluster potential)."
    else:
        interpretation = "DIFFERENCE DETECTED: Field Binary and Isolated MSPs differ significantly. This suggests intrinsic evolutionary differences, complicating the GC environmental interpretation."

    stats_res["interpretation"] = interpretation

    # Save JSON
    with open(OUT_JSON, "w") as f:
        json.dump(stats_res, f, indent=2)
        
    # Save Markdown
    md = [
        "# Field Binary vs Isolated MSP Analysis",
        "",
        "## Purpose",
        "To test if the Binary vs Isolated Pdot difference observed in Globular Clusters exists in the Field.",
        "If the difference vanishes in the field, the GC signal is likely environmental.",
        "",
        "## Sample Selection",
        f"- Source: ATNF Pulsar Catalog",
        f"- Criteria: P < {MSP_PERIOD_MAX_S*1000:.1f} ms",
        "- Exclusion: Associated with Globular Clusters",
        "- Requirement: Measured positive Pdot",
        "",
        "## Results",
        "",
        "| Metric | Binary MSPs | Isolated MSPs |",
        "|---|---|---|",
        f"| Count | {len(binaries)} | {len(isolated)} |",
        f"| Mean log10(Pdot_1e20) | {mean_bin:.3f} | {mean_iso:.3f} |",
        f"| Std Dev | {std_bin:.3f} | {std_iso:.3f} |",
        "",
        f"**Difference (Binary - Isolated):** {diff:.3f} dex",
        "",
        "### Statistical Tests",
        f"- Welch's t-test p-value: **{p_val:.4g}**",
        f"- Mann-Whitney U p-value: **{u_p:.4g}**",
        "",
        "## Interpretation",
        f"{interpretation}"
    ]
    
    with open(OUT_MD, "w") as f:
        f.write("\n".join(md))

    print(json.dumps(stats_res, indent=2))
    print(f"Results saved to {OUT_MD}")

if __name__ == "__main__":
    analyze_field_binaries()
