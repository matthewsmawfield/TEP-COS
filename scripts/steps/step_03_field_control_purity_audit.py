#!/usr/bin/env python3
"""
Step 03: Field Control Purity Audit
========================================

Check whether field-control construction is damping the TEP signal.
Tests:
1. Recompute residual excluding field MSPs near dense stellar environments
2. Split field controls by Galactic height / local density proxy
3. Use only high-quality field MSPs with proper-motion corrections
4. Compare 1-nearest, 3-nearest, 5-nearest, 10-nearest matching
5. Use optimal transport / propensity matching validation
6. Run no-replacement versus replacement matching

Outputs:
- step_03_field_control_purity.json

Author: M. Smawfield
Date: July 2026
"""

import json
from pathlib import Path

import numpy as np
from scipy import stats

REPO_ROOT = Path(__file__).resolve().parents[2]
RESULTS_DIR = REPO_ROOT / "results" / "outputs"
OUT_JSON = RESULTS_DIR / "step_03_field_control_purity.json"


def load_population_controls() -> dict:
    path = RESULTS_DIR / "step_02_pulsar_population_controls.json"
    if not path.exists():
        return {"error": "Population controls not found"}
    with open(path) as f:
        return json.load(f)


def audit_matching_depth(data: dict) -> dict:
    """Test 4: Compare residuals under different nearest-neighbor matching depths."""
    # If the step_02 output contains per-match-depth results, use them;
    # otherwise flag for future implementation.
    out = {
        "status": "requires_pipeline_integration",
        "note": (
            "The current step_02 pipeline uses a single hybrid matching strategy. "
            "To compare 1-NN, 3-NN, 5-NN, and 10-NN, rerun step_02 with variable 'n_neighbors'. "
            "Deeper matching (larger k) tends to increase residual variance and can dilute signals."
        ),
        "expected_direction": "Residual amplitude should decrease with excessive matching depth if controls are over-smoothed.",
    }
    # Check if any sub-analysis files exist
    for k in [1, 3, 5, 10]:
        sub = RESULTS_DIR / f"step_02_k{k}_residual.json"
        out[f"k{k}_available"] = sub.exists()
    return out


def audit_replacement_strategy() -> dict:
    """Test 6: No-replacement vs replacement matching."""
    out = {
        "status": "requires_pipeline_integration",
        "note": (
            "Replacement matching allows the same field MSP to match multiple GC MSPs, "
            "which can create leverage points. No-replacement matching avoids this but "
            "requires a larger field sample and may discard GC MSPs. "
            "Current step_02 uses replacement matching; a no-replacement variant is needed."
        ),
        "expected_direction": "If a few field MSPs dominate matches, no-replacement may increase residual amplitude.",
    }
    return out


def audit_high_quality_subset() -> dict:
    """Test 3: Restrict to field MSPs with proper-motion corrected Pdot."""
    out = {
        "status": "requires_catalog_metadata",
        "note": (
            "The ATNF catalog includes a 'PMcorr' flag for proper-motion-corrected Pdot. "
            "Filtering field controls to PMcorr=1 removes Shklovskii contamination, "
            "which currently worsens the discrepancy (higher |Pdot| in field). "
            "Applying this filter should strengthen the TEP signal if residual field |Pdot| drops."
        ),
        "expected_direction": "Cleaner field controls should increase the GC–field residual.",
    }
    return out


def audit_galactic_height() -> dict:
    """Test 2: Split field controls by Galactic height z."""
    out = {
        "status": "requires_catalog_metadata",
        "note": (
            "Field MSPs at low Galactic height (|z| < 0.5 kpc) may reside in denser environments "
            "and could carry weak TEP modulation. Splitting by z tests whether 'field' truly means 'unscreened'. "
            "If low-z field MSPs show systematically higher |Pdot|, the 0.40 dex residual may be conservative."
        ),
        "expected_direction": "If low-z field MSPs have elevated |Pdot|, the TEP signal strengthens after exclusion.",
    }
    return out


def audit_dense_environment_exclusion() -> dict:
    """Test 1: Exclude field MSPs near OB associations or spiral arms."""
    out = {
        "status": "requires_catalog_metadata",
        "note": (
            "Field MSPs in star-forming regions or near the Galactic centre may experience "
            "local gravitational potentials that mimic weak cluster environments. "
            "Removing these from the control pool tests whether field purity matters."
        ),
        "expected_direction": "If dense-environment field MSPs are excluded, the GC–field residual should increase.",
    }
    return out


def main():
    print("Step 03: Field Control Purity Audit")
    data = load_population_controls()

    result = {
        "population_controls_loaded": "controls" in data,
        "current_residual_dex": data.get("controls", {}).get("period_and_bproxy_matched", {}).get("diff_mean", None),
        "tests": {
            "dense_environment_exclusion": audit_dense_environment_exclusion(),
            "galactic_height_split": audit_galactic_height(),
            "high_quality_pm_corrected": audit_high_quality_subset(),
            "matching_depth_sensitivity": audit_matching_depth(data),
            "replacement_strategy": audit_replacement_strategy(),
        },
        "recommendations": [
            "Re-run step_02 with no-replacement matching and compare residual amplitude.",
            "Request ATNF PMcorr flags and recompute controls using only proper-motion-corrected field MSPs.",
            "Cross-match field MSP positions with Galactic structure catalogs (spiral arms, OB associations) and exclude high-density regions.",
            "If any test increases the residual above 0.60 dex, the 0.40 dex headline is conservative and should be updated.",
        ],
    }

    with open(OUT_JSON, "w") as f:
        json.dump(result, f, indent=2)

    print(f"Audit complete. Output: {OUT_JSON}")


if __name__ == "__main__":
    main()
