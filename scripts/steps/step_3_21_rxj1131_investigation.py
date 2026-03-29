#!/usr/bin/env python3
"""
Step 3.21: RXJ1131 Anomaly Investigation & Lensing Data Audit
=============================================================

Purpose: Investigate the RXJ1131 chromaticity anomaly and audit
the lensing analysis for methodological weaknesses.

Key Questions:
1. Is the RXJ1131 chromaticity result from real data or simulation?
2. What is the actual multi-band data availability?
3. Are there systematic errors in the temporal shear analysis?
4. Can TEP be distinguished from microlensing with current data?

Critical Finding from Initial Review:
- The step_3_20_lensing_chromaticity.py uses SIMULATED data
- RXJ1131 "chromaticity" is a random simulation outcome, not real
- This is a major methodological weakness in the lensing channel
"""

import numpy as np
import pandas as pd
from scipy import stats
import json
from pathlib import Path
import os

REPO_ROOT = Path(__file__).resolve().parents[2]
RESULTS_DIR = REPO_ROOT / "results" / "outputs"
DATA_DIR = REPO_ROOT / "data" / "cosmograil"
OUTPUT_JSON = RESULTS_DIR / "step_3_21_rxj1131_investigation.json"
OUTPUT_MD = RESULTS_DIR / "step_3_21_rxj1131_investigation.md"

# Random seed for reproducibility
RANDOM_SEED = 42
np.random.seed(RANDOM_SEED)

def check_multiband_data_availability():
    """
    Check what actual multi-band data files exist.
    """
    print("=" * 70)
    print("DATA AVAILABILITY AUDIT")
    print("=" * 70)
    
    if not DATA_DIR.exists():
        print(f"Data directory not found: {DATA_DIR}")
        return {}
    
    # List all files in cosmograil directory
    all_files = list(DATA_DIR.glob("**/*"))
    
    # Look for multi-band indicators
    multiband_systems = {}
    
    for f in all_files:
        fname = f.name.lower()
        
        # Check for RXJ1131 files
        if "rxj1131" in fname or "rxj" in fname:
            if "rxj1131" not in multiband_systems:
                multiband_systems["rxj1131"] = {"files": [], "bands": set()}
            multiband_systems["rxj1131"]["files"].append(str(f.relative_to(REPO_ROOT)))
            
            # Extract band from filename
            for band in ["g", "r", "i", "z", "R", "V", "I", "B"]:
                if f"_{band.lower()}." in fname or f"_{band.upper()}." in fname:
                    multiband_systems["rxj1131"]["bands"].add(band.upper())
        
        # Check for other systems
        for system in ["desj0408", "he0435", "q2237", "pg1115"]:
            if system in fname:
                if system not in multiband_systems:
                    multiband_systems[system] = {"files": [], "bands": set()}
                multiband_systems[system]["files"].append(str(f.relative_to(REPO_ROOT)))
                
                # Extract band
                for band in ["g", "r", "i", "z", "R", "V", "I", "B"]:
                    if f"_{band.lower()}." in fname or f"_{band.upper()}." in fname:
                        multiband_systems[system]["bands"].add(band.upper())
    
    print(f"\nSystems with data files found:")
    for system, info in multiband_systems.items():
        bands = sorted(info["bands"]) if info["bands"] else ["unknown/single-band"]
        print(f"  {system.upper()}: {len(info['files'])} files, bands: {bands}")
    
    return multiband_systems

def analyze_chromaticity_methodology():
    """
    Analyze the methodology used in chromaticity testing.
    """
    print("\n" + "=" * 70)
    print("METHODOLOGY AUDIT")
    print("=" * 70)
    
    issues = []
    strengths = []
    
    # Check if analysis is simulation-based
    chromaticity_script = REPO_ROOT / "scripts/steps/step_3_20_lensing_chromaticity.py"
    
    if chromaticity_script.exists():
        with open(chromaticity_script, 'r') as f:
            content = f.read()
        
        if "simulate" in content.lower() or "random.normal" in content:
            issues.append({
                "severity": "CRITICAL",
                "issue": "Chromaticity test uses SIMULATED data, not real multi-band analysis",
                "details": "The RXJ1131 'chromaticity' is from random noise simulation, not actual data",
                "impact": "Lensing channel cannot distinguish TEP from microlensing - no real multi-band test performed"
            })
        
        if "has_multiband" in content:
            # Check what systems actually have multi-band flagged
            if "RXJ1131" in content and "has_multiband" in content:
                issues.append({
                    "severity": "HIGH",
                    "issue": "RXJ1131 flagged as 'has_multiband' but no real analysis performed",
                    "details": "Simulation parameters used instead of actual data",
                    "impact": "False confidence in multi-band coverage"
                })
    
    # Check temporal shear results
    temporal_shear_file = RESULTS_DIR / "step_3_0_cosmograil_temporal_shear.json"
    
    if temporal_shear_file.exists():
        with open(temporal_shear_file, 'r') as f:
            temporal_data = json.load(f)
        
        # Analyze significance levels
        significances = []
        for system, data in temporal_data.get("systems", {}).items():
            for pair, pair_data in data.get("image_pairs", {}).items():
                if "gamma" in pair_data and "gamma_sigma" in pair_data:
                    significances.append({
                        "system": system,
                        "pair": pair,
                        "sigma": pair_data["gamma_sigma"]
                    })
        
        if significances:
            print(f"\nTemporal shear significance levels:")
            for sig in significances:
                status = "SIGNIFICANT" if sig["sigma"] > 2 else "MARGINAL" if sig["sigma"] > 1 else "NULL"
                print(f"  {sig['system']} {sig['pair']}: {sig['sigma']:.2f}σ ({status})")
            
            # Count significant detections
            n_significant = sum(1 for s in significances if s["sigma"] > 2)
            n_total = len(significances)
            
            if n_significant == 0:
                issues.append({
                    "severity": "HIGH",
                    "issue": f"No significant ( >2σ) temporal shear detections in {n_total} measurements",
                    "details": "All gamma measurements are consistent with null (Γ=0)",
                    "impact": "No evidence for TEP-GL temporal shear effect"
                })
            elif n_significant < n_total / 2:
                issues.append({
                    "severity": "MODERATE",
                    "issue": f"Only {n_significant}/{n_total} measurements show >2σ significance",
                    "details": "Mixed results - some systems show signal, others don't",
                    "impact": "Inconsistent evidence for TEP-GL"
                })
    
    # Check gamma values for consistency
    chromaticity_file = RESULTS_DIR / "step_3_20_lensing_chromaticity.json"
    
    if chromaticity_file.exists():
        with open(chromaticity_file, 'r') as f:
            chromaticity_data = json.load(f)
        
        chromatic_count = chromaticity_data.get("summary", {}).get("chromatic_count", 0)
        achromatic_count = chromaticity_data.get("summary", {}).get("achromatic_count", 0)
        
        print(f"\nChromaticity test results:")
        print(f"  Achromatic systems: {achromatic_count}")
        print(f"  Chromatic systems: {chromatic_count}")
        
        # Note that these are simulations
        issues.append({
            "severity": "CRITICAL",
            "issue": "Chromaticity results are from SIMULATION, not data",
            "details": f"Reported {achromatic_count} achromatic, {chromatic_count} chromatic from random noise",
            "impact": "Cannot use lensing as evidence until real multi-band analysis performed"
        })
    
    return issues, strengths

def evaluate_tep_discrimination_capability():
    """
    Evaluate whether current data can distinguish TEP from GR/microlensing.
    """
    print("\n" + "=" * 70)
    print("TEP DISCRIMINATION ASSESSMENT")
    print("=" * 70)
    
    assessment = {}
    
    # Key test: Can we distinguish TEP-GL from standard GR lensing?
    # TEP-GL predicts: scale-dependent delays (Γ ≠ 0)
    # Standard GR predicts: constant delays (Γ = 0)
    # Microlensing predicts: chromatic effects (different Γ per band)
    
    print("\nRequired discriminating tests:")
    print("  1. Γ ≠ 0 (temporal shear detection) - NOT ACHIEVED")
    print("     Current: Marginal detections (1-2σ), no >3σ confirmation")
    print()
    print("  2. Achromatic (ΔΓ = 0 across bands) - NOT TESTED")
    print("     Current: Simulation only, no real multi-band analysis")
    print()
    print("  3. Microlensing independence - UNCLEAR")
    print("     Current: Detrending applied but effectiveness uncertain")
    
    assessment["gamma_detection"] = {
        "status": "NOT_ACHIEVED",
        "current_significance": "1-2σ marginal",
        "required": ">3σ for confirmation",
        "verdict": "Insufficient evidence for TEP-GL temporal shear"
    }
    
    assessment["chromaticity"] = {
        "status": "NOT_TESTED",
        "current": "Simulation only",
        "required": "Real multi-band data analysis",
        "verdict": "Cannot distinguish TEP from microlensing"
    }
    
    assessment["microlensing_separation"] = {
        "status": "UNCLEAR",
        "method": "Gaussian detrending (200-day window)",
        "effectiveness": "Unknown - no validation performed",
        "verdict": "Systematic uncertainty unquantified"
    }
    
    return assessment

def generate_recommendations(issues, assessment):
    """
    Generate specific recommendations for strengthening the lensing channel.
    """
    print("\n" + "=" * 70)
    print("RECOMMENDATIONS")
    print("=" * 70)
    
    recommendations = []
    
    # Critical recommendations
    critical_issues = [i for i in issues if i["severity"] == "CRITICAL"]
    if critical_issues:
        recommendations.append({
            "priority": "URGENT",
            "action": "Halt lensing claims until real multi-band analysis performed",
            "rationale": "Current results are simulation-based, not data-driven",
            "implementation": "Remove or re-label lensing as 'preliminary simulation'"
        })
    
    # Data acquisition
    recommendations.append({
        "priority": "HIGH",
        "action": "Acquire and analyze archival COSMOGRAIL multi-band data",
        "rationale": "Need real g, r, i, R, V band light curves for chromaticity test",
        "implementation": "Contact COSMOGRAIL collaboration for data access"
    })
    
    # Methodology improvements
    recommendations.append({
        "priority": "HIGH",
        "action": "Validate microlensing detrending effectiveness",
        "rationale": "200-day Gaussian smoothing may not fully remove microlensing",
        "implementation": "Test on simulated microlensing + TEP signals"
    })
    
    # Statistical requirements
    recommendations.append({
        "priority": "MEDIUM",
        "action": "Increase sample size or observation epochs",
        "rationale": "Current significance levels (1-2σ) insufficient for discovery",
        "implementation": "Target >3σ per system or combine multiple systems"
    })
    
    # Alternative approaches
    recommendations.append({
        "priority": "MEDIUM",
        "action": "Consider alternative TEP-GL observables",
        "rationale": "Temporal shear may be too subtle with current precision",
        "implementation": "Explore frequency-dependent magnification or time-delay ratios"
    })
    
    for rec in recommendations:
        print(f"\n[{rec['priority']}] {rec['action']}")
        print(f"    Rationale: {rec['rationale']}")
        print(f"    Implementation: {rec['implementation']}")
    
    return recommendations

def main():
    """
    Main investigation pipeline.
    """
    print("=" * 70)
    print("RXJ1131 ANOMALY INVESTIGATION & LENSING CHANNEL AUDIT")
    print("=" * 70)
    print()
    
    # 1. Check data availability
    multiband_data = check_multiband_data_availability()
    
    # 2. Analyze methodology
    issues, strengths = analyze_chromaticity_methodology()
    
    # 3. Evaluate discrimination capability
    assessment = evaluate_tep_discrimination_capability()
    
    # 4. Generate recommendations
    recommendations = generate_recommendations(issues, assessment)
    
    # 5. Overall verdict
    print("\n" + "=" * 70)
    print("OVERALL ASSESSMENT")
    print("=" * 70)
    
    critical_count = sum(1 for i in issues if i["severity"] == "CRITICAL")
    high_count = sum(1 for i in issues if i["severity"] == "HIGH")
    
    print(f"\nIssues found:")
    print(f"  CRITICAL: {critical_count}")
    print(f"  HIGH: {high_count}")
    print(f"  MODERATE: {sum(1 for i in issues if i['severity'] == 'MODERATE')}")
    
    if critical_count > 0:
        overall_verdict = "NOT_PUBLICATION_READY"
        overall_assessment = "Lensing channel relies on simulations, not data. Cannot be presented as evidence."
    elif high_count > 0:
        overall_verdict = "WEAK"
        overall_assessment = "Significant methodological issues limit confidence in results."
    else:
        overall_verdict = "ACCEPTABLE"
        overall_assessment = "Methodology sound but significance marginal."
    
    print(f"\nOverall verdict: {overall_verdict}")
    print(f"Assessment: {overall_assessment}")
    
    # Save results
    os.makedirs(RESULTS_DIR, exist_ok=True)
    
    results = {
        "timestamp": pd.Timestamp.now().isoformat(),
        "multiband_data_availability": {k: {"files": len(v["files"]), "bands": list(v["bands"])} 
                                          for k, v in multiband_data.items()},
        "methodology_issues": issues,
        "methodology_strengths": strengths,
        "tep_discrimination": assessment,
        "recommendations": recommendations,
        "overall_verdict": overall_verdict,
        "overall_assessment": overall_assessment,
        "key_finding": "RXJ1131 chromaticity result is from SIMULATION, not real data analysis"
    }
    
    with open(OUTPUT_JSON, 'w') as f:
        json.dump(results, f, indent=2)
    
    # Write markdown report
    with open(OUTPUT_MD, 'w') as f:
        f.write("# RXJ1131 Anomaly Investigation & Lensing Channel Audit\n\n")
        f.write(f"**Investigation Date:** {results['timestamp']}\n\n")
        f.write("## Executive Summary\n\n")
        f.write(f"**Overall Verdict:** {overall_verdict}\n\n")
        f.write(f"{overall_assessment}\n\n")
        f.write(f"**Key Finding:** {results['key_finding']}\n\n")
        
        f.write("## Data Availability\n\n")
        if multiband_data:
            for system, info in multiband_data.items():
                f.write(f"- **{system.upper()}**: {len(info['files'])} files, bands: {list(info['bands'])}\n")
        else:
            f.write("No multi-band data files found in cosmograil directory.\n")
        
        f.write("\n## Methodology Issues\n\n")
        for issue in issues:
            f.write(f"### [{issue['severity']}] {issue['issue']}\n\n")
            f.write(f"- **Details:** {issue['details']}\n")
            f.write(f"- **Impact:** {issue['impact']}\n\n")
        
        f.write("## TEP Discrimination Assessment\n\n")
        for test, data in assessment.items():
            f.write(f"### {test.replace('_', ' ').title()}\n\n")
            for key, value in data.items():
                f.write(f"- **{key}:** {value}\n")
            f.write("\n")
        
        f.write("## Recommendations\n\n")
        for rec in recommendations:
            f.write(f"### [{rec['priority']}] {rec['action']}\n\n")
            f.write(f"- **Rationale:** {rec['rationale']}\n")
            f.write(f"- **Implementation:** {rec['implementation']}\n\n")
    
    print(f"\n{'='*70}")
    print(f"Results saved to:")
    print(f"  JSON: {OUTPUT_JSON}")
    print(f"  Markdown: {OUTPUT_MD}")
    print(f"{'='*70}")

if __name__ == "__main__":
    main()
