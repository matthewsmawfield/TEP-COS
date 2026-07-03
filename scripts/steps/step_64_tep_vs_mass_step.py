#!/usr/bin/env python3
"""
Step 7.1: SN Ia TEP vs. Mass Step Discrimination Analysis
==========================================================

This script performs a critical reanalysis to distinguish TEP effects
from the standard SN Ia "mass step" systematic.

The Problem:
- SN Ia mB shows correlation with host σ (see Step 7.0 output for exact values)
- Partial correlation controlling for host mass shows null result
- Standard interpretation questions whether this is the "mass step" effect

Key Question:
Can we discriminate TEP from the mass step? Both predict correlations
with host potential depth (σ), but through different mechanisms:

1. Mass Step (Astrophysical):
   - Mechanism: Metallicity-driven progenitor differences
   - More massive galaxies → higher metallicity → brighter SNe
   - Amplitude: ΔmB ≈ 0.05-0.10 mag across mass range
   - Observable: Correlation with stellar mass, NOT σ specifically

2. TEP (Gravitational Time Dilation):
   - Mechanism: Clock rate variation in gravitational potential
   - Deeper potential → slower time → brighter apparent magnitude
   - Amplitude: ΔmB ≈ 0.05-0.15 mag depending on potential depth
   - Observable: Correlation with σ (velocity dispersion traces potential)

Discrimination Tests:
1. σ vs. Mass Residual Test: If TEP, residuals after mass correction
   should correlate with σ (time dilation independent of metallicity)
   
2. Continuous Gradient vs. Linear: TEP predicts continuous suppression
   of correlation with screening strength (Temporal Shear suppression modulating
   the Temporal Topology); mass step predicts smooth metallicity gradient
   
3. Evolution with Redshift: TEP amplitude redshift-independent;
   mass step may evolve with galaxy formation history

4. x1 (Stretch) Test: TEP predicts no x1-σ correlation (fossil);
   Mass step might predict weak correlation through progenitor age

Expected Outcomes:
- If pure mass step: r(mB,σ|mass) ≈ 0, smooth σ dependence
- If pure TEP: r(mB,σ|mass) > 0, gradual suppression with screening strength
- If combined: Need to fit both components simultaneously
"""

import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import pearsonr, spearmanr
import json
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Constants
RESULTS_DIR = Path('results/outputs')
OUTPUT_JSON = RESULTS_DIR / 'step_64_tep_vs_mass_step.json'
OUTPUT_MD = RESULTS_DIR / 'step_64_tep_vs_mass_step.md'

# Load previous results
def load_sn_ia_results():
    """Load the SN Ia analysis results."""
    mB_path = RESULTS_DIR / 'step_62_sn_ia_mB_sigma.json'
    stretch_path = RESULTS_DIR / 'step_62_sn_ia_stretch_sigma.json'
    
    with open(mB_path, 'r') as f:
        mB_data = json.load(f)
    with open(stretch_path, 'r') as f:
        stretch_data = json.load(f)
    
    return mB_data, stretch_data

def analyze_mass_step_vs_tep(mB_data, stretch_data):
    """
    Discriminate mass step from TEP effects.
    """
    print("=" * 70)
    print("TEP vs. MASS STEP DISCRIMINATION ANALYSIS")
    print("=" * 70)
    
    results = {
        'analysis_type': 'TEP vs Mass Step Discrimination',
        'timestamp': pd.Timestamp.now().isoformat(),
        'tests': []
    }
    
    # Test 1: Partial Correlation Significance
    print("\n[TEST 1] Partial Correlation Analysis")
    print("-" * 50)
    
    r_raw = mB_data['partial_correlation']['r_raw']
    r_partial = mB_data['partial_correlation']['r_partial_mass_controlled']
    p_partial = mB_data['partial_correlation']['p_partial']
    r_mass_mB = mB_data['partial_correlation']['r_mass_mB']
    r_mass_sigma = mB_data['partial_correlation']['r_mass_sigma']
    
    print(f"Raw correlation r(mB, σ):        {r_raw:.3f}")
    print(f"Partial correlation r(mB, σ|M):    {r_partial:.3f} (p={p_partial:.3f})")
    print(f"Correlation r(mB, M):              {r_mass_mB:.3f}")
    print(f"Correlation r(σ, M):               {r_mass_sigma:.3f}")
    
    # Key discrimination: If TEP were the primary driver, we'd expect
    # partial correlation to remain significant (time dilation independent of mass)
    partial_significant = p_partial < 0.05
    
    if not partial_significant and abs(r_partial) < 0.1:
        discrimination = "MASS_STEP_DOMINATED"
        interpretation = "Correlation fully explained by host mass; no residual TEP signal"
    elif partial_significant and r_partial > 0:
        discrimination = "TEP_DETECTED"
        interpretation = "Residual correlation with σ after mass correction suggests TEP"
    else:
        discrimination = "AMBIGUOUS"
        interpretation = "Unable to discriminate; mixed or null signal"
    
    print(f"\nDiscrimination verdict: {discrimination}")
    print(f"Interpretation: {interpretation}")
    
    test1 = {
        'name': 'Partial Correlation Test',
        'r_raw': r_raw,
        'r_partial': r_partial,
        'p_partial': p_partial,
        'r_mass_mB': r_mass_mB,
        'discrimination': discrimination,
        'interpretation': interpretation
    }
    results['tests'].append(test1)
    
    # Test 2: Collinearity Assessment
    print("\n[TEST 2] σ-Mass Collinearity Assessment")
    print("-" * 50)
    
    # σ and mass are highly correlated, making discrimination difficult
    r_sigma_mass = r_mass_sigma
    vif = 1 / (1 - r_sigma_mass**2)  # Variance Inflation Factor
    
    print(f"Correlation r(σ, M):               {r_sigma_mass:.3f}")
    print(f"Variance Inflation Factor (VIF):   {vif:.2f}")
    
    if vif > 5:
        collinearity = "HIGH"
        collinearity_note = "σ and mass are highly collinear; discrimination is challenging"
    elif vif > 2:
        collinearity = "MODERATE"
        collinearity_note = "Some collinearity; marginal discrimination possible"
    else:
        collinearity = "LOW"
        collinearity_note = "σ and mass are separable; good discrimination potential"
    
    print(f"Collinearity level: {collinearity}")
    print(f"Note: {collinearity_note}")
    
    test2 = {
        'name': 'Collinearity Assessment',
        'r_sigma_mass': r_sigma_mass,
        'vif': vif,
        'collinearity_level': collinearity,
        'note': collinearity_note
    }
    results['tests'].append(test2)
    
    # Test 3: Step Function vs. Linear (Shape Discrimination)
    print("\n[TEST 3] Functional Form Discrimination")
    print("-" * 50)
    
    # TEP predicts step-function at screening threshold
    # Mass step predicts smoother metallicity gradient
    
    unscreened_r = mB_data['screening_test']['unscreened_correlation_r']
    screened_r = mB_data['screening_test']['screened_correlation_r']
    
    print(f"Unscreened regime (σ < 165) r:     {unscreened_r:.3f}")
    print(f"Screened regime (σ ≥ 165) r:       {screened_r:.3f}")
    
    # If TEP: unscreened should show correlation, screened should not
    # If mass step: both regimes should show similar correlations
    
    r_diff = unscreened_r - screened_r
    
    if unscreened_r > 0.15 and abs(screened_r) < 0.15:
        step_verdict = "TEP_CONSISTENT"
        step_note = "Correlation present in unscreened, absent in screened - TEP signature"
    elif abs(r_diff) < 0.1:
        step_verdict = "MASS_STEP_CONSISTENT"
        step_note = "Similar correlations in both regimes - consistent with mass step"
    else:
        step_verdict = "UNCLEAR"
        step_note = "Inconsistent pattern; needs further investigation"
    
    print(f"Correlation difference:            {r_diff:.3f}")
    print(f"Step-function verdict: {step_verdict}")
    print(f"Note: {step_note}")
    
    test3 = {
        'name': 'Functional Form Test',
        'unscreened_r': unscreened_r,
        'screened_r': screened_r,
        'r_difference': r_diff,
        'verdict': step_verdict,
        'note': step_note
    }
    results['tests'].append(test3)
    
    # Test 4: x1 (Stretch) Comparison
    print("\n[TEST 4] x1 (Stretch) Observable Comparison")
    print("-" * 50)
    
    x1_r = stretch_data['stretch_sigma']['r_pearson']
    x1_p = stretch_data['stretch_sigma']['p_pearson']
    x1_r_partial = stretch_data['stretch_sigma']['r_partial']
    x1_p_partial = stretch_data['stretch_sigma']['p_partial']
    
    print(f"x1 vs σ correlation:               r = {x1_r:.3f} (p = {x1_p:.4f})")
    
    # TEP predicts: mB should correlate (rate), x1 should not (fossil)
    # Mass step might predict: both could correlate through metallicity
    
    mB_significant = mB_data['pearson']['p_value'] < 0.05
    x1_significant = x1_p < 0.05
    
    if mB_significant and not x1_significant:
        observable_verdict = "TEP_CONSISTENT"
        observable_note = "mB correlates but x1 does not - matches TEP RATE vs FOSSIL framework"
    elif mB_significant and x1_significant:
        observable_verdict = "MASS_STEP_CONSISTENT"
        observable_note = "Both observables correlate - suggests common driver (mass/metallicity)"
    else:
        observable_verdict = "UNCLEAR"
        observable_note = "Unexpected pattern"
    
    print(f"Observable verdict: {observable_verdict}")
    print(f"Note: {observable_note}")
    
    test4 = {
        'name': 'Observable Comparison Test',
        'mB_correlation_significant': mB_significant,
        'x1_correlation_significant': x1_significant,
        'x1_r': x1_r,
        'x1_p': x1_p,
        'verdict': observable_verdict,
        'note': observable_note
    }
    results['tests'].append(test4)
    
    # Overall Assessment
    print("\n" + "=" * 70)
    print("OVERALL DISCRIMINATION ASSESSMENT")
    print("=" * 70)
    
    # Count evidence for each hypothesis
    tep_evidence = 0
    mass_step_evidence = 0
    
    if discrimination == "TEP_DETECTED":
        tep_evidence += 2
    elif discrimination == "MASS_STEP_DOMINATED":
        mass_step_evidence += 2
    
    if step_verdict == "TEP_CONSISTENT":
        tep_evidence += 1
    elif step_verdict == "MASS_STEP_CONSISTENT":
        mass_step_evidence += 1
    
    if observable_verdict == "TEP_CONSISTENT":
        tep_evidence += 1
    elif observable_verdict == "MASS_STEP_CONSISTENT":
        mass_step_evidence += 1
    
    print(f"\nEvidence tally:")
    print(f"  TEP hypothesis:         {tep_evidence} points")
    print(f"  Mass step hypothesis:   {mass_step_evidence} points")
    
    if mass_step_evidence > tep_evidence:
        overall_verdict = "MASS_STEP_DOMINATED"
        overall_note = "SN Ia mB-σ correlation is primarily the standard mass step effect"
        tep_contribution = "None detected above mass step baseline"
    elif tep_evidence > mass_step_evidence:
        overall_verdict = "TEP_SUPPORTED"
        overall_note = "Evidence favors TEP interpretation over pure mass step"
        tep_contribution = f"Residual correlation r = {r_partial:.3f}"
    else:
        overall_verdict = "INCONCLUSIVE"
        overall_note = "Cannot reliably discriminate TEP from mass step with current data"
        tep_contribution = "Ambiguous"
    
    print(f"\nOverall verdict: {overall_verdict}")
    print(f"Assessment: {overall_note}")
    print(f"TEP contribution: {tep_contribution}")
    
    results['overall'] = {
        'tep_evidence_score': tep_evidence,
        'mass_step_evidence_score': mass_step_evidence,
        'verdict': overall_verdict,
        'assessment': overall_note,
        'tep_contribution': tep_contribution
    }
    
    # Recommendations
    recommendations = []
    if collinearity == "HIGH":
        recommendations.append(
            "HIGH PRIORITY: σ and mass are highly collinear (VIF={:.1f}). ".format(vif) +
            "Need larger sample with better σ-M separation for discrimination."
        )
    
    if discrimination == "MASS_STEP_DOMINATED":
        recommendations.append(
            "SN Ia channel cannot distinguish TEP from mass step. " +
            "Consider this exploratory only; do not present as primary evidence."
        )
    
    if step_verdict == "UNCLEAR":
        recommendations.append(
            "Screening threshold analysis shows unexpected pattern. " +
            "Verify σ_screen = 165 km/s calibration or check for systematic errors."
        )
    
    recommendations.append(
        "Future work: Use independent mass estimates (e.g., SED fitting) " +
        "to break σ-M degeneracy."
    )
    
    results['recommendations'] = recommendations
    
    print("\n" + "=" * 70)
    print("RECOMMENDATIONS")
    print("=" * 70)
    for i, rec in enumerate(recommendations, 1):
        print(f"\n{i}. {rec}")
    
    return results

def save_results(results):
    """Save results to JSON and Markdown."""
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    
    # Save JSON
    with open(OUTPUT_JSON, 'w') as f:
        json.dump(results, f, indent=2)
    
    # Save Markdown
    with open(OUTPUT_MD, 'w') as f:
        f.write("# SN Ia TEP vs. Mass Step Discrimination Analysis\n\n")
        f.write(f"**Analysis Date:** {results['timestamp']}\n\n")
        f.write("## Summary\n\n")
        
        overall = results['overall']
        f.write(f"**Overall Verdict:** {overall['verdict']}\n\n")
        f.write(f"**Assessment:** {overall['assessment']}\n\n")
        f.write(f"**TEP Contribution:** {overall['tep_contribution']}\n\n")
        
        f.write("## Individual Tests\n\n")
        for test in results['tests']:
            f.write(f"### {test['name']}\n\n")
            for key, value in test.items():
                if key != 'name':
                    f.write(f"- **{key}:** {value}\n")
            f.write("\n")
        
        f.write("## Recommendations\n\n")
        for i, rec in enumerate(results['recommendations'], 1):
            f.write(f"{i}. {rec}\n")
    
    print(f"\n\nResults saved to:")
    print(f"  JSON: {OUTPUT_JSON}")
    print(f"  Markdown: {OUTPUT_MD}")

def main():
    print("STEP 7.1: TEP vs. Mass Step Discrimination")
    print("=" * 70)
    print()
    
    # Load previous results
    mB_data, stretch_data = load_sn_ia_results()
    
    # Run discrimination analysis
    results = analyze_mass_step_vs_tep(mB_data, stretch_data)
    
    # Save results
    save_results(results)
    
    print("\n" + "=" * 70)
    print("Analysis complete.")
    print("=" * 70)

if __name__ == "__main__":
    main()
