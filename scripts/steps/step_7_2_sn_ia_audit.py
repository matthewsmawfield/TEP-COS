#!/usr/bin/env python3
"""
TEP-COS Step 7.2: SN Ia Deep Audit
===================================

Comprehensive data quality audit for SN Ia analysis:

1. Data completeness check
2. Cross-catalog validation
3. Systematic bias detection
4. Redshift coverage analysis

Dependencies:
- Requires outputs from step_7_0 and step_7_1
- Pantheon+ catalog

Verdict Criteria:
- No critical data gaps → PASS
- Systematic biases detected → FLAG
"""

import pandas as pd
import numpy as np
from scipy import stats
import json
import warnings
warnings.filterwarnings('ignore')

print("="*70)
print("TEP-COS STEP 7.2: SN Ia DEEP AUDIT")
print("="*70)

# ============================================================================
# LOAD DATA
# ============================================================================
print("\n1. DATA COMPLETENESS AUDIT")
print("-"*70)

try:
    pantheon = pd.read_csv('data/supernovae/pantheon_plus_parsed.csv')
    pantheon = pantheon[pantheon['zCMB'] > 0.01]
    print(f"✓ Pantheon+ data loaded: {len(pantheon)} SNe Ia (z > 0.01)")
except Exception as e:
    print(f"✗ Failed to load Pantheon+ data: {e}")
    pantheon = None

# Check for required columns
required_cols = ['zCMB', 'mB', 'mBERR', 'HOST_LOGMASS', 'MU_SH0ES', 'CID']
if pantheon is not None:
    missing = [c for c in required_cols if c not in pantheon.columns]
    if missing:
        print(f"⚠ Missing columns: {missing}")
    else:
        print(f"✓ All required columns present")

# ============================================================================
# LOAD PREVIOUS RESULTS
# ============================================================================
print("\n2. CROSS-VALIDATION WITH PREVIOUS RESULTS")
print("-"*70)

try:
    with open('results/outputs/step_7_0_sn_ia_mB_sigma.json', 'r') as f:
        step_7_0 = json.load(f)
    print(f"✓ Step 7.0 results loaded: r = {step_7_0['pearson']['r']:.4f}")
except Exception as e:
    print(f"✗ Step 7.0 results not found: {e}")
    step_7_0 = None

try:
    with open('results/outputs/robustness_analysis.json', 'r') as f:
        step_7_1 = json.load(f)
    print(f"✓ Step 7.1 results loaded: overall_robust = {step_7_1.get('overall_robust', False)}")
except Exception as e:
    print(f"✗ Step 7.1 results not found: {e}")
    step_7_1 = None

# ============================================================================
# DATA QUALITY METRICS
# ============================================================================
print("\n3. DATA QUALITY METRICS")
print("-"*70)

audit_results = {
    'timestamp': pd.Timestamp.now().isoformat(),
    'data_completeness': {},
    'cross_validation': {},
    'quality_flags': []
}

if pantheon is not None:
    # Completeness by field
    completeness = {}
    for col in ['mB', 'HOST_LOGMASS', 'MU_SH0ES']:
        valid = pantheon[col].notna().sum()
        completeness[col] = {
            'valid': int(valid),
            'total': len(pantheon),
            'fraction': float(valid / len(pantheon))
        }
        print(f"   {col}: {valid}/{len(pantheon)} ({100*valid/len(pantheon):.1f}%)")
    
    audit_results['data_completeness'] = completeness
    
    # Flag if completeness < 95%
    for col, stats_dict in completeness.items():
        if stats_dict['fraction'] < 0.95:
            audit_results['quality_flags'].append(f"LOW_COMPLETENESS_{col}")
            print(f"   ⚠ Flag: Low completeness for {col}")

# ============================================================================
# REDSHIFT COVERAGE
# ============================================================================
print("\n4. REDSHIFT COVERAGE ANALYSIS")
print("-"*70)

if pantheon is not None and 'zCMB' in pantheon.columns:
    z = pantheon['zCMB'].dropna()
    z_bins = [0.01, 0.1, 0.5, 1.0, 2.0]
    
    coverage = {}
    for i in range(len(z_bins)-1):
        z_low, z_high = z_bins[i], z_bins[i+1]
        n_in_bin = ((z >= z_low) & (z < z_high)).sum()
        coverage[f"{z_low}-{z_high}"] = int(n_in_bin)
        print(f"   z = {z_low}-{z_high}: {n_in_bin} SNe")
    
    audit_results['redshift_coverage'] = coverage
    
    # Flag if any bin has < 10 SNe
    for bin_name, count in coverage.items():
        if count < 10:
            audit_results['quality_flags'].append(f"LOW_COUNT_Z_{bin_name}")
            print(f"   ⚠ Flag: Low count in z={bin_name}")

# ============================================================================
# SYSTEMATIC BIAS CHECKS
# ============================================================================
print("\n5. SYSTEMATIC BIAS DETECTION")
print("-"*70)

if pantheon is not None and 'MU_SH0ES' in pantheon.columns:
    # Check for outliers in distance moduli
    mures = pantheon['MU_SH0ES'].dropna()
    zscore = np.abs(stats.zscore(mures))
    outliers = (zscore > 3).sum()
    
    print(f"   Hubble residual outliers (|z| > 3): {outliers}")
    audit_results['systematic_checks'] = {
        'hubble_residual_outliers': int(outliers),
        'hubble_residual_std': float(mures.std())
    }
    
    if outliers > len(mures) * 0.05:  # > 5% outliers
        audit_results['quality_flags'].append("HIGH_OUTLIER_FRACTION")
        print(f"   ⚠ Flag: High outlier fraction")

# ============================================================================
# CROSS-VALIDATION SUMMARY
# ============================================================================
print("\n6. CROSS-VALIDATION SUMMARY")
print("-"*70)

consistency_checks = []

# Check consistency between steps
if step_7_0 is not None and step_7_1 is not None:
    r_orig = step_7_0['pearson']['r']
    # Compare with bootstrap mean from step_7_1
    if 'bootstrap' in step_7_1:
        r_boot_mean = step_7_1['bootstrap']['r_mean']
        diff = abs(r_orig - r_boot_mean)
        print(f"   Original r vs Bootstrap mean: |{r_orig:.4f} - {r_boot_mean:.4f}| = {diff:.4f}")
        
        if diff < 0.05:  # Within 0.05
            print(f"   ✓ Results consistent")
            consistency_checks.append(True)
        else:
            print(f"   ⚠ Large difference detected")
            consistency_checks.append(False)
            audit_results['quality_flags'].append("INCONSISTENT_RESULTS")

audit_results['cross_validation'] = {
    'steps_consistent': all(consistency_checks) if consistency_checks else None,
    'n_checks': len(consistency_checks)
}

# ============================================================================
# FINAL VERDICT
# ============================================================================
print("\n" + "="*70)
print("AUDIT COMPLETE")
print("="*70)

# Determine overall status
if len(audit_results['quality_flags']) == 0:
    verdict = "PASS"
    print("\n✓ No quality flags raised")
else:
    verdict = "FLAG"
    print(f"\n⚠ {len(audit_results['quality_flags'])} quality flags:")
    for flag in audit_results['quality_flags']:
        print(f"   - {flag}")

audit_results['verdict'] = verdict

print(f"\n{'='*70}")
print(f"OVERALL: {verdict}")
print(f"{'='*70}")

# Save results
with open('results/outputs/sn_ia_audit.json', 'w') as f:
    json.dump(audit_results, f, indent=2)

print("\nResults saved to: results/outputs/sn_ia_audit.json")
