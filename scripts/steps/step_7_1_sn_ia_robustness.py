#!/usr/bin/env python3
"""
TEP-COS Step 7.1: SN Ia σ-mB Robustness Validation
==================================================

Academic Standard Robustness Analysis
-------------------------------------
This script performs comprehensive robustness tests on the SN Ia peak 
magnitude vs host velocity dispersion correlation (Step 7.0):

1. Outlier removal (z-score > 3)
2. Subsample analysis (cluster exclusion, quality cuts)
3. Bootstrap confidence intervals (10,000 samples)

Dependencies:
- Requires output from step_7_0_sn_ia_stretch_test.py
- Pantheon+ and SDSS specObj catalogs

Verdict Criteria:
- All subsample tests significant → ROBUST
- Any test null → CONDITIONAL
"""

import pandas as pd
import numpy as np
from scipy import stats
from scipy.stats import pearsonr, spearmanr, linregress, ttest_ind
import json
import warnings
warnings.filterwarnings('ignore')

# Statistical thresholds
SIGNIFICANCE_THRESHOLD = 0.05
BOOTSTRAP_SAMPLES_LARGE = 10000

np.random.seed(42)

print("="*70)
print("TEP-COS STEP 7.1: SN Ia σ-mB ROBUSTNESS VALIDATION")
print("="*70)

# ============================================================================
# LOAD DATA
# ============================================================================
print("\n1. LOADING DATA")
print("-"*70)

pantheon = pd.read_csv('data/supernovae/pantheon_plus_parsed.csv')
# Apply z > 0.01 filter to match step_7_0
pantheon = pantheon[pantheon['zCMB'] > 0.01]
specobj = pd.read_csv('data/supernovae/sdss_sigma_specobj_matches.csv')
specobj_valid = specobj[specobj['sigma_host'] > 0]
merged = pantheon.merge(specobj_valid[['CID', 'sigma_host']], on='CID', how='inner')
valid_full = merged[merged['mB'].notna()].copy()

# Apply quality cuts (matches step_7_0: σ > 30, σ <= 400)
valid_full = valid_full[(valid_full['sigma_host'] >= 30) & (valid_full['sigma_host'] <= 400)]

print(f"Quality-filtered sample: {len(valid_full)} SNe (z > 0.01, 30 <= σ <= 400 km/s)")

results = {
    'timestamp': pd.Timestamp.now().isoformat(),
    'full_sample_size': len(valid_full),
    'tests': {}
}

# ============================================================================
# TEST 1: OUTLIER REMOVAL
# ============================================================================
print("\n2. TEST 1: OUTLIER REMOVAL (Z-score > 3)")
print("-"*70)

z_sigma = np.abs(stats.zscore(valid_full['sigma_host']))
z_mb = np.abs(stats.zscore(valid_full['mB']))

# Remove outliers
outlier_mask = (z_sigma > 3) | (z_mb > 3)
valid_no_outliers = valid_full[~outlier_mask].copy()

print(f"Removed {outlier_mask.sum()} outliers")
print(f"Clean sample: {len(valid_no_outliers)} SNe")

log_sigma_clean = np.log10(valid_no_outliers['sigma_host'])
mB_clean = valid_no_outliers['mB']

r_clean, p_clean = pearsonr(log_sigma_clean, mB_clean)
slope_clean, intercept_clean, r_val_clean, p_val_clean, std_err_clean = linregress(log_sigma_clean, mB_clean)
# Correct significance: convert two-tailed p-value to Gaussian sigma
sig_clean = abs(stats.norm.ppf(p_clean / 2))

print(f"  r = {r_clean:.6f}")
print(f"  p = {p_clean:.6e}")
print(f"  Significance = {sig_clean:.2f}σ")
print(f"  Slope = {slope_clean:.4f} ± {std_err_clean:.4f}")

results['tests']['outlier_removal'] = {
    'n_sample': len(valid_no_outliers),
    'n_removed': int(outlier_mask.sum()),
    'pearson_r': float(r_clean),
    'pearson_p': float(p_clean),
    'significance_sigma': float(sig_clean),
    'slope': float(slope_clean),
    'slope_err': float(std_err_clean),
    'verdict': 'TEP-CONSISTENT' if (p_clean < 0.05 and r_clean > 0.05) else 'NULL'
}

# ============================================================================
# TEST 2: SUBSAMPLE ANALYSIS
# ============================================================================
print("\n3. TEST 2: SUBSAMPLE ANALYSIS")
print("-"*70)

# Define subsamples
subsamples = {
    'no_clusters': valid_full[valid_full['sigma_host'] < 300],  # Exclude cluster members
    'no_low_quality': valid_full[valid_full['sigma_host'] > 40],  # Exclude low-σ
    'clean_range': valid_full[(valid_full['sigma_host'] >= 40) & (valid_full['sigma_host'] < 300)],
    'low_z': valid_full[valid_full['zCMB'] < 0.3],  # Nearby SNe
    'mid_sigma': valid_full[(valid_full['sigma_host'] >= 60) & (valid_full['sigma_host'] < 200)]
}

subsample_results = {}
for name, data in subsamples.items():
    if len(data) < 20:
        print(f"  {name}: Insufficient data (n={len(data)})")
        continue
    
    log_sigma = np.log10(data['sigma_host'])
    mB = data['mB']
    
    r, p = pearsonr(log_sigma, mB)
    slope, intercept, r_val, p_val, std_err = linregress(log_sigma, mB)
    # Correct significance: convert two-tailed p-value to Gaussian sigma
    sig = abs(stats.norm.ppf(p / 2))
    
    print(f"  {name}: n={len(data)}, r={r:.4f}, p={p:.4e}, σ={sig:.2f}σ")
    
    subsample_results[name] = {
        'n_sample': len(data),
        'pearson_r': float(r),
        'pearson_p': float(p),
        'significance_sigma': float(sig),
        'slope': float(slope),
        'slope_err': float(std_err),
        'verdict': 'TEP-CONSISTENT' if (p < SIGNIFICANCE_THRESHOLD and r > 0.05) else 'NULL'
    }

results['tests']['subsamples'] = subsample_results

# ============================================================================
# TEST 3: BOOTSTRAP CONFIDENCE INTERVALS
# ============================================================================
print("\n4. TEST 3: BOOTSTRAP CONFIDENCE INTERVALS")
print("-"*70)

n_bootstrap = BOOTSTRAP_SAMPLES_LARGE
boot_slopes = []
boot_intercepts = []
boot_r = []

n_sample = len(valid_full)
log_sigma_full = np.log10(valid_full['sigma_host'])
mB_full = valid_full['mB']

for i in range(n_bootstrap):
    # Resample with replacement
    idx = np.random.choice(n_sample, size=n_sample, replace=True)
    log_s_boot = log_sigma_full.iloc[idx]
    mB_boot = mB_full.iloc[idx]
    
    slope, intercept, r_val, p_val, std_err = linregress(log_s_boot, mB_boot)
    boot_slopes.append(slope)
    boot_intercepts.append(intercept)
    boot_r.append(r_val)

boot_slopes = np.array(boot_slopes)
boot_r = np.array(boot_r)

# Calculate confidence intervals
slope_ci_95 = (np.percentile(boot_slopes, 2.5), np.percentile(boot_slopes, 97.5))
slope_ci_99 = (np.percentile(boot_slopes, 0.5), np.percentile(boot_slopes, 99.5))
r_ci_95 = (np.percentile(boot_r, 2.5), np.percentile(boot_r, 97.5))

print(f"Bootstrap samples: {n_bootstrap}")
print(f"Slope distribution:")
print(f"  Mean: {boot_slopes.mean():.4f}")
print(f"  Std: {boot_slopes.std():.4f}")
print(f"  95% CI: [{slope_ci_95[0]:.4f}, {slope_ci_95[1]:.4f}]")
print(f"  99% CI: [{slope_ci_99[0]:.4f}, {slope_ci_99[1]:.4f}]")
print(f"  P(slope > 0): {(boot_slopes > 0).mean():.4f}")

print(f"\nr distribution:")
print(f"  Mean: {boot_r.mean():.4f}")
print(f"  95% CI: [{r_ci_95[0]:.4f}, {r_ci_95[1]:.4f}]")

results['tests']['bootstrap'] = {
    'n_bootstrap': n_bootstrap,
    'slope_mean': float(boot_slopes.mean()),
    'slope_std': float(boot_slopes.std()),
    'slope_ci_95': [float(slope_ci_95[0]), float(slope_ci_95[1])],
    'slope_ci_99': [float(slope_ci_99[0]), float(slope_ci_99[1])],
    'p_slope_positive': float((boot_slopes > 0).mean()),
    'r_ci_95': [float(r_ci_95[0]), float(r_ci_95[1])]
}

# ============================================================================
# COMPARISON WITH ORIGINAL
# ============================================================================
print("\n5. COMPARISON WITH ORIGINAL ANALYSIS")
print("-"*70)

# Load original results
with open('results/outputs/step_7_0_sn_ia_mB_sigma.json', 'r') as f:
    orig = json.load(f)

print(f"Original: r={orig['pearson']['r']:.4f}, σ={orig['pearson']['significance_sigma']:.2f}σ")
print(f"No outliers: r={r_clean:.4f}, σ={sig_clean:.2f}σ")
print(f"Bootstrap mean: r={boot_r.mean():.4f}")

# All significant?
all_sig = all([
    orig['pearson']['p_value'] < 0.05,
    p_clean < 0.05,
    slope_ci_95[0] > 0  # 95% CI excludes zero
])

print(f"\nAll tests significant: {all_sig}")

# ============================================================================
# FINAL SUMMARY
# ============================================================================
print("\n" + "="*70)
print("ROBUSTNESS ASSESSMENT COMPLETE")
print("="*70)

print("\n1. OUTLIER REMOVAL:")
print(f"   Signal persists: {results['tests']['outlier_removal']['verdict']}")

print("\n2. SUBSAMPLE TESTS:")
for name, res in subsample_results.items():
    print(f"   {name}: {res['verdict']} ({res['significance_sigma']:.2f}σ)")

print("\n3. BOOTSTRAP:")
print(f"   Slope 95% CI: [{slope_ci_95[0]:.3f}, {slope_ci_95[1]:.3f}]")
print(f"   P(slope > 0) = {(boot_slopes > 0).mean():.4f}")

overall_robust = all_sig and all(r['verdict'] == 'TEP-CONSISTENT' for r in subsample_results.values())
results['overall_robust'] = overall_robust
results['all_significant'] = all_sig

print(f"\n{'='*70}")
print(f"OVERALL: {'Consistent with TEP prediction' if overall_robust else 'Conditional - requires further validation'}")
print(f"{'='*70}")

# Save results
with open('results/outputs/robustness_analysis.json', 'w') as f:
    json.dump(results, f, indent=2)

print("\nResults saved to: results/outputs/robustness_analysis.json")
