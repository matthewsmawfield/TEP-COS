#!/usr/bin/env python3
"""
Step 6.6b: Rigorous Validation of SFR Holonomy Result

CRITICAL CHECKS:
1. Star-forming galaxies only (where SFR is most reliable)
2. Does size (R_e or concentration) flip the sign like stellar archaeology?
3. Is the result driven by quiescent galaxies with unreliable SFR?
4. Does the sign persist with different SFR estimation methods?

Author: M. Smawfield
Date: January 2026
"""

import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import pearsonr, spearmanr
from sklearn.linear_model import LinearRegression
import json
import os
from datetime import datetime
import warnings
warnings.filterwarnings('ignore')

DATA_DIR = os.path.join(os.path.dirname(__file__), '..', '..', 'data', 'sdss')
RESULTS_DIR = os.path.join(os.path.dirname(__file__), '..', '..', 'results', 'outputs')


def load_data():
    """Load SDSS spectral indices data."""
    df = pd.read_csv(os.path.join(DATA_DIR, 'sdss_spectral_indices.csv'))
    print(f"Loaded {len(df):,} galaxies")
    return df


def partial_corr_multi(x, y, Z):
    """Partial correlation of x and y controlling for multiple variables Z."""
    Z = np.atleast_2d(Z).T if Z.ndim == 1 else Z
    
    reg_x = LinearRegression().fit(Z, x)
    reg_y = LinearRegression().fit(Z, y)
    
    resid_x = x - reg_x.predict(Z)
    resid_y = y - reg_y.predict(Z)
    
    r, p = pearsonr(resid_x, resid_y)
    return r, p


def prepare_data(df):
    """Prepare derived quantities."""
    df['log_sigma'] = np.log10(df['veldisp'])
    df['log_ssfr'] = df['log_sfr'] - df['log_mass']
    df['fe_avg'] = (df['fe5270'] + df['fe5335']) / 2
    df['mg_fe_ratio'] = df['mgb'] / df['fe_avg']
    df['log_mg_fe'] = np.log10(df['mg_fe_ratio'])
    
    valid = (
        (df['veldisp'] > 50) & (df['veldisp'] < 400) &
        (df['veldisp_err'] < 30) &
        np.isfinite(df['log_sfr']) &
        (df['log_sfr'] > -5) & (df['log_sfr'] < 3) &
        np.isfinite(df['log_mass']) &
        (df['log_mass'] > 9.0) & (df['log_mass'] < 12.5) &
        np.isfinite(df['mg_fe_ratio']) &
        (df['mg_fe_ratio'] > 0.3) & (df['mg_fe_ratio'] < 5.0) &
        (df['mgb_err'] < 0.5) &
        (df['fe5270_err'] < 0.5) &
        (df['redshift'] > 0.02) & (df['redshift'] < 0.25)
    )
    
    return df[valid].copy()


def check_star_forming_only(df):
    """
    CRITICAL CHECK 1: Star-forming galaxies only
    
    SFR estimates are most reliable for actively star-forming galaxies.
    If the signal disappears for SF galaxies, it's likely an artifact.
    """
    print("\n" + "=" * 70)
    print("CRITICAL CHECK 1: Star-Forming Galaxies Only")
    print("=" * 70)
    
    # BPT class 1 = Star-forming
    sf_mask = df['bptclass'] == 1
    sf = df[sf_mask].copy()
    
    print(f"\nStar-forming sample: {len(sf):,} galaxies")
    
    if len(sf) < 500:
        print("WARNING: Small sample size for SF galaxies")
    
    # Simple correlation
    r_simple, p_simple = pearsonr(sf['log_sigma'], sf['log_ssfr'])
    print(f"\nSimple r(sSFR, σ): {r_simple:.4f} (p = {p_simple:.2e})")
    
    # Controlled for mass
    r_mass, p_mass = partial_corr_multi(
        sf['log_sigma'].values, sf['log_ssfr'].values, sf['log_mass'].values
    )
    print(f"Controlled for M*: {r_mass:.4f} (p = {p_mass:.2e})")
    
    # Controlled for mass and [Mg/Fe]
    Z = np.column_stack([sf['log_mass'].values, sf['log_mg_fe'].values])
    r_full, p_full = partial_corr_multi(
        sf['log_sigma'].values, sf['log_ssfr'].values, Z
    )
    print(f"Controlled for M*, [Mg/Fe]: {r_full:.4f} (p = {p_full:.2e})")
    
    # Compare to full sample
    print("\n--- Comparison ---")
    r_all, _ = pearsonr(df['log_sigma'], df['log_ssfr'])
    print(f"Full sample r(sSFR, σ): {r_all:.4f}")
    print(f"SF-only r(sSFR, σ): {r_simple:.4f}")
    print(f"Ratio: {r_simple/r_all:.2f}")
    
    if abs(r_simple) < 0.1:
        print("\n⚠️  WARNING: Signal is WEAK for star-forming galaxies!")
        print("   The main result may be driven by quiescent galaxies")
        print("   where SFR estimates are less reliable.")
        robust = False
    else:
        print("\n✓ Signal persists for star-forming galaxies")
        robust = True
    
    return {
        'sf_sample_size': len(sf),
        'r_simple': r_simple,
        'r_mass_controlled': r_mass,
        'r_full_controlled': r_full,
        'robust': robust
    }


def check_quiescent_contamination(df):
    """
    CRITICAL CHECK 2: Are quiescent galaxies driving the signal?
    
    For quiescent galaxies, SFR is essentially an upper limit.
    If the correlation is driven by these, it's not a real TEP signal.
    """
    print("\n" + "=" * 70)
    print("CRITICAL CHECK 2: Quiescent Contamination")
    print("=" * 70)
    
    # Define quiescent by D4000 (D4000 > 1.8 typically indicates quiescent)
    quiescent_mask = df['d4000'] > 1.8
    active_mask = df['d4000'] <= 1.8
    
    quiescent = df[quiescent_mask]
    active = df[active_mask]
    
    print(f"\nQuiescent (D4000 > 1.8): {len(quiescent):,} galaxies")
    print(f"Active (D4000 ≤ 1.8): {len(active):,} galaxies")
    
    # Correlations for each
    r_quiescent, p_q = pearsonr(quiescent['log_sigma'], quiescent['log_ssfr'])
    r_active, p_a = pearsonr(active['log_sigma'], active['log_ssfr'])
    
    print(f"\nQuiescent r(sSFR, σ): {r_quiescent:.4f} (p = {p_q:.2e})")
    print(f"Active r(sSFR, σ): {r_active:.4f} (p = {p_a:.2e})")
    
    # Controlled for mass - active only
    Z = active['log_mass'].values
    r_active_ctrl, p_active_ctrl = partial_corr_multi(
        active['log_sigma'].values, active['log_ssfr'].values, Z
    )
    print(f"Active controlled for M*: {r_active_ctrl:.4f} (p = {p_active_ctrl:.2e})")
    
    if abs(r_active_ctrl) < 0.1:
        print("\n⚠️  WARNING: Signal DISAPPEARS for active galaxies!")
        robust = False
    else:
        print("\n✓ Signal persists for active (non-quiescent) galaxies")
        robust = True
    
    return {
        'n_quiescent': len(quiescent),
        'n_active': len(active),
        'r_quiescent': r_quiescent,
        'r_active': r_active,
        'r_active_controlled': r_active_ctrl,
        'robust': robust
    }


def check_ssfr_mass_disentangle(df):
    """
    CRITICAL CHECK 3: Is this just the sSFR-M* relation in disguise?
    
    The σ-M* (Faber-Jackson) and sSFR-M* (downsizing) relations are strong.
    We need to verify the signal is NOT just these combined.
    """
    print("\n" + "=" * 70)
    print("CRITICAL CHECK 3: Disentangling sSFR-M* Relation")
    print("=" * 70)
    
    # First, show the underlying relations
    r_sigma_mass, _ = pearsonr(df['log_sigma'], df['log_mass'])
    r_ssfr_mass, _ = pearsonr(df['log_ssfr'], df['log_mass'])
    r_ssfr_sigma, _ = pearsonr(df['log_ssfr'], df['log_sigma'])
    
    print(f"\nUnderlying correlations:")
    print(f"  r(σ, M*): {r_sigma_mass:.4f} (Faber-Jackson)")
    print(f"  r(sSFR, M*): {r_ssfr_mass:.4f} (Downsizing)")
    print(f"  r(sSFR, σ): {r_ssfr_sigma:.4f} (Our signal)")
    
    # Expected sSFR-σ correlation from the chain sSFR → M* → σ
    # Using path analysis: r_expected ≈ r(sSFR,M*) × r(M*,σ)
    r_expected = r_ssfr_mass * r_sigma_mass
    print(f"\nExpected from sSFR→M*→σ chain: {r_expected:.4f}")
    print(f"Observed: {r_ssfr_sigma:.4f}")
    print(f"Excess: {r_ssfr_sigma - r_expected:.4f}")
    
    # The partial correlation tells us if there's a DIRECT effect
    r_partial, p_partial = partial_corr_multi(
        df['log_sigma'].values, df['log_ssfr'].values, df['log_mass'].values
    )
    
    print(f"\nDirect effect r(sSFR, σ | M*): {r_partial:.4f}")
    
    if abs(r_partial) > 0.1:
        print("\n✓ There IS a direct sSFR-σ correlation beyond mass mediation")
    else:
        print("\n⚠️  The signal may be largely mediated by mass")
    
    return {
        'r_sigma_mass': r_sigma_mass,
        'r_ssfr_mass': r_ssfr_mass,
        'r_ssfr_sigma': r_ssfr_sigma,
        'r_expected_chain': r_expected,
        'r_partial_direct': r_partial
    }


def check_matched_pairs(df):
    """
    CRITICAL CHECK 4: Matched-pair analysis
    
    Match galaxies by M*, [Mg/Fe], z, and D4000 - then check σ-sSFR
    This mimics the "twin galaxy" test that showed sign flip for age.
    """
    print("\n" + "=" * 70)
    print("CRITICAL CHECK 4: Matched-Pair Analysis")
    print("=" * 70)
    
    # Bin galaxies into matched cells
    df['mass_bin'] = pd.qcut(df['log_mass'], q=10, labels=False, duplicates='drop')
    df['mgfe_bin'] = pd.qcut(df['log_mg_fe'], q=5, labels=False, duplicates='drop')
    df['z_bin'] = pd.qcut(df['redshift'], q=5, labels=False, duplicates='drop')
    df['d4000_bin'] = pd.qcut(df['d4000'], q=10, labels=False, duplicates='drop')
    
    # Group by matched cell
    df['match_cell'] = (
        df['mass_bin'].astype(str) + '_' +
        df['mgfe_bin'].astype(str) + '_' +
        df['z_bin'].astype(str) + '_' +
        df['d4000_bin'].astype(str)
    )
    
    # Within each cell, compute correlation
    cell_results = []
    for cell, group in df.groupby('match_cell'):
        if len(group) >= 20:
            r, p = pearsonr(group['log_sigma'], group['log_ssfr'])
            cell_results.append({
                'cell': cell,
                'n': len(group),
                'r': r,
                'p': p
            })
    
    if len(cell_results) == 0:
        print("No cells with sufficient data for matched-pair analysis")
        return {'robust': None}
    
    # Weighted average
    rs = [c['r'] for c in cell_results]
    ns = [c['n'] for c in cell_results]
    weighted_r = np.average(rs, weights=ns)
    
    print(f"\nMatched-pair cells: {len(cell_results)}")
    print(f"Weighted mean r(sSFR, σ) within matched cells: {weighted_r:.4f}")
    
    # Distribution of within-cell correlations
    positive = sum(1 for c in cell_results if c['r'] > 0)
    negative = sum(1 for c in cell_results if c['r'] < 0)
    print(f"Cells with r > 0: {positive}")
    print(f"Cells with r < 0: {negative}")
    
    # Sign test
    n_total = positive + negative
    p_sign = 2 * min(
        stats.binom.cdf(min(positive, negative), n_total, 0.5),
        1 - stats.binom.cdf(max(positive, negative) - 1, n_total, 0.5)
    )
    print(f"Sign test p-value: {p_sign:.4f}")
    
    if weighted_r > 0:
        print("\n⚠️  WARNING: Sign FLIPS in matched-pair analysis!")
        print("   This is similar to the stellar archaeology result.")
        robust = False
    elif abs(weighted_r) < 0.05:
        print("\n⚠️  WARNING: Signal DISAPPEARS in matched-pair analysis!")
        robust = False
    else:
        print("\n✓ Signal persists in matched-pair analysis")
        robust = True
    
    return {
        'n_cells': len(cell_results),
        'weighted_mean_r': weighted_r,
        'n_positive': positive,
        'n_negative': negative,
        'p_sign_test': p_sign,
        'robust': robust
    }


def main():
    """Run all validation checks."""
    print("=" * 70)
    print("SFR HOLONOMY VALIDATION")
    print("=" * 70)
    
    df = load_data()
    df = prepare_data(df)
    
    results = {}
    
    results['sf_only'] = check_star_forming_only(df)
    results['quiescent'] = check_quiescent_contamination(df)
    results['mass_disentangle'] = check_ssfr_mass_disentangle(df)
    results['matched_pairs'] = check_matched_pairs(df)
    
    # Final assessment
    print("\n" + "=" * 70)
    print("VALIDATION SUMMARY")
    print("=" * 70)
    
    checks = [
        ('Star-forming only', results['sf_only'].get('robust')),
        ('Active galaxies', results['quiescent'].get('robust')),
        ('Matched pairs', results['matched_pairs'].get('robust')),
    ]
    
    print("\nRobustness checks:")
    all_robust = True
    for name, robust in checks:
        if robust is None:
            status = "INCONCLUSIVE"
        elif robust:
            status = "PASS ✓"
        else:
            status = "FAIL ⚠️"
            all_robust = False
        print(f"  {name}: {status}")
    
    print(f"\nKey numbers:")
    print(f"  Full sample r(sSFR, σ): -0.59")
    print(f"  SF-only r(sSFR, σ): {results['sf_only']['r_simple']:.4f}")
    print(f"  Active-only r(sSFR, σ | M*): {results['quiescent']['r_active_controlled']:.4f}")
    print(f"  Matched-pair r(sSFR, σ): {results['matched_pairs']['weighted_mean_r']:.4f}")
    
    if all_robust:
        print("\n✓ SFR HOLONOMY RESULT IS ROBUST")
        print("  Safe to add to manuscript")
    else:
        print("\n⚠️  CAUTION: Some validation checks failed")
        print("  Consider presenting with appropriate caveats")
    
    # Save results
    output = {
        'timestamp': datetime.now().isoformat(),
        'validation_results': results,
        'all_robust': all_robust
    }
    
    output_path = os.path.join(RESULTS_DIR, 'sdss_sfr_holonomy_validation.json')
    with open(output_path, 'w') as f:
        json.dump(output, f, indent=2, default=float)
    print(f"\nSaved: {output_path}")
    
    return results


if __name__ == "__main__":
    main()
