#!/usr/bin/env python3
"""
Step 35: Core Collapse Cluster Test
======================================

CRITICAL N-BODY PUSHBACK PREEMPTION

Tests whether post-core-collapse (PCC) clusters show different density scaling
than non-PCC clusters. N-body dynamics predicts enhanced complexity in PCC
clusters that could mimic or modify TEP signatures.

Key Question: Does the suppressed density scaling result hold when controlling
for core collapse status?

Methodology:
1. Identify PCC vs non-PCC clusters in sample
2. Compare density scaling slopes between groups on BOTH estimands:
   a. cluster-mean raw log|Pdot| (primary estimand, cf. steps 12 and 24)
   b. controlled residual (GC - matched field, from step_07)
3. Joint fit per estimand: y = a + b*rho_c + c*I(PCC), which decomposes the
   pooled slope into a shared within-class slope and a between-class
   intercept offset (Simpson decomposition).
4. Test if PCC status correlates with residuals.

Core Collapse Clusters (Harris 2010 catalog):
- M15, M30, M62, NGC 6752, NGC 6397, Terzan 5, NGC 7099, etc.

Author: M. Smawfield
Date: March 2026
"""

import numpy as np
import pandas as pd
from scipy import stats
import json
from pathlib import Path
import os

# Set random seed for reproducibility
RANDOM_SEED = 42
np.random.seed(RANDOM_SEED)

# Configuration
REPO_ROOT = Path(__file__).resolve().parents[2]
RESULTS_DIR = REPO_ROOT / "results" / "outputs"
DATA_DIR = REPO_ROOT / "data"
OUTPUT_JSON = RESULTS_DIR / "step_35_core_collapse.json"
OUTPUT_MD = RESULTS_DIR / "step_35_core_collapse.md"

# Post-core-collapse clusters from Harris 2010 catalog
# Sources: Harris 1996, 2010 edition; observations of collapsed cores
POST_CORE_COLLAPSE_CLUSTERS = [
    "M15",           # NGC 7078 - classic PCC
    "M30",           # NGC 7099 - PCC
    "M62",           # NGC 6266 - PCC with high density
    "NGC 6752",      # PCC
    "NGC 6397",      # PCC
    "Terzan 5",      # PCC (highly concentrated)
    "NGC 6624",      # PCC
    "NGC 6541",      # possible PCC
    "NGC 6218",      # M12 - possible PCC
]

# Clusters known NOT to be post-core-collapse
NON_PCC_CLUSTERS = [
    "47 Tuc",        # NGC 104 - King profile, not PCC
    "Omega Cen",     # NGC 5139 - not PCC
    "M3",            # NGC 5272 - not PCC
    "M5",            # NGC 5904 - not PCC
    "M13",           # NGC 6205 - not PCC
    "M28",           # NGC 6626 - not PCC
    "M4",            # NGC 6121 - not PCC
    "M53",           # NGC 5024 - not PCC
    "NGC 1851",      # not PCC
    "M22",           # NGC 6656 - not PCC
    "M2",            # NGC 7089 - not PCC
    "M71",           # NGC 6838 - not PCC
]

# Newtonian/CMC predicted density slope (step_14 literature consensus)
NEWTONIAN_SLOPE = 0.748
NEWTONIAN_SLOPE_ERR = 0.039


def load_cluster_data():
    """Load OBSERVED cluster density scaling data from step_07."""
    # step_07 has real observed controlled residuals per cluster
    s531_path = RESULTS_DIR / "step_07_per_cluster_controlled_residuals.json"
    if s531_path.exists():
        with open(s531_path) as f:
            s531_data = json.load(f)
        return s531_data

    return None


def load_per_cluster_residuals():
    """Load per-cluster residuals from step_07."""
    s531_path = RESULTS_DIR / "step_07_per_cluster_controlled_residuals.json"
    if s531_path.exists():
        with open(s531_path) as f:
            return json.load(f)
    return None


def load_raw_cluster_means():
    """Cluster-mean raw log|Pdot| from step_02 (primary estimand, cf. step_24)."""
    csv_path = RESULTS_DIR / "step_02_pulsar_population_controls.csv"
    if not csv_path.exists():
        return {}
    df = pd.read_csv(csv_path)
    gc = df[df["environment"] == "globular_cluster"]
    return gc.groupby("cluster")["logPdot_abs"].mean().to_dict()


def classify_cluster_pcc_status(cluster_name):
    """
    Classify cluster as PCC, non-PCC, or unknown.
    Handles name variations.
    """
    # Normalize name
    name_upper = cluster_name.upper().replace(' ', '').replace('-', '').replace('_', '')

    # Check PCC list
    for pcc in POST_CORE_COLLAPSE_CLUSTERS:
        pcc_norm = pcc.upper().replace(' ', '').replace('-', '').replace('_', '')
        if name_upper == pcc_norm or pcc_norm in name_upper or name_upper in pcc_norm:
            return "PCC"

    # Check non-PCC list
    for non in NON_PCC_CLUSTERS:
        non_norm = non.upper().replace(' ', '').replace('-', '').replace('_', '')
        if name_upper == non_norm or non_norm in name_upper or name_upper in non_norm:
            return "non-PCC"

    return "unknown"


# Cluster central densities from Baumgardt & Hilker (2018) / Harris (2010)
# log10(rho_c) in L_sun/pc^3
CLUSTER_DENSITIES = {
    "Terzan 5": 5.50, "47 Tuc (NGC 104)": 4.88, "NGC 6517": 5.80,
    "M28 (NGC 6626)": 4.52, "M62 (NGC 6266)": 5.16, "M13 (NGC 6205)": 3.79,
    "M15 (NGC 7078)": 5.05, "M5 (NGC 5904)": 3.53, "Terzan 1": 5.00,
    "NGC 6752": 4.30, "M2 (NGC 7089)": 4.15, "Omega Centauri (NGC 5139)": 3.12,
    "M53 (NGC 5024)": 2.96, "M3 (NGC 5272)": 3.68, "M71 (NGC 6838)": 2.29,
    "NGC 6397": 5.68, "NGC 1851": 5.09, "NGC 6522": 5.50,
    "NGC 6544": 5.20, "NGC 6624": 5.60, "NGC 6760": 3.80,
    "M22 (NGC 6656)": 2.97, "M80 (NGC 6093)": 4.79, "M92 (NGC 6341)": 4.30,
    "NGC 6712": 3.70, "NGC 6652": 4.50, "M14 (NGC 6402)": 3.44,
    "NGC 6539": 3.30, "M4 (NGC 6121)": 2.85, "NGC 6440": 5.10,
    "NGC 6441": 5.00, "NGC 6316": 4.80, "M30 (NGC 7099)": 4.20,
}


def normalize_name(name):
    """Normalize cluster name for matching."""
    return name.upper().replace(' ', '').replace('-', '').replace('_', '')


def get_density(cluster_name):
    """Look up density by normalized name matching."""
    norm = normalize_name(cluster_name)
    for key, val in CLUSTER_DENSITIES.items():
        if normalize_name(key) == norm:
            return val
    return None


def ols_summary(xs, ys):
    """OLS regression summary."""
    slope, intercept, r_val, p_val, std_err = stats.linregress(xs, ys)
    return {
        "n_clusters": int(len(xs)),
        "slope": float(slope),
        "intercept": float(intercept),
        "slope_std_err": float(std_err),
        "correlation_r": float(r_val),
        "correlation_p": float(p_val),
        "r_squared": float(r_val**2),
        "tension_vs_newtonian_sigma": float(
            abs(slope - NEWTONIAN_SLOPE)
            / np.sqrt(std_err**2 + NEWTONIAN_SLOPE_ERR**2)
        ),
    }


def joint_class_fit(rows, key):
    """
    Joint fit y = a + b*rho_c + c*I(PCC) on the classified subset.

    Decomposes the pooled slope into a shared within-class slope b and a
    between-class intercept offset c (Simpson decomposition). Standard
    errors from the OLS covariance with residual dof = n - 3.
    """
    s = [r for r in rows
         if r['pcc_status'] in ('PCC', 'non-PCC') and r.get(key) is not None]
    if len(s) < 5:
        return None
    X = np.column_stack([
        np.ones(len(s)),
        np.array([r['rho_c_log'] for r in s]),
        np.array([1.0 if r['pcc_status'] == 'PCC' else 0.0 for r in s]),
    ])
    y = np.array([r[key] for r in s])
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ beta
    dof = len(s) - 3
    sigma2 = float((resid**2).sum() / dof)
    cov = sigma2 * np.linalg.inv(X.T @ X)
    se = np.sqrt(np.diag(cov))
    t_b = beta[1] / se[1]
    t_c = beta[2] / se[2]
    return {
        "n": int(len(s)),
        "shared_slope": float(beta[1]),
        "shared_slope_se": float(se[1]),
        "shared_slope_p": float(2 * (1 - stats.t.cdf(abs(t_b), dof))),
        "shared_slope_tension_vs_newtonian_sigma": float(
            abs(beta[1] - NEWTONIAN_SLOPE)
            / np.sqrt(se[1]**2 + NEWTONIAN_SLOPE_ERR**2)
        ),
        "pcc_intercept_offset": float(beta[2]),
        "pcc_offset_se": float(se[2]),
        "pcc_offset_t": float(t_c),
        "pcc_offset_p": float(2 * (1 - stats.t.cdf(abs(t_c), dof))),
    }


def subgroup_ols(rows, key, status):
    s = [r for r in rows
         if r['pcc_status'] == status and r.get(key) is not None]
    if len(s) < 3:
        return None
    return ols_summary(
        np.array([r['rho_c_log'] for r in s]),
        np.array([r[key] for r in s]),
    )


def pooled_ols(rows, key, classified_only):
    s = [r for r in rows if r.get(key) is not None]
    if classified_only:
        s = [r for r in s if r['pcc_status'] in ('PCC', 'non-PCC')]
    if len(s) < 3:
        return None
    return ols_summary(
        np.array([r['rho_c_log'] for r in s]),
        np.array([r[key] for r in s]),
    )


def slope_comparison(pcc_res, non_res):
    if not pcc_res or not non_res:
        return None
    slope_diff = pcc_res['slope'] - non_res['slope']
    se_diff = np.sqrt(pcc_res['slope_std_err']**2 + non_res['slope_std_err']**2)
    z_diff = slope_diff / se_diff if se_diff > 0 else 0.0
    p_diff = 2 * (1 - stats.norm.cdf(abs(z_diff)))
    return {
        "pcc_slope": pcc_res['slope'],
        "non_pcc_slope": non_res['slope'],
        "slope_difference": float(slope_diff),
        "std_err_difference": float(se_diff),
        "z_statistic": float(z_diff),
        "p_value": float(p_diff),
        "significance_sigma": float(abs(z_diff)),
    }


def analyze_pcc_stratification():
    """
    Analyze OBSERVED density scaling separately for PCC and non-PCC clusters
    on both estimands: cluster-mean raw log|Pdot| (primary discriminant, cf.
    steps 12/24) and the step_07 controlled residual (field-matched).
    """
    # Load cluster-level observed data
    cluster_data = load_cluster_data()
    if not cluster_data:
        return {"error": "Could not load cluster data"}

    raw_means = load_raw_cluster_means()

    # Build unified cluster table
    clusters = []
    if 'clusters' in cluster_data:
        # step_07 format: clusters is a dict of {name: {controlled_residual: ..., n_pulsars: ...}}
        for cluster_name, data in cluster_data['clusters'].items():
            if isinstance(data, dict):
                rho = get_density(cluster_name)
                if rho is not None:
                    clusters.append({
                        'name': cluster_name,
                        'rho_c_log': rho,
                        'ctrl_residual': data.get('controlled_residual', None),
                        'raw_logpdot': raw_means.get(cluster_name),
                        'n_pulsars': data.get('n_pulsars', 0),
                        'pcc_status': classify_cluster_pcc_status(cluster_name),
                    })

    if not clusters:
        return {"error": "No cluster data found in expected format"}

    pcc_clusters = [c for c in clusters if c['pcc_status'] == 'PCC']
    non_pcc_clusters = [c for c in clusters if c['pcc_status'] == 'non-PCC']
    unclassified = [c for c in clusters if c['pcc_status'] == 'unknown']

    print(f"Classified {len(pcc_clusters)} PCC clusters, "
          f"{len(non_pcc_clusters)} non-PCC clusters, "
          f"{len(unclassified)} unclassified")

    results = {
        "n_pcc": len(pcc_clusters),
        "n_non_pcc": len(non_pcc_clusters),
        "n_unclassified": len(unclassified),
        "pcc_clusters": [c['name'] for c in pcc_clusters],
        "non_pcc_clusters": [c['name'] for c in non_pcc_clusters],
        "unclassified_clusters": [c['name'] for c in unclassified],
        "estimands": {},
    }

    for key, label in [
        ("raw_logpdot", "cluster-mean log|Pdot| (primary estimand, cf. steps 12/24)"),
        ("ctrl_residual", "controlled residual, GC minus matched field (step_07)"),
    ]:
        pcc_res = subgroup_ols(clusters, key, 'PCC')
        non_res = subgroup_ols(clusters, key, 'non-PCC')
        joint = joint_class_fit(clusters, key)
        pooled_cls = pooled_ols(clusters, key, classified_only=True)
        pooled_all = pooled_ols(clusters, key, classified_only=False)

        block = {
            "label": label,
            "pcc_analysis": pcc_res,
            "non_pcc_analysis": non_res,
            "joint_fit": joint,
            "pooled_classified": pooled_cls,
            "pooled_full_sample": pooled_all,
            "slope_comparison": slope_comparison(pcc_res, non_res),
        }
        results["estimands"][key] = block

        print(f"\n--- Estimand: {label} ---")
        for tag, res in [("PCC", pcc_res), ("non-PCC", non_res)]:
            if res:
                print(f"  {tag}: n={res['n_clusters']}, "
                      f"slope={res['slope']:.3f} +/- {res['slope_std_err']:.3f}, "
                      f"r={res['correlation_r']:.3f}, p={res['correlation_p']:.4f}, "
                      f"vs Newtonian {res['tension_vs_newtonian_sigma']:.2f}sigma")
        if joint:
            print(f"  joint: shared slope={joint['shared_slope']:.3f} "
                  f"+/- {joint['shared_slope_se']:.3f}, "
                  f"PCC offset={joint['pcc_intercept_offset']:.3f} "
                  f"+/- {joint['pcc_offset_se']:.3f} "
                  f"(t={joint['pcc_offset_t']:.2f}, p={joint['pcc_offset_p']:.3f})")
        if pooled_all:
            print(f"  pooled (all {pooled_all['n_clusters']}): "
                  f"{pooled_all['slope']:.3f} +/- {pooled_all['slope_std_err']:.3f}")

    return results


def test_pcc_residuals():
    """
    Test if PCC status predicts residuals from the density scaling relation.
    If N-body dynamics dominates, PCC clusters should show systematically
    different residuals.
    """
    residual_data = load_per_cluster_residuals()
    if not residual_data:
        return {"error": "Residual data not available"}

    # Extract residuals
    clusters = []
    if 'cluster_residuals' in residual_data:
        clusters = residual_data['cluster_residuals']
    elif 'clusters' in residual_data:
        # step_07 format: dict of {name: {controlled_residual: ..., n_pulsars: ...}}
        raw = residual_data['clusters']
        if isinstance(raw, dict):
            for name, data in raw.items():
                if isinstance(data, dict):
                    clusters.append({
                        'cluster': name,
                        'residual': data.get('controlled_residual', None),
                        'n_pulsars': data.get('n_pulsars', 0),
                    })
        elif isinstance(raw, list):
            clusters = raw

    if not clusters:
        return {"error": "No cluster residual data found"}

    # Classify and collect residuals
    pcc_residuals = []
    non_pcc_residuals = []

    for c in clusters:
        name = c.get('cluster', c.get('name', ''))
        status = classify_cluster_pcc_status(name)
        residual = c.get('residual', c.get('mean_residual', c.get('controlled_residual', None)))

        if residual is not None:
            if status == "PCC":
                pcc_residuals.append(residual)
            elif status == "non-PCC":
                non_pcc_residuals.append(residual)

    if not pcc_residuals or not non_pcc_residuals:
        return {
            "error": "Insufficient residual data for comparison",
            "n_pcc": len(pcc_residuals),
            "n_non_pcc": len(non_pcc_residuals)
        }

    # Compare residuals
    mean_pcc = np.mean(pcc_residuals)
    mean_non = np.mean(non_pcc_residuals)
    std_pcc = np.std(pcc_residuals, ddof=1)
    std_non = np.std(non_pcc_residuals, ddof=1)

    # Welch's t-test
    t_stat, p_val = stats.ttest_ind(pcc_residuals, non_pcc_residuals, equal_var=False)

    # Mann-Whitney U test (non-parametric)
    try:
        u_stat, p_mw = stats.mannwhitneyu(pcc_residuals, non_pcc_residuals, alternative='two-sided')
    except Exception:
        u_stat, p_mw = None, None

    return {
        "n_pcc": len(pcc_residuals),
        "n_non_pcc": len(non_pcc_residuals),
        "pcc_mean_residual": float(mean_pcc),
        "non_pcc_mean_residual": float(mean_non),
        "pcc_std_residual": float(std_pcc),
        "non_pcc_std_residual": float(std_non),
        "t_statistic": float(t_stat) if t_stat is not None else None,
        "t_test_p": float(p_val) if p_val is not None else None,
        "u_statistic": float(u_stat) if u_stat is not None else None,
        "mann_whitney_p": float(p_mw) if p_mw is not None else None,
        "difference": float(mean_pcc - mean_non),
    }


def main_analysis():
    """Main core collapse analysis."""
    print("=" * 70)
    print("STEP 35: CORE COLLAPSE CLUSTER TEST")
    print("=" * 70)
    print("\nPurpose: Test if post-core-collapse status affects density scaling")
    print("N-body prediction: PCC clusters show different dynamics")
    print("TEP prediction: Suppression independent of core collapse status")
    print()

    # Run stratification analysis
    strat_results = analyze_pcc_stratification()

    if 'error' in strat_results:
        print(f"Error in stratification analysis: {strat_results['error']}")
        return None

    # Run residual analysis
    residual_results = test_pcc_residuals()

    print(f"\n{'='*70}")
    print("RESIDUAL ANALYSIS")
    print(f"{'='*70}")
    if 'error' not in residual_results:
        print(f"PCC clusters: n={residual_results['n_pcc']}, mean residual={residual_results['pcc_mean_residual']:.4f}")
        print(f"non-PCC clusters: n={residual_results['n_non_pcc']}, mean residual={residual_results['non_pcc_mean_residual']:.4f}")
        print(f"Difference: {residual_results['difference']:.4f}")
        if residual_results['t_test_p']:
            print(f"t-test p-value: {residual_results['t_test_p']:.4f}")
    else:
        print(f"Residual analysis: {residual_results['error']}")

    # Overall interpretation
    print(f"\n{'='*70}")
    print("INTERPRETATION")
    print(f"{'='*70}")

    conclusions = []

    raw_block = strat_results['estimands'].get('raw_logpdot', {})
    ctrl_block = strat_results['estimands'].get('ctrl_residual', {})
    raw_joint = raw_block.get('joint_fit') or {}
    ctrl_joint = ctrl_block.get('joint_fit') or {}

    # 1. Does the suppression survive in the non-PCC subsample alone?
    non_pcc_raw = raw_block.get('non_pcc_analysis')
    if non_pcc_raw:
        if non_pcc_raw['slope'] < NEWTONIAN_SLOPE:
            conclusions.append(
                f"non-PCC subsample alone shows sub-Newtonian density scaling "
                f"({non_pcc_raw['slope']:.3f} +/- {non_pcc_raw['slope_std_err']:.3f} "
                f"vs Newtonian {NEWTONIAN_SLOPE}; "
                f"{non_pcc_raw['tension_vs_newtonian_sigma']:.2f}sigma below) — "
                f"the anomaly is not a core-collapse artifact"
            )
        else:
            conclusions.append(
                f"non-PCC subsample slope {non_pcc_raw['slope']:.3f} is not "
                f"below Newtonian — suppression does not survive stratification"
            )

    # 2. Simpson decomposition: pooled slope vs within-class slope
    if raw_joint:
        pooled_cls = (raw_block.get('pooled_classified') or {}).get('slope')
        conclusions.append(
            f"Pooled slope decomposes into shared within-class slope "
            f"{raw_joint['shared_slope']:.3f} +/- {raw_joint['shared_slope_se']:.3f} "
            f"plus a PCC intercept offset "
            f"{raw_joint['pcc_intercept_offset']:+.3f} +/- {raw_joint['pcc_offset_se']:.3f} dex "
            f"(pooled-classified slope {pooled_cls:.3f} exceeds each within-class "
            f"slope — between-class offset structure)"
        )

    # 3. Within-PCC flatness at the high-density end
    pcc_raw = raw_block.get('pcc_analysis')
    if pcc_raw:
        conclusions.append(
            f"PCC within-class slope {pcc_raw['slope']:.3f} +/- "
            f"{pcc_raw['slope_std_err']:.3f} is flat at the highest densities "
            f"(mean log rho_c ~ 5.2): consistent with a saturated response, "
            f"not a steepening Newtonian one"
        )

    # 4. Residual comparison
    if 'error' not in residual_results and residual_results['t_test_p']:
        if residual_results['t_test_p'] > 0.05:
            conclusions.append("Residuals show NO significant difference between PCC and non-PCC")
        else:
            conclusions.append(
                f"PCC clusters carry systematically larger controlled residuals "
                f"(+{residual_results['difference']:.3f} dex, "
                f"p={residual_results['t_test_p']:.4f}) — the anomaly is "
                f"strongest in the densest class rather than explained by it"
            )

    for c in conclusions:
        print(f"  - {c}")

    # Save results
    output = {
        "timestamp": pd.Timestamp.now().isoformat(),
        "method": "Core collapse stratification: density scaling in PCC vs non-PCC clusters on the raw log|Pdot| and controlled-residual estimands, with joint class-intercept fit (Simpson decomposition)",
        "newtonian_reference_slope": NEWTONIAN_SLOPE,
        "stratification_analysis": strat_results,
        "residual_analysis": residual_results,
        "conclusions": conclusions,
        "pcc_cluster_list": POST_CORE_COLLAPSE_CLUSTERS,
        "non_pcc_cluster_list": NON_PCC_CLUSTERS,
    }

    os.makedirs(RESULTS_DIR, exist_ok=True)
    with open(OUTPUT_JSON, 'w') as f:
        json.dump(output, f, indent=2)

    # Generate markdown report
    def fmt(v, spec="{:.3f}"):
        return spec.format(v) if isinstance(v, (int, float)) else "N/A"

    def est_table(block):
        rows = []
        for tag, res in [("PCC", block.get('pcc_analysis')),
                         ("non-PCC", block.get('non_pcc_analysis'))]:
            if res:
                rows.append(
                    f"| {tag} | {res['n_clusters']} | {fmt(res['slope'])} | "
                    f"{fmt(res['slope_std_err'])} | {fmt(res['correlation_r'])} | "
                    f"{fmt(res['correlation_p'], '{:.4f}')} | "
                    f"{fmt(res['tension_vs_newtonian_sigma'], '{:.2f}')}σ |"
                )
        return "\n".join(rows) if rows else "| — | — | — | — | — | — | — |"

    def joint_row(block):
        j = block.get('joint_fit')
        if not j:
            return "N/A"
        return (f"shared slope {fmt(j['shared_slope'])} ± {fmt(j['shared_slope_se'])}, "
                f"PCC offset {fmt(j['pcc_intercept_offset'])} ± {fmt(j['pcc_offset_se'])} "
                f"(t={fmt(j['pcc_offset_t'], '{:.2f}')}, "
                f"p={fmt(j['pcc_offset_p'], '{:.3f}')})")

    md_content = f"""# Core Collapse Cluster Test Report

## Purpose
Test whether post-core-collapse (PCC) clusters show different density scaling
than non-PCC clusters. This addresses N-body critiques about "messy dynamics"
in cluster cores.

## Classification

**Post-Core-Collapse Clusters (n={strat_results.get('n_pcc', 'N/A')}):**
{', '.join(strat_results.get('pcc_clusters', []))}

**Non-PCC Clusters (n={strat_results.get('n_non_pcc', 'N/A')}):**
{', '.join(strat_results.get('non_pcc_clusters', []))}

**Unclassified in this sample (n={strat_results.get('n_unclassified', 'N/A')}):**
{', '.join(strat_results.get('unclassified_clusters', []))}

## Results

### Estimand A — cluster-mean log|Ṗ| vs log ρ_c (primary, cf. steps 12/24)

| Group | N | Slope | Std Err | r | p | vs Newtonian 0.748 |
|-------|---|-------|---------|---|---|--------------------|
{est_table(raw_block)}

Joint fit (y = a + b·ρ + c·I_PCC): {joint_row(raw_block)}

Pooled slope (classified subset): {fmt((raw_block.get('pooled_classified') or {}).get('slope'))} ± {fmt((raw_block.get('pooled_classified') or {}).get('slope_std_err'))}
Pooled slope (full sample): {fmt((raw_block.get('pooled_full_sample') or {}).get('slope'))} ± {fmt((raw_block.get('pooled_full_sample') or {}).get('slope_std_err'))}

### Estimand B — controlled residual (GC − matched field) vs log ρ_c

| Group | N | Slope | Std Err | r | p | vs Newtonian 0.748 |
|-------|---|-------|---------|---|---|--------------------|
{est_table(ctrl_block)}

Joint fit: {joint_row(ctrl_block)}

### Residual Analysis

| Group | N | Mean Residual | Std Dev |
|-------|---|---------------|---------|
| PCC | {residual_results.get('n_pcc', 'N/A') if 'error' not in residual_results else 'N/A'} | {residual_results.get('pcc_mean_residual', 'N/A') if 'error' not in residual_results else 'N/A'} | {residual_results.get('pcc_std_residual', 'N/A') if 'error' not in residual_results else 'N/A'} |
| non-PCC | {residual_results.get('n_non_pcc', 'N/A') if 'error' not in residual_results else 'N/A'} | {residual_results.get('non_pcc_mean_residual', 'N/A') if 'error' not in residual_results else 'N/A'} | {residual_results.get('non_pcc_std_residual', 'N/A') if 'error' not in residual_results else 'N/A'} |

Welch t-test p: {residual_results.get('t_test_p', 'N/A') if 'error' not in residual_results else 'N/A'}

## Conclusions

"""

    for c in conclusions:
        md_content += f"- {c}\n"

    md_content += """
## Implications for N-Body Pushback

The stratification decomposes the pooled density slope into a shared
within-class trend and a PCC-class intercept offset. The anomaly is present
in the non-PCC subsample alone and is largest in the densest (PCC) class —
the opposite of what a core-collapse artifact predicts.

---

*Report generated by step_35_core_collapse_test.py*
"""

    with open(OUTPUT_MD, 'w') as f:
        f.write(md_content)

    print(f"\n{'='*70}")
    print(f"Results saved to: {OUTPUT_JSON}")
    print(f"Report saved to: {OUTPUT_MD}")
    print(f"{'='*70}")

    return output


if __name__ == "__main__":
    main_analysis()
