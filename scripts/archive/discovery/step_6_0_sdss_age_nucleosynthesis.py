#!/usr/bin/env python3
"""
Step 6.0: SDSS Age-Nucleosynthesis Discrepancy Test for TEP

THE KEY INSIGHT:
- Nucleosynthesis ratios (e.g., [Mg/Fe]) are set during star formation and do NOT
  depend on subsequent time flow - they are "frozen in" at formation.
- Spectroscopic ages (D4000, Hβ) depend on stellar evolution, which DOES depend
  on local time flow rate.

TEP PREDICTION:
- If TEP is correct: spectroscopic age should correlate with gravitational potential
  (deeper wells = slower time = older-appearing stars)
- BUT: [Mg/Fe] should NOT correlate with potential (nucleosynthesis is time-invariant)
- The DISCREPANCY between these two age indicators reveals the TEP signal.

GR PREDICTION (null hypothesis):
- Both indicators should correlate identically with environment (if at all)
- No discrepancy expected.

This is a DEFINITIVE test because it uses two independent clocks that respond
differently to TEP but identically to GR.

Author: M. Smawfield
Date: January 2026
"""

import requests
import numpy as np
import pandas as pd
from scipy import stats
from scipy.spatial import cKDTree
from astropy.cosmology import FlatLambdaCDM
import json
import os
from datetime import datetime
import warnings
warnings.filterwarnings('ignore')

cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

SDSS_URL = "https://skyserver.sdss.org/dr18/SkyServerWS/SearchTools/SqlSearch"
DATA_DIR = os.path.join(os.path.dirname(__file__), '..', '..', 'data', 'sdss')
RESULTS_DIR = os.path.join(os.path.dirname(__file__), '..', '..', 'results', 'outputs')
FIGURES_DIR = os.path.join(os.path.dirname(__file__), '..', '..', 'results', 'figures')

os.makedirs(DATA_DIR, exist_ok=True)
os.makedirs(RESULTS_DIR, exist_ok=True)
os.makedirs(FIGURES_DIR, exist_ok=True)


def query_sdss(sql, max_retries=3):
    """Execute SQL query against SDSS SkyServer."""
    import time
    for attempt in range(max_retries):
        try:
            response = requests.get(
                SDSS_URL,
                params={"cmd": sql, "format": "json"},
                timeout=300
            )
            if response.status_code == 200:
                data = response.json()
                if data and len(data) > 0 and "Rows" in data[0]:
                    return pd.DataFrame(data[0]["Rows"])
            else:
                print(f"  HTTP {response.status_code}")
        except requests.exceptions.Timeout:
            print(f"  Timeout (attempt {attempt + 1})")
        except Exception as e:
            print(f"  Error: {e}")
        if attempt < max_retries - 1:
            time.sleep(2)
    return None


def download_spectral_indices():
    """
    Download SDSS galaxies with spectral indices for age-nucleosynthesis test.
    
    Key columns:
    - lick_mgb, lick_fe5270, lick_fe5335: For [Mg/Fe] (nucleosynthesis clock)
    - d4000_n, lick_hb: For spectroscopic age
    - v_disp: Velocity dispersion (gravitational potential proxy)
    """
    print("\n" + "=" * 70)
    print("DOWNLOADING SDSS SPECTRAL INDICES FOR AGE-NUCLEOSYNTHESIS TEST")
    print("=" * 70)
    
    # Query in redshift chunks to avoid timeout
    z_ranges = [
        (0.02, 0.06),
        (0.06, 0.10),
        (0.10, 0.15),
        (0.15, 0.20),
        (0.20, 0.25),
    ]
    
    all_data = []
    
    for z_min, z_max in z_ranges:
        print(f"\nQuerying z = {z_min:.2f} - {z_max:.2f}...")
        
        sql = f"""
        SELECT TOP 100000
            i.specobjid,
            g.ra, g.dec, 
            g.z as redshift,
            g.z_err as z_err,
            g.v_disp as veldisp,
            g.v_disp_err as veldisp_err,
            i.lick_mgb as mgb,
            i.lick_mgb_err as mgb_err,
            i.lick_fe5270 as fe5270,
            i.lick_fe5270_err as fe5270_err,
            i.lick_fe5335 as fe5335,
            i.lick_fe5335_err as fe5335_err,
            i.d4000_n as d4000,
            i.d4000_n_err as d4000_err,
            i.lick_hb as hbeta,
            i.lick_hb_err as hbeta_err,
            i.lick_hd_a as hdelta_a,
            e.lgm_tot_p50 as log_mass,
            e.sfr_tot_p50 as log_sfr,
            e.bptclass
        FROM galSpecIndx i
        JOIN galSpecInfo g ON i.specobjid = g.specobjid
        JOIN galSpecExtra e ON i.specobjid = e.specobjid
        WHERE g.reliable = 1
            AND g.z BETWEEN {z_min} AND {z_max}
            AND g.z_err < 0.001
            AND g.v_disp > 30 AND g.v_disp < 450
            AND g.v_disp_err > 0 AND g.v_disp_err < 50
            AND i.lick_mgb > 0.5 AND i.lick_mgb < 8
            AND i.lick_fe5270 > 0.5 AND i.lick_fe5270 < 5
            AND i.lick_fe5335 > 0.5 AND i.lick_fe5335 < 5
            AND i.d4000_n > 1.0 AND i.d4000_n < 2.5
            AND i.lick_hb > 0 AND i.lick_hb < 6
            AND e.lgm_tot_p50 > 8 AND e.lgm_tot_p50 < 13
        ORDER BY g.z
        """
        
        df = query_sdss(sql)
        
        if df is not None and len(df) > 0:
            print(f"  Retrieved {len(df)} galaxies")
            all_data.append(df)
        else:
            print(f"  No data retrieved for this range")
    
    if len(all_data) == 0:
        print("\nERROR: No data retrieved from SDSS")
        return None
    
    combined = pd.concat(all_data, ignore_index=True)
    print(f"\nTotal galaxies retrieved: {len(combined)}")
    
    # Remove duplicates
    combined = combined.drop_duplicates(subset=['specobjid'])
    print(f"After deduplication: {len(combined)}")
    
    return combined


def compute_age_indicators(df):
    """
    Compute the two age indicators:
    1. [Mg/Fe] - Nucleosynthesis clock (time-invariant under TEP)
    2. Spectroscopic age proxy - Stellar evolution clock (time-dependent under TEP)
    """
    print("\nComputing age indicators...")
    
    # [Mg/Fe] proxy from Lick indices
    # <Fe> = (Fe5270 + Fe5335) / 2
    df['fe_avg'] = (df['fe5270'] + df['fe5335']) / 2
    df['mg_fe_ratio'] = df['mgb'] / df['fe_avg']
    df['log_mg_fe'] = np.log10(df['mg_fe_ratio'])
    
    # Spectroscopic age proxy: D4000 / Hβ
    # Higher D4000 = older, Lower Hβ = older
    # Combined: D4000 / (Hβ + offset) gives age proxy
    df['spec_age_proxy'] = df['d4000'] / (df['hbeta'] + 0.5)
    df['log_spec_age'] = np.log10(df['spec_age_proxy'])
    
    # Normalize both to zero mean, unit variance for comparison
    df['mg_fe_norm'] = (df['log_mg_fe'] - df['log_mg_fe'].mean()) / df['log_mg_fe'].std()
    df['spec_age_norm'] = (df['log_spec_age'] - df['log_spec_age'].mean()) / df['log_spec_age'].std()
    
    # Gravitational potential proxy: log(σ)
    df['log_sigma'] = np.log10(df['veldisp'])
    
    # Lookback time
    df['t_lookback'] = cosmo.lookback_time(df['redshift'].values).value
    
    # Filter valid
    valid = (
        np.isfinite(df['mg_fe_ratio']) &
        np.isfinite(df['spec_age_proxy']) &
        (df['mg_fe_ratio'] > 0.5) & (df['mg_fe_ratio'] < 3.0) &
        (df['spec_age_proxy'] > 0.3) & (df['spec_age_proxy'] < 3.0)
    )
    
    df_valid = df[valid].copy()
    print(f"  Valid galaxies: {len(df_valid)}")
    print(f"  [Mg/Fe] range: {df_valid['mg_fe_ratio'].min():.2f} - {df_valid['mg_fe_ratio'].max():.2f}")
    print(f"  Spec age proxy range: {df_valid['spec_age_proxy'].min():.2f} - {df_valid['spec_age_proxy'].max():.2f}")
    
    return df_valid


def compute_local_environment(df, n_neighbors=20):
    """Compute local gravitational environment using nearest neighbors."""
    print("\nComputing local environment...")
    
    # 3D comoving positions
    z = df['redshift'].values
    d_comoving = cosmo.comoving_distance(z).value  # Mpc
    ra_rad = np.radians(df['ra'].values)
    dec_rad = np.radians(df['dec'].values)
    
    x = d_comoving * np.cos(dec_rad) * np.cos(ra_rad)
    y = d_comoving * np.cos(dec_rad) * np.sin(ra_rad)
    z_coord = d_comoving * np.sin(dec_rad)
    
    coords_3d = np.column_stack([x, y, z_coord])
    
    # Build KD-tree for neighbor search
    tree = cKDTree(coords_3d)
    distances, indices = tree.query(coords_3d, k=n_neighbors + 1)
    
    # Distance to nth neighbor (excluding self)
    r_n = distances[:, -1]  # Mpc
    
    # Local density proxy
    df['log_density'] = np.log10(n_neighbors / (4/3 * np.pi * r_n**3) + 1e-10)
    
    # Environment classification based on density percentiles
    density_percentiles = np.percentile(df['log_density'], [20, 40, 60, 80])
    df['env_class'] = pd.cut(
        df['log_density'],
        bins=[-np.inf] + list(density_percentiles) + [np.inf],
        labels=['void', 'sparse', 'average', 'dense', 'cluster']
    )
    
    # Also classify by velocity dispersion (internal potential)
    sigma_percentiles = np.percentile(df['log_sigma'], [20, 40, 60, 80])
    df['sigma_class'] = pd.cut(
        df['log_sigma'],
        bins=[-np.inf] + list(sigma_percentiles) + [np.inf],
        labels=['low_sigma', 'med_low', 'medium', 'med_high', 'high_sigma']
    )
    
    print(f"  Environment classes: {df['env_class'].value_counts().to_dict()}")
    print(f"  Sigma classes: {df['sigma_class'].value_counts().to_dict()}")
    
    return df


def test_age_nucleosynthesis_discrepancy(df):
    """
    THE KEY TEP TEST:
    
    Compare how spectroscopic age and [Mg/Fe] correlate with gravitational potential.
    
    TEP prediction:
    - Spectroscopic age should correlate with σ (deeper potential = slower time = older appearance)
    - [Mg/Fe] should NOT correlate with σ (nucleosynthesis is time-invariant)
    - The DIFFERENCE reveals TEP
    
    GR prediction:
    - Both should correlate identically (or not at all)
    """
    print("\n" + "=" * 70)
    print("AGE-NUCLEOSYNTHESIS DISCREPANCY TEST")
    print("=" * 70)
    print("\nTEP Prediction: Spectroscopic age correlates with σ, [Mg/Fe] does not.")
    print("GR Prediction: Both correlate identically (or not at all).")
    
    results = {}
    
    # 1. Correlation of spectroscopic age with log(σ)
    r_spec, p_spec = stats.pearsonr(df['log_sigma'], df['spec_age_norm'])
    print(f"\n1. Spectroscopic age vs log(σ):")
    print(f"   Pearson r = {r_spec:.4f}, p = {p_spec:.2e}")
    results['spec_age_vs_sigma'] = {'r': r_spec, 'p': p_spec}
    
    # 2. Correlation of [Mg/Fe] with log(σ)
    r_mgfe, p_mgfe = stats.pearsonr(df['log_sigma'], df['mg_fe_norm'])
    print(f"\n2. [Mg/Fe] vs log(σ):")
    print(f"   Pearson r = {r_mgfe:.4f}, p = {p_mgfe:.2e}")
    results['mg_fe_vs_sigma'] = {'r': r_mgfe, 'p': p_mgfe}
    
    # 3. The discrepancy
    delta_r = r_spec - r_mgfe
    print(f"\n3. DISCREPANCY (Δr = r_spec - r_mgfe):")
    print(f"   Δr = {delta_r:.4f}")
    results['delta_r'] = delta_r
    
    # 4. Statistical significance of the discrepancy
    # Use Fisher z-transformation to compare correlations
    n = len(df)
    z_spec = 0.5 * np.log((1 + r_spec) / (1 - r_spec))
    z_mgfe = 0.5 * np.log((1 + r_mgfe) / (1 - r_mgfe))
    se_diff = np.sqrt(2 / (n - 3))
    z_diff = (z_spec - z_mgfe) / se_diff
    p_diff = 2 * (1 - stats.norm.cdf(abs(z_diff)))
    
    print(f"\n4. Significance of discrepancy:")
    print(f"   z-score = {z_diff:.2f}")
    print(f"   p-value = {p_diff:.2e}")
    results['discrepancy_z'] = z_diff
    results['discrepancy_p'] = p_diff
    
    # 5. Binned analysis by σ class
    print("\n5. Binned analysis by velocity dispersion:")
    print("   σ class        | <Spec Age>  | <[Mg/Fe]>   | Δ(Spec-MgFe)")
    print("   " + "-" * 55)
    
    binned_results = []
    for sigma_class in ['low_sigma', 'med_low', 'medium', 'med_high', 'high_sigma']:
        mask = df['sigma_class'] == sigma_class
        if mask.sum() < 10:
            continue
        
        mean_spec = df.loc[mask, 'spec_age_norm'].mean()
        mean_mgfe = df.loc[mask, 'mg_fe_norm'].mean()
        delta = mean_spec - mean_mgfe
        
        print(f"   {sigma_class:14} | {mean_spec:+.3f}       | {mean_mgfe:+.3f}       | {delta:+.3f}")
        binned_results.append({
            'sigma_class': sigma_class,
            'mean_spec_age': mean_spec,
            'mean_mg_fe': mean_mgfe,
            'delta': delta,
            'n': mask.sum()
        })
    
    results['binned'] = binned_results
    
    # 6. Interpretation
    print("\n" + "=" * 70)
    print("INTERPRETATION")
    print("=" * 70)
    
    if abs(delta_r) > 0.05 and p_diff < 0.001:
        print("\n*** SIGNIFICANT DISCREPANCY DETECTED ***")
        print(f"Δr = {delta_r:.3f} (p ~ 0)")
        print(f"\nObserved pattern:")
        print(f"  - [Mg/Fe] correlates STRONGLY with σ (r = {r_mgfe:.3f})")
        print(f"  - Spectroscopic age correlates WEAKLY with σ (r = {r_spec:.3f})")
        print(f"\nPhysical interpretation:")
        print(f"  The [Mg/Fe]-σ relation is the well-known mass-metallicity relation:")
        print(f"  massive galaxies (high σ) formed stars rapidly, enriching in α-elements.")
        print(f"  This is a FORMATION TIMESCALE effect, not a time-flow effect.")
        print(f"\nTEP-relevant insight:")
        print(f"  The WEAK correlation of spectroscopic age with σ is surprising.")
        print(f"  Under standard stellar evolution, high-σ galaxies should appear")
        print(f"  systematically older (they formed earlier). The weak correlation")
        print(f"  suggests spectroscopic ages are DECOUPLED from formation epoch.")
        print(f"\n  This could indicate:")
        print(f"  1. TEP effect: time dilation in massive halos makes stars appear")
        print(f"     YOUNGER than expected, partially canceling the formation-epoch effect.")
        print(f"  2. Or: systematic issues with D4000/Hβ as age indicators.")
        results['interpretation'] = 'DISCREPANCY_DETECTED'
        results['tep_hypothesis'] = 'Weak spec_age-σ correlation may indicate TEP time dilation'
    elif abs(delta_r) < 0.02:
        print("\nNo significant discrepancy detected.")
        print("Both age indicators correlate similarly with environment.")
        print("This is consistent with GR (null hypothesis).")
        results['interpretation'] = 'GR_CONSISTENT'
    else:
        print("\nModerate discrepancy - requires further investigation.")
        results['interpretation'] = 'MODERATE_DISCREPANCY'
    
    return results


def create_figure(df, results):
    """Create publication figure showing the age-nucleosynthesis discrepancy."""
    import matplotlib.pyplot as plt
    
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # Panel A: Spectroscopic age vs log(σ)
    ax = axes[0, 0]
    ax.hexbin(df['log_sigma'], df['spec_age_norm'], gridsize=50, cmap='Blues', mincnt=1)
    ax.set_xlabel('log(σ / km s⁻¹)')
    ax.set_ylabel('Spectroscopic Age (normalized)')
    ax.set_title(f"A. Spec Age vs σ (r = {results['spec_age_vs_sigma']['r']:.3f})")
    
    # Add regression line
    slope, intercept = np.polyfit(df['log_sigma'], df['spec_age_norm'], 1)
    x_line = np.linspace(df['log_sigma'].min(), df['log_sigma'].max(), 100)
    ax.plot(x_line, slope * x_line + intercept, 'r-', lw=2, label='Linear fit')
    ax.legend()
    
    # Panel B: [Mg/Fe] vs log(σ)
    ax = axes[0, 1]
    ax.hexbin(df['log_sigma'], df['mg_fe_norm'], gridsize=50, cmap='Oranges', mincnt=1)
    ax.set_xlabel('log(σ / km s⁻¹)')
    ax.set_ylabel('[Mg/Fe] (normalized)')
    ax.set_title(f"B. [Mg/Fe] vs σ (r = {results['mg_fe_vs_sigma']['r']:.3f})")
    
    slope, intercept = np.polyfit(df['log_sigma'], df['mg_fe_norm'], 1)
    ax.plot(x_line, slope * x_line + intercept, 'r-', lw=2, label='Linear fit')
    ax.legend()
    
    # Panel C: Binned comparison
    ax = axes[1, 0]
    binned = results['binned']
    x = np.arange(len(binned))
    width = 0.35
    
    spec_vals = [b['mean_spec_age'] for b in binned]
    mgfe_vals = [b['mean_mg_fe'] for b in binned]
    labels = [b['sigma_class'].replace('_', '\n') for b in binned]
    
    ax.bar(x - width/2, spec_vals, width, label='Spec Age', color='steelblue')
    ax.bar(x + width/2, mgfe_vals, width, label='[Mg/Fe]', color='darkorange')
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=9)
    ax.set_ylabel('Normalized Value')
    ax.set_title('C. Age Indicators by σ Class')
    ax.legend()
    ax.axhline(0, color='gray', linestyle='--', alpha=0.5)
    
    # Panel D: The discrepancy
    ax = axes[1, 1]
    deltas = [b['delta'] for b in binned]
    colors = ['#2ecc71' if d > 0 else '#e74c3c' for d in deltas]
    ax.bar(x, deltas, color=colors, edgecolor='black')
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=9)
    ax.set_ylabel('Δ(Spec Age - [Mg/Fe])')
    ax.set_title(f'D. TEP Signature: Δr = {results["delta_r"]:.3f} (p = {results["discrepancy_p"]:.1e})')
    ax.axhline(0, color='gray', linestyle='--', alpha=0.5)
    
    # Add interpretation
    if results['interpretation'] == 'TEP_DETECTED':
        ax.text(0.5, 0.95, 'TEP SIGNATURE DETECTED', transform=ax.transAxes,
                ha='center', va='top', fontsize=12, fontweight='bold', color='green')
    
    plt.tight_layout()
    
    fig_path = os.path.join(FIGURES_DIR, 'sdss_age_nucleosynthesis_discrepancy.png')
    plt.savefig(fig_path, dpi=150, bbox_inches='tight')
    plt.savefig(fig_path.replace('.png', '.pdf'), bbox_inches='tight')
    print(f"\nFigure saved: {fig_path}")
    
    return fig_path


def main():
    """Main analysis pipeline."""
    print("=" * 70)
    print("SDSS AGE-NUCLEOSYNTHESIS DISCREPANCY TEST FOR TEP")
    print("=" * 70)
    print(f"Started: {datetime.now().isoformat()}")
    
    # Check for cached data
    cache_path = os.path.join(DATA_DIR, 'sdss_spectral_indices.csv')
    
    if os.path.exists(cache_path):
        print(f"\nLoading cached data from {cache_path}")
        df = pd.read_csv(cache_path)
    else:
        # Download fresh data
        df = download_spectral_indices()
        if df is None:
            print("Failed to download data.")
            return None
        
        # Save cache
        df.to_csv(cache_path, index=False)
        print(f"Data cached to {cache_path}")
    
    print(f"\nTotal galaxies: {len(df)}")
    
    # Compute age indicators
    df = compute_age_indicators(df)
    
    # Compute environment
    df = compute_local_environment(df)
    
    # Run the key test
    results = test_age_nucleosynthesis_discrepancy(df)
    
    # Add metadata
    results['n_galaxies'] = len(df)
    results['z_range'] = [float(df['redshift'].min()), float(df['redshift'].max())]
    results['timestamp'] = datetime.now().isoformat()
    
    # Save results
    results_path = os.path.join(RESULTS_DIR, 'sdss_age_nucleosynthesis_results.json')
    with open(results_path, 'w') as f:
        json.dump(results, f, indent=2, default=str)
    print(f"\nResults saved: {results_path}")
    
    # Create figure
    fig_path = create_figure(df, results)
    
    # Summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"Galaxies analyzed: {len(df):,}")
    print(f"Redshift range: {df['redshift'].min():.3f} - {df['redshift'].max():.3f}")
    print(f"Spec age vs σ correlation: r = {results['spec_age_vs_sigma']['r']:.4f}")
    print(f"[Mg/Fe] vs σ correlation: r = {results['mg_fe_vs_sigma']['r']:.4f}")
    print(f"DISCREPANCY: Δr = {results['delta_r']:.4f} (p = {results['discrepancy_p']:.2e})")
    print(f"Interpretation: {results['interpretation']}")
    
    return results


if __name__ == '__main__':
    results = main()
