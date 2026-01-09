#!/usr/bin/env python3
"""
Test K: Galaxy Size-Age Degeneracy Breaking

Hypothesis:
Under TEP, compact galaxies at fixed [Mg/Fe] should appear YOUNGER (lower D4000)
than extended galaxies. This is OPPOSITE to the standard expectation where
compact galaxies are thought to be older (formed at high redshift).

TEP Prediction:
At fixed [Mg/Fe] (formation timescale):
    r(D4000, Compactness) < 0     (compact galaxies appear younger)

This is because compact → deeper central potential → slower time flow → younger appearance
"""

import os
import sys
import json
import requests
import numpy as np
import pandas as pd
from scipy import stats
from astropy.cosmology import Planck18 as cosmo
import matplotlib.pyplot as plt

# Paths
BASE_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATA_DIR = os.path.join(BASE_DIR, 'data', 'sdss')
OUTPUT_DIR = os.path.join(BASE_DIR, 'results', 'outputs')
FIGURE_DIR = os.path.join(BASE_DIR, 'results', 'figures')

os.makedirs(DATA_DIR, exist_ok=True)
os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(FIGURE_DIR, exist_ok=True)

SDSS_URL = "https://skyserver.sdss.org/dr18/SkyServerWS/SearchTools/SqlSearch"

def query_sdss(sql, max_rows=500000):
    """Query SDSS SkyServer."""
    params = {
        'cmd': sql,
        'format': 'csv'
    }
    try:
        response = requests.get(SDSS_URL, params=params, timeout=300)
        response.raise_for_status()
        
        from io import StringIO
        df = pd.read_csv(StringIO(response.text))
        
        if len(df) == 0 or 'error' in df.columns[0].lower():
            print(f"Query returned error or empty: {df.head()}")
            return None
        
        return df
    except Exception as e:
        print(f"Query failed: {e}")
        return None


def download_data():
    """Download size-age data from SDSS or use cached data."""
    # First check for existing cached data with size measurements
    twin_cache = os.path.join(DATA_DIR, 'sdss_twin_base_sample_with_size.csv')
    
    if os.path.exists(twin_cache):
        print(f"Loading cached data from {twin_cache}")
        df = pd.read_csv(twin_cache)
        
        # Rename columns to match expected format
        column_map = {
            'specobjid': 'specObjID',
            'veldisp': 'sigma_stars',
            'veldisp_err': 'sigma_stars_err',
            'mgb': 'Mgb',
            'fe5270': 'Fe5270',
            'fe5335': 'Fe5335',
            'd4000': 'D4000',
            'hbeta': 'Hbeta',
            'log_mass': 'logMass',
            'petroR50_r_arcsec': 'Re_arcsec'
        }
        df = df.rename(columns=column_map)
        
        # Filter for quality
        mask = (
            (df['sigma_stars'] > 50) & (df['sigma_stars'] < 400) &
            (df['D4000'] > 1.3) &  # Include broader range for comparison
            (df['Mgb'] > 0) &
            (df['Fe5270'] > 0) &
            (df['Fe5335'] > 0) &
            (df['Re_arcsec'] > 0.5) & (df['Re_arcsec'] < 50) &
            (df['logMass'] > 9.0) & (df['logMass'] < 12.5)
        )
        df = df[mask].copy()
        print(f"Using {len(df)} galaxies after quality cuts")
        return df
    
    # Fallback: try SDSS query
    cache_file = os.path.join(DATA_DIR, 'test_k_size_age.csv')
    if os.path.exists(cache_file):
        print(f"Loading cached data from {cache_file}")
        return pd.read_csv(cache_file)
    
    # Simplified query
    sql = """
    SELECT TOP 50000
        g.specObjID,
        g.z AS redshift,
        g.velDisp AS sigma_stars,
        g.velDispErr AS sigma_stars_err,
        i.d4000_n AS D4000,
        i.lick_mgb AS Mgb,
        i.lick_fe5270 AS Fe5270,
        i.lick_fe5335 AS Fe5335,
        s.logMass
        
    FROM galSpecInfo g
    JOIN galSpecIndx i ON g.specObjID = i.specObjID
    JOIN stellarMassFSPSGranWideDust s ON g.specObjID = s.specObjID
    
    WHERE 
        g.z BETWEEN 0.02 AND 0.15
        AND g.velDisp > 80 AND g.velDisp < 400
        AND i.d4000_n > 1.5
        AND i.lick_mgb > 0
        AND s.logMass > 9.5
    """
    
    print("Querying SDSS for size-age data...")
    df = query_sdss(sql)
    
    if df is not None and len(df) > 100:
        df.to_csv(cache_file, index=False)
        print(f"Saved {len(df)} galaxies to {cache_file}")
        return df
    else:
        print("Query failed or returned insufficient data")
        return None


def compute_physical_quantities(df):
    """Compute physical sizes and derived quantities."""
    # Angular diameter distance (Mpc)
    z = df['redshift'].values
    D_A = cosmo.angular_diameter_distance(z).value  # Mpc
    
    # Physical half-light radius (kpc)
    # R_e (arcsec) → R_e (kpc) = R_e (arcsec) × D_A (Mpc) × 1000 / 206265
    arcsec_to_kpc = D_A * 1000 / 206265
    df['Re_kpc'] = df['Re_arcsec'] * arcsec_to_kpc
    
    # Compactness = log(M*) - 2×log(R_e) [log surface mass density]
    df['Compactness'] = df['logMass'] - 2 * np.log10(df['Re_kpc'])
    
    # [Mg/Fe] proxy: log10(Mgb / <Fe>)
    # <Fe> = 0.5 × (Fe5270 + Fe5335)
    df['Fe_avg'] = 0.5 * (df['Fe5270'] + df['Fe5335'])
    # Avoid log of negative values
    valid_idx = (df['Mgb'] > 0) & (df['Fe_avg'] > 0)
    df['MgFe'] = np.nan
    df.loc[valid_idx, 'MgFe'] = np.log10(df.loc[valid_idx, 'Mgb'] / df.loc[valid_idx, 'Fe_avg'])
    
    # Combined index [MgFe]' = sqrt(Mgb × <Fe>)
    df['MgFe_prime'] = np.sqrt(df['Mgb'] * df['Fe_avg'])
    
    return df


def analyze_size_age_correlation(df):
    """Analyze the size-age correlation at fixed [Mg/Fe]."""
    results = {}
    
    # Clean data
    mask = (
        df['D4000'].notna() & 
        df['Compactness'].notna() & 
        df['MgFe'].notna() &
        df['logMass'].notna() &
        np.isfinite(df['D4000']) &
        np.isfinite(df['Compactness']) &
        np.isfinite(df['MgFe'])
    )
    clean = df[mask].copy()
    
    print(f"\nAnalyzing {len(clean)} galaxies with valid data")
    results['n_galaxies'] = len(clean)
    
    # 1. Raw correlation: D4000 vs Compactness
    r_raw, p_raw = stats.pearsonr(clean['D4000'], clean['Compactness'])
    print(f"\n1. Raw correlation D4000 vs Compactness:")
    print(f"   r = {r_raw:.4f}, p = {p_raw:.2e}")
    print(f"   Standard expectation: r > 0 (compact = old)")
    print(f"   TEP expectation: r < 0 (compact = young appearance)")
    
    results['raw_correlation'] = {
        'r': float(r_raw),
        'p': float(p_raw),
        'interpretation': 'TEP-consistent' if r_raw < 0 else 'Standard'
    }
    
    # 2. Partial correlation controlling for [Mg/Fe]
    # D4000 vs Compactness | [Mg/Fe]
    from scipy.stats import pearsonr
    
    # Residualize D4000 against [Mg/Fe]
    slope_age, intercept_age, _, _, _ = stats.linregress(clean['MgFe'], clean['D4000'])
    clean['D4000_resid'] = clean['D4000'] - (slope_age * clean['MgFe'] + intercept_age)
    
    # Residualize Compactness against [Mg/Fe]
    slope_comp, intercept_comp, _, _, _ = stats.linregress(clean['MgFe'], clean['Compactness'])
    clean['Compactness_resid'] = clean['Compactness'] - (slope_comp * clean['MgFe'] + intercept_comp)
    
    r_partial, p_partial = pearsonr(clean['D4000_resid'], clean['Compactness_resid'])
    print(f"\n2. Partial correlation D4000 vs Compactness | [Mg/Fe]:")
    print(f"   r_partial = {r_partial:.4f}, p = {p_partial:.2e}")
    
    results['partial_correlation_MgFe'] = {
        'r': float(r_partial),
        'p': float(p_partial),
        'controlled_for': '[Mg/Fe]'
    }
    
    # 3. Control for mass AND [Mg/Fe]
    from sklearn.linear_model import LinearRegression
    
    X_controls = clean[['MgFe', 'logMass']].values
    
    # Residualize D4000
    reg_d4000 = LinearRegression().fit(X_controls, clean['D4000'])
    clean['D4000_resid_full'] = clean['D4000'] - reg_d4000.predict(X_controls)
    
    # Residualize Compactness
    reg_comp = LinearRegression().fit(X_controls, clean['Compactness'])
    clean['Compactness_resid_full'] = clean['Compactness'] - reg_comp.predict(X_controls)
    
    r_full, p_full = pearsonr(clean['D4000_resid_full'], clean['Compactness_resid_full'])
    print(f"\n3. Partial correlation D4000 vs Compactness | [Mg/Fe], M*:")
    print(f"   r_partial = {r_full:.4f}, p = {p_full:.2e}")
    
    results['partial_correlation_full'] = {
        'r': float(r_full),
        'p': float(p_full),
        'controlled_for': '[Mg/Fe], logMass'
    }
    
    # 4. Binned analysis: D4000 vs Compactness in [Mg/Fe] bins
    print("\n4. Binned analysis by [Mg/Fe]:")
    mgfe_bins = np.percentile(clean['MgFe'], [0, 33, 67, 100])
    mgfe_labels = ['Low [Mg/Fe]', 'Mid [Mg/Fe]', 'High [Mg/Fe]']
    clean['MgFe_bin'] = pd.cut(clean['MgFe'], bins=mgfe_bins, labels=mgfe_labels)
    
    binned_results = []
    for label in mgfe_labels:
        subset = clean[clean['MgFe_bin'] == label]
        if len(subset) > 30:
            r_bin, p_bin = pearsonr(subset['D4000'], subset['Compactness'])
            print(f"   {label}: r = {r_bin:.4f}, p = {p_bin:.2e}, n = {len(subset)}")
            binned_results.append({
                'bin': label,
                'r': float(r_bin),
                'p': float(p_bin),
                'n': len(subset)
            })
    
    results['binned_by_MgFe'] = binned_results
    
    # 5. Compare compact vs extended populations
    compactness_median = clean['Compactness'].median()
    compact = clean[clean['Compactness'] > compactness_median]
    extended = clean[clean['Compactness'] <= compactness_median]
    
    # At fixed [Mg/Fe], compare D4000
    mgfe_overlap = (
        (clean['MgFe'] > np.percentile(compact['MgFe'], 25)) &
        (clean['MgFe'] < np.percentile(extended['MgFe'], 75))
    )
    
    compact_overlap = compact[compact['MgFe'].between(
        np.percentile(compact['MgFe'], 25),
        np.percentile(compact['MgFe'], 75)
    )]
    extended_overlap = extended[extended['MgFe'].between(
        np.percentile(compact['MgFe'], 25),
        np.percentile(compact['MgFe'], 75)
    )]
    
    if len(compact_overlap) > 30 and len(extended_overlap) > 30:
        d4000_diff = compact_overlap['D4000'].mean() - extended_overlap['D4000'].mean()
        t_stat, p_ttest = stats.ttest_ind(compact_overlap['D4000'], extended_overlap['D4000'])
        
        print(f"\n5. Compact vs Extended at fixed [Mg/Fe]:")
        print(f"   ΔD4000 (compact - extended) = {d4000_diff:.4f}")
        print(f"   t = {t_stat:.2f}, p = {p_ttest:.2e}")
        print(f"   Standard expects: ΔD4000 > 0 (compact older)")
        print(f"   TEP expects: ΔD4000 < 0 (compact appears younger)")
        
        results['compact_vs_extended'] = {
            'delta_D4000': float(d4000_diff),
            't_stat': float(t_stat),
            'p_value': float(p_ttest),
            'n_compact': len(compact_overlap),
            'n_extended': len(extended_overlap),
            'interpretation': 'TEP-consistent' if d4000_diff < 0 else 'Standard'
        }
    
    # 6. Sigma bins for TEP gradient
    print("\n6. D4000 residual vs velocity dispersion:")
    sigma_bins = [(80, 150), (150, 200), (200, 250), (250, 400)]
    sigma_results = []
    for lo, hi in sigma_bins:
        subset = clean[(clean['sigma_stars'] >= lo) & (clean['sigma_stars'] < hi)]
        if len(subset) > 50:
            mean_resid = subset['D4000_resid'].mean()
            std_resid = subset['D4000_resid'].std() / np.sqrt(len(subset))
            print(f"   σ = {lo}-{hi}: ΔD4000_resid = {mean_resid:.4f} ± {std_resid:.4f}, n = {len(subset)}")
            sigma_results.append({
                'sigma_range': f"{lo}-{hi}",
                'mean_D4000_resid': float(mean_resid),
                'stderr': float(std_resid),
                'n': len(subset)
            })
    
    results['by_sigma'] = sigma_results
    
    # Overall interpretation
    if r_partial < 0 and p_partial < 0.05:
        results['verdict'] = 'Signal'
        results['interpretation'] = ('Compact galaxies appear YOUNGER at fixed [Mg/Fe]. '
                                    'This is TEP-consistent: deeper potential = slower time flow = younger appearance.')
    elif r_partial > 0 and p_partial < 0.05:
        results['verdict'] = 'Contradicted'
        results['interpretation'] = ('Compact galaxies appear OLDER at fixed [Mg/Fe]. '
                                    'This follows standard expectation, contradicting TEP.')
    else:
        results['verdict'] = 'Null'
        results['interpretation'] = 'No significant correlation between compactness and age at fixed [Mg/Fe].'
    
    print(f"\n=== VERDICT: {results['verdict']} ===")
    print(results['interpretation'])
    
    return results, clean


def create_figure(df, results):
    """Create visualization of size-age analysis."""
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # Clean data
    mask = df['D4000'].notna() & df['Compactness'].notna() & df['MgFe'].notna()
    clean = df[mask]
    
    # 1. Raw D4000 vs Compactness colored by [Mg/Fe]
    ax1 = axes[0, 0]
    scatter = ax1.scatter(clean['Compactness'], clean['D4000'], 
                         c=clean['MgFe'], cmap='viridis', alpha=0.3, s=5)
    plt.colorbar(scatter, ax=ax1, label='[Mg/Fe]')
    ax1.set_xlabel('Compactness (log Σ*)')
    ax1.set_ylabel('D4000 (Age indicator)')
    ax1.set_title(f"Raw: r = {results['raw_correlation']['r']:.3f}")
    
    # 2. D4000 residual vs Compactness residual (after controlling for [Mg/Fe])
    ax2 = axes[0, 1]
    if 'D4000_resid' in clean.columns:
        ax2.scatter(clean['Compactness_resid'], clean['D4000_resid'], 
                   c=clean['sigma_stars'], cmap='plasma', alpha=0.3, s=5)
        ax2.axhline(0, color='gray', linestyle='--', alpha=0.5)
        ax2.axvline(0, color='gray', linestyle='--', alpha=0.5)
        ax2.set_xlabel('Compactness residual | [Mg/Fe]')
        ax2.set_ylabel('D4000 residual | [Mg/Fe]')
        r_part = results['partial_correlation_MgFe']['r']
        ax2.set_title(f"Partial (|[Mg/Fe]): r = {r_part:.3f}")
    
    # 3. Binned by [Mg/Fe]: D4000 vs Compactness
    ax3 = axes[1, 0]
    colors = ['blue', 'green', 'red']
    for i, (bin_result, color) in enumerate(zip(results.get('binned_by_MgFe', []), colors)):
        label = bin_result['bin']
        subset = clean[clean['MgFe_bin'] == label] if 'MgFe_bin' in clean.columns else clean
        if len(subset) > 0:
            # Bin means
            comp_bins = np.percentile(subset['Compactness'], np.linspace(10, 90, 6))
            comp_centers = 0.5 * (comp_bins[:-1] + comp_bins[1:])
            d4000_means = []
            d4000_errs = []
            for j in range(len(comp_bins) - 1):
                in_bin = subset[(subset['Compactness'] >= comp_bins[j]) & 
                               (subset['Compactness'] < comp_bins[j+1])]
                if len(in_bin) > 5:
                    d4000_means.append(in_bin['D4000'].mean())
                    d4000_errs.append(in_bin['D4000'].std() / np.sqrt(len(in_bin)))
                else:
                    d4000_means.append(np.nan)
                    d4000_errs.append(np.nan)
            
            ax3.errorbar(comp_centers, d4000_means, yerr=d4000_errs, 
                        fmt='o-', color=color, label=f"{label} (r={bin_result['r']:.2f})",
                        capsize=3)
    
    ax3.set_xlabel('Compactness (log Σ*)')
    ax3.set_ylabel('D4000')
    ax3.legend(fontsize=8)
    ax3.set_title('D4000 vs Compactness by [Mg/Fe] bin')
    
    # 4. Summary panel
    ax4 = axes[1, 1]
    ax4.axis('off')
    
    summary_text = f"""
TEST K: Size-Age Degeneracy Breaking

HYPOTHESIS:
Under TEP, compact galaxies should appear
YOUNGER at fixed [Mg/Fe] (opposite to standard).

RESULTS:
• Sample: {results['n_galaxies']:,} early-type galaxies
• Raw r(D4000, Compactness): {results['raw_correlation']['r']:.4f}
• Partial r | [Mg/Fe]: {results['partial_correlation_MgFe']['r']:.4f}
  p = {results['partial_correlation_MgFe']['p']:.2e}
• Partial r | [Mg/Fe], M*: {results['partial_correlation_full']['r']:.4f}
  p = {results['partial_correlation_full']['p']:.2e}

INTERPRETATION:
Standard: r > 0 (compact = old)
TEP: r < 0 (compact = young appearance)

VERDICT: {results['verdict']}
{results['interpretation'][:100]}...
"""
    ax4.text(0.05, 0.95, summary_text, transform=ax4.transAxes, 
             fontsize=10, verticalalignment='top', fontfamily='monospace',
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    plt.tight_layout()
    fig_path = os.path.join(FIGURE_DIR, 'test_k_size_age_degeneracy.png')
    plt.savefig(fig_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"\nFigure saved to {fig_path}")
    
    return fig_path


def main():
    print("=" * 60)
    print("TEST K: Galaxy Size-Age Degeneracy Breaking")
    print("=" * 60)
    
    # Download data
    df = download_data()
    
    if df is None or len(df) < 100:
        print("ERROR: Insufficient data for analysis")
        return None
    
    print(f"\nLoaded {len(df)} galaxies")
    
    # Compute physical quantities
    df = compute_physical_quantities(df)
    
    # Analyze correlations
    results, clean_df = analyze_size_age_correlation(df)
    
    # Create figure
    fig_path = create_figure(clean_df, results)
    results['figure_path'] = fig_path
    
    # Save results
    output_file = os.path.join(OUTPUT_DIR, 'sdss_test_k_size_age_results.json')
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\nResults saved to {output_file}")
    
    return results


if __name__ == '__main__':
    results = main()
