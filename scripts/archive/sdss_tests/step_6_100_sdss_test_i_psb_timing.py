#!/usr/bin/env python3
"""
Test I: Post-Starburst Timing Anomaly

Hypothesis:
Post-starburst (E+A/K+A) galaxies show strong Balmer absorption (recent A-stars)
and weak emission (no current SF). Under TEP, if time flows slower in deeper 
potentials, PSB galaxies with high σ should appear younger post-burst.

TEP Prediction:
  r(HδA, σ) > 0   (stronger Balmer absorption at high σ = younger-appearing)
  
Standard Prediction:
  No correlation expected between HδA and σ at fixed time since burst.
"""

import os
import sys
import json
import requests
import numpy as np
import pandas as pd
from scipy import stats
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

def query_sdss(sql):
    """Query SDSS SkyServer."""
    params = {'cmd': sql, 'format': 'csv'}
    try:
        response = requests.get(SDSS_URL, params=params, timeout=300)
        response.raise_for_status()
        from io import StringIO
        df = pd.read_csv(StringIO(response.text))
        if len(df) == 0 or 'error' in df.columns[0].lower():
            return None
        return df
    except Exception as e:
        print(f"Query failed: {e}")
        return None


def download_psb_data():
    """Download post-starburst galaxy data."""
    cache_file = os.path.join(DATA_DIR, 'test_i_psb.csv')
    
    if os.path.exists(cache_file):
        print(f"Loading cached PSB data from {cache_file}")
        return pd.read_csv(cache_file)
    
    # Query for post-starburst galaxies
    # Simplified query to avoid HTTP 500
    sql = """
    SELECT TOP 50000
        g.specObjID,
        g.z AS redshift,
        g.velDisp AS sigma_stars,
        g.velDispErr AS sigma_stars_err,
        
        -- Balmer absorption indices
        i.lick_hd_a AS HdeltaA,
        i.lick_hd_a_err AS HdeltaA_err,
        i.lick_hg_a AS HgammaA,
        i.lick_hg_a_err AS HgammaA_err,
        i.lick_hb AS Hbeta,
        
        -- Age indicators
        i.d4000_n AS D4000,
        i.d4000_n_err AS D4000_err,
        
        -- Metallicity
        i.lick_mgb AS Mgb,
        i.lick_fe5270 AS Fe5270,
        
        -- Stellar mass
        s.logMass
        
    FROM galSpecInfo g
    JOIN galSpecIndx i ON g.specObjID = i.specObjID
    JOIN stellarMassFSPSGranWideDust s ON g.specObjID = s.specObjID
    
    WHERE 
        g.z BETWEEN 0.02 AND 0.20
        AND g.velDisp > 50 AND g.velDisp < 400
        AND g.velDispErr > 0 AND g.velDispErr < 30
        
        -- Strong Balmer absorption (PSB selection)
        AND i.lick_hd_a > 3.0
        AND i.lick_hd_a_err > 0 AND i.lick_hd_a_err < 2.0
        
        -- Quality
        AND i.d4000_n > 0 AND i.d4000_n_err < 0.1
        AND s.logMass > 9.0
    """
    
    print("Querying SDSS for post-starburst candidates...")
    df = query_sdss(sql)
    
    if df is not None and len(df) > 50:
        df.to_csv(cache_file, index=False)
        print(f"Saved {len(df)} PSB candidates to {cache_file}")
        return df
    
    print("Query failed, trying to construct from existing data...")
    return None


def load_from_existing_data():
    """Try to construct PSB sample from existing galaxy data."""
    twin_file = os.path.join(DATA_DIR, 'sdss_twin_base_sample_with_size.csv')
    
    if os.path.exists(twin_file):
        print(f"Loading from {twin_file}")
        df = pd.read_csv(twin_file)
        
        # Rename columns
        df = df.rename(columns={
            'specobjid': 'specObjID',
            'veldisp': 'sigma_stars',
            'veldisp_err': 'sigma_stars_err',
            'd4000': 'D4000',
            'hbeta': 'Hbeta',
            'log_mass': 'logMass'
        })
        
        # This dataset doesn't have HdeltaA, so we can't do PSB selection
        # But we can look at Hbeta as a proxy for recent star formation
        # Higher Hbeta = younger population
        
        print("Note: Using Hbeta as young population proxy (no HδA available)")
        return df, 'Hbeta'
    
    return None, None


def analyze_psb_timing(df, young_indicator='HdeltaA'):
    """Analyze post-starburst timing vs velocity dispersion."""
    results = {}
    
    if young_indicator not in df.columns:
        print(f"Column {young_indicator} not found")
        return {'verdict': 'Skipped', 'interpretation': f'No {young_indicator} data available'}
    
    # Clean data
    mask = (
        df[young_indicator].notna() &
        df['sigma_stars'].notna() &
        np.isfinite(df[young_indicator]) &
        np.isfinite(df['sigma_stars'])
    )
    clean = df[mask].copy()
    
    print(f"\nAnalyzing {len(clean)} galaxies with {young_indicator} data")
    results['n_galaxies'] = len(clean)
    results['indicator'] = young_indicator
    
    # 1. Raw correlation: young indicator vs σ
    r_raw, p_raw = stats.pearsonr(clean['sigma_stars'], clean[young_indicator])
    print(f"\n1. Raw correlation {young_indicator} vs σ:")
    print(f"   r = {r_raw:.4f}, p = {p_raw:.2e}")
    print(f"   TEP expects: r > 0 (high σ = younger-appearing)")
    
    results['raw_correlation'] = {
        'r': float(r_raw),
        'p': float(p_raw)
    }
    
    # 2. If we have D4000, control for it (D4000 is overall age proxy)
    if 'D4000' in clean.columns:
        from sklearn.linear_model import LinearRegression
        
        mask_d4000 = clean['D4000'].notna()
        subset = clean[mask_d4000]
        
        if len(subset) > 100:
            # Residualize young indicator against D4000
            X = subset[['D4000']].values
            y = subset[young_indicator].values
            
            reg = LinearRegression().fit(X, y)
            resid = y - reg.predict(X)
            
            # Correlation of residual with σ
            r_partial, p_partial = stats.pearsonr(subset['sigma_stars'], resid)
            print(f"\n2. Partial correlation {young_indicator} vs σ | D4000:")
            print(f"   r_partial = {r_partial:.4f}, p = {p_partial:.2e}")
            
            results['partial_correlation'] = {
                'r': float(r_partial),
                'p': float(p_partial),
                'controlled_for': 'D4000'
            }
    
    # 3. If we have mass, control for mass and D4000
    if 'logMass' in clean.columns and 'D4000' in clean.columns:
        mask_full = clean['logMass'].notna() & clean['D4000'].notna()
        subset = clean[mask_full]
        
        if len(subset) > 100:
            X = subset[['D4000', 'logMass']].values
            y = subset[young_indicator].values
            
            reg = LinearRegression().fit(X, y)
            resid = y - reg.predict(X)
            
            r_full, p_full = stats.pearsonr(subset['sigma_stars'], resid)
            print(f"\n3. Partial correlation {young_indicator} vs σ | D4000, M*:")
            print(f"   r_partial = {r_full:.4f}, p = {p_full:.2e}")
            
            results['partial_full'] = {
                'r': float(r_full),
                'p': float(p_full),
                'controlled_for': 'D4000, logMass'
            }
    
    # 4. Binned analysis by σ
    print(f"\n4. Mean {young_indicator} by σ bin:")
    sigma_bins = [(50, 100), (100, 150), (150, 200), (200, 300), (300, 400)]
    binned_results = []
    
    for lo, hi in sigma_bins:
        subset = clean[(clean['sigma_stars'] >= lo) & (clean['sigma_stars'] < hi)]
        if len(subset) > 20:
            mean_val = subset[young_indicator].mean()
            std_val = subset[young_indicator].std()
            sem_val = std_val / np.sqrt(len(subset))
            print(f"   σ = {lo}-{hi}: {young_indicator} = {mean_val:.3f} ± {sem_val:.3f}, n = {len(subset)}")
            binned_results.append({
                'sigma_range': f"{lo}-{hi}",
                'mean': float(mean_val),
                'stderr': float(sem_val),
                'n': len(subset)
            })
    
    results['binned_by_sigma'] = binned_results
    
    # 5. Linear regression slope
    slope, intercept, r_val, p_slope, stderr = stats.linregress(
        np.log10(clean['sigma_stars']), clean[young_indicator]
    )
    print(f"\n5. Linear fit: {young_indicator} = {slope:.3f} × log(σ) + {intercept:.3f}")
    print(f"   slope = {slope:.4f} ± {stderr:.4f}, p = {p_slope:.2e}")
    
    results['linear_fit'] = {
        'slope': float(slope),
        'intercept': float(intercept),
        'slope_err': float(stderr),
        'p_value': float(p_slope)
    }
    
    # Verdict
    # For HδA: positive correlation with σ is TEP-consistent
    # For Hβ: negative correlation (higher Hβ = younger) - need to think about sign
    if young_indicator == 'HdeltaA':
        if r_raw > 0 and p_raw < 0.05:
            results['verdict'] = 'Signal'
            results['interpretation'] = ('Stronger Balmer absorption at high σ. '
                                        'TEP-consistent: deeper potential = younger appearance.')
        elif r_raw < 0 and p_raw < 0.05:
            results['verdict'] = 'Contradicted'
            results['interpretation'] = ('Weaker Balmer absorption at high σ. '
                                        'Contradicts TEP prediction.')
        else:
            results['verdict'] = 'Null'
            results['interpretation'] = 'No significant correlation between Balmer absorption and σ.'
    else:  # Hbeta proxy
        # Higher Hbeta = younger stars, so negative correlation with σ would be TEP-consistent
        # (high σ = younger = higher Hbeta... wait, we expect them to LOOK younger)
        # Actually positive correlation: high σ → younger appearance → higher Hbeta
        if r_raw > 0 and p_raw < 0.05:
            results['verdict'] = 'Signal'
            results['interpretation'] = (f'Higher {young_indicator} at high σ. '
                                        'TEP-consistent: deeper potential = younger appearance.')
        elif r_raw < 0 and p_raw < 0.05:
            results['verdict'] = 'Contradicted'
            results['interpretation'] = (f'Lower {young_indicator} at high σ. '
                                        'Standard aging, contradicts TEP.')
        else:
            results['verdict'] = 'Null'
            results['interpretation'] = f'No significant correlation between {young_indicator} and σ.'
    
    print(f"\n=== VERDICT: {results['verdict']} ===")
    print(results['interpretation'])
    
    return results, clean


def create_figure(df, results, young_indicator):
    """Create visualization."""
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    mask = df[young_indicator].notna() & df['sigma_stars'].notna()
    clean = df[mask]
    
    # 1. Scatter plot
    ax1 = axes[0, 0]
    if 'D4000' in clean.columns:
        scatter = ax1.scatter(clean['sigma_stars'], clean[young_indicator],
                             c=clean['D4000'], cmap='viridis', alpha=0.3, s=10)
        plt.colorbar(scatter, ax=ax1, label='D4000')
    else:
        ax1.scatter(clean['sigma_stars'], clean[young_indicator], alpha=0.3, s=10)
    
    # Add regression line
    x_fit = np.linspace(clean['sigma_stars'].min(), clean['sigma_stars'].max(), 100)
    slope = results['linear_fit']['slope']
    intercept = results['linear_fit']['intercept']
    y_fit = slope * np.log10(x_fit) + intercept
    ax1.plot(x_fit, y_fit, 'r-', linewidth=2, label=f'Slope = {slope:.3f}')
    
    ax1.set_xlabel('σ (km/s)')
    ax1.set_ylabel(young_indicator)
    ax1.set_title(f'{young_indicator} vs Velocity Dispersion')
    ax1.legend()
    
    # 2. Binned means
    ax2 = axes[0, 1]
    if results.get('binned_by_sigma'):
        x_centers = []
        y_means = []
        y_errs = []
        for b in results['binned_by_sigma']:
            lo, hi = map(int, b['sigma_range'].split('-'))
            x_centers.append((lo + hi) / 2)
            y_means.append(b['mean'])
            y_errs.append(b['stderr'])
        
        ax2.errorbar(x_centers, y_means, yerr=y_errs, fmt='o-', capsize=5, markersize=8)
        ax2.set_xlabel('σ (km/s)')
        ax2.set_ylabel(f'Mean {young_indicator}')
        ax2.set_title(f'Binned {young_indicator} vs σ')
    
    # 3. Distribution of young indicator
    ax3 = axes[1, 0]
    ax3.hist(clean[young_indicator], bins=50, edgecolor='black', alpha=0.7)
    ax3.axvline(clean[young_indicator].mean(), color='red', linestyle='--',
               label=f'Mean = {clean[young_indicator].mean():.2f}')
    ax3.set_xlabel(young_indicator)
    ax3.set_ylabel('Count')
    ax3.set_title(f'Distribution of {young_indicator}')
    ax3.legend()
    
    # 4. Summary
    ax4 = axes[1, 1]
    ax4.axis('off')
    
    summary = f"""
TEST I: Post-Starburst Timing Anomaly

HYPOTHESIS:
Under TEP, post-starburst galaxies with high σ
should appear younger (stronger Balmer absorption).

TEP Prediction: r({young_indicator}, σ) > 0

RESULTS:
• N = {results['n_galaxies']:,} galaxies
• Raw r = {results['raw_correlation']['r']:.4f}
  p = {results['raw_correlation']['p']:.2e}
• Linear slope = {results['linear_fit']['slope']:.4f} ± {results['linear_fit']['slope_err']:.4f}

VERDICT: {results['verdict']}
{results['interpretation'][:80]}...
"""
    ax4.text(0.05, 0.95, summary, transform=ax4.transAxes,
            fontsize=10, verticalalignment='top', fontfamily='monospace',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    plt.tight_layout()
    fig_path = os.path.join(FIGURE_DIR, 'test_i_psb_timing.png')
    plt.savefig(fig_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"\nFigure saved to {fig_path}")
    
    return fig_path


def main():
    print("=" * 60)
    print("TEST I: Post-Starburst Timing Anomaly")
    print("=" * 60)
    
    # Try to download PSB data
    df = download_psb_data()
    young_indicator = 'HdeltaA'
    
    # Fallback to existing data
    if df is None or len(df) < 50:
        print("PSB query failed, using existing data with Hbeta proxy...")
        df, young_indicator = load_from_existing_data()
    
    if df is None or len(df) < 50:
        print("ERROR: Insufficient data for analysis")
        results = {
            'verdict': 'Skipped',
            'interpretation': 'Could not retrieve post-starburst data.'
        }
        output_file = os.path.join(OUTPUT_DIR, 'sdss_test_i_psb_timing_results.json')
        with open(output_file, 'w') as f:
            json.dump(results, f, indent=2)
        return results
    
    print(f"\nLoaded {len(df)} galaxies, using {young_indicator} as young population indicator")
    
    # Analyze
    results, clean_df = analyze_psb_timing(df, young_indicator)
    
    # Create figure
    fig_path = create_figure(clean_df, results, young_indicator)
    results['figure_path'] = fig_path
    
    # Save results
    output_file = os.path.join(OUTPUT_DIR, 'sdss_test_i_psb_timing_results.json')
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\nResults saved to {output_file}")
    
    return results


if __name__ == '__main__':
    results = main()
