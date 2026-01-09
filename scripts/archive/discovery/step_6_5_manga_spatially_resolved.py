#!/usr/bin/env python3
"""
Step 6.5: MaNGA Spatially-Resolved Age Gradient Test for TEP

MaNGA provides IFU spectroscopy of ~10,000 galaxies, allowing us to measure
age and abundance GRADIENTS within individual galaxies.

TEP PREDICTION:
If time flows slower in deeper gravitational potentials, then:
- Galaxy centers (deeper Φ) should appear YOUNGER than outskirts
- But [α/Fe] gradients should be INDEPENDENT of local potential
- The age gradient should correlate with the potential gradient

This is a powerful test because it controls for galaxy-to-galaxy variations
by looking at gradients WITHIN each galaxy.

DATA:
- mangaPipe3D: Stellar population fits from Pipe3D
- Key columns: Age_LW_Re_fit, alpha_Age_LW_Re_fit (age gradient)
- vel_sigma_Re: Central velocity dispersion

Author: M. Smawfield
Date: January 2026
"""

import requests
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
import json
import os
from datetime import datetime

SDSS_URL = "https://skyserver.sdss.org/dr18/SkyServerWS/SearchTools/SqlSearch"
DATA_DIR = os.path.join(os.path.dirname(__file__), '..', '..', 'data', 'manga')
RESULTS_DIR = os.path.join(os.path.dirname(__file__), '..', '..', 'results', 'outputs')
FIGURES_DIR = os.path.join(os.path.dirname(__file__), '..', '..', 'results', 'figures')

os.makedirs(DATA_DIR, exist_ok=True)


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
                if response.status_code == 500:
                    print(f"  Error: {response.text[:200]}")
        except Exception as e:
            print(f"  Error: {e}")
        if attempt < max_retries - 1:
            time.sleep(2)
    return None


def download_manga_data():
    """Download MaNGA Pipe3D data with age gradients."""
    print("\n" + "=" * 70)
    print("DOWNLOADING MaNGA PIPE3D DATA")
    print("=" * 70)
    
    # Simplified query - no WHERE clause
    sql = """
    SELECT TOP 10000
        plateifu,
        mangaid,
        objra, objdec,
        log_Mass,
        Age_LW_Re_fit,
        alpha_Age_LW_Re_fit,
        Age_MW_Re_fit,
        alpha_Age_MW_Re_fit,
        ZH_LW_Re_fit,
        alpha_ZH_LW_Re_fit,
        vel_sigma_Re,
        Re_kpc
    FROM mangaPipe3D
    """
    
    print("Querying mangaPipe3D...")
    df = query_sdss(sql)
    
    if df is None or len(df) == 0:
        print("ERROR: No data from mangaPipe3D")
        return None
    
    print(f"Retrieved {len(df)} galaxies")
    return df


def prepare_data(df):
    """Prepare MaNGA data for analysis."""
    print("\nPreparing data...")
    
    # Convert to numeric
    numeric_cols = ['log_Mass', 'Age_LW_Re_fit', 'alpha_Age_LW_Re_fit',
                    'Age_MW_Re_fit', 'alpha_Age_MW_Re_fit',
                    'ZH_LW_Re_fit', 'alpha_ZH_LW_Re_fit',
                    'vel_sigma_Re', 'Re_kpc']
    
    for col in numeric_cols:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors='coerce')
    
    # Replace sentinel values (-9999) with NaN
    df = df.replace(-9999, np.nan)
    df = df.replace(-9999.0, np.nan)
    
    # vel_sigma_Re appears to be in log units (values ~0.3-0.5)
    # Convert to linear: σ = 10^vel_sigma_Re * 100 km/s (approximate)
    # Actually looking at range 0.2-1.0, this is likely log10(σ/100)
    # So σ = 10^vel_sigma_Re * 100
    df['sigma_linear'] = 10**df['vel_sigma_Re'] * 100  # km/s
    
    print(f"σ range after conversion: {df['sigma_linear'].min():.0f} - {df['sigma_linear'].max():.0f} km/s")
    
    # Filter valid - Age values ~9 are log(age/yr), so 10^9 = 1 Gyr
    # Age_LW_Re_fit appears to be log10(age/yr)
    valid = (
        np.isfinite(df['Age_LW_Re_fit']) &
        np.isfinite(df['alpha_Age_LW_Re_fit']) &
        np.isfinite(df['sigma_linear']) &
        (df['Age_LW_Re_fit'] > 8) &  # > 100 Myr
        (df['Age_LW_Re_fit'] < 10.2) &  # < 15 Gyr
        (df['sigma_linear'] > 30) &
        (df['sigma_linear'] < 400)
    )
    
    df_valid = df[valid].copy()
    print(f"Valid galaxies: {len(df_valid)}")
    
    # Convert age to Gyr
    df_valid['Age_Gyr'] = 10**(df_valid['Age_LW_Re_fit'] - 9)  # Gyr
    
    # Use sigma_linear for analysis
    df_valid['vel_sigma_Re'] = df_valid['sigma_linear']
    df_valid['log_sigma'] = np.log10(df_valid['vel_sigma_Re'])
    
    # Potential proxy
    c = 299792.458  # km/s
    df_valid['phi_c2'] = (df_valid['vel_sigma_Re'])**2 / c**2
    
    return df_valid


def test_age_gradients(df):
    """
    THE KEY TEP TEST:
    
    Age gradients within galaxies should correlate with potential gradients.
    
    TEP prediction:
    - Negative age gradient (younger center) at higher σ
    - Metallicity gradient should NOT correlate with σ
    
    The gradient α is defined as: property(R) = property(Re) + α × log(R/Re)
    Negative α means the property DECREASES toward the center.
    """
    print("\n" + "=" * 70)
    print("MaNGA AGE GRADIENT TEST")
    print("=" * 70)
    
    results = {}
    
    # 1. Age gradient vs σ
    print("\n1. Age gradient (α_Age) vs log(σ):")
    r_age_grad, p_age = stats.pearsonr(df['log_sigma'], df['alpha_Age_LW_Re_fit'])
    print(f"   r = {r_age_grad:.4f}, p = {p_age:.2e}")
    results['age_gradient_vs_sigma'] = {'r': r_age_grad, 'p': p_age}
    
    # 2. Metallicity gradient vs σ
    print("\n2. Metallicity gradient (α_ZH) vs log(σ):")
    valid_zh = np.isfinite(df['alpha_ZH_LW_Re_fit'])
    if valid_zh.sum() > 100:
        r_zh_grad, p_zh = stats.pearsonr(df.loc[valid_zh, 'log_sigma'], 
                                          df.loc[valid_zh, 'alpha_ZH_LW_Re_fit'])
        print(f"   r = {r_zh_grad:.4f}, p = {p_zh:.2e}")
        results['zh_gradient_vs_sigma'] = {'r': r_zh_grad, 'p': p_zh}
    else:
        print("   Insufficient data")
        r_zh_grad = 0
    
    # 3. Central age vs σ
    print("\n3. Central age (Age_Re) vs log(σ):")
    r_age_central, p_age_c = stats.pearsonr(df['log_sigma'], df['Age_LW_Re_fit'])
    print(f"   r = {r_age_central:.4f}, p = {p_age_c:.2e}")
    results['central_age_vs_sigma'] = {'r': r_age_central, 'p': p_age_c}
    
    # 4. Discrepancy
    delta_r = r_age_grad - r_zh_grad
    print(f"\n4. DISCREPANCY (Δr = r_age_grad - r_zh_grad):")
    print(f"   Δr = {delta_r:.4f}")
    results['delta_r'] = delta_r
    
    # 5. Binned analysis
    print("\n5. Binned analysis by σ:")
    print("-" * 70)
    print(f"{'σ bin':<15} {'N':>6} {'<α_Age>':>12} {'<α_ZH>':>12} {'<Age_Re>':>12}")
    print("-" * 70)
    
    sigma_bins = pd.qcut(df['vel_sigma_Re'], q=5, labels=['low', 'med_low', 'medium', 'med_high', 'high'])
    df['sigma_bin'] = sigma_bins
    
    binned_results = []
    for sbin in ['low', 'med_low', 'medium', 'med_high', 'high']:
        mask = df['sigma_bin'] == sbin
        if mask.sum() < 10:
            continue
        
        mean_alpha_age = df.loc[mask, 'alpha_Age_LW_Re_fit'].mean()
        mean_alpha_zh = df.loc[mask & np.isfinite(df['alpha_ZH_LW_Re_fit']), 'alpha_ZH_LW_Re_fit'].mean()
        mean_age = df.loc[mask, 'Age_LW_Re_fit'].mean()
        mean_sigma = df.loc[mask, 'vel_sigma_Re'].mean()
        
        print(f"{sbin:<15} {mask.sum():>6} {mean_alpha_age:>+12.4f} {mean_alpha_zh:>+12.4f} {mean_age:>12.2f}")
        
        binned_results.append({
            'sigma_bin': sbin,
            'mean_sigma': mean_sigma,
            'n': mask.sum(),
            'alpha_age': mean_alpha_age,
            'alpha_zh': mean_alpha_zh,
            'age_re': mean_age
        })
    
    print("-" * 70)
    results['binned'] = binned_results
    
    # 6. Interpretation
    print("\n" + "=" * 70)
    print("INTERPRETATION")
    print("=" * 70)
    
    if r_age_grad < -0.05:
        print("\n*** TEP-CONSISTENT PATTERN ***")
        print("  Age gradients become MORE NEGATIVE at higher σ.")
        print("  This means galaxy centers appear YOUNGER relative to outskirts")
        print("  in more massive galaxies - consistent with TEP time dilation.")
        results['interpretation'] = 'TEP_CONSISTENT'
    elif r_age_grad > 0.05:
        print("\n  Age gradients become MORE POSITIVE at higher σ.")
        print("  Galaxy centers appear OLDER relative to outskirts in massive galaxies.")
        print("  This is OPPOSITE to TEP prediction.")
        results['interpretation'] = 'TEP_INCONSISTENT'
    else:
        print("\n  No significant correlation between age gradient and σ.")
        print("  Inconclusive for TEP.")
        results['interpretation'] = 'INCONCLUSIVE'
    
    return results


def test_mass_weighted_ages(df):
    """Test with mass-weighted ages (more robust to recent SF)."""
    print("\n" + "=" * 70)
    print("MASS-WEIGHTED AGE ANALYSIS")
    print("=" * 70)
    
    results = {}
    
    valid = np.isfinite(df['Age_MW_Re_fit']) & np.isfinite(df['alpha_Age_MW_Re_fit'])
    df_mw = df[valid].copy()
    
    if len(df_mw) < 100:
        print("Insufficient data for mass-weighted analysis")
        return results
    
    print(f"\nUsing {len(df_mw)} galaxies with mass-weighted ages")
    
    # Age gradient vs σ
    r_age_mw, p = stats.pearsonr(df_mw['log_sigma'], df_mw['alpha_Age_MW_Re_fit'])
    print(f"\nMass-weighted age gradient vs σ: r = {r_age_mw:.4f}, p = {p:.2e}")
    results['mw_age_gradient_vs_sigma'] = {'r': r_age_mw, 'p': p}
    
    # Central age vs σ
    r_central_mw, p = stats.pearsonr(df_mw['log_sigma'], df_mw['Age_MW_Re_fit'])
    print(f"Mass-weighted central age vs σ: r = {r_central_mw:.4f}, p = {p:.2e}")
    results['mw_central_age_vs_sigma'] = {'r': r_central_mw, 'p': p}
    
    return results


def create_manga_figure(df, results):
    """Create publication figure for MaNGA analysis."""
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # Panel A: Age gradient vs σ
    ax = axes[0, 0]
    ax.scatter(df['log_sigma'], df['alpha_Age_LW_Re_fit'], alpha=0.3, s=10, c='steelblue')
    ax.axhline(0, color='gray', linestyle='--')
    
    # Fit line
    m, b = np.polyfit(df['log_sigma'], df['alpha_Age_LW_Re_fit'], 1)
    x_line = np.linspace(df['log_sigma'].min(), df['log_sigma'].max(), 100)
    ax.plot(x_line, m * x_line + b, 'r-', lw=2)
    
    ax.set_xlabel('log(σ / km s⁻¹)')
    ax.set_ylabel('Age Gradient α_Age')
    ax.set_title(f"A. Age Gradient vs σ (r = {results['age_gradient_vs_sigma']['r']:.3f})")
    
    # Panel B: Central age vs σ
    ax = axes[0, 1]
    ax.scatter(df['log_sigma'], df['Age_LW_Re_fit'], alpha=0.3, s=10, c='darkorange')
    ax.set_xlabel('log(σ / km s⁻¹)')
    ax.set_ylabel('Central Age (Gyr)')
    ax.set_title(f"B. Central Age vs σ (r = {results['central_age_vs_sigma']['r']:.3f})")
    
    # Panel C: Binned comparison
    ax = axes[1, 0]
    binned = results['binned']
    sigmas = [b['mean_sigma'] for b in binned]
    alpha_ages = [b['alpha_age'] for b in binned]
    alpha_zhs = [b['alpha_zh'] for b in binned]
    
    ax.plot(sigmas, alpha_ages, 'o-', color='steelblue', markersize=10, lw=2, label='Age gradient')
    ax.plot(sigmas, alpha_zhs, 's-', color='darkorange', markersize=10, lw=2, label='[Z/H] gradient')
    ax.axhline(0, color='gray', linestyle='--')
    ax.set_xlabel('σ (km/s)')
    ax.set_ylabel('Gradient α')
    ax.set_title('C. Gradients by σ Bin')
    ax.legend()
    
    # Panel D: Summary
    ax = axes[1, 1]
    ax.axis('off')
    
    summary = f"""
MaNGA SPATIALLY-RESOLVED TEST

Sample: {len(df):,} galaxies with IFU data

CORRELATIONS WITH log(σ):
  Age gradient: r = {results['age_gradient_vs_sigma']['r']:.4f}
  [Z/H] gradient: r = {results.get('zh_gradient_vs_sigma', {}).get('r', 0):.4f}
  Central age: r = {results['central_age_vs_sigma']['r']:.4f}

INTERPRETATION:
  {results['interpretation']}

TEP PREDICTION:
  Negative age gradient (younger centers)
  should correlate with higher σ.
  
  α_Age < 0 means center is YOUNGER than outskirts.
  If TEP causes time dilation in galaxy centers,
  we expect α_Age to become more negative at higher σ.
"""
    ax.text(0.05, 0.95, summary, transform=ax.transAxes,
            fontsize=10, verticalalignment='top', fontfamily='monospace',
            bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))
    
    plt.tight_layout()
    
    fig_path = os.path.join(FIGURES_DIR, 'manga_age_gradients.png')
    plt.savefig(fig_path, dpi=150, bbox_inches='tight')
    plt.savefig(fig_path.replace('.png', '.pdf'), bbox_inches='tight')
    print(f"\nFigure saved: {fig_path}")
    
    return fig_path


def main():
    """Main analysis."""
    print("=" * 70)
    print("MaNGA SPATIALLY-RESOLVED AGE GRADIENT TEST FOR TEP")
    print("=" * 70)
    print(f"Started: {datetime.now().isoformat()}")
    
    # Check for cached data
    cache_path = os.path.join(DATA_DIR, 'manga_pipe3d.csv')
    
    if os.path.exists(cache_path):
        print(f"\nLoading cached data from {cache_path}")
        df = pd.read_csv(cache_path)
    else:
        df = download_manga_data()
        if df is None:
            print("Failed to download data.")
            return None
        df.to_csv(cache_path, index=False)
        print(f"Data cached to {cache_path}")
    
    # Prepare data
    df = prepare_data(df)
    
    # Run tests
    results = test_age_gradients(df)
    results['mass_weighted'] = test_mass_weighted_ages(df)
    
    # Add metadata
    results['n_galaxies'] = len(df)
    results['timestamp'] = datetime.now().isoformat()
    
    # Create figure
    fig_path = create_manga_figure(df, results)
    
    # Save results
    results_path = os.path.join(RESULTS_DIR, 'manga_age_gradients_results.json')
    
    def convert_types(obj):
        if isinstance(obj, (np.integer, np.floating)):
            return float(obj)
        elif isinstance(obj, dict):
            return {k: convert_types(v) for k, v in obj.items()}
        elif isinstance(obj, list):
            return [convert_types(i) for i in obj]
        return obj
    
    with open(results_path, 'w') as f:
        json.dump(convert_types(results), f, indent=2)
    print(f"\nResults saved: {results_path}")
    
    # Summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"Galaxies analyzed: {len(df):,}")
    print(f"Age gradient vs σ: r = {results['age_gradient_vs_sigma']['r']:.4f}")
    print(f"Interpretation: {results['interpretation']}")
    
    return results


if __name__ == '__main__':
    results = main()
