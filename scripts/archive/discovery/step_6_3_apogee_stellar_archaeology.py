#!/usr/bin/env python3
"""
Step 6.3: APOGEE Stellar Archaeology Test for TEP

This script tests TEP predictions using Milky Way stars from APOGEE.
The key insight: stars at different Galactocentric radii experience
different gravitational potentials.

TEP PREDICTION:
- Stars in the inner Galaxy (deeper potential) should appear YOUNGER
  spectroscopically than their nucleosynthesis ages suggest.
- Stars in the outer Galaxy (shallower potential) should show no offset.
- The [α/Fe] ratio (nucleosynthesis clock) should be INDEPENDENT of
  Galactic radius at fixed metallicity.

This is an independent validation of the SDSS galaxy result at stellar scales.

DATA:
- APOGEE DR17 via SDSS SkyServer
- aspcapStar: stellar parameters and abundances
- apogee_starhorse: distances and ages

Author: M. Smawfield
Date: January 2026
"""

import requests
import numpy as np
import pandas as pd
from scipy import stats
from astropy.coordinates import Galactocentric
import astropy.units as u
import matplotlib.pyplot as plt
import json
import os
from datetime import datetime
import warnings
warnings.filterwarnings('ignore')

SDSS_URL = "https://skyserver.sdss.org/dr18/SkyServerWS/SearchTools/SqlSearch"
DATA_DIR = os.path.join(os.path.dirname(__file__), '..', '..', 'data', 'apogee')
RESULTS_DIR = os.path.join(os.path.dirname(__file__), '..', '..', 'results', 'outputs')
FIGURES_DIR = os.path.join(os.path.dirname(__file__), '..', '..', 'results', 'figures')

os.makedirs(DATA_DIR, exist_ok=True)
os.makedirs(RESULTS_DIR, exist_ok=True)
os.makedirs(FIGURES_DIR, exist_ok=True)

# Galactic constants
R_sun = 8.178  # kpc, Galactocentric distance of Sun (GRAVITY Collaboration 2019)
Z_sun = 0.025  # kpc, height above Galactic plane


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


def download_apogee_data():
    """
    Download APOGEE stars with abundances and distances.
    Query tables separately and merge locally to avoid JOIN timeout.
    """
    print("\n" + "=" * 70)
    print("DOWNLOADING APOGEE DATA")
    print("=" * 70)
    
    # Query aspcapStar (abundances) - simplified query
    print("\n1. Querying aspcapStar (abundances)...")
    sql_aspcap = """
    SELECT TOP 50000
        apogee_id,
        teff, logg,
        m_h, alpha_m,
        fe_h,
        c_fe, n_fe
    FROM aspcapStar
    WHERE aspcapflag = 0
    """
    
    df_aspcap = query_sdss(sql_aspcap)
    if df_aspcap is None or len(df_aspcap) == 0:
        print("ERROR: No data from aspcapStar")
        return None
    print(f"  Retrieved {len(df_aspcap)} stars from aspcapStar")
    
    # Query apogee_starhorse (distances)
    print("\n2. Querying apogee_starhorse (distances)...")
    sql_starhorse = """
    SELECT TOP 200000
        apogee_id,
        glon, glat,
        dist50, dist16, dist84,
        mass50
    FROM apogee_starhorse
    WHERE dist50 > 0 AND dist50 < 20
    """
    
    df_starhorse = query_sdss(sql_starhorse)
    if df_starhorse is None or len(df_starhorse) == 0:
        print("ERROR: No data from apogee_starhorse")
        return None
    print(f"  Retrieved {len(df_starhorse)} stars from apogee_starhorse")
    
    # Merge on apogee_id
    print("\n3. Merging tables...")
    df = pd.merge(df_aspcap, df_starhorse, on='apogee_id', how='inner')
    print(f"  Merged: {len(df)} stars with both abundances and distances")
    
    if len(df) == 0:
        print("ERROR: No matching stars after merge")
        return None
    
    return df


def compute_galactocentric_coordinates(df):
    """
    Compute Galactocentric radius and height for each star.
    
    Uses astropy's Galactocentric frame with standard parameters.
    """
    print("\nComputing Galactocentric coordinates...")
    
    # Convert to Galactocentric coordinates
    # Distance is in kpc from StarHorse
    glon = df['glon'].values * u.deg
    glat = df['glat'].values * u.deg
    dist = df['dist50'].values * u.kpc
    
    # Heliocentric Galactic coordinates
    x_helio = dist * np.cos(glat) * np.cos(glon)
    y_helio = dist * np.cos(glat) * np.sin(glon)
    z_helio = dist * np.sin(glat)
    
    # Convert to Galactocentric (Sun at R_sun, Z_sun)
    x_gc = R_sun * u.kpc - x_helio
    y_gc = -y_helio  # Sign convention
    z_gc = z_helio + Z_sun * u.kpc
    
    # Galactocentric radius
    R_gc = np.sqrt(x_gc**2 + y_gc**2)
    
    df['R_gc'] = R_gc.to(u.kpc).value
    df['z_gc'] = z_gc.to(u.kpc).value
    df['abs_z_gc'] = np.abs(df['z_gc'])
    
    print(f"  R_gc range: {df['R_gc'].min():.1f} - {df['R_gc'].max():.1f} kpc")
    print(f"  |z_gc| range: {df['abs_z_gc'].min():.2f} - {df['abs_z_gc'].max():.2f} kpc")
    
    return df


def compute_galactic_potential(df):
    """
    Compute the Galactic gravitational potential at each star's position.
    
    Uses a simple model: Φ = -v_c² × ln(R/R_0)
    where v_c ≈ 220 km/s is the circular velocity.
    
    For TEP, we care about Φ/c².
    """
    print("\nComputing Galactic potential...")
    
    v_c = 220  # km/s, circular velocity
    c = 299792.458  # km/s
    
    # Potential relative to solar position
    # Φ - Φ_sun = v_c² × ln(R_sun/R)
    # Positive for R < R_sun (deeper potential)
    df['delta_phi'] = v_c**2 * np.log(R_sun / df['R_gc'])
    df['delta_phi_c2'] = df['delta_phi'] / c**2
    
    # Also include vertical component (approximate)
    # Φ_z ≈ ν_z² × z² / 2 where ν_z ≈ 70 km/s/kpc
    nu_z = 70  # km/s/kpc
    df['phi_z'] = 0.5 * (nu_z * df['abs_z_gc'])**2
    df['phi_z_c2'] = df['phi_z'] / c**2
    
    # Total potential difference
    df['total_phi_c2'] = df['delta_phi_c2'] + df['phi_z_c2']
    
    print(f"  ΔΦ/c² range: {df['delta_phi_c2'].min():.2e} to {df['delta_phi_c2'].max():.2e}")
    
    return df


def compute_spectroscopic_age_proxy(df):
    """
    Compute spectroscopic age proxies for APOGEE stars.
    
    For red giants, the key age indicators are:
    1. log(g) at fixed Teff and [M/H] - lower log(g) = more evolved = older
    2. [C/N] ratio - decreases with age due to dredge-up
    3. Mass from asteroseismology (if available) - lower mass = older
    
    We use log(g) residual as the primary age proxy.
    """
    print("\nComputing spectroscopic age proxies...")
    
    # Filter to red giants (most reliable ages)
    giants = (df['logg'] < 3.5) & (df['teff'] < 5500)
    df_giants = df[giants].copy()
    print(f"  Red giants: {len(df_giants)} stars")
    
    # Fit log(g) vs Teff and [M/H] to get age proxy
    # Younger stars have higher log(g) at fixed Teff/[M/H]
    from sklearn.linear_model import LinearRegression
    
    X = df_giants[['teff', 'm_h']].values
    y = df_giants['logg'].values
    
    model = LinearRegression()
    model.fit(X, y)
    
    logg_predicted = model.predict(X)
    df_giants['logg_resid'] = df_giants['logg'] - logg_predicted
    
    # Positive residual = higher log(g) than expected = YOUNGER
    # Negative residual = lower log(g) than expected = OLDER
    
    print(f"  log(g) residual range: {df_giants['logg_resid'].min():.3f} to {df_giants['logg_resid'].max():.3f}")
    
    # Also compute [C/N] if available
    if 'c_fe' in df_giants.columns and 'n_fe' in df_giants.columns:
        df_giants['c_n'] = df_giants['c_fe'] - df_giants['n_fe']
        # Higher [C/N] = younger
        print(f"  [C/N] range: {df_giants['c_n'].min():.2f} to {df_giants['c_n'].max():.2f}")
    
    return df_giants


def test_age_vs_galactic_radius(df):
    """
    THE KEY TEP TEST:
    
    Compare spectroscopic age proxy vs Galactocentric radius at fixed [M/H].
    
    TEP prediction:
    - Inner Galaxy (R < R_sun): stars appear YOUNGER (positive log(g) residual)
    - Outer Galaxy (R > R_sun): stars appear normal
    - [α/Fe] should NOT depend on R at fixed [M/H]
    """
    print("\n" + "=" * 70)
    print("APOGEE STELLAR ARCHAEOLOGY TEST")
    print("=" * 70)
    print("\nTEP Prediction: Inner Galaxy stars appear YOUNGER at fixed [M/H].")
    print("GR Prediction: No systematic age offset with Galactic radius.")
    
    results = {}
    
    # Bin by metallicity to control for chemical evolution
    df['mh_bin'] = pd.cut(df['m_h'], bins=5, labels=['very_low', 'low', 'solar', 'high', 'very_high'])
    
    # 1. Age proxy vs R_gc
    print("\n1. log(g) residual (age proxy) vs R_gc:")
    r_age, p_age = stats.pearsonr(df['R_gc'], df['logg_resid'])
    print(f"   Pearson r = {r_age:.4f}, p = {p_age:.2e}")
    results['age_vs_R'] = {'r': r_age, 'p': p_age}
    
    # 2. [α/M] vs R_gc
    print("\n2. [α/M] vs R_gc:")
    r_alpha, p_alpha = stats.pearsonr(df['R_gc'], df['alpha_m'])
    print(f"   Pearson r = {r_alpha:.4f}, p = {p_alpha:.2e}")
    results['alpha_vs_R'] = {'r': r_alpha, 'p': p_alpha}
    
    # 3. Discrepancy
    delta_r = r_age - r_alpha
    print(f"\n3. DISCREPANCY (Δr = r_age - r_alpha):")
    print(f"   Δr = {delta_r:.4f}")
    results['delta_r'] = delta_r
    
    # 4. Binned analysis by R_gc
    print("\n4. Binned analysis by Galactocentric radius:")
    print("-" * 70)
    print(f"{'R_gc (kpc)':<15} {'N':>8} {'<log(g)_resid>':>15} {'<[α/M]>':>12} {'<[M/H]>':>10}")
    print("-" * 70)
    
    R_bins = [0, 4, 6, 8, 10, 12, 20]
    binned_results = []
    
    for i in range(len(R_bins) - 1):
        R_min, R_max = R_bins[i], R_bins[i+1]
        mask = (df['R_gc'] >= R_min) & (df['R_gc'] < R_max)
        
        if mask.sum() < 50:
            continue
        
        sub = df[mask]
        mean_logg_resid = sub['logg_resid'].mean()
        mean_alpha = sub['alpha_m'].mean()
        mean_mh = sub['m_h'].mean()
        
        print(f"{R_min:>3} - {R_max:<3} kpc    {mask.sum():>8} {mean_logg_resid:>+15.4f} {mean_alpha:>12.3f} {mean_mh:>10.3f}")
        
        binned_results.append({
            'R_min': R_min,
            'R_max': R_max,
            'R_mid': (R_min + R_max) / 2,
            'n': mask.sum(),
            'logg_resid': mean_logg_resid,
            'alpha_m': mean_alpha,
            'm_h': mean_mh
        })
    
    print("-" * 70)
    results['binned'] = binned_results
    
    # 5. Control for metallicity
    print("\n5. Partial correlation controlling for [M/H]:")
    
    # Residualize R_gc and age proxy against [M/H]
    from scipy.stats import pearsonr
    
    def partial_corr(x, y, z):
        r_xy, _ = pearsonr(x, y)
        r_xz, _ = pearsonr(x, z)
        r_yz, _ = pearsonr(y, z)
        num = r_xy - r_xz * r_yz
        den = np.sqrt((1 - r_xz**2) * (1 - r_yz**2))
        return num / den if den > 0 else 0
    
    r_age_partial = partial_corr(df['R_gc'], df['logg_resid'], df['m_h'])
    r_alpha_partial = partial_corr(df['R_gc'], df['alpha_m'], df['m_h'])
    
    print(f"   r(Age, R_gc | [M/H]) = {r_age_partial:.4f}")
    print(f"   r([α/M], R_gc | [M/H]) = {r_alpha_partial:.4f}")
    results['partial_age'] = r_age_partial
    results['partial_alpha'] = r_alpha_partial
    
    # 6. Interpretation
    print("\n" + "=" * 70)
    print("INTERPRETATION")
    print("=" * 70)
    
    if r_age > 0 and r_alpha < 0:
        print("\n*** PATTERN CONSISTENT WITH TEP ***")
        print("  - Age proxy INCREASES with R_gc (inner stars appear YOUNGER)")
        print("  - [α/M] DECREASES with R_gc (standard chemical evolution)")
        print("  - The opposite signs suggest different physical origins.")
        results['interpretation'] = 'TEP_CONSISTENT'
    elif r_age > 0 and r_alpha > 0:
        print("\n  Both indicators increase with R_gc.")
        print("  This could be chemical evolution effects.")
        results['interpretation'] = 'CHEMICAL_EVOLUTION'
    else:
        print("\n  Pattern requires further investigation.")
        results['interpretation'] = 'AMBIGUOUS'
    
    return results


def test_age_vs_potential(df):
    """
    Direct test: age proxy vs gravitational potential.
    """
    print("\n" + "=" * 70)
    print("AGE VS GRAVITATIONAL POTENTIAL")
    print("=" * 70)
    
    results = {}
    
    # Correlation with potential
    r_age_phi, p_age_phi = stats.pearsonr(df['delta_phi_c2'], df['logg_resid'])
    r_alpha_phi, p_alpha_phi = stats.pearsonr(df['delta_phi_c2'], df['alpha_m'])
    
    print(f"\nCorrelation with ΔΦ/c² (positive = deeper potential):")
    print(f"  Age proxy: r = {r_age_phi:.4f}, p = {p_age_phi:.2e}")
    print(f"  [α/M]: r = {r_alpha_phi:.4f}, p = {p_alpha_phi:.2e}")
    
    results['age_vs_phi'] = {'r': r_age_phi, 'p': p_age_phi}
    results['alpha_vs_phi'] = {'r': r_alpha_phi, 'p': p_alpha_phi}
    
    # TEP prediction: positive r_age_phi (younger at deeper potential)
    if r_age_phi > 0:
        print("\n  *** Stars in deeper potentials appear YOUNGER ***")
        print("  This is consistent with TEP time dilation.")
    else:
        print("\n  Stars in deeper potentials appear OLDER.")
        print("  This is opposite to TEP prediction.")
    
    return results


def create_apogee_figure(df, results, phi_results):
    """Create publication figure for APOGEE analysis."""
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # Panel A: Age proxy vs R_gc
    ax = axes[0, 0]
    ax.hexbin(df['R_gc'], df['logg_resid'], gridsize=40, cmap='Blues', mincnt=1)
    ax.axhline(0, color='red', linestyle='--', lw=2)
    ax.axvline(R_sun, color='orange', linestyle='--', lw=2, label=f'R☉ = {R_sun} kpc')
    ax.set_xlabel('Galactocentric Radius R (kpc)')
    ax.set_ylabel('log(g) residual (age proxy)')
    ax.set_title(f"A. Age Proxy vs R (r = {results['age_vs_R']['r']:.3f})")
    ax.legend()
    
    # Panel B: [α/M] vs R_gc
    ax = axes[0, 1]
    ax.hexbin(df['R_gc'], df['alpha_m'], gridsize=40, cmap='Oranges', mincnt=1)
    ax.axvline(R_sun, color='blue', linestyle='--', lw=2, label=f'R☉ = {R_sun} kpc')
    ax.set_xlabel('Galactocentric Radius R (kpc)')
    ax.set_ylabel('[α/M]')
    ax.set_title(f"B. [α/M] vs R (r = {results['alpha_vs_R']['r']:.3f})")
    ax.legend()
    
    # Panel C: Binned comparison
    ax = axes[1, 0]
    binned = results['binned']
    R_mids = [b['R_mid'] for b in binned]
    logg_resids = [b['logg_resid'] for b in binned]
    alphas = [b['alpha_m'] for b in binned]
    
    ax2 = ax.twinx()
    ax.plot(R_mids, logg_resids, 'o-', color='steelblue', markersize=10, lw=2, label='Age proxy')
    ax2.plot(R_mids, alphas, 's-', color='darkorange', markersize=10, lw=2, label='[α/M]')
    
    ax.axhline(0, color='steelblue', linestyle=':', alpha=0.5)
    ax.axvline(R_sun, color='gray', linestyle='--', alpha=0.5)
    
    ax.set_xlabel('Galactocentric Radius R (kpc)')
    ax.set_ylabel('log(g) residual', color='steelblue')
    ax2.set_ylabel('[α/M]', color='darkorange')
    ax.set_title('C. Binned Analysis')
    
    lines1, labels1 = ax.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax.legend(lines1 + lines2, labels1 + labels2, loc='upper right')
    
    # Panel D: Summary
    ax = axes[1, 1]
    ax.axis('off')
    
    summary = f"""
APOGEE STELLAR ARCHAEOLOGY TEST

Sample: {len(df):,} red giant stars
R_gc range: {df['R_gc'].min():.1f} - {df['R_gc'].max():.1f} kpc

CORRELATIONS WITH R_gc:
  Age proxy: r = {results['age_vs_R']['r']:.4f}
  [α/M]: r = {results['alpha_vs_R']['r']:.4f}
  Discrepancy: Δr = {results['delta_r']:.4f}

PARTIAL CORRELATIONS (controlling for [M/H]):
  Age proxy: r = {results['partial_age']:.4f}
  [α/M]: r = {results['partial_alpha']:.4f}

CORRELATIONS WITH ΔΦ/c²:
  Age proxy: r = {phi_results['age_vs_phi']['r']:.4f}
  [α/M]: r = {phi_results['alpha_vs_phi']['r']:.4f}

INTERPRETATION:
{results['interpretation']}

TEP PREDICTION:
  Inner Galaxy stars (deeper Φ) should appear
  YOUNGER due to time dilation.
"""
    ax.text(0.05, 0.95, summary, transform=ax.transAxes,
            fontsize=10, verticalalignment='top', fontfamily='monospace',
            bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))
    
    plt.tight_layout()
    
    fig_path = os.path.join(FIGURES_DIR, 'apogee_stellar_archaeology.png')
    plt.savefig(fig_path, dpi=150, bbox_inches='tight')
    plt.savefig(fig_path.replace('.png', '.pdf'), bbox_inches='tight')
    print(f"\nFigure saved: {fig_path}")
    
    return fig_path


def main():
    """Main analysis pipeline."""
    print("=" * 70)
    print("APOGEE STELLAR ARCHAEOLOGY TEST FOR TEP")
    print("=" * 70)
    print(f"Started: {datetime.now().isoformat()}")
    
    # Check for cached data
    cache_path = os.path.join(DATA_DIR, 'apogee_starhorse.csv')
    
    if os.path.exists(cache_path):
        print(f"\nLoading cached data from {cache_path}")
        df = pd.read_csv(cache_path)
    else:
        # Download fresh data
        df = download_apogee_data()
        if df is None:
            print("Failed to download data.")
            return None
        
        # Save cache
        df.to_csv(cache_path, index=False)
        print(f"Data cached to {cache_path}")
    
    print(f"\nTotal stars: {len(df)}")
    
    # Compute Galactocentric coordinates
    df = compute_galactocentric_coordinates(df)
    
    # Compute Galactic potential
    df = compute_galactic_potential(df)
    
    # Compute spectroscopic age proxy
    df = compute_spectroscopic_age_proxy(df)
    
    # Run the key test
    results = test_age_vs_galactic_radius(df)
    
    # Test vs potential directly
    phi_results = test_age_vs_potential(df)
    
    # Combine results
    results['phi_test'] = phi_results
    results['n_stars'] = len(df)
    results['R_gc_range'] = [float(df['R_gc'].min()), float(df['R_gc'].max())]
    results['timestamp'] = datetime.now().isoformat()
    
    # Save results
    results_path = os.path.join(RESULTS_DIR, 'apogee_stellar_archaeology_results.json')
    
    def convert_types(obj):
        if isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, dict):
            return {k: convert_types(v) for k, v in obj.items()}
        elif isinstance(obj, list):
            return [convert_types(i) for i in obj]
        return obj
    
    with open(results_path, 'w') as f:
        json.dump(convert_types(results), f, indent=2)
    print(f"\nResults saved: {results_path}")
    
    # Create figure
    fig_path = create_apogee_figure(df, results, phi_results)
    
    # Summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"Stars analyzed: {len(df):,}")
    print(f"R_gc range: {df['R_gc'].min():.1f} - {df['R_gc'].max():.1f} kpc")
    print(f"Age proxy vs R_gc: r = {results['age_vs_R']['r']:.4f}")
    print(f"[α/M] vs R_gc: r = {results['alpha_vs_R']['r']:.4f}")
    print(f"Interpretation: {results['interpretation']}")
    
    return results


if __name__ == '__main__':
    results = main()
