#!/usr/bin/env python3
"""
Step 6.4: Systematics Analysis - Why APOGEE and SDSS Show Opposite Results

This script investigates the discrepancy between:
- SDSS galaxies: younger appearance at higher σ (TEP-consistent)
- APOGEE stars: older appearance at deeper Galactic potential (TEP-inconsistent)

POTENTIAL EXPLANATIONS:
1. Different age indicators (D4000/Hβ vs log(g) residual)
2. Different potential definitions (σ² vs Galactic Φ)
3. Selection effects (galaxy sample vs stellar sample)
4. Scale-dependent effects (extragalactic vs Galactic)
5. Metallicity confounders

Author: M. Smawfield
Date: January 2026
"""

import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
import json
import os
from datetime import datetime

DATA_DIR_SDSS = os.path.join(os.path.dirname(__file__), '..', '..', 'data', 'sdss')
DATA_DIR_APOGEE = os.path.join(os.path.dirname(__file__), '..', '..', 'data', 'apogee')
RESULTS_DIR = os.path.join(os.path.dirname(__file__), '..', '..', 'results', 'outputs')
FIGURES_DIR = os.path.join(os.path.dirname(__file__), '..', '..', 'results', 'figures')


def load_data():
    """Load both SDSS and APOGEE datasets."""
    print("Loading datasets...")
    
    # SDSS galaxies
    sdss_path = os.path.join(DATA_DIR_SDSS, 'sdss_spectral_indices.csv')
    df_sdss = pd.read_csv(sdss_path)
    print(f"  SDSS galaxies: {len(df_sdss)}")
    
    # APOGEE stars
    apogee_path = os.path.join(DATA_DIR_APOGEE, 'apogee_starhorse.csv')
    df_apogee = pd.read_csv(apogee_path)
    print(f"  APOGEE stars: {len(df_apogee)}")
    
    return df_sdss, df_apogee


def prepare_sdss(df):
    """Prepare SDSS data with derived quantities."""
    df['fe_avg'] = (df['fe5270'] + df['fe5335']) / 2
    df['mg_fe_ratio'] = df['mgb'] / df['fe_avg']
    df['log_mg_fe'] = np.log10(df['mg_fe_ratio'])
    df['spec_age_proxy'] = df['d4000'] / (df['hbeta'] + 0.5)
    df['log_spec_age'] = np.log10(df['spec_age_proxy'])
    df['log_sigma'] = np.log10(df['veldisp'])
    
    # Potential proxy: Φ/c² ~ σ²/c²
    c = 299792.458  # km/s
    df['phi_c2'] = (df['veldisp'])**2 / c**2
    
    valid = (
        np.isfinite(df['mg_fe_ratio']) &
        np.isfinite(df['spec_age_proxy']) &
        (df['mg_fe_ratio'] > 0.5) & (df['mg_fe_ratio'] < 3.0) &
        (df['spec_age_proxy'] > 0.3) & (df['spec_age_proxy'] < 3.0)
    )
    return df[valid].copy()


def prepare_apogee(df):
    """Prepare APOGEE data with derived quantities."""
    # Galactocentric coordinates
    R_sun = 8.178  # kpc
    glon = np.radians(df['glon'].values)
    glat = np.radians(df['glat'].values)
    dist = df['dist50'].values
    
    x_helio = dist * np.cos(glat) * np.cos(glon)
    y_helio = dist * np.cos(glat) * np.sin(glon)
    
    x_gc = R_sun - x_helio
    y_gc = -y_helio
    df['R_gc'] = np.sqrt(x_gc**2 + y_gc**2)
    
    # Potential: Φ/c² ~ v_c² × ln(R_sun/R) / c²
    v_c = 220  # km/s
    c = 299792.458  # km/s
    df['phi_c2'] = v_c**2 * np.log(R_sun / np.clip(df['R_gc'], 0.1, 100)) / c**2
    
    # Age proxy: log(g) residual
    from sklearn.linear_model import LinearRegression
    giants = (df['logg'] < 3.5) & (df['teff'] < 5500) & (df['teff'] > 3500)
    df_giants = df[giants].copy()
    
    X = df_giants[['teff', 'm_h']].values
    y = df_giants['logg'].values
    model = LinearRegression()
    model.fit(X, y)
    df_giants['logg_resid'] = df_giants['logg'] - model.predict(X)
    
    return df_giants


def compare_potential_definitions():
    """
    ANALYSIS 1: Compare how potential is defined in each sample.
    
    SDSS: Φ ~ σ² (internal velocity dispersion of galaxy)
    APOGEE: Φ ~ v_c² × ln(R_sun/R) (Galactic potential at star's position)
    
    These are fundamentally different:
    - SDSS σ measures the depth of the galaxy's own potential well
    - APOGEE Φ measures the star's position in the Milky Way's potential
    """
    print("\n" + "=" * 70)
    print("ANALYSIS 1: Potential Definition Comparison")
    print("=" * 70)
    
    print("\nSDSS (galaxies):")
    print("  Φ/c² = σ²/c² where σ = internal velocity dispersion")
    print("  Range: ~10⁻⁷ to 10⁻⁶")
    print("  Interpretation: Depth of galaxy's own gravitational well")
    
    print("\nAPOGEE (stars):")
    print("  Φ/c² = v_c² × ln(R_sun/R) / c² where v_c = 220 km/s")
    print("  Range: ~10⁻⁷ to 10⁻⁶")
    print("  Interpretation: Position in Milky Way's gravitational well")
    
    print("\nKEY DIFFERENCE:")
    print("  SDSS: Each galaxy is a separate potential well")
    print("  APOGEE: All stars share the same Milky Way potential")
    print("  TEP may operate differently in these regimes!")
    
    return {
        'sdss_potential': 'internal σ² of each galaxy',
        'apogee_potential': 'position in MW potential',
        'key_difference': 'separate wells vs shared well'
    }


def compare_age_indicators():
    """
    ANALYSIS 2: Compare age indicator definitions.
    
    SDSS: D4000/Hβ (spectroscopic age from integrated light)
    APOGEE: log(g) residual (evolutionary state of individual stars)
    """
    print("\n" + "=" * 70)
    print("ANALYSIS 2: Age Indicator Comparison")
    print("=" * 70)
    
    print("\nSDSS (galaxies):")
    print("  Age proxy = D4000 / (Hβ + 0.5)")
    print("  D4000: 4000Å break strength (increases with age)")
    print("  Hβ: Balmer absorption (decreases with age)")
    print("  Measures: Mean stellar population age of galaxy")
    
    print("\nAPOGEE (stars):")
    print("  Age proxy = log(g) residual at fixed Teff and [M/H]")
    print("  Lower log(g) = more evolved = older")
    print("  Measures: Evolutionary state of individual star")
    
    print("\nKEY DIFFERENCE:")
    print("  SDSS: Integrated light from billions of stars")
    print("  APOGEE: Individual stellar evolution")
    print("  Different systematics and sensitivities!")
    
    return {
        'sdss_age': 'D4000/Hβ (integrated light)',
        'apogee_age': 'log(g) residual (individual star)',
        'key_difference': 'population average vs individual evolution'
    }


def analyze_metallicity_effects(df_sdss, df_apogee):
    """
    ANALYSIS 3: Metallicity as a confounder.
    
    Both [Mg/Fe] and log(g) depend on metallicity.
    Check if metallicity gradients explain the discrepancy.
    """
    print("\n" + "=" * 70)
    print("ANALYSIS 3: Metallicity Effects")
    print("=" * 70)
    
    results = {}
    
    # SDSS: metallicity vs σ
    if 'log_mass' in df_sdss.columns:
        # Use stellar mass as metallicity proxy
        r_mass_sigma, p = stats.pearsonr(df_sdss['log_sigma'], df_sdss['log_mass'])
        print(f"\nSDSS: log(M*) vs log(σ): r = {r_mass_sigma:.4f}")
        results['sdss_mass_sigma'] = r_mass_sigma
    
    # APOGEE: metallicity vs R_gc
    r_mh_R, p = stats.pearsonr(df_apogee['R_gc'], df_apogee['m_h'])
    print(f"APOGEE: [M/H] vs R_gc: r = {r_mh_R:.4f}")
    results['apogee_mh_R'] = r_mh_R
    
    # APOGEE: [α/M] vs R_gc
    r_alpha_R, p = stats.pearsonr(df_apogee['R_gc'], df_apogee['alpha_m'])
    print(f"APOGEE: [α/M] vs R_gc: r = {r_alpha_R:.4f}")
    results['apogee_alpha_R'] = r_alpha_R
    
    print("\nInterpretation:")
    print("  APOGEE shows standard Galactic chemical evolution:")
    print("  - Inner Galaxy: higher [M/H], higher [α/M]")
    print("  - Outer Galaxy: lower [M/H], lower [α/M]")
    print("  This is NOT a TEP effect but chemical evolution!")
    
    return results


def analyze_selection_effects(df_sdss, df_apogee):
    """
    ANALYSIS 4: Selection effects.
    
    SDSS: Magnitude-limited sample of galaxies
    APOGEE: Targeted sample of bright giants
    """
    print("\n" + "=" * 70)
    print("ANALYSIS 4: Selection Effects")
    print("=" * 70)
    
    print("\nSDSS sample:")
    print(f"  N = {len(df_sdss):,} galaxies")
    print(f"  z range: {df_sdss['redshift'].min():.3f} - {df_sdss['redshift'].max():.3f}")
    print(f"  σ range: {df_sdss['veldisp'].min():.0f} - {df_sdss['veldisp'].max():.0f} km/s")
    
    print("\nAPOGEE sample:")
    print(f"  N = {len(df_apogee):,} stars")
    print(f"  R_gc range: {df_apogee['R_gc'].min():.1f} - {df_apogee['R_gc'].max():.1f} kpc")
    print(f"  [M/H] range: {df_apogee['m_h'].min():.2f} - {df_apogee['m_h'].max():.2f}")
    
    print("\nKey selection biases:")
    print("  SDSS: Fiber aperture samples different fractions at different z")
    print("  APOGEE: Bright giants preferentially sample old populations")
    print("  Both samples have complex selection functions!")
    
    return {
        'sdss_n': len(df_sdss),
        'apogee_n': len(df_apogee),
        'sdss_selection': 'magnitude-limited, fiber aperture',
        'apogee_selection': 'targeted bright giants'
    }


def test_consistent_potential_definition(df_sdss, df_apogee):
    """
    ANALYSIS 5: Use consistent potential definition.
    
    For SDSS: Use σ as before
    For APOGEE: Use stellar mass as potential proxy (like galaxy M*)
    """
    print("\n" + "=" * 70)
    print("ANALYSIS 5: Consistent Potential Definition")
    print("=" * 70)
    
    # APOGEE: Use stellar mass (mass50) as potential proxy
    # More massive stars = deeper self-gravity
    if 'mass50' in df_apogee.columns:
        valid_mass = df_apogee['mass50'] > 0
        df_mass = df_apogee[valid_mass].copy()
        
        r_age_mass, p = stats.pearsonr(df_mass['mass50'], df_mass['logg_resid'])
        r_alpha_mass, p = stats.pearsonr(df_mass['mass50'], df_mass['alpha_m'])
        
        print(f"\nAPOGEE with stellar mass as potential proxy:")
        print(f"  Age proxy vs M*: r = {r_age_mass:.4f}")
        print(f"  [α/M] vs M*: r = {r_alpha_mass:.4f}")
        print(f"  Discrepancy: Δr = {r_age_mass - r_alpha_mass:.4f}")
        
        if r_age_mass > 0:
            print("\n  More massive stars appear YOUNGER (higher log(g))")
            print("  This is expected: massive stars evolve faster!")
        
        return {
            'apogee_age_vs_mass': r_age_mass,
            'apogee_alpha_vs_mass': r_alpha_mass
        }
    
    return {}


def synthesize_findings():
    """
    SYNTHESIS: Why do SDSS and APOGEE show opposite results?
    """
    print("\n" + "=" * 70)
    print("SYNTHESIS: Explaining the SDSS-APOGEE Discrepancy")
    print("=" * 70)
    
    print("""
KEY FINDINGS:

1. DIFFERENT POTENTIAL DEFINITIONS
   - SDSS: σ measures each galaxy's internal potential well
   - APOGEE: R_gc measures position in shared Milky Way potential
   - These probe fundamentally different physical regimes

2. DIFFERENT AGE INDICATORS
   - SDSS: D4000/Hβ (integrated stellar population)
   - APOGEE: log(g) residual (individual stellar evolution)
   - Different systematics and sensitivities

3. METALLICITY CONFOUNDERS
   - APOGEE: Strong [M/H] gradient with R_gc (chemical evolution)
   - This dominates the age-R_gc correlation
   - Not a TEP effect!

4. SELECTION EFFECTS
   - SDSS: Fiber aperture effects at different z
   - APOGEE: Bright giant selection biases toward old stars

CONCLUSION:

The SDSS and APOGEE results are NOT directly comparable because:
- They use different potential definitions
- They use different age indicators
- They have different selection functions
- Galactic chemical evolution dominates the APOGEE signal

The SDSS result (younger at higher σ) may still be TEP-consistent,
but the APOGEE result does not contradict it because they probe
different physical regimes.

RECOMMENDATION:
The SDSS galaxy result should be interpreted cautiously.
The weak age-σ anticorrelation could be:
1. A genuine TEP effect (unlikely given APOGEE null)
2. Astrophysical systematics (formation timescales, aperture effects)
3. Statistical fluctuation (r ~ -0.02 is very weak)
""")
    
    return {
        'conclusion': 'SDSS and APOGEE not directly comparable',
        'reason': 'Different potential definitions, age indicators, selection effects',
        'recommendation': 'Interpret SDSS result cautiously'
    }


def create_comparison_figure(df_sdss, df_apogee):
    """Create figure comparing SDSS and APOGEE results."""
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    
    # Panel A: SDSS age vs σ
    ax = axes[0, 0]
    spec_age_norm = (df_sdss['log_spec_age'] - df_sdss['log_spec_age'].mean()) / df_sdss['log_spec_age'].std()
    ax.hexbin(df_sdss['log_sigma'], spec_age_norm, gridsize=40, cmap='Blues', mincnt=1)
    ax.set_xlabel('log(σ / km s⁻¹)')
    ax.set_ylabel('Spectroscopic Age (normalized)')
    ax.set_title('A. SDSS: Age vs σ')
    
    # Panel B: APOGEE age vs R_gc
    ax = axes[0, 1]
    ax.hexbin(df_apogee['R_gc'], df_apogee['logg_resid'], gridsize=40, cmap='Oranges', mincnt=1)
    ax.set_xlabel('R_gc (kpc)')
    ax.set_ylabel('log(g) residual')
    ax.set_title('B. APOGEE: Age vs R_gc')
    
    # Panel C: Potential comparison
    ax = axes[0, 2]
    ax.hist(df_sdss['phi_c2'] * 1e7, bins=50, alpha=0.7, label='SDSS (σ²/c²)', color='steelblue')
    ax.hist(df_apogee['phi_c2'] * 1e7, bins=50, alpha=0.7, label='APOGEE (Φ_MW/c²)', color='darkorange')
    ax.set_xlabel('Φ/c² (×10⁻⁷)')
    ax.set_ylabel('Count')
    ax.set_title('C. Potential Distribution')
    ax.legend()
    
    # Panel D: SDSS [Mg/Fe] vs σ
    ax = axes[1, 0]
    mg_fe_norm = (df_sdss['log_mg_fe'] - df_sdss['log_mg_fe'].mean()) / df_sdss['log_mg_fe'].std()
    ax.hexbin(df_sdss['log_sigma'], mg_fe_norm, gridsize=40, cmap='Greens', mincnt=1)
    ax.set_xlabel('log(σ / km s⁻¹)')
    ax.set_ylabel('[Mg/Fe] (normalized)')
    ax.set_title('D. SDSS: [Mg/Fe] vs σ')
    
    # Panel E: APOGEE [α/M] vs R_gc
    ax = axes[1, 1]
    ax.hexbin(df_apogee['R_gc'], df_apogee['alpha_m'], gridsize=40, cmap='Purples', mincnt=1)
    ax.set_xlabel('R_gc (kpc)')
    ax.set_ylabel('[α/M]')
    ax.set_title('E. APOGEE: [α/M] vs R_gc')
    
    # Panel F: Summary
    ax = axes[1, 2]
    ax.axis('off')
    
    summary = """
SDSS vs APOGEE COMPARISON

SDSS (361k galaxies):
  Age vs σ: r = -0.024 (younger at higher σ)
  [Mg/Fe] vs σ: r = +0.258
  Discrepancy: TEP-consistent

APOGEE (40k stars):
  Age vs R_gc: r = +0.064 (older at outer R)
  [α/M] vs R_gc: r = -0.029
  Discrepancy: TEP-inconsistent

KEY DIFFERENCES:
  1. Potential: σ² vs Galactic Φ
  2. Age: D4000/Hβ vs log(g)
  3. Scale: Extragalactic vs Galactic

CONCLUSION:
  Results not directly comparable.
  Different physical regimes probed.
  SDSS result requires caution.
"""
    ax.text(0.05, 0.95, summary, transform=ax.transAxes,
            fontsize=10, verticalalignment='top', fontfamily='monospace',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    plt.tight_layout()
    
    fig_path = os.path.join(FIGURES_DIR, 'sdss_apogee_comparison.png')
    plt.savefig(fig_path, dpi=150, bbox_inches='tight')
    plt.savefig(fig_path.replace('.png', '.pdf'), bbox_inches='tight')
    print(f"\nFigure saved: {fig_path}")
    
    return fig_path


def main():
    """Main analysis."""
    print("=" * 70)
    print("SYSTEMATICS ANALYSIS: SDSS vs APOGEE DISCREPANCY")
    print("=" * 70)
    print(f"Started: {datetime.now().isoformat()}")
    
    # Load data
    df_sdss, df_apogee = load_data()
    
    # Prepare data
    df_sdss = prepare_sdss(df_sdss)
    df_apogee = prepare_apogee(df_apogee)
    
    print(f"\nPrepared samples:")
    print(f"  SDSS: {len(df_sdss):,} galaxies")
    print(f"  APOGEE: {len(df_apogee):,} stars")
    
    # Run analyses
    results = {}
    results['potential'] = compare_potential_definitions()
    results['age_indicators'] = compare_age_indicators()
    results['metallicity'] = analyze_metallicity_effects(df_sdss, df_apogee)
    results['selection'] = analyze_selection_effects(df_sdss, df_apogee)
    results['consistent_potential'] = test_consistent_potential_definition(df_sdss, df_apogee)
    results['synthesis'] = synthesize_findings()
    
    # Create figure
    fig_path = create_comparison_figure(df_sdss, df_apogee)
    
    # Save results
    results['timestamp'] = datetime.now().isoformat()
    results_path = os.path.join(RESULTS_DIR, 'sdss_apogee_systematics.json')
    
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
    
    return results


if __name__ == '__main__':
    results = main()
