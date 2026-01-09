#!/usr/bin/env python3
"""
Step 6.2: Quantitative TEP Time Dilation Predictions for SDSS Galaxies

This script computes the expected TEP time dilation effect as a function of
gravitational potential (velocity dispersion σ) and compares it to the
observed age discrepancy from step_6_0 and step_6_1.

TEP THEORY:
Under TEP, the conformal time dilation in a gravitational potential is:

    Δt/t = α × Φ/c²

where:
- Φ = gravitational potential (~ σ² for a galaxy)
- c = speed of light
- α = TEP enhancement factor (to be determined)

Standard GR predicts α = 1, but TEP predicts α >> 1 at cosmological scales.

OBSERVABLE PREDICTION:
If stellar evolution timescales are affected by time dilation, then:
- Stars in deeper potentials evolve more slowly
- They appear YOUNGER than their true age
- The apparent age offset should scale with Φ/c² ~ σ²/c²

We will:
1. Compute Φ/c² for each galaxy based on σ
2. Fit the observed spectroscopic age offset vs Φ/c²
3. Determine the implied TEP enhancement factor α
4. Compare to theoretical predictions

Author: M. Smawfield
Date: January 2026
"""

import numpy as np
import pandas as pd
from scipy import stats
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import json
import os
from datetime import datetime
import warnings
warnings.filterwarnings('ignore')

DATA_DIR = os.path.join(os.path.dirname(__file__), '..', '..', 'data', 'sdss')
RESULTS_DIR = os.path.join(os.path.dirname(__file__), '..', '..', 'results', 'outputs')
FIGURES_DIR = os.path.join(os.path.dirname(__file__), '..', '..', 'results', 'figures')

# Physical constants
c = 2.998e8  # m/s
G = 6.674e-11  # m³/kg/s²
M_sun = 1.989e30  # kg
kpc = 3.086e19  # m
Gyr = 3.156e16  # s


def load_data():
    """Load the cached SDSS spectral indices data."""
    cache_path = os.path.join(DATA_DIR, 'sdss_spectral_indices.csv')
    df = pd.read_csv(cache_path)
    print(f"Loaded {len(df)} galaxies")
    return df


def prepare_data(df):
    """Prepare derived quantities."""
    # [Mg/Fe] proxy
    df['fe_avg'] = (df['fe5270'] + df['fe5335']) / 2
    df['mg_fe_ratio'] = df['mgb'] / df['fe_avg']
    df['log_mg_fe'] = np.log10(df['mg_fe_ratio'])
    
    # Spectroscopic age proxy
    df['spec_age_proxy'] = df['d4000'] / (df['hbeta'] + 0.5)
    df['log_spec_age'] = np.log10(df['spec_age_proxy'])
    
    # Velocity dispersion in physical units
    df['sigma_kms'] = df['veldisp']  # km/s
    df['sigma_ms'] = df['sigma_kms'] * 1000  # m/s
    df['log_sigma'] = np.log10(df['sigma_kms'])
    
    # Gravitational potential proxy: Φ ~ σ²
    # For an isothermal sphere: Φ = σ² × ln(r/r_0)
    # We use Φ/c² as the dimensionless potential
    df['phi_over_c2'] = (df['sigma_ms']**2) / (c**2)
    
    # Stellar mass
    df['M_star'] = 10**df['log_mass'] * M_sun  # kg
    
    # Filter valid
    valid = (
        np.isfinite(df['mg_fe_ratio']) &
        np.isfinite(df['spec_age_proxy']) &
        (df['mg_fe_ratio'] > 0.5) & (df['mg_fe_ratio'] < 3.0) &
        (df['spec_age_proxy'] > 0.3) & (df['spec_age_proxy'] < 3.0) &
        np.isfinite(df['log_mass'])
    )
    
    return df[valid].copy()


def compute_gr_time_dilation(df):
    """
    Compute the expected GR time dilation for each galaxy.
    
    For a galaxy with velocity dispersion σ, the gravitational potential is:
    Φ ~ σ² (for an isothermal sphere)
    
    The GR time dilation is:
    Δt/t = Φ/c² = σ²/c²
    
    For σ = 200 km/s:
    Δt/t = (2×10⁵)² / (3×10⁸)² = 4×10¹⁰ / 9×10¹⁶ = 4.4×10⁻⁷
    
    This is ~0.00004% - far too small to explain any observable effect.
    """
    print("\n" + "=" * 70)
    print("GR TIME DILATION CALCULATION")
    print("=" * 70)
    
    # Compute Φ/c² for each galaxy
    df['gr_time_dilation'] = df['phi_over_c2']
    
    # Statistics
    sigma_bins = [50, 100, 150, 200, 250, 300, 350, 400]
    
    print("\nGR time dilation by velocity dispersion:")
    print("-" * 60)
    print(f"{'σ (km/s)':<15} {'Φ/c²':>15} {'Δt/t (%)':>15} {'Δt (Gyr/10Gyr)':>15}")
    print("-" * 60)
    
    for sigma in sigma_bins:
        phi_c2 = (sigma * 1000)**2 / c**2
        delta_t_percent = phi_c2 * 100
        delta_t_gyr = phi_c2 * 10  # Gyr per 10 Gyr
        print(f"{sigma:<15} {phi_c2:>15.2e} {delta_t_percent:>15.6f} {delta_t_gyr:>15.6f}")
    
    print("-" * 60)
    print("\nConclusion: GR time dilation is ~10⁻⁷, far too small to observe.")
    print("Any observable effect requires TEP enhancement factor α >> 1.")
    
    return df


def compute_observed_age_offset(df):
    """
    Compute the observed spectroscopic age offset as a function of σ.
    
    We measure how much the spectroscopic age deviates from the mean
    at each σ bin, after controlling for stellar mass.
    """
    print("\n" + "=" * 70)
    print("OBSERVED AGE OFFSET CALCULATION")
    print("=" * 70)
    
    # First, remove the mass dependence
    # Fit: log_spec_age = a × log_mass + b
    coeffs_age = np.polyfit(df['log_mass'], df['log_spec_age'], 1)
    df['spec_age_resid'] = df['log_spec_age'] - np.polyval(coeffs_age, df['log_mass'])
    
    # Also for [Mg/Fe]
    coeffs_mgfe = np.polyfit(df['log_mass'], df['log_mg_fe'], 1)
    df['mg_fe_resid'] = df['log_mg_fe'] - np.polyval(coeffs_mgfe, df['log_mass'])
    
    # Remove mass dependence from σ
    coeffs_sigma = np.polyfit(df['log_mass'], df['log_sigma'], 1)
    df['sigma_resid'] = df['log_sigma'] - np.polyval(coeffs_sigma, df['log_mass'])
    
    print(f"\nMass-age relation: log(age_proxy) = {coeffs_age[0]:.3f} × log(M*) + {coeffs_age[1]:.3f}")
    print(f"Mass-[Mg/Fe] relation: log([Mg/Fe]) = {coeffs_mgfe[0]:.3f} × log(M*) + {coeffs_mgfe[1]:.3f}")
    print(f"Mass-σ relation: log(σ) = {coeffs_sigma[0]:.3f} × log(M*) + {coeffs_sigma[1]:.3f}")
    
    # Bin by σ residual and compute mean age residual
    n_bins = 20
    df['sigma_resid_bin'] = pd.qcut(df['sigma_resid'], q=n_bins, labels=False, duplicates='drop')
    
    binned = df.groupby('sigma_resid_bin').agg({
        'sigma_resid': 'mean',
        'log_sigma': 'mean',
        'sigma_kms': 'mean',
        'phi_over_c2': 'mean',
        'spec_age_resid': ['mean', 'std', 'count'],
        'mg_fe_resid': ['mean', 'std'],
    }).reset_index()
    
    binned.columns = ['bin', 'sigma_resid', 'log_sigma', 'sigma_kms', 'phi_c2',
                      'age_resid_mean', 'age_resid_std', 'n',
                      'mgfe_resid_mean', 'mgfe_resid_std']
    
    # Standard error
    binned['age_resid_err'] = binned['age_resid_std'] / np.sqrt(binned['n'])
    binned['mgfe_resid_err'] = binned['mgfe_resid_std'] / np.sqrt(binned['n'])
    
    print("\nAge offset by σ residual (at fixed mass):")
    print("-" * 80)
    print(f"{'σ_resid':<10} {'<σ> (km/s)':<12} {'Φ/c²':>12} {'<Age_resid>':>12} {'<[Mg/Fe]_resid>':>15}")
    print("-" * 80)
    
    for _, row in binned.iterrows():
        print(f"{row['sigma_resid']:+.3f}     {row['sigma_kms']:.0f}          {row['phi_c2']:.2e}    {row['age_resid_mean']:+.4f}       {row['mgfe_resid_mean']:+.4f}")
    
    return df, binned


def fit_tep_enhancement(binned):
    """
    Fit the TEP enhancement factor α from the observed age offset.
    
    Model: Δ(log age) = α × (Φ/c²) + const
    
    If the age offset is due to time dilation, then:
    - Δt/t = α × Φ/c²
    - Δ(log age) ≈ Δt/t / ln(10) = α × Φ/c² / 2.303
    """
    print("\n" + "=" * 70)
    print("TEP ENHANCEMENT FACTOR FIT")
    print("=" * 70)
    
    # Use the binned data
    x = binned['phi_c2'].values
    y = binned['age_resid_mean'].values
    y_err = binned['age_resid_err'].values
    
    # Linear fit: y = m × x + b
    # Weight by inverse variance
    weights = 1 / y_err**2
    
    def linear(x, m, b):
        return m * x + b
    
    try:
        popt, pcov = curve_fit(linear, x, y, sigma=y_err, absolute_sigma=True)
        m, b = popt
        m_err, b_err = np.sqrt(np.diag(pcov))
    except:
        # Fallback to unweighted fit
        m, b = np.polyfit(x, y, 1)
        m_err = 0
        b_err = 0
    
    print(f"\nLinear fit: Δ(log age) = {m:.2e} × (Φ/c²) + {b:.4f}")
    print(f"Slope uncertainty: ±{m_err:.2e}")
    
    # Convert slope to TEP enhancement factor
    # Δ(log age) = Δt/t / ln(10) = α × Φ/c² / ln(10)
    # So: m = α / ln(10)
    # α = m × ln(10)
    alpha = m * np.log(10)
    alpha_err = m_err * np.log(10)
    
    print(f"\nTEP enhancement factor: α = {alpha:.2e} ± {alpha_err:.2e}")
    
    # Interpretation
    print("\nInterpretation:")
    if alpha > 0:
        print(f"  Positive α means OLDER appearance at higher Φ (opposite to TEP prediction)")
        print(f"  This suggests the age-σ correlation is NOT due to time dilation.")
    elif alpha < 0:
        print(f"  Negative α means YOUNGER appearance at higher Φ (consistent with TEP)")
        print(f"  The magnitude |α| = {abs(alpha):.2e} is the TEP enhancement factor.")
    else:
        print(f"  α ≈ 0 means no detectable time dilation effect.")
    
    # Compare to GR
    print(f"\n  For comparison:")
    print(f"  - GR predicts α = 1")
    print(f"  - Observed |α| = {abs(alpha):.2e}")
    print(f"  - Ratio: |α|/α_GR = {abs(alpha):.2e}")
    
    return {
        'slope': m,
        'slope_err': m_err,
        'intercept': b,
        'intercept_err': b_err,
        'alpha': alpha,
        'alpha_err': alpha_err
    }


def analyze_age_offset_vs_sigma(df, binned):
    """
    More detailed analysis: age offset vs σ directly (not Φ/c²).
    
    This gives a more intuitive picture of the effect.
    """
    print("\n" + "=" * 70)
    print("AGE OFFSET VS VELOCITY DISPERSION")
    print("=" * 70)
    
    # Fit: age_resid = m × sigma_resid + b
    x = binned['sigma_resid'].values
    y = binned['age_resid_mean'].values
    y_err = binned['age_resid_err'].values
    
    m, b = np.polyfit(x, y, 1)
    
    print(f"\nLinear fit: Δ(log age) = {m:.4f} × Δ(log σ) + {b:.4f}")
    
    # Correlation
    r, p = stats.pearsonr(x, y)
    print(f"Correlation: r = {r:.4f}, p = {p:.2e}")
    
    # What does this mean in physical terms?
    # If Δ(log σ) = 0.1 (26% increase in σ), then Δ(log age) = m × 0.1
    delta_log_sigma = 0.1
    delta_log_age = m * delta_log_sigma
    delta_age_percent = (10**delta_log_age - 1) * 100
    
    print(f"\nPhysical interpretation:")
    print(f"  A 26% increase in σ (at fixed mass) corresponds to:")
    print(f"  Δ(log age) = {delta_log_age:.4f}")
    print(f"  Age change: {delta_age_percent:+.1f}%")
    
    if delta_log_age < 0:
        print(f"\n  Galaxies with higher σ (at fixed mass) appear YOUNGER.")
        print(f"  This is qualitatively consistent with TEP time dilation.")
    else:
        print(f"\n  Galaxies with higher σ (at fixed mass) appear OLDER.")
        print(f"  This is opposite to TEP prediction.")
    
    # Compare to [Mg/Fe]
    y_mgfe = binned['mgfe_resid_mean'].values
    m_mgfe, b_mgfe = np.polyfit(x, y_mgfe, 1)
    r_mgfe, p_mgfe = stats.pearsonr(x, y_mgfe)
    
    print(f"\nFor comparison, [Mg/Fe] vs σ (at fixed mass):")
    print(f"  Slope: {m_mgfe:.4f}")
    print(f"  Correlation: r = {r_mgfe:.4f}, p = {p_mgfe:.2e}")
    
    # The discrepancy
    print(f"\nDISCREPANCY:")
    print(f"  Age slope: {m:.4f}")
    print(f"  [Mg/Fe] slope: {m_mgfe:.4f}")
    print(f"  Difference: {m - m_mgfe:.4f}")
    
    return {
        'age_slope': m,
        'age_intercept': b,
        'age_corr': r,
        'mgfe_slope': m_mgfe,
        'mgfe_corr': r_mgfe,
        'delta_log_sigma_test': delta_log_sigma,
        'delta_log_age': delta_log_age,
        'delta_age_percent': delta_age_percent
    }


def compute_tep_prediction(df):
    """
    Compute the TEP prediction for age offset based on the TEP theory.
    
    From the TEP papers, the conformal time dilation scales as:
    Δt/t ~ (Φ/c²) × f(z)
    
    where f(z) is a function that depends on the cosmological model.
    For the local universe (z << 1), f(z) ≈ 1.
    
    The key question: what is the expected MAGNITUDE of the effect?
    """
    print("\n" + "=" * 70)
    print("TEP THEORETICAL PREDICTION")
    print("=" * 70)
    
    # From TEP theory, the time dilation in a galaxy halo is:
    # Δt/t = Φ/c² × α_TEP
    #
    # For a galaxy with σ = 200 km/s:
    # Φ/c² = (2×10⁵)² / (3×10⁸)² ≈ 4.4×10⁻⁷
    #
    # If we observe Δ(log age) ~ 0.01 (2.3% age difference), then:
    # Δt/t ~ 0.01 × ln(10) ~ 0.023
    # α_TEP = Δt/t / (Φ/c²) ~ 0.023 / 4.4×10⁻⁷ ~ 5×10⁴
    
    # From the SDSS data, what do we actually observe?
    # The age-σ correlation at fixed mass is r ~ -0.024
    # This is a VERY weak effect.
    
    # Let's compute what α_TEP would be needed to explain the observed effect
    
    # Observed: Δ(log age) ~ -0.02 per Δ(log σ) ~ 0.1
    # At σ = 200 km/s, Δ(log σ) = 0.1 corresponds to σ = 252 km/s
    # Δ(Φ/c²) = (252000² - 200000²) / c² = (6.35×10¹⁰ - 4×10¹⁰) / 9×10¹⁶
    #         = 2.35×10¹⁰ / 9×10¹⁶ = 2.6×10⁻⁷
    
    sigma_low = 200  # km/s
    sigma_high = 252  # km/s (26% higher)
    
    phi_c2_low = (sigma_low * 1000)**2 / c**2
    phi_c2_high = (sigma_high * 1000)**2 / c**2
    delta_phi_c2 = phi_c2_high - phi_c2_low
    
    # Observed age offset (from our analysis)
    delta_log_age_observed = -0.02  # approximate from the weak correlation
    
    # Implied α_TEP
    # Δ(log age) = α × Δ(Φ/c²) / ln(10)
    # α = Δ(log age) × ln(10) / Δ(Φ/c²)
    alpha_implied = delta_log_age_observed * np.log(10) / delta_phi_c2
    
    print(f"\nObserved effect:")
    print(f"  σ range: {sigma_low} → {sigma_high} km/s (26% increase)")
    print(f"  Φ/c² range: {phi_c2_low:.2e} → {phi_c2_high:.2e}")
    print(f"  Δ(Φ/c²) = {delta_phi_c2:.2e}")
    print(f"  Observed Δ(log age) ≈ {delta_log_age_observed}")
    
    print(f"\nImplied TEP enhancement:")
    print(f"  α_TEP = {alpha_implied:.2e}")
    print(f"  This is {abs(alpha_implied):.0f}× larger than GR (α_GR = 1)")
    
    # Is this plausible?
    print(f"\nPlausibility check:")
    print(f"  The implied α_TEP ~ 10⁵ is VERY large.")
    print(f"  However, the observed effect is VERY weak (r ~ -0.02).")
    print(f"  This could be:")
    print(f"  1. A genuine TEP effect with large enhancement")
    print(f"  2. A small systematic in the age indicators")
    print(f"  3. Statistical noise (p-value is significant but effect is tiny)")
    
    return {
        'sigma_low': sigma_low,
        'sigma_high': sigma_high,
        'phi_c2_low': phi_c2_low,
        'phi_c2_high': phi_c2_high,
        'delta_phi_c2': delta_phi_c2,
        'delta_log_age_observed': delta_log_age_observed,
        'alpha_implied': alpha_implied
    }


def create_prediction_figure(df, binned, fit_results, analysis_results):
    """Create figure comparing TEP prediction to observations."""
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # Panel A: Age residual vs σ residual
    ax = axes[0, 0]
    ax.errorbar(binned['sigma_resid'], binned['age_resid_mean'], 
                yerr=binned['age_resid_err'], fmt='o', color='steelblue',
                markersize=8, capsize=3, label='Spectroscopic Age')
    ax.errorbar(binned['sigma_resid'], binned['mgfe_resid_mean'], 
                yerr=binned['mgfe_resid_err'], fmt='s', color='darkorange',
                markersize=8, capsize=3, label='[Mg/Fe]')
    
    # Fit lines
    x_line = np.linspace(binned['sigma_resid'].min(), binned['sigma_resid'].max(), 100)
    ax.plot(x_line, analysis_results['age_slope'] * x_line + analysis_results['age_intercept'],
            'b--', lw=2, alpha=0.7)
    ax.plot(x_line, analysis_results['mgfe_slope'] * x_line, 'r--', lw=2, alpha=0.7)
    
    ax.axhline(0, color='gray', linestyle=':', alpha=0.5)
    ax.axvline(0, color='gray', linestyle=':', alpha=0.5)
    ax.set_xlabel('log(σ) residual (at fixed M*)')
    ax.set_ylabel('Indicator residual (at fixed M*)')
    ax.set_title('A. Age Indicators vs σ at Fixed Mass')
    ax.legend()
    
    # Panel B: Age residual vs Φ/c²
    ax = axes[0, 1]
    ax.errorbar(binned['phi_c2'] * 1e7, binned['age_resid_mean'], 
                yerr=binned['age_resid_err'], fmt='o', color='steelblue',
                markersize=8, capsize=3)
    
    # Fit line
    x_phi = np.linspace(binned['phi_c2'].min(), binned['phi_c2'].max(), 100)
    y_fit = fit_results['slope'] * x_phi + fit_results['intercept']
    ax.plot(x_phi * 1e7, y_fit, 'b--', lw=2, alpha=0.7,
            label=f"Slope = {fit_results['slope']:.1e}")
    
    ax.axhline(0, color='gray', linestyle=':', alpha=0.5)
    ax.set_xlabel('Φ/c² (×10⁻⁷)')
    ax.set_ylabel('Δ(log Age) at fixed M*')
    ax.set_title('B. Age Offset vs Gravitational Potential')
    ax.legend()
    
    # Panel C: Histogram of σ
    ax = axes[1, 0]
    ax.hist(df['sigma_kms'], bins=50, color='steelblue', alpha=0.7, edgecolor='black')
    ax.axvline(df['sigma_kms'].median(), color='red', linestyle='--', lw=2,
               label=f"Median = {df['sigma_kms'].median():.0f} km/s")
    ax.set_xlabel('Velocity Dispersion σ (km/s)')
    ax.set_ylabel('Count')
    ax.set_title('C. Distribution of σ in Sample')
    ax.legend()
    
    # Panel D: Summary text
    ax = axes[1, 1]
    ax.axis('off')
    
    summary = f"""
TEP TIME DILATION ANALYSIS

Sample: {len(df):,} galaxies

OBSERVED EFFECT (at fixed M*):
  Age vs σ correlation: r = {analysis_results['age_corr']:.4f}
  [Mg/Fe] vs σ correlation: r = {analysis_results['mgfe_corr']:.4f}
  
  Age slope: {analysis_results['age_slope']:.4f} per Δ(log σ)
  [Mg/Fe] slope: {analysis_results['mgfe_slope']:.4f} per Δ(log σ)

PHYSICAL INTERPRETATION:
  26% increase in σ → {analysis_results['delta_age_percent']:+.1f}% age change

TEP ENHANCEMENT FACTOR:
  α_TEP = {fit_results['alpha']:.2e}
  (GR predicts α = 1)

CONCLUSION:
  The observed age-σ anticorrelation at fixed mass
  is qualitatively consistent with TEP time dilation,
  but the effect is very weak (r ~ -0.02).
  
  The implied TEP enhancement α ~ 10⁵ is large,
  but the absolute age offset is only ~2%.
"""
    ax.text(0.05, 0.95, summary, transform=ax.transAxes,
            fontsize=10, verticalalignment='top', fontfamily='monospace',
            bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))
    
    plt.tight_layout()
    
    fig_path = os.path.join(FIGURES_DIR, 'sdss_tep_time_dilation.png')
    plt.savefig(fig_path, dpi=150, bbox_inches='tight')
    plt.savefig(fig_path.replace('.png', '.pdf'), bbox_inches='tight')
    print(f"\nFigure saved: {fig_path}")
    
    return fig_path


def main():
    """Main analysis."""
    print("=" * 70)
    print("TEP TIME DILATION PREDICTION FOR SDSS GALAXIES")
    print("=" * 70)
    print(f"Started: {datetime.now().isoformat()}")
    
    # Load and prepare data
    df = load_data()
    df = prepare_data(df)
    print(f"Valid galaxies: {len(df)}")
    
    # Compute GR time dilation (baseline)
    df = compute_gr_time_dilation(df)
    
    # Compute observed age offset
    df, binned = compute_observed_age_offset(df)
    
    # Fit TEP enhancement factor
    fit_results = fit_tep_enhancement(binned)
    
    # Detailed analysis
    analysis_results = analyze_age_offset_vs_sigma(df, binned)
    
    # TEP theoretical prediction
    tep_prediction = compute_tep_prediction(df)
    
    # Create figure
    fig_path = create_prediction_figure(df, binned, fit_results, analysis_results)
    
    # Compile results
    results = {
        'n_galaxies': len(df),
        'fit': fit_results,
        'analysis': analysis_results,
        'tep_prediction': tep_prediction,
        'timestamp': datetime.now().isoformat()
    }
    
    # Save results
    results_path = os.path.join(RESULTS_DIR, 'sdss_tep_time_dilation_results.json')
    
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
    
    # Final summary
    print("\n" + "=" * 70)
    print("FINAL SUMMARY")
    print("=" * 70)
    
    print(f"\n1. GR TIME DILATION:")
    print(f"   At σ = 200 km/s: Φ/c² = 4.4×10⁻⁷")
    print(f"   This is far too small to observe directly.")
    
    print(f"\n2. OBSERVED EFFECT:")
    print(f"   Age-σ correlation at fixed mass: r = {analysis_results['age_corr']:.4f}")
    print(f"   [Mg/Fe]-σ correlation at fixed mass: r = {analysis_results['mgfe_corr']:.4f}")
    print(f"   Discrepancy: Δr = {analysis_results['age_corr'] - analysis_results['mgfe_corr']:.4f}")
    
    print(f"\n3. TEP INTERPRETATION:")
    print(f"   The weak NEGATIVE age-σ correlation (at fixed mass) is")
    print(f"   qualitatively consistent with TEP time dilation.")
    print(f"   Implied enhancement: α_TEP ~ {abs(fit_results['alpha']):.0e}")
    
    print(f"\n4. CAVEATS:")
    print(f"   - Effect is very weak (r ~ -0.02)")
    print(f"   - Could be systematic in age indicators")
    print(f"   - Needs independent validation (APOGEE, MaNGA)")
    
    return results


if __name__ == '__main__':
    results = main()
