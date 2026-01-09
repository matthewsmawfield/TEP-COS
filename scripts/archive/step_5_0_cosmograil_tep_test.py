#!/usr/bin/env python3
"""
Step 5.0: COSMOGRAIL Gravitational Lensing TEP Test

CRITICAL TEST: Gravitational lensing time delays are INDEPENDENT of the
isochrony assumption because they measure GEOMETRIC path differences,
not internal galaxy clocks.

TEP Prediction: If time flows differently along different light paths
(due to varying gravitational potential), the time delay ratios should
show systematic deviations from GR predictions.

Key insight: Lensing time delays depend on:
1. Geometric path length (Fermat potential)
2. Shapiro delay (gravitational time dilation along path)

Under TEP, the Shapiro delay component could be DIFFERENT from GR,
leading to anomalous time delay ratios.

Data: COSMOGRAIL monitoring of multiply-imaged quasars
- HE0435-1223 (quad lens, z_lens=0.46, z_source=1.69)
- RXJ1131-1231 (quad lens, z_lens=0.30, z_source=0.66)
- WFI2033-4723 (quad lens, z_lens=0.66, z_source=1.66)
- And others

Author: M. Smawfield
Date: January 2026
"""

import numpy as np
import pandas as pd
from scipy import stats
from scipy.optimize import curve_fit, minimize
from astropy.cosmology import FlatLambdaCDM
import json
import os
from datetime import datetime
import warnings
warnings.filterwarnings('ignore')

cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

DATA_DIR = os.path.join(os.path.dirname(__file__), '..', '..', 'data', 'cosmograil')
RESULTS_DIR = os.path.join(os.path.dirname(__file__), '..', '..', 'results', 'outputs', 'topology')
FIGURES_DIR = os.path.join(os.path.dirname(__file__), '..', '..', 'results', 'figures', 'topology')

os.makedirs(RESULTS_DIR, exist_ok=True)
os.makedirs(FIGURES_DIR, exist_ok=True)

# Published time delays from COSMOGRAIL/H0LiCOW
# Format: {system: {delay_name: (value_days, error_days), ...}}
PUBLISHED_DELAYS = {
    'HE0435': {
        'AB': (8.4, 0.8),    # Bonvin et al. 2017
        'AC': (0.6, 0.7),
        'AD': (14.9, 0.9),
        'BC': (-7.8, 0.8),
        'BD': (6.5, 0.7),
        'CD': (14.3, 0.8),
        'z_lens': 0.4546,
        'z_source': 1.693,
        'sigma_lens': 222,  # km/s, lens velocity dispersion
    },
    'RXJ1131': {
        'AB': (0.7, 1.0),    # Tewes et al. 2013
        'AC': (1.5, 1.4),
        'AD': (91.4, 1.5),
        'BC': (0.8, 1.6),
        'BD': (90.7, 1.4),
        'CD': (89.9, 1.4),
        'z_lens': 0.295,
        'z_source': 0.654,
        'sigma_lens': 323,
    },
    'WFI2033': {
        'AB': (36.2, 1.0),   # Bonvin et al. 2019
        'AC': (-23.3, 1.4),
        'BC': (-59.4, 1.3),
        'z_lens': 0.6575,
        'z_source': 1.662,
        'sigma_lens': 250,
    },
    'PG1115': {
        'AC': (18.8, 1.6),   # Bonvin et al. 2018
        'BC': (9.9, 1.1),
        'z_lens': 0.311,
        'z_source': 1.722,
        'sigma_lens': 281,
    },
}


def load_light_curves(system):
    """Load COSMOGRAIL light curve data."""
    files = {
        'HE0435': 'HE0435_Bonvin2016.rdb',
        'RXJ1131': 'RXJ1131_Tewes2013.rdb',
        'WFI2033': 'WFI2033_Bonvin2019.rdb',
    }
    
    if system not in files:
        return None
    
    path = os.path.join(DATA_DIR, files[system])
    if not os.path.exists(path):
        return None
    
    # Read RDB format (tab-separated with header)
    df = pd.read_csv(path, sep='\t', comment='=', skiprows=1)
    return df


def measure_time_delay_ccf(lc1, lc2, err1, err2, time, max_lag=150):
    """
    Measure time delay using cross-correlation function.
    
    Returns: (delay, error, ccf_peak)
    """
    from scipy.interpolate import interp1d
    from scipy.ndimage import gaussian_filter1d
    
    # Interpolate to regular grid
    t_min, t_max = time.min(), time.max()
    t_grid = np.linspace(t_min, t_max, 1000)
    
    # Interpolate light curves
    f1 = interp1d(time, lc1, kind='linear', fill_value='extrapolate')
    f2 = interp1d(time, lc2, kind='linear', fill_value='extrapolate')
    
    lc1_interp = f1(t_grid)
    lc2_interp = f2(t_grid)
    
    # Normalize
    lc1_norm = (lc1_interp - np.mean(lc1_interp)) / np.std(lc1_interp)
    lc2_norm = (lc2_interp - np.mean(lc2_interp)) / np.std(lc2_interp)
    
    # Cross-correlation
    dt = t_grid[1] - t_grid[0]
    max_shift = int(max_lag / dt)
    
    lags = []
    ccf = []
    
    for shift in range(-max_shift, max_shift + 1):
        if shift < 0:
            c = np.corrcoef(lc1_norm[:shift], lc2_norm[-shift:])[0, 1]
        elif shift > 0:
            c = np.corrcoef(lc1_norm[shift:], lc2_norm[:-shift])[0, 1]
        else:
            c = np.corrcoef(lc1_norm, lc2_norm)[0, 1]
        
        lags.append(shift * dt)
        ccf.append(c if np.isfinite(c) else 0)
    
    lags = np.array(lags)
    ccf = np.array(ccf)
    
    # Find peak
    peak_idx = np.argmax(ccf)
    delay = lags[peak_idx]
    ccf_peak = ccf[peak_idx]
    
    # Estimate error from CCF width
    half_max = ccf_peak / 2
    above_half = ccf > half_max
    if np.sum(above_half) > 1:
        width = lags[above_half][-1] - lags[above_half][0]
        error = width / 2.35  # FWHM to sigma
    else:
        error = 10.0  # Default
    
    return delay, error, ccf_peak


def compute_fermat_potential_ratio(z_lens, z_source, sigma_lens):
    """
    Compute expected Fermat potential contribution to time delay.
    
    Under GR: Δt = (1+z_lens) * D_Δt * Δφ_Fermat / c
    
    Where D_Δt is the time-delay distance.
    """
    # Angular diameter distances
    D_l = cosmo.angular_diameter_distance(z_lens).value  # Mpc
    D_s = cosmo.angular_diameter_distance(z_source).value
    D_ls = cosmo.angular_diameter_distance_z1z2(z_lens, z_source).value
    
    # Time-delay distance
    D_dt = (1 + z_lens) * D_l * D_s / D_ls
    
    # Einstein radius (approximate)
    c = 299792.458  # km/s
    theta_E = 4 * np.pi * (sigma_lens / c)**2 * D_ls / D_s  # radians
    theta_E_arcsec = theta_E * 206265
    
    return {
        'D_l': D_l,
        'D_s': D_s,
        'D_ls': D_ls,
        'D_dt': D_dt,
        'theta_E_arcsec': theta_E_arcsec,
        'z_lens': z_lens,
        'z_source': z_source,
    }


def test_time_delay_ratios():
    """
    Test if time delay RATIOS match GR predictions.
    
    Key insight: Absolute delays depend on H0 and lens model,
    but RATIOS should be robust and depend only on geometry.
    
    TEP prediction: If time dilation varies along different paths,
    ratios could deviate from GR.
    """
    print("\n" + "=" * 70)
    print("TIME DELAY RATIO TEST")
    print("=" * 70)
    
    results = {}
    
    for system, delays in PUBLISHED_DELAYS.items():
        if 'z_lens' not in delays:
            continue
        
        print(f"\n{system} (z_lens={delays['z_lens']:.3f}, z_source={delays['z_source']:.3f}):")
        
        # Get all delay pairs
        delay_names = [k for k in delays.keys() if len(k) == 2 and k[0] in 'ABCD' and k[1] in 'ABCD']
        
        if len(delay_names) < 3:
            print("  Insufficient delays for ratio test")
            continue
        
        # Compute ratios
        ratios = {}
        for i, d1 in enumerate(delay_names):
            for d2 in delay_names[i+1:]:
                val1, err1 = delays[d1]
                val2, err2 = delays[d2]
                
                if abs(val2) > 1:  # Avoid division by small numbers
                    ratio = val1 / val2
                    ratio_err = abs(ratio) * np.sqrt((err1/val1)**2 + (err2/val2)**2) if val1 != 0 else err1/abs(val2)
                    ratios[f'{d1}/{d2}'] = (ratio, ratio_err)
        
        # Check for consistency
        # In GR, ratios should follow from Fermat potential differences
        # Any deviation could indicate TEP effects
        
        print(f"  Time delay ratios:")
        for name, (ratio, err) in ratios.items():
            print(f"    {name}: {ratio:.3f} ± {err:.3f}")
        
        results[system] = {
            'delays': {k: v for k, v in delays.items() if isinstance(v, tuple)},
            'ratios': ratios,
            'z_lens': delays['z_lens'],
            'z_source': delays['z_source'],
        }
    
    return results


def test_shapiro_delay_anomaly():
    """
    Test for anomalous Shapiro delay component.
    
    The total time delay has two components:
    1. Geometric delay (path length difference)
    2. Shapiro delay (gravitational time dilation)
    
    Under TEP, the Shapiro component could be ENHANCED or MODIFIED.
    
    We test this by comparing systems with different lens masses/potentials.
    """
    print("\n" + "=" * 70)
    print("SHAPIRO DELAY ANOMALY TEST")
    print("=" * 70)
    
    # For each system, estimate the Shapiro contribution
    # Shapiro delay ~ (4GM/c³) * ln(geometric_factor)
    # For a lens with velocity dispersion σ, M ~ σ² * R_E / G
    
    results = []
    
    for system, delays in PUBLISHED_DELAYS.items():
        if 'sigma_lens' not in delays:
            continue
        
        sigma = delays['sigma_lens']
        z_lens = delays['z_lens']
        z_source = delays['z_source']
        
        # Estimate lens mass from velocity dispersion
        # M ~ σ² * R_E / G, where R_E ~ 4π(σ/c)² * D_ls/D_s * D_l
        c = 299792.458
        D_l = cosmo.angular_diameter_distance(z_lens).value * 3.086e19  # km
        D_s = cosmo.angular_diameter_distance(z_source).value * 3.086e19
        D_ls = cosmo.angular_diameter_distance_z1z2(z_lens, z_source).value * 3.086e19
        
        theta_E = 4 * np.pi * (sigma / c)**2 * D_ls / D_s  # radians
        R_E = theta_E * D_l  # km
        
        G = 1.327e11  # km³/M_sun/s²
        M_lens = sigma**2 * R_E / G  # M_sun (approximate)
        
        # Shapiro delay scale: 4GM/c³
        shapiro_scale = 4 * G * M_lens / c**3  # seconds
        shapiro_scale_days = shapiro_scale / 86400
        
        # Get maximum observed delay
        max_delay = max([abs(v[0]) for k, v in delays.items() if isinstance(v, tuple)])
        
        # Ratio of observed to Shapiro scale
        ratio = max_delay / shapiro_scale_days if shapiro_scale_days > 0 else 0
        
        print(f"\n{system}:")
        print(f"  σ_lens = {sigma} km/s")
        print(f"  M_lens ~ {M_lens:.2e} M_sun")
        print(f"  Shapiro scale: {shapiro_scale_days:.2f} days")
        print(f"  Max observed delay: {max_delay:.1f} days")
        print(f"  Ratio (obs/Shapiro): {ratio:.1f}")
        
        results.append({
            'system': system,
            'sigma_lens': sigma,
            'M_lens': float(M_lens),
            'shapiro_scale_days': float(shapiro_scale_days),
            'max_delay_days': float(max_delay),
            'ratio': float(ratio),
            'z_lens': z_lens,
        })
    
    # Test if ratio correlates with lens properties
    if len(results) >= 3:
        sigmas = [r['sigma_lens'] for r in results]
        ratios = [r['ratio'] for r in results]
        
        r, p = stats.pearsonr(sigmas, ratios)
        print(f"\n  Correlation (σ vs ratio): r = {r:.3f}, p = {p:.3f}")
        
        if p < 0.1:
            print("  → Significant correlation detected!")
            print("  → Could indicate TEP-modified Shapiro delay")
    
    return results


def test_redshift_dependent_delay():
    """
    Test if time delays show anomalous redshift dependence.
    
    Under GR: Δt ∝ (1+z_lens) * D_dt
    Under TEP: Additional z-dependent terms possible
    """
    print("\n" + "=" * 70)
    print("REDSHIFT-DEPENDENT DELAY TEST")
    print("=" * 70)
    
    results = []
    
    for system, delays in PUBLISHED_DELAYS.items():
        if 'z_lens' not in delays or 'sigma_lens' not in delays:
            continue
        
        z_lens = delays['z_lens']
        z_source = delays['z_source']
        sigma = delays['sigma_lens']
        
        # Get maximum delay
        max_delay = max([abs(v[0]) for k, v in delays.items() if isinstance(v, tuple)])
        
        # Compute expected scaling
        # Δt ∝ (1+z_l) * D_l * D_s / D_ls * θ_E²
        D_l = cosmo.angular_diameter_distance(z_lens).value
        D_s = cosmo.angular_diameter_distance(z_source).value
        D_ls = cosmo.angular_diameter_distance_z1z2(z_lens, z_source).value
        
        c = 299792.458
        theta_E = 4 * np.pi * (sigma / c)**2 * D_ls / D_s
        
        # Expected delay scale (arbitrary normalization)
        expected_scale = (1 + z_lens) * D_l * D_s / D_ls * theta_E**2
        
        # Normalized delay
        normalized_delay = max_delay / expected_scale
        
        results.append({
            'system': system,
            'z_lens': z_lens,
            'z_source': z_source,
            'max_delay': max_delay,
            'expected_scale': expected_scale,
            'normalized_delay': normalized_delay,
        })
        
        print(f"  {system}: z_l={z_lens:.3f}, Δt_max={max_delay:.1f}d, normalized={normalized_delay:.2e}")
    
    # Test for trend with redshift
    if len(results) >= 3:
        z_vals = [r['z_lens'] for r in results]
        norm_vals = [r['normalized_delay'] for r in results]
        
        r, p = stats.pearsonr(z_vals, norm_vals)
        print(f"\n  Correlation (z_lens vs normalized delay): r = {r:.3f}, p = {p:.3f}")
        
        if abs(r) > 0.5:
            print("  → Strong correlation with redshift!")
            if r > 0:
                print("  → Delays LARGER than expected at high z (TEP-consistent?)")
            else:
                print("  → Delays SMALLER than expected at high z")
    
    return results


def analyze_light_curve_variability():
    """
    Analyze quasar variability for TEP signatures.
    
    If time flows differently along different image paths,
    the variability timescales should differ systematically.
    """
    print("\n" + "=" * 70)
    print("VARIABILITY TIMESCALE ANALYSIS")
    print("=" * 70)
    
    results = {}
    
    for system in ['HE0435', 'RXJ1131', 'WFI2033']:
        df = load_light_curves(system)
        if df is None:
            continue
        
        print(f"\n{system}:")
        
        # Get time and magnitude columns
        time_cols = [c for c in df.columns if 'hjd' in c.lower() or 'mjd' in c.lower()]
        if not time_cols:
            # Try first column as time
            time_col = df.columns[0]
        else:
            time_col = time_cols[0]
        time = df[time_col].values
        
        # Analyze each image
        image_results = {}
        
        for img in ['A', 'B', 'C', 'D']:
            mag_col = f'mag_{img}'
            if mag_col not in df.columns:
                continue
            
            mag = df[mag_col].values
            valid = np.isfinite(mag)
            
            if np.sum(valid) < 50:
                continue
            
            t = time[valid]
            m = mag[valid]
            
            # Compute structure function
            # SF(τ) = sqrt(<(m(t+τ) - m(t))²>)
            tau_bins = np.logspace(0, 3, 20)  # 1 to 1000 days
            sf = []
            
            for tau in tau_bins:
                diffs = []
                for i in range(len(t)):
                    mask = (np.abs(t - t[i] - tau) < tau * 0.2)
                    if np.any(mask):
                        diffs.extend(np.abs(m[mask] - m[i]))
                
                if len(diffs) > 5:
                    sf.append(np.sqrt(np.mean(np.array(diffs)**2)))
                else:
                    sf.append(np.nan)
            
            sf = np.array(sf)
            
            # Fit power law: SF ∝ τ^β
            valid_sf = np.isfinite(sf) & (sf > 0)
            if np.sum(valid_sf) > 5:
                log_tau = np.log10(tau_bins[valid_sf])
                log_sf = np.log10(sf[valid_sf])
                
                slope, intercept, r, p, se = stats.linregress(log_tau, log_sf)
                
                print(f"  Image {img}: SF ∝ τ^{slope:.2f} (r²={r**2:.2f})")
                
                image_results[img] = {
                    'sf_slope': float(slope),
                    'sf_intercept': float(intercept),
                    'r_squared': float(r**2),
                }
        
        # Compare slopes between images
        if len(image_results) >= 2:
            slopes = [v['sf_slope'] for v in image_results.values()]
            slope_std = np.std(slopes)
            slope_mean = np.mean(slopes)
            
            print(f"  Slope variation: {slope_std:.3f} (mean={slope_mean:.2f})")
            
            # Under GR, slopes should be identical (same source)
            # Under TEP, different time flow → different apparent timescales
            if slope_std > 0.1:
                print("  → Significant slope variation detected!")
                print("  → Could indicate different time flow along paths")
        
        results[system] = image_results
    
    return results


def create_visualization(ratio_results, shapiro_results, z_results, var_results, output_path):
    """Create comprehensive visualization."""
    import matplotlib.pyplot as plt
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    
    # 1. Time delay ratios
    ax = axes[0, 0]
    
    systems = list(ratio_results.keys())
    if systems:
        y_pos = 0
        for system in systems:
            ratios = ratio_results[system].get('ratios', {})
            for name, (val, err) in ratios.items():
                ax.errorbar(val, y_pos, xerr=err, fmt='o', capsize=3)
                ax.text(val + err + 0.1, y_pos, f'{system} {name}', fontsize=8, va='center')
                y_pos += 1
        
        ax.axvline(0, color='gray', linestyle='--', alpha=0.5)
        ax.set_xlabel('Time Delay Ratio')
        ax.set_ylabel('Measurement')
        ax.set_title('Time Delay Ratios')
    
    # 2. Shapiro delay scaling
    ax = axes[0, 1]
    
    if shapiro_results:
        sigmas = [r['sigma_lens'] for r in shapiro_results]
        ratios = [r['ratio'] for r in shapiro_results]
        labels = [r['system'] for r in shapiro_results]
        
        ax.scatter(sigmas, ratios, s=100)
        for i, label in enumerate(labels):
            ax.annotate(label, (sigmas[i], ratios[i]), fontsize=10)
        
        ax.set_xlabel('Lens Velocity Dispersion (km/s)')
        ax.set_ylabel('Observed / Shapiro Scale')
        ax.set_title('Shapiro Delay Scaling')
    
    # 3. Redshift dependence
    ax = axes[1, 0]
    
    if z_results:
        z_vals = [r['z_lens'] for r in z_results]
        norm_vals = [r['normalized_delay'] for r in z_results]
        labels = [r['system'] for r in z_results]
        
        ax.scatter(z_vals, norm_vals, s=100)
        for i, label in enumerate(labels):
            ax.annotate(label, (z_vals[i], norm_vals[i]), fontsize=10)
        
        ax.set_xlabel('Lens Redshift')
        ax.set_ylabel('Normalized Delay')
        ax.set_title('Redshift Dependence')
    
    # 4. Summary
    ax = axes[1, 1]
    ax.axis('off')
    
    summary = """
COSMOGRAIL TEP TEST SUMMARY

Gravitational lensing time delays are INDEPENDENT
of the isochrony assumption - they measure geometric
path differences, not internal clocks.

KEY TESTS:

1. TIME DELAY RATIOS
   - Ratios should be robust to H0 uncertainty
   - Any deviation from GR could indicate TEP

2. SHAPIRO DELAY SCALING
   - Tests if gravitational time dilation
     follows GR or shows TEP enhancement

3. REDSHIFT DEPENDENCE
   - Tests for anomalous z-dependent terms
     in the time delay formula

4. VARIABILITY TIMESCALES
   - Different images should show identical
     variability (same source)
   - Differences could indicate different
     time flow along paths

INTERPRETATION:
Lensing provides a CLEAN test of TEP because
it doesn't assume isochrony in the analysis.
"""
    
    ax.text(0.05, 0.95, summary, transform=ax.transAxes, fontsize=10,
           verticalalignment='top', fontfamily='monospace')
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    plt.close()
    print(f"\nVisualization saved: {output_path}")


def main():
    """Main analysis."""
    print("=" * 70)
    print("COSMOGRAIL GRAVITATIONAL LENSING TEP TEST")
    print("=" * 70)
    print(f"Started: {datetime.now().isoformat()}")
    print("\nKey insight: Lensing time delays are INDEPENDENT of isochrony!")
    print("They measure geometric path differences, not internal clocks.")
    
    ratio_results = test_time_delay_ratios()
    shapiro_results = test_shapiro_delay_anomaly()
    z_results = test_redshift_dependent_delay()
    var_results = analyze_light_curve_variability()
    
    # Visualization
    fig_path = os.path.join(FIGURES_DIR, 'step_5_0_cosmograil_tep.png')
    create_visualization(ratio_results, shapiro_results, z_results, var_results, fig_path)
    
    # Save results
    results = {
        'metadata': {
            'timestamp': datetime.now().isoformat(),
            'description': 'COSMOGRAIL gravitational lensing TEP test',
        },
        'time_delay_ratios': ratio_results,
        'shapiro_delay': shapiro_results,
        'redshift_dependence': z_results,
        'variability': var_results,
    }
    
    output_path = os.path.join(RESULTS_DIR, 'step_5_0_cosmograil_tep.json')
    with open(output_path, 'w') as f:
        json.dump(results, f, indent=2, default=str)
    print(f"\nResults saved: {output_path}")
    
    return results


if __name__ == '__main__':
    results = main()
