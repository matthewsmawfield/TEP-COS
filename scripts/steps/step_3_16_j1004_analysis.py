#!/usr/bin/env python3
"""
Step 3.16: SDSS J1004+4112 Temporal Shear Analysis

Analyzes the 14.5-year r-band light curves from Munoz et al. (2022)
for temporal shear. This is a cluster-lensed quasar at z_source = 1.734,
providing a test of the high-z scaling prediction.

System properties:
- z_lens = 0.68 (cluster)
- z_source = 1.734
- 4 bright images (A, B, C, D)
- Einstein radius ~ 14.62 arcsec (largest known)
- Known delays: Δt_DC = 2458 days, Δt_BC = 782 days, Δt_AC = 825 days
"""

import numpy as np
import pandas as pd
from scipy import signal
from scipy.interpolate import interp1d
from scipy.stats import pearsonr, linregress
import json
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

DATA_DIR = Path("/Users/matthewsmawfield/www/TEP-COS/data/cosmograil")
OUTPUT_DIR = Path("/Users/matthewsmawfield/www/TEP-COS/results/outputs")

def load_j1004_data():
    """Load J1004+4112 light curve data."""
    filepath = DATA_DIR / "j1004_JApJ937_34.dat"
    
    # Format: JD-2450000, magA, errA, magB, errB, magC, errC, magD, errD
    df = pd.read_csv(filepath, sep=r'\s+', header=None,
                     names=['JD', 'magA', 'errA', 'magB', 'errB', 
                            'magC', 'errC', 'magD', 'errD'])
    
    # Convert to MJD
    df['MJD'] = df['JD'] + 2450000 - 2400000.5
    
    return df

def gaussian_smooth(t, y, sigma):
    """Apply Gaussian smoothing."""
    if len(t) < 10:
        return np.zeros_like(y)
    
    t_min, t_max = t.min(), t.max()
    t_grid = np.linspace(t_min, t_max, 1000)
    
    try:
        f = interp1d(t, y, kind='linear', fill_value='extrapolate')
        y_grid = f(t_grid)
        
        dt = t_grid[1] - t_grid[0]
        kernel_size = int(3 * sigma / dt)
        if kernel_size < 3:
            kernel_size = 3
        kernel = signal.windows.gaussian(kernel_size * 2 + 1, sigma / dt)
        kernel /= kernel.sum()
        
        y_smooth = np.convolve(y_grid, kernel, mode='same')
        
        f_smooth = interp1d(t_grid, y_smooth, kind='linear', fill_value='extrapolate')
        return f_smooth(t)
    except:
        return np.zeros_like(y)

def bandpass_filter(t, y, tau):
    """Apply bandpass filter using difference of Gaussians."""
    sigma_low = tau * 0.5
    sigma_high = tau * 2.0
    
    smooth_low = gaussian_smooth(t, y, sigma_low)
    smooth_high = gaussian_smooth(t, y, sigma_high)
    
    return smooth_low - smooth_high

def estimate_delay_iccf(t1, y1, t2, y2, search_center, search_width=100):
    """Estimate time delay using ICCF with mode-lock."""
    search_range = (search_center - search_width, search_center + search_width)
    
    t_min = max(t1.min(), t2.min())
    t_max = min(t1.max(), t2.max())
    
    if t_max <= t_min:
        return np.nan, np.nan
    
    t_common = np.linspace(t_min, t_max, 500)
    
    try:
        f1 = interp1d(t1, y1, kind='linear', bounds_error=False, fill_value=np.nan)
        f2 = interp1d(t2, y2, kind='linear', bounds_error=False, fill_value=np.nan)
    except:
        return np.nan, np.nan
    
    y1_interp = f1(t_common)
    y2_interp = f2(t_common)
    
    valid = ~(np.isnan(y1_interp) | np.isnan(y2_interp))
    if valid.sum() < 20:
        return np.nan, np.nan
    
    y1_valid = y1_interp[valid]
    y2_valid = y2_interp[valid]
    t_valid = t_common[valid]
    
    y1_norm = (y1_valid - y1_valid.mean()) / (y1_valid.std() + 1e-10)
    y2_norm = (y2_valid - y2_valid.mean()) / (y2_valid.std() + 1e-10)
    
    delays = np.linspace(search_range[0], search_range[1], 201)
    correlations = []
    
    for delay in delays:
        t_shifted = t_valid + delay
        try:
            f2_shift = interp1d(t2, y2, kind='linear', bounds_error=False, fill_value=np.nan)
            y2_shifted = f2_shift(t_shifted)
            valid_shift = ~np.isnan(y2_shifted)
            if valid_shift.sum() < 10:
                correlations.append(-1)
                continue
            r, _ = pearsonr(y1_norm[valid_shift], y2_shifted[valid_shift])
            correlations.append(r if not np.isnan(r) else -1)
        except:
            correlations.append(-1)
    
    correlations = np.array(correlations)
    best_idx = np.argmax(correlations)
    best_delay = delays[best_idx]
    best_corr = correlations[best_idx]
    
    return best_delay, best_corr

def compute_temporal_shear(t1, y1, t2, y2, broadband_delay, tau_values=[20, 40, 80, 160, 320], 
                           detrend_sigma=400, mode_lock_window=100):
    """Compute temporal shear Γ for an image pair."""
    # Detrend with longer window for this long baseline
    trend1 = gaussian_smooth(t1, y1, detrend_sigma)
    trend2 = gaussian_smooth(t2, y2, detrend_sigma)
    
    y1_detrend = y1 - trend1
    y2_detrend = y2 - trend2
    
    delays = []
    log_taus = []
    correlations = []
    
    for tau in tau_values:
        y1_filt = bandpass_filter(t1, y1_detrend, tau)
        y2_filt = bandpass_filter(t2, y2_detrend, tau)
        
        # Use mode-lock around broadband delay
        delay, corr = estimate_delay_iccf(t1, y1_filt, t2, y2_filt, 
                                          broadband_delay, mode_lock_window)
        
        if not np.isnan(delay) and corr > 0.2:
            delays.append(delay)
            log_taus.append(np.log10(tau))
            correlations.append(corr)
    
    if len(delays) < 3:
        return np.nan, np.nan, []
    
    delays = np.array(delays)
    log_taus = np.array(log_taus)
    
    # Linear fit
    slope, intercept, r_value, p_value, std_err = linregress(log_taus, delays)
    
    return slope, std_err, list(zip(tau_values[:len(delays)], delays, correlations))

def main():
    print("="*70)
    print("SDSS J1004+4112 TEMPORAL SHEAR ANALYSIS")
    print("14.5-year r-band monitoring (Munoz et al. 2022)")
    print("="*70)
    
    # System parameters
    z_lens = 0.68
    z_source = 1.734
    theta_E = 14.62  # arcsec
    
    # Known time delays (days)
    known_delays = {
        'A-C': 825.23,
        'B-C': 782.20,
        'D-C': 2458.47
    }
    
    print(f"\nSystem: z_lens = {z_lens}, z_source = {z_source}")
    print(f"Einstein radius: {theta_E}\"")
    print(f"Known delays: A-C = {known_delays['A-C']:.1f}d, B-C = {known_delays['B-C']:.1f}d, D-C = {known_delays['D-C']:.1f}d")
    
    # TEP prediction for this system
    geom_factor = (1 + z_source) / (1 + z_lens)
    print(f"Geometric factor (1+z_S)/(1+z_L) = {geom_factor:.2f}")
    
    # Load data
    df = load_j1004_data()
    print(f"\nLoaded {len(df)} epochs spanning {(df['MJD'].max() - df['MJD'].min())/365.25:.1f} years")
    
    t = df['MJD'].values
    images = {
        'A': df['magA'].values,
        'B': df['magB'].values,
        'C': df['magC'].values,
        'D': df['magD'].values
    }
    
    # Analyze each pair
    pairs = [('A', 'C'), ('B', 'C'), ('A', 'B')]
    results = {}
    
    print("\n" + "-"*70)
    print("TEMPORAL SHEAR RESULTS")
    print("-"*70)
    
    for img1, img2 in pairs:
        pair_name = f"{img1}-{img2}"
        
        # Get broadband delay
        if f"{img1}-C" in known_delays:
            bd = known_delays[f"{img1}-C"]
        elif f"{img2}-C" in known_delays:
            bd = -known_delays[f"{img2}-C"]
        else:
            bd = known_delays.get(f"{img1}-C", 0) - known_delays.get(f"{img2}-C", 0)
        
        # For A-B pair
        if pair_name == 'A-B':
            bd = known_delays['A-C'] - known_delays['B-C']
        
        gamma, sigma, details = compute_temporal_shear(
            t, images[img1], t, images[img2], 
            broadband_delay=bd,
            tau_values=[40, 80, 160, 320, 640]
        )
        
        results[pair_name] = {
            'gamma': float(gamma) if not np.isnan(gamma) else None,
            'sigma': float(sigma) if not np.isnan(sigma) else None,
            'broadband_delay': bd,
            'details': details
        }
        
        if not np.isnan(gamma):
            sig = abs(gamma) / sigma if sigma > 0 else 0
            print(f"\n{pair_name}:")
            print(f"  Broadband delay: {bd:.1f} days")
            print(f"  Temporal shear Γ = {gamma:.1f} ± {sigma:.1f} days/decade ({sig:.1f}σ)")
            if details:
                print(f"  Scale-dependent delays:")
                for tau, delay, corr in details:
                    print(f"    τ={tau:4d}d: Δt = {delay:7.1f}d (r={corr:.2f})")
        else:
            print(f"\n{pair_name}: Insufficient data for fit")
    
    # Summary
    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)
    
    valid_results = {k: v for k, v in results.items() if v['gamma'] is not None}
    
    if valid_results:
        gammas = [v['gamma'] for v in valid_results.values()]
        mean_gamma = np.mean(gammas)
        
        print(f"\nMean |Γ| across pairs: {np.mean(np.abs(gammas)):.1f} days/decade")
        print(f"Individual Γ values: {[f'{g:.1f}' for g in gammas]}")
        
        # Compare to TEP prediction
        # From DESJ0408 at z_S=2.375: Γ ~ 33 days/decade
        # Scaling: Γ ∝ (1+z_S)/(1+z_L)
        desj0408_gamma = 33
        desj0408_factor = (1 + 2.375) / (1 + 0.597)
        j1004_factor = geom_factor
        predicted_gamma = desj0408_gamma * (j1004_factor / desj0408_factor)
        
        print(f"\nTEP prediction (scaled from DESJ0408): |Γ| ~ {predicted_gamma:.1f} days/decade")
        print(f"Observed mean |Γ|: {np.mean(np.abs(gammas)):.1f} days/decade")
    
    # Save results
    output = {
        'system': 'SDSS J1004+4112',
        'z_lens': z_lens,
        'z_source': z_source,
        'theta_E_arcsec': theta_E,
        'geometric_factor': geom_factor,
        'n_epochs': len(df),
        'baseline_years': (df['MJD'].max() - df['MJD'].min()) / 365.25,
        'known_delays': known_delays,
        'temporal_shear_results': results,
        'source': 'Munoz et al. 2022 (J/ApJ/937/34)'
    }
    
    output_path = OUTPUT_DIR / "j1004_temporal_shear_results.json"
    with open(output_path, 'w') as f:
        json.dump(output, f, indent=2, default=str)
    
    print(f"\nResults saved to: {output_path}")
    
    return output

if __name__ == "__main__":
    main()
