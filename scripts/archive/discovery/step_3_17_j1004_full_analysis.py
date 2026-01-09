#!/usr/bin/env python3
"""
Step 3.17: Full J1004+4112 Analysis with All Image Pairs

Analyzes all 6 image pairs (A-B, A-C, A-D, B-C, B-D, C-D) for temporal shear.
Uses longer tau scales appropriate for this 14.5-year baseline with very long delays.
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
    df = pd.read_csv(filepath, sep=r'\s+', header=None,
                     names=['JD', 'magA', 'errA', 'magB', 'errB', 
                            'magC', 'errC', 'magD', 'errD'])
    df['MJD'] = df['JD'] + 2450000 - 2400000.5
    return df

def gaussian_smooth(t, y, sigma):
    """Apply Gaussian smoothing."""
    if len(t) < 10:
        return np.zeros_like(y)
    
    t_min, t_max = t.min(), t.max()
    n_points = min(2000, int((t_max - t_min) / 1))
    t_grid = np.linspace(t_min, t_max, n_points)
    
    try:
        f = interp1d(t, y, kind='linear', fill_value='extrapolate')
        y_grid = f(t_grid)
        
        dt = t_grid[1] - t_grid[0]
        kernel_size = int(3 * sigma / dt)
        if kernel_size < 3:
            kernel_size = 3
        if kernel_size > len(t_grid) // 2:
            kernel_size = len(t_grid) // 2
        kernel = signal.windows.gaussian(kernel_size * 2 + 1, sigma / dt)
        kernel /= kernel.sum()
        
        y_smooth = np.convolve(y_grid, kernel, mode='same')
        
        f_smooth = interp1d(t_grid, y_smooth, kind='linear', fill_value='extrapolate')
        return f_smooth(t)
    except:
        return np.zeros_like(y)

def bandpass_filter(t, y, tau):
    """Apply bandpass filter."""
    sigma_low = tau * 0.5
    sigma_high = tau * 2.0
    
    smooth_low = gaussian_smooth(t, y, sigma_low)
    smooth_high = gaussian_smooth(t, y, sigma_high)
    
    return smooth_low - smooth_high

def estimate_delay_iccf(t1, y1, t2, y2, search_center, search_width=150):
    """Estimate time delay using ICCF."""
    search_range = (search_center - search_width, search_center + search_width)
    
    t_min = max(t1.min(), t2.min())
    t_max = min(t1.max(), t2.max())
    
    if t_max <= t_min + 100:
        return np.nan, np.nan
    
    t_common = np.linspace(t_min, t_max, 800)
    
    try:
        f1 = interp1d(t1, y1, kind='linear', bounds_error=False, fill_value=np.nan)
        f2 = interp1d(t2, y2, kind='linear', bounds_error=False, fill_value=np.nan)
    except:
        return np.nan, np.nan
    
    y1_interp = f1(t_common)
    y2_interp = f2(t_common)
    
    valid = ~(np.isnan(y1_interp) | np.isnan(y2_interp))
    if valid.sum() < 50:
        return np.nan, np.nan
    
    y1_valid = y1_interp[valid]
    t_valid = t_common[valid]
    
    y1_norm = (y1_valid - y1_valid.mean()) / (y1_valid.std() + 1e-10)
    
    delays = np.linspace(search_range[0], search_range[1], 301)
    correlations = []
    
    for delay in delays:
        t_shifted = t_valid + delay
        mask = (t_shifted >= t2.min()) & (t_shifted <= t2.max())
        if mask.sum() < 30:
            correlations.append(-1)
            continue
        try:
            f2_eval = interp1d(t2, y2, kind='linear', bounds_error=False, fill_value=np.nan)
            y2_at_shifted = f2_eval(t_shifted[mask])
            y2_norm = (y2_at_shifted - np.nanmean(y2_at_shifted)) / (np.nanstd(y2_at_shifted) + 1e-10)
            
            valid_both = ~np.isnan(y2_norm)
            if valid_both.sum() < 20:
                correlations.append(-1)
                continue
            
            r, _ = pearsonr(y1_norm[mask][valid_both], y2_norm[valid_both])
            correlations.append(r if not np.isnan(r) else -1)
        except:
            correlations.append(-1)
    
    correlations = np.array(correlations)
    best_idx = np.argmax(correlations)
    best_delay = delays[best_idx]
    best_corr = correlations[best_idx]
    
    return best_delay, best_corr

def compute_temporal_shear(t1, y1, t2, y2, broadband_delay, 
                           tau_values=[50, 100, 200, 400, 800], 
                           detrend_sigma=500, mode_lock_window=150):
    """Compute temporal shear Γ."""
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
        
        # Check variance
        if np.std(y1_filt) < 0.001 or np.std(y2_filt) < 0.001:
            continue
        
        delay, corr = estimate_delay_iccf(t1, y1_filt, t2, y2_filt, 
                                          broadband_delay, mode_lock_window)
        
        if not np.isnan(delay) and corr > 0.15:
            delays.append(delay)
            log_taus.append(np.log10(tau))
            correlations.append(corr)
    
    if len(delays) < 3:
        return np.nan, np.nan, []
    
    delays = np.array(delays)
    log_taus = np.array(log_taus)
    
    slope, intercept, r_value, p_value, std_err = linregress(log_taus, delays)
    
    return slope, std_err, list(zip(tau_values[:len(delays)], delays, correlations))

def main():
    print("="*70)
    print("SDSS J1004+4112 FULL TEMPORAL SHEAR ANALYSIS")
    print("All 6 image pairs with extended tau range")
    print("="*70)
    
    # System parameters
    z_lens = 0.68
    z_source = 1.734
    
    # Known time delays (relative to C)
    delays_to_C = {
        'A': 825.23,   # A leads C by 825 days
        'B': 782.20,   # B leads C by 782 days
        'C': 0,
        'D': -2458.47  # D trails C by 2458 days
    }
    
    df = load_j1004_data()
    print(f"\nLoaded {len(df)} epochs spanning {(df['MJD'].max() - df['MJD'].min())/365.25:.1f} years")
    
    t = df['MJD'].values
    images = {
        'A': df['magA'].values,
        'B': df['magB'].values,
        'C': df['magC'].values,
        'D': df['magD'].values
    }
    
    # All 6 pairs
    all_pairs = [('A', 'B'), ('A', 'C'), ('A', 'D'), ('B', 'C'), ('B', 'D'), ('C', 'D')]
    results = {}
    
    print("\n" + "-"*70)
    print("TEMPORAL SHEAR RESULTS (ALL PAIRS)")
    print("-"*70)
    
    for img1, img2 in all_pairs:
        pair_name = f"{img1}-{img2}"
        
        # Compute broadband delay
        bd = delays_to_C[img1] - delays_to_C[img2]
        
        gamma, sigma, details = compute_temporal_shear(
            t, images[img1], t, images[img2], 
            broadband_delay=bd,
            tau_values=[50, 100, 200, 400, 800, 1200]
        )
        
        results[pair_name] = {
            'gamma': float(gamma) if not np.isnan(gamma) else None,
            'sigma': float(sigma) if not np.isnan(sigma) else None,
            'broadband_delay': bd,
            'n_tau_points': len(details),
            'details': [(int(t), float(d), float(c)) for t,d,c in details] if details else []
        }
        
        print(f"\n{pair_name}:")
        print(f"  Broadband delay: {bd:.1f} days")
        
        if not np.isnan(gamma):
            sig = abs(gamma) / sigma if sigma > 0 else 0
            print(f"  Temporal shear Γ = {gamma:.1f} ± {sigma:.1f} days/decade ({sig:.1f}σ)")
            if details:
                print(f"  Scale-dependent delays ({len(details)} points):")
                for tau, delay, corr in details:
                    print(f"    τ={tau:5d}d: Δt = {delay:8.1f}d (r={corr:.3f})")
        else:
            print(f"  Insufficient data for fit ({len(details)} valid tau points)")
    
    # Summary statistics
    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)
    
    valid_results = {k: v for k, v in results.items() if v['gamma'] is not None}
    
    print(f"\nValid pairs: {len(valid_results)}/{len(all_pairs)}")
    
    if valid_results:
        gammas = [v['gamma'] for v in valid_results.values()]
        sigmas = [v['sigma'] for v in valid_results.values()]
        
        print(f"\nIndividual Γ values:")
        for pair, res in valid_results.items():
            sig = abs(res['gamma']) / res['sigma'] if res['sigma'] > 0 else 0
            print(f"  {pair}: Γ = {res['gamma']:+.1f} ± {res['sigma']:.1f} days/dec ({sig:.1f}σ)")
        
        # Weighted mean
        weights = [1/s**2 if s > 0 else 0 for s in sigmas]
        if sum(weights) > 0:
            weighted_mean = sum(g*w for g,w in zip(gammas, weights)) / sum(weights)
            weighted_err = 1 / np.sqrt(sum(weights))
            print(f"\nWeighted mean Γ = {weighted_mean:.1f} ± {weighted_err:.1f} days/decade")
        
        # Sign consistency
        n_positive = sum(1 for g in gammas if g > 0)
        n_negative = sum(1 for g in gammas if g < 0)
        print(f"\nSign distribution: {n_positive} positive, {n_negative} negative")
        
        # TEP prediction comparison
        geom_factor = (1 + z_source) / (1 + z_lens)
        desj0408_gamma = 33  # From DESJ0408
        desj0408_factor = (1 + 2.375) / (1 + 0.597)
        predicted = desj0408_gamma * (geom_factor / desj0408_factor)
        
        print(f"\nTEP Prediction (scaled from DESJ0408):")
        print(f"  Geometric factor: {geom_factor:.2f}")
        print(f"  Predicted |Γ|: ~{predicted:.0f} days/decade")
        print(f"  Observed mean |Γ|: {np.mean(np.abs(gammas)):.1f} days/decade")
        
        ratio = np.mean(np.abs(gammas)) / predicted if predicted > 0 else 0
        print(f"  Ratio (observed/predicted): {ratio:.2f}")
    
    # Save results
    output = {
        'system': 'SDSS J1004+4112',
        'z_lens': z_lens,
        'z_source': z_source,
        'n_epochs': len(df),
        'baseline_years': float((df['MJD'].max() - df['MJD'].min()) / 365.25),
        'delays_to_C': delays_to_C,
        'temporal_shear_results': results,
        'summary': {
            'n_valid_pairs': len(valid_results),
            'gammas': [v['gamma'] for v in valid_results.values()] if valid_results else [],
            'mean_abs_gamma': float(np.mean(np.abs([v['gamma'] for v in valid_results.values()]))) if valid_results else None
        }
    }
    
    output_path = OUTPUT_DIR / "j1004_full_analysis_results.json"
    with open(output_path, 'w') as f:
        json.dump(output, f, indent=2, default=str)
    
    print(f"\nResults saved to: {output_path}")
    
    return output

if __name__ == "__main__":
    main()
