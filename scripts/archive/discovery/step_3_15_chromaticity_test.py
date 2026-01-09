#!/usr/bin/env python3
"""
Step 3.15: Multi-band Chromaticity Test for Temporal Shear

Tests the achromaticity prediction of TEP: if temporal shear is gravitational,
it should be wavelength-independent (Γ_blue ≈ Γ_red).

Systems with multi-band data:
- Q2237+0305: g, r, V, I bands
- HE1104-1805: B, R, I, J bands
- HE0435-1223: V, R bands
- Q2237 (Vakulik): V, R, I bands
"""

import numpy as np
import pandas as pd
from scipy import signal
from scipy.interpolate import interp1d
from scipy.stats import pearsonr
import json
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

DATA_DIR = Path("/Users/matthewsmawfield/www/TEP-COS/data/cosmograil")
OUTPUT_DIR = Path("/Users/matthewsmawfield/www/TEP-COS/results/outputs")

def load_q2237_band(band):
    """Load Q2237 data for a specific band."""
    filepath = DATA_DIR / f"q2237_JAA637A89_{band}.csv"
    if not filepath.exists():
        return None
    df = pd.read_csv(filepath)
    return df

def load_he1104_band(band):
    """Load HE1104 data for a specific band."""
    filepath = DATA_DIR / f"he1104_JApJ798_95_{band}.csv"
    if not filepath.exists():
        return None
    df = pd.read_csv(filepath)
    return df

def load_he0435_band(band):
    """Load HE0435 data for a specific band."""
    filepath = DATA_DIR / f"he0435_JAA703A250_{band}.csv"
    if not filepath.exists():
        return None
    df = pd.read_csv(filepath)
    return df

def load_q2237_vakulik_band(band):
    """Load Q2237 Vakulik data for a specific band."""
    filepath = DATA_DIR / f"q2237_vakulik_JAA420_447_{band}.csv"
    if not filepath.exists():
        return None
    df = pd.read_csv(filepath)
    return df

def gaussian_smooth(t, y, sigma):
    """Apply Gaussian smoothing to remove slow trends."""
    if len(t) < 5:
        return np.zeros_like(y)
    
    # Create interpolated function
    t_min, t_max = t.min(), t.max()
    t_grid = np.linspace(t_min, t_max, 500)
    
    try:
        f = interp1d(t, y, kind='linear', fill_value='extrapolate')
        y_grid = f(t_grid)
        
        # Gaussian kernel
        dt = t_grid[1] - t_grid[0]
        kernel_size = int(3 * sigma / dt)
        if kernel_size < 3:
            kernel_size = 3
        kernel = signal.windows.gaussian(kernel_size * 2 + 1, sigma / dt)
        kernel /= kernel.sum()
        
        # Convolve
        y_smooth = np.convolve(y_grid, kernel, mode='same')
        
        # Interpolate back
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

def estimate_delay_iccf(t1, y1, t2, y2, search_range=(-50, 50)):
    """Estimate time delay using interpolated cross-correlation."""
    # Create common time grid
    t_min = max(t1.min(), t2.min())
    t_max = min(t1.max(), t2.max())
    
    if t_max <= t_min:
        return np.nan, np.nan
    
    t_common = np.linspace(t_min, t_max, 200)
    
    try:
        f1 = interp1d(t1, y1, kind='linear', bounds_error=False, fill_value=np.nan)
        f2 = interp1d(t2, y2, kind='linear', bounds_error=False, fill_value=np.nan)
    except:
        return np.nan, np.nan
    
    y1_interp = f1(t_common)
    y2_interp = f2(t_common)
    
    # Remove NaNs
    valid = ~(np.isnan(y1_interp) | np.isnan(y2_interp))
    if valid.sum() < 10:
        return np.nan, np.nan
    
    y1_valid = y1_interp[valid]
    y2_valid = y2_interp[valid]
    t_valid = t_common[valid]
    
    # Normalize
    y1_norm = (y1_valid - y1_valid.mean()) / (y1_valid.std() + 1e-10)
    y2_norm = (y2_valid - y2_valid.mean()) / (y2_valid.std() + 1e-10)
    
    # Search for best delay
    delays = np.linspace(search_range[0], search_range[1], 101)
    correlations = []
    
    for delay in delays:
        t_shifted = t_valid + delay
        try:
            f2_shift = interp1d(t2, y2, kind='linear', bounds_error=False, fill_value=np.nan)
            y2_shifted = f2_shift(t_shifted)
            valid_shift = ~np.isnan(y2_shifted)
            if valid_shift.sum() < 5:
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

def compute_temporal_shear(t1, y1, t2, y2, tau_values=[20, 40, 80, 160], detrend_sigma=200):
    """Compute temporal shear Γ for an image pair."""
    # Detrend
    trend1 = gaussian_smooth(t1, y1, detrend_sigma)
    trend2 = gaussian_smooth(t2, y2, detrend_sigma)
    
    y1_detrend = y1 - trend1
    y2_detrend = y2 - trend2
    
    delays = []
    log_taus = []
    
    for tau in tau_values:
        # Bandpass filter
        y1_filt = bandpass_filter(t1, y1_detrend, tau)
        y2_filt = bandpass_filter(t2, y2_detrend, tau)
        
        # Estimate delay
        delay, corr = estimate_delay_iccf(t1, y1_filt, t2, y2_filt)
        
        if not np.isnan(delay) and corr > 0.3:
            delays.append(delay)
            log_taus.append(np.log10(tau))
    
    if len(delays) < 2:
        return np.nan, np.nan
    
    # Linear fit: delay = Γ * log10(tau) + offset
    delays = np.array(delays)
    log_taus = np.array(log_taus)
    
    A = np.vstack([log_taus, np.ones(len(log_taus))]).T
    try:
        result = np.linalg.lstsq(A, delays, rcond=None)
        gamma = result[0][0]
        
        # Estimate uncertainty from residuals
        residuals = delays - (gamma * log_taus + result[0][1])
        if len(residuals) > 2:
            sigma_gamma = np.std(residuals) / np.sqrt(len(residuals))
        else:
            sigma_gamma = np.nan
        
        return gamma, sigma_gamma
    except:
        return np.nan, np.nan

def analyze_q2237():
    """Analyze Q2237 chromaticity across 4 bands."""
    print("\n" + "="*60)
    print("Q2237+0305 CHROMATICITY ANALYSIS (g, r, V, I bands)")
    print("="*60)
    
    bands = ['g', 'r', 'V', 'I']
    results = {}
    
    for band in bands:
        df = load_q2237_band(band)
        if df is None:
            print(f"  {band}-band: No data")
            continue
        
        t = df['MJD'].values
        
        # Analyze A-B pair
        if 'mA' in df.columns and 'mB' in df.columns:
            yA = df['mA'].values
            yB = df['mB'].values
            
            gamma, sigma = compute_temporal_shear(t, yA, t, yB)
            results[f'{band}_AB'] = {'gamma': gamma, 'sigma': sigma}
            print(f"  {band}-band A-B: Γ = {gamma:.1f} ± {sigma:.1f} days/decade")
    
    # Compute chromaticity
    if 'g_AB' in results and 'r_AB' in results:
        dGamma_gr = results['g_AB']['gamma'] - results['r_AB']['gamma']
        sigma_gr = np.sqrt(results['g_AB']['sigma']**2 + results['r_AB']['sigma']**2)
        print(f"\n  ΔΓ(g-r) = {dGamma_gr:.1f} ± {sigma_gr:.1f} days/decade")
        results['chromaticity_gr'] = {'delta_gamma': dGamma_gr, 'sigma': sigma_gr}
    
    if 'V_AB' in results and 'I_AB' in results:
        dGamma_VI = results['V_AB']['gamma'] - results['I_AB']['gamma']
        sigma_VI = np.sqrt(results['V_AB']['sigma']**2 + results['I_AB']['sigma']**2)
        print(f"  ΔΓ(V-I) = {dGamma_VI:.1f} ± {sigma_VI:.1f} days/decade")
        results['chromaticity_VI'] = {'delta_gamma': dGamma_VI, 'sigma': sigma_VI}
    
    return results

def analyze_he1104():
    """Analyze HE1104 chromaticity across 4 bands."""
    print("\n" + "="*60)
    print("HE1104-1805 CHROMATICITY ANALYSIS (B, R, I, J bands)")
    print("="*60)
    
    bands = ['B', 'R', 'I', 'J']
    results = {}
    
    for band in bands:
        df = load_he1104_band(band)
        if df is None:
            print(f"  {band}-band: No data")
            continue
        
        # HE1104 format: HJD, A, e_A, B, e_B, Filt
        t = df['HJD'].values - 2450000  # Convert to MJD-like
        yA = df['A'].values
        yB = df['B'].values
        
        gamma, sigma = compute_temporal_shear(t, yA, t, yB)
        results[f'{band}_AB'] = {'gamma': gamma, 'sigma': sigma}
        print(f"  {band}-band A-B: Γ = {gamma:.1f} ± {sigma:.1f} days/decade")
    
    # Compute chromaticity relative to R-band
    if 'R_AB' in results:
        for band in ['B', 'I', 'J']:
            if f'{band}_AB' in results:
                dGamma = results[f'{band}_AB']['gamma'] - results['R_AB']['gamma']
                sigma = np.sqrt(results[f'{band}_AB']['sigma']**2 + results['R_AB']['sigma']**2)
                print(f"\n  ΔΓ({band}-R) = {dGamma:.1f} ± {sigma:.1f} days/decade")
                results[f'chromaticity_{band}R'] = {'delta_gamma': dGamma, 'sigma': sigma}
    
    return results

def analyze_he0435():
    """Analyze HE0435 chromaticity (V, R bands)."""
    print("\n" + "="*60)
    print("HE0435-1223 CHROMATICITY ANALYSIS (V, R bands)")
    print("="*60)
    
    bands = ['V', 'R']
    results = {}
    
    for band in bands:
        df = load_he0435_band(band)
        if df is None:
            print(f"  {band}-band: No data")
            continue
        
        # Check columns
        print(f"  {band}-band columns: {df.columns.tolist()}")
        
        if len(df.columns) >= 3:
            t = df.iloc[:, 0].values
            yA = df.iloc[:, 1].values
            yB = df.iloc[:, 3].values if len(df.columns) > 3 else df.iloc[:, 2].values
            
            gamma, sigma = compute_temporal_shear(t, yA, t, yB)
            results[f'{band}_AB'] = {'gamma': gamma, 'sigma': sigma}
            print(f"  {band}-band A-B: Γ = {gamma:.1f} ± {sigma:.1f} days/decade")
    
    # Compute chromaticity
    if 'V_AB' in results and 'R_AB' in results:
        dGamma_VR = results['V_AB']['gamma'] - results['R_AB']['gamma']
        sigma_VR = np.sqrt(results['V_AB']['sigma']**2 + results['R_AB']['sigma']**2)
        print(f"\n  ΔΓ(V-R) = {dGamma_VR:.1f} ± {sigma_VR:.1f} days/decade")
        results['chromaticity_VR'] = {'delta_gamma': dGamma_VR, 'sigma': sigma_VR}
    
    return results

def analyze_q2237_vakulik():
    """Analyze Q2237 Vakulik chromaticity (V, R, I bands)."""
    print("\n" + "="*60)
    print("Q2237+0305 VAKULIK CHROMATICITY ANALYSIS (V, R, I bands)")
    print("="*60)
    print("  Note: This system has ~zero physical delay (null control)")
    
    bands = ['V', 'R', 'I']
    results = {}
    
    for band in bands:
        df = load_q2237_vakulik_band(band)
        if df is None:
            print(f"  {band}-band: No data")
            continue
        
        print(f"  {band}-band: {len(df)} data points")
        print(f"    Columns: {df.columns.tolist()}")
    
    return results

def main():
    print("="*60)
    print("MULTI-BAND CHROMATICITY TEST FOR TEMPORAL SHEAR")
    print("TEP Prediction: Achromatic (ΔΓ ≈ 0)")
    print("Microlensing: Chromatic (ΔΓ ≠ 0)")
    print("="*60)
    
    all_results = {}
    
    # Q2237 (4 bands)
    all_results['Q2237'] = analyze_q2237()
    
    # HE1104 (4 bands)
    all_results['HE1104'] = analyze_he1104()
    
    # HE0435 (2 bands)
    all_results['HE0435'] = analyze_he0435()
    
    # Q2237 Vakulik (null control)
    all_results['Q2237_Vakulik'] = analyze_q2237_vakulik()
    
    # Summary
    print("\n" + "="*60)
    print("CHROMATICITY SUMMARY")
    print("="*60)
    
    chromaticity_tests = []
    
    for system, results in all_results.items():
        for key, val in results.items():
            if 'chromaticity' in key:
                dG = val['delta_gamma']
                sig = val['sigma']
                significance = abs(dG) / sig if sig > 0 else 0
                status = "ACHROMATIC" if significance < 2 else "CHROMATIC"
                print(f"  {system} {key}: ΔΓ = {dG:.1f} ± {sig:.1f} ({significance:.1f}σ) → {status}")
                chromaticity_tests.append({
                    'system': system,
                    'test': key,
                    'delta_gamma': dG,
                    'sigma': sig,
                    'significance': significance,
                    'status': status
                })
    
    # Save results
    output = {
        'all_results': {},
        'chromaticity_tests': chromaticity_tests,
        'summary': {
            'n_tests': len(chromaticity_tests),
            'n_achromatic': sum(1 for t in chromaticity_tests if t['status'] == 'ACHROMATIC'),
            'n_chromatic': sum(1 for t in chromaticity_tests if t['status'] == 'CHROMATIC')
        }
    }
    
    # Convert numpy types for JSON
    for system, results in all_results.items():
        output['all_results'][system] = {}
        for key, val in results.items():
            if isinstance(val, dict):
                output['all_results'][system][key] = {
                    k: float(v) if isinstance(v, (np.floating, np.integer)) else v 
                    for k, v in val.items()
                }
    
    output_path = OUTPUT_DIR / "chromaticity_test_results.json"
    with open(output_path, 'w') as f:
        json.dump(output, f, indent=2, default=str)
    
    print(f"\nResults saved to: {output_path}")
    print(f"\nConclusion: {output['summary']['n_achromatic']}/{output['summary']['n_tests']} tests are achromatic (< 2σ)")
    
    return output

if __name__ == "__main__":
    main()
