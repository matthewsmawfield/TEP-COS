#!/usr/bin/env python3
"""
Step 3.0 (Stable): COSMOGRAIL Temporal Shear Analysis with Mode Locking

This script re-runs the temporal shear analysis but enforces a "Mode Lock":
Multiscale time delays are constrained to be within ±50 days of the 
broadband (global) delay. This prevents the estimator from jumping to 
alias modes (e.g. seasonal gaps) at small smoothing scales.

If the "Temporal Shear" signal is real, it should persist (perhaps slightly weaker).
If the signal is an artifact of mode jumping, it will vanish.
"""

import argparse
import json
import logging
import sys
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
from scipy import interpolate, stats
from scipy.ndimage import gaussian_filter1d

# Reuse core classes/functions from step 3.0 to ensure consistency
from step_3_0_cosmograil_temporal_shear import (
    LightCurve, 
    LensSystem, 
    parse_rdb_file, 
    detrend_lightcurve, 
    estimate_delay_correlation,
    fit_gamma,
    bandpass_filter,
    print_status
)

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="[%(asctime)s] [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)

def compute_stable_multiscale_delays(
    lc1: LightCurve,
    lc2: LightCurve,
    broadband_delay: float,
    tau_values: List[float],
    lock_window: float = 50.0,
    min_correlation: float = 0.3,
) -> Dict[float, Tuple[float, float, float]]:
    """
    Compute time delay at multiple variability timescales with Mode Locking.
    The search range is restricted to [broadband_delay - lock_window, broadband_delay + lock_window].
    """
    results = {}
    
    # Define restricted lag range centered on broadband delay
    lag_range = (broadband_delay - lock_window, broadband_delay + lock_window)
    
    for tau in tau_values:
        # Bandpass filter both curves
        lc1_bp = bandpass_filter(lc1, tau)
        lc2_bp = bandpass_filter(lc2, tau)
        
        # Estimate delay in restricted range
        delay, corr, err = estimate_delay_correlation(lc1_bp, lc2_bp, lag_range)
        
        # Reject low-correlation estimates as unreliable
        if not np.isfinite(corr) or corr < min_correlation:
            delay = np.nan
            err = np.nan
        
        results[tau] = (delay, corr, err)
    
    return results

def analyze_system_stable(
    system: LensSystem,
    detrend_window: float = 200.0,
    tau_values: Optional[List[float]] = None,
    global_lag_range: Tuple[float, float] = (-200, 200),
    lock_window: float = 50.0,
) -> Dict:
    """
    Full temporal shear analysis with Mode Locking.
    """
    if tau_values is None:
        tau_values = [5, 10, 20, 40, 80, 160]
    
    results = {
        "system_id": system.system_id,
        "mode": "stable_mode_locked",
        "lock_window_days": lock_window,
        "pairs": {},
    }
    
    # Detrend all light curves
    detrended = {
        label: detrend_lightcurve(lc, detrend_window)
        for label, lc in system.light_curves.items()
    }
    
    # Analyze each image pair
    for l1, l2 in system.get_image_pairs():
        pair_key = f"{l1}-{l2}"
        
        lc1 = detrended[l1]
        lc2 = detrended[l2]
        
        # 1. Broadband delay (Global search)
        delay_bb, corr_bb, err_bb = estimate_delay_correlation(lc1, lc2, global_lag_range)
        
        if not np.isfinite(delay_bb):
            print_status(f"  Pair {pair_key}: Failed broadband delay", "WARNING")
            continue
            
        # 2. Stable Multi-scale delays (Restricted search)
        multiscale = compute_stable_multiscale_delays(
            lc1, lc2, 
            broadband_delay=delay_bb, 
            tau_values=tau_values, 
            lock_window=lock_window
        )
        
        # Extract for Gamma fit
        delays = [multiscale[tau][0] for tau in tau_values]
        corrs = [multiscale[tau][1] for tau in tau_values]
        errs = [multiscale[tau][2] for tau in tau_values]
        
        # Fit Gamma
        gamma, gamma_err, intercept, r_sq = fit_gamma(
            tau_values, delays, errs, correlations=corrs, min_valid_points=4
        )
        
        # Significance
        if np.isfinite(gamma) and np.isfinite(gamma_err) and gamma_err > 0:
            gamma_sigma = abs(gamma) / gamma_err
        else:
            gamma_sigma = np.nan
        
        results["pairs"][pair_key] = {
            "broadband": {
                "delay_days": delay_bb,
                "correlation": corr_bb,
            },
            "multiscale": {
                str(tau): {
                    "delay_days": multiscale[tau][0],
                    "uncertainty_days": multiscale[tau][2],
                }
                for tau in tau_values
            },
            "gamma": {
                "value": gamma,
                "uncertainty": gamma_err,
                "sigma": gamma_sigma,
                "r_squared": r_sq,
            },
        }
    
    return results

def main():
    data_dir = Path("data/cosmograil")
    output_dir = Path("results/outputs")
    output_file = output_dir / "step_3_0_stable_gamma_analysis.json"
    
    # Systems to check (Focus on the 'detections' first)
    # Plus Q0957 (needs separate loader, but we can do COSMOGRAIL first)
    
    rdb_files = sorted(data_dir.glob("*.rdb"))
    tau_values = [5, 10, 20, 40, 80, 160]
    
    all_results = {}
    
    print_status("Starting STABLE (Mode-Locked) Analysis...")
    
    for rdb_file in rdb_files:
        try:
            system = parse_rdb_file(rdb_file)
            if system.n_images < 2: continue
            
            print_status(f"Processing {system.system_id}")
            res = analyze_system_stable(system, tau_values=tau_values)
            all_results[system.system_id] = res
            
        except Exception as e:
            print_status(f"Error {rdb_file.name}: {e}", "ERROR")
            
    # Save
    output_dir.mkdir(parents=True, exist_ok=True)
    with open(output_file, "w") as f:
        json.dump(all_results, f, indent=2, default=lambda x: None if not np.isfinite(x) else x)
        
    # Print comparison for key pairs
    targets = [
        ('DESJ0408', 'A-D'),
        ('DESJ0408', 'B-D'),
        ('PG1115', 'B-C'),
        ('PG1115', 'A-B'),
        ('J1206', 'A-B')
    ]
    
    print("\n" + "="*60)
    print(f"{'Pair':<15} {'Gamma (Stable)':<20} {'Verdict'}")
    print("-" * 60)
    
    for sys_id, pair_id in targets:
        if sys_id in all_results and pair_id in all_results[sys_id]['pairs']:
            g = all_results[sys_id]['pairs'][pair_id]['gamma']
            val = g['value']
            sig = g['sigma']
            
            if val is not None:
                verdict = "VANISHED" if sig < 2 else "PERSISTS"
                print(f"{sys_id} {pair_id:<5} {val:>6.1f} (sig={sig:.1f})      {verdict}")
            else:
                print(f"{sys_id} {pair_id:<5}    NaN")
    print("="*60)

if __name__ == "__main__":
    main()
