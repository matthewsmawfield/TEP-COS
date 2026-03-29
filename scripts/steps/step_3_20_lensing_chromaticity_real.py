#!/usr/bin/env python3
"""
Step 3.20: Real Multi-Band Lensing Chromaticity Analysis
========================================================

Performs actual temporal shear measurements per band on available multi-band data,
then tests achromaticity (TEP prediction: ΔΓ = 0) vs chromaticity (microlensing).

Systems analyzed:
- Q2237 (Goicoechea 2020): g, r, V, I bands
- HE0435 (Sorgenfrei 2025): R, V bands

Author: TEP Collaboration
Date: March 2026
"""

import json
import sys
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import numpy as np
from scipy import stats
import pandas as pd

# Add repo root to path for imports
REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "scripts" / "steps"))

from step_3_0_cosmograil_temporal_shear import (
    LightCurve, LensSystem, detrend_lightcurve,
    estimate_delay_iccf, compute_multiscale_delays, fit_gamma,
    parse_multiband_csv
)

# Configuration
DATA_DIR = REPO_ROOT / "data" / "cosmograil"
RESULTS_DIR = REPO_ROOT / "results" / "outputs"
OUTPUT_JSON = RESULTS_DIR / "step_3_20_lensing_chromaticity.json"

RANDOM_SEED = 42
np.random.seed(RANDOM_SEED)

# Analysis parameters
TAU_VALUES = [40.0, 80.0, 160.0]
DETREND_WINDOW = 200.0
MODE_LOCK_WINDOW = 50.0
MIN_VARIANCE_FRAC = 0.02


def analyze_system_band(system_name: str, band: str, data_file: Path) -> Optional[Dict]:
    """Analyze temporal shear for a single system-band combination."""
    
    # Parse the data
    lens_system = parse_multiband_csv(data_file)
    if lens_system is None or len(lens_system.light_curves) < 2:
        return None
    
    results = {
        "system": system_name,
        "band": band,
        "n_images": len(lens_system.light_curves),
        "image_labels": lens_system.image_labels,
        "pairs": {}
    }
    
    # Analyze each image pair
    for img1, img2 in lens_system.get_image_pairs():
        pair_key = f"{img1}-{img2}"
        lc1 = lens_system.light_curves[img1]
        lc2 = lens_system.light_curves[img2]
        
        # Detrend to remove microlensing
        lc1_dt = detrend_lightcurve(lc1, DETREND_WINDOW)
        lc2_dt = detrend_lightcurve(lc2, DETREND_WINDOW)
        
        # Measure broadband delay
        bb_delay, bb_corr, bb_err = estimate_delay_iccf(
            lc1_dt, lc2_dt, lag_range=(-200, 200), lag_step=1.0
        )
        
        if not np.isfinite(bb_delay) or not np.isfinite(bb_corr):
            results["pairs"][pair_key] = {"status": "FAILED", "reason": "No broadband correlation"}
            continue
        
        # Multi-scale analysis
        ms_results = compute_multiscale_delays(
            lc1_dt, lc2_dt,
            tau_values=TAU_VALUES,
            broadband_delay=bb_delay,
            mode_lock_window=MODE_LOCK_WINDOW,
            min_variance_fraction=MIN_VARIANCE_FRAC,
            estimator="iccf"
        )
        
        # Extract delays and fit gamma
        tau_list = []
        delay_list = []
        err_list = []
        corr_list = []
        
        for tau in TAU_VALUES:
            delay, corr, err = ms_results.get(tau, (np.nan, np.nan, np.nan))
            if np.isfinite(delay) and np.isfinite(corr) and corr > 0.3:
                tau_list.append(tau)
                delay_list.append(delay)
                err_list.append(err if np.isfinite(err) else 10.0)
                corr_list.append(corr)
        
        if len(tau_list) < 2:
            results["pairs"][pair_key] = {
                "status": "INSUFFICIENT_POINTS",
                "n_valid": len(tau_list),
                "broadband_delay": float(bb_delay),
                "broadband_correlation": float(bb_corr)
            }
            continue
        
        # Fit gamma
        gamma, gamma_err, intercept, r_sq = fit_gamma(tau_list, delay_list, err_list, corr_list)
        
        results["pairs"][pair_key] = {
            "status": "OK",
            "broadband_delay": float(bb_delay),
            "broadband_correlation": float(bb_corr),
            "gamma": {
                "value": float(gamma),
                "uncertainty": float(gamma_err),
                "sigma": abs(gamma/gamma_err) if gamma_err > 0 else np.nan,
                "r_squared": float(r_sq),
                "n_points": len(tau_list)
            },
            "multiscale": {str(t): {"delay": float(d), "correlation": float(c)} 
                          for t, d, c in zip(tau_list, delay_list, corr_list)}
        }
    
    return results


def test_achromaticity(band_results: List[Dict]) -> Dict:
    """Test if gamma values are consistent across bands (achromatic)."""
    
    # Collect gamma values per band
    gamma_by_band = {}
    for result in band_results:
        band = result["band"]
        gammas = []
        uncertainties = []
        
        for pair_key, pair_data in result.get("pairs", {}).items():
            if pair_data.get("status") == "OK":
                gamma = pair_data["gamma"]["value"]
                unc = pair_data["gamma"]["uncertainty"]
                if np.isfinite(gamma) and np.isfinite(unc) and unc > 0:
                    gammas.append(gamma)
                    uncertainties.append(unc)
        
        if gammas:
            # Weighted average per band
            weights = [1/u**2 for u in uncertainties]
            weighted_mean = sum(g*w for g,w in zip(gammas, weights)) / sum(weights)
            weighted_err = np.sqrt(1/sum(weights))
            gamma_by_band[band] = {
                "gamma_mean": float(weighted_mean),
                "gamma_err": float(weighted_err),
                "n_pairs": len(gammas),
                "individual_gammas": gammas
            }
    
    if len(gamma_by_band) < 2:
        return {
            "status": "INSUFFICIENT_BANDS",
            "n_bands": len(gamma_by_band),
            "message": "Need at least 2 bands for achromaticity test"
        }
    
    # Test consistency: chi-square for agreement
    bands = list(gamma_by_band.keys())
    gammas = [gamma_by_band[b]["gamma_mean"] for b in bands]
    errs = [gamma_by_band[b]["gamma_err"] for b in bands]
    
    # Weighted global mean
    weights = [1/e**2 for e in errs]
    global_mean = sum(g*w for g,w in zip(gammas, weights)) / sum(weights)
    
    # Chi-square test for consistency
    chi2 = sum(((g - global_mean)/e)**2 for g, e in zip(gammas, errs))
    ndof = len(bands) - 1
    p_consistency = 1 - stats.chi2.cdf(chi2, ndof) if ndof > 0 else np.nan
    
    # Band differences
    band_diffs = {}
    for i, b1 in enumerate(bands):
        for b2 in bands[i+1:]:
            g1, e1 = gamma_by_band[b1]["gamma_mean"], gamma_by_band[b1]["gamma_err"]
            g2, e2 = gamma_by_band[b2]["gamma_mean"], gamma_by_band[b2]["gamma_err"]
            diff = g1 - g2
            diff_err = np.sqrt(e1**2 + e2**2)
            sigma = diff / diff_err if diff_err > 0 else np.nan
            band_diffs[f"{b1}-{b2}"] = {
                "delta_gamma": float(diff),
                "uncertainty": float(diff_err),
                "sigma": float(sigma),
                "consistent_at_2sigma": abs(sigma) < 2 if np.isfinite(sigma) else None
            }
    
    return {
        "status": "OK",
        "n_bands": len(bands),
        "bands": bands,
        "gamma_by_band": gamma_by_band,
        "global_mean": float(global_mean),
        "chi2": float(chi2),
        "ndof": ndof,
        "p_consistency": float(p_consistency),
        "achromatic_at_2sigma": chi2 < 4.0,  # Rough 2-sigma for 1 dof
        "band_differences": band_diffs,
        "interpretation": "achromatic" if chi2 < 4.0 else "inconsistent_or_chromatic"
    }


def analyze_all_systems():
    """Main analysis pipeline."""
    print("=" * 70)
    print("STEP 3.20: REAL MULTI-BAND CHROMATICITY ANALYSIS")
    print("=" * 70)
    
    # Define multi-band datasets
    multiband_systems = {
        "Q2237": {
            "bands": {
                "g": DATA_DIR / "q2237_JAA637A89_g.csv",
                "r": DATA_DIR / "q2237_JAA637A89_r.csv",
                "V": DATA_DIR / "q2237_JAA637A89_V.csv",
                "I": DATA_DIR / "q2237_JAA637A89_I.csv",
            },
            "z_lens": 0.039,
            "z_source": 1.695
        },
        "HE0435": {
            "bands": {
                "R": DATA_DIR / "he0435_JAA703A250_R.csv",
                "V": DATA_DIR / "he0435_JAA703A250_V.csv",
            },
            "z_lens": 0.454,
            "z_source": 1.693
        }
    }
    
    all_results = {}
    
    for system_name, system_info in multiband_systems.items():
        print(f"\n{'='*70}")
        print(f"System: {system_name}")
        print(f"z_lens = {system_info['z_lens']}, z_source = {system_info['z_source']}")
        print(f"{'='*70}")
        
        band_results = []
        
        for band, filepath in system_info["bands"].items():
            if not filepath.exists():
                print(f"  [{band}] File not found: {filepath}")
                continue
            
            print(f"\n  [{band}] Analyzing {filepath.name}...")
            result = analyze_system_band(system_name, band, filepath)
            
            if result is None:
                print(f"  [{band}] FAILED: Could not parse data")
                continue
            
            # Print summary
            n_ok = sum(1 for p in result["pairs"].values() if p.get("status") == "OK")
            n_total = len(result["pairs"])
            print(f"  [{band}] Successful pairs: {n_ok}/{n_total}")
            
            for pair_key, pair_data in result["pairs"].items():
                if pair_data.get("status") == "OK":
                    g = pair_data["gamma"]["value"]
                    sig = pair_data["gamma"]["sigma"]
                    print(f"    {pair_key}: Γ = {g:+.1f} ± {pair_data['gamma']['uncertainty']:.1f} days/decade ({sig:.1f}σ)")
            
            band_results.append(result)
        
        # Achromaticity test
        if len(band_results) >= 2:
            print(f"\n  Achromaticity test:")
            achro_test = test_achromaticity(band_results)
            
            if achro_test["status"] == "OK":
                print(f"    χ² = {achro_test['chi2']:.2f} (ndof={achro_test['ndof']})")
                print(f"    p(consistency) = {achro_test['p_consistency']:.3f}")
                print(f"    Interpretation: {achro_test['interpretation']}")
                
                for diff_key, diff_data in achro_test.get("band_differences", {}).items():
                    print(f"    ΔΓ({diff_key}) = {diff_data['delta_gamma']:+.1f} ± {diff_data['uncertainty']:.1f} days/decade")
            else:
                print(f"    Status: {achro_test['status']}")
            
            all_results[system_name] = {
                "band_results": band_results,
                "achromaticity_test": achro_test,
                "redshifts": {"lens": system_info["z_lens"], "source": system_info["z_source"]}
            }
        else:
            all_results[system_name] = {
                "band_results": band_results,
                "achromaticity_test": {"status": "INSUFFICIENT_BANDS"},
                "redshifts": {"lens": system_info["z_lens"], "source": system_info["z_source"]}
            }
    
    # Save results
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    output = {
        "timestamp": pd.Timestamp.now().isoformat(),
        "method": "Real multi-band temporal shear analysis",
        "systems": all_results,
        "summary": {
            "n_systems_analyzed": len(all_results),
            "systems": list(all_results.keys())
        }
    }
    
    with open(OUTPUT_JSON, 'w') as f:
        json.dump(output, f, indent=2, default=str)
    
    print(f"\n{'='*70}")
    print(f"Results saved to: {OUTPUT_JSON}")
    print(f"{'='*70}")
    
    return output


if __name__ == "__main__":
    analyze_all_systems()
