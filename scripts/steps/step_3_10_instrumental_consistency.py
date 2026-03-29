#!/usr/bin/env python3
"""
Step 3.10: Instrumental Consistency Check (WFI2033)

Tests if the Temporal Shear signal is consistent across different telescopes/instruments.
TEP Prediction: Signal is physical -> Consistent across instruments.
Systematic Error Prediction: Signal is due to instrument artifacts -> Inconsistent.

Data:
- WFI2033: EulerCAM (ecam) vs SMARTS
- Period: Overlapping monitoring campaigns

Methodology:
1. Load WFI2033 light curves from both instruments.
2. Compute Gamma independently for each instrument.
3. Check for consistency (Gamma_ecam vs Gamma_smarts).
"""

import numpy as np
from pathlib import Path
import json
from scipy import stats, interpolate
from dataclasses import dataclass
from typing import Dict, List, Tuple
from scipy.ndimage import gaussian_filter1d

@dataclass
class LightCurve:
    label: str
    t: np.ndarray
    mag: np.ndarray
    magerr: np.ndarray

@dataclass
class LensSystem:
    system_id: str
    light_curves: Dict[str, LightCurve]

DATA_DIR = Path("data/cosmograil")
RESULTS_DIR = Path("results/outputs")

def load_wfi2033_instrument(filename: str, label_prefix: str):
    fpath = DATA_DIR / filename
    if not fpath.exists():
        print(f"Missing {fpath}")
        return None
        
    try:
        # Format: MHJD A errA B errB C errC [D errD]
        raw = np.loadtxt(fpath)
        
        cols = raw.shape[1]
        lcs = {}
        t = raw[:, 0]
        
        # 7 columns -> A, B, C
        if cols == 7:
            images = ['A', 'B', 'C']
            for i, img in enumerate(images):
                idx = 1 + i*2
                lcs[img] = LightCurve(img, t, raw[:, idx], raw[:, idx+1])
        # 9 columns -> A, B, C, D
        elif cols >= 9:
            images = ['A', 'B', 'C', 'D']
            for i, img in enumerate(images):
                idx = 1 + i*2
                lcs[img] = LightCurve(img, t, raw[:, idx], raw[:, idx+1])
        else:
            print(f"Unknown format for {filename} with {cols} columns")
            return None
            
        return LensSystem(label_prefix, lcs)
        
    except Exception as e:
        print(f"Error loading {filename}: {e}")
        return None

def analyze_gamma(system: LensSystem):
    # Standard Gamma estimation
    taus = [10, 20, 40, 80] # WFI2033 allows slightly longer scales?
    lcs = system.light_curves
    labels = sorted(lcs.keys())
    
    gammas = []
    
    def get_filtered(lc, tau):
        if len(lc.t) < 10: return None
        
        # Fill gaps
        dt = 1.0
        tg = np.arange(lc.t.min(), lc.t.max(), dt)
        if len(tg) < 50: return None
        
        y = interpolate.interp1d(lc.t, lc.mag, bounds_error=False, fill_value=np.nan)(tg)
        mask = np.isnan(y)
        if np.sum(mask) / len(mask) > 0.5: return None
        if np.all(mask): return None
        
        y[mask] = np.interp(np.flatnonzero(mask), np.flatnonzero(~mask), y[~mask])
        
        s = tau / 2.35482 / dt  # Precise FWHM to sigma (2*sqrt(2*ln(2)))
        f = gaussian_filter1d(y, s)
        return f

    for i, l1 in enumerate(labels):
        for l2 in labels[i+1:]:
            delays = []
            valid_taus = []
            
            c1 = lcs[l1]
            c2 = lcs[l2]
            
            for tau in taus:
                y1 = get_filtered(c1, tau)
                y2 = get_filtered(c2, tau)
                
                if y1 is None or y2 is None: continue
                
                n_ov = min(len(y1), len(y2))
                if n_ov < 2 * tau: continue
                
                lags = np.arange(-30, 30, 0.5)
                corrs = []
                for l in lags:
                    shift = int(l)
                    if shift > 0:
                        v1 = y1[shift:]
                        v2 = y2[:-shift]
                    elif shift < 0:
                        v1 = y1[:shift]
                        v2 = y2[-shift:]
                    else:
                        v1 = y1
                        v2 = y2
                    
                    if len(v1) < 20: 
                        corrs.append(0)
                        continue
                        
                    r = np.corrcoef(v1, v2)[0,1]
                    corrs.append(r)
                
                if np.max(corrs) > 0.4:
                    best_lag = lags[np.argmax(corrs)]
                    delays.append(best_lag)
                    valid_taus.append(tau)
            
            if len(valid_taus) >= 3:
                res = stats.linregress(np.log10(valid_taus), delays)
                gammas.append(res.slope)

    if not gammas:
        return None
        
    return {
        'gamma_mean': np.mean(gammas),
        'gamma_std': np.std(gammas),
        'n_pairs': len(gammas)
    }

def main():
    print("Running WFI2033 Instrumental Consistency Analysis...")
    
    ecam = load_wfi2033_instrument("WFI2033_ecam.dat", "WFI2033_ECAM")
    smarts = load_wfi2033_instrument("WFI2033_smarts.dat", "WFI2033_SMARTS")
    
    res_ecam = analyze_gamma(ecam) if ecam else None
    res_smarts = analyze_gamma(smarts) if smarts else None
    
    print("\nResults:")
    if res_ecam:
        print(f"EulerCAM: Gamma = {res_ecam['gamma_mean']:.1f} +/- {res_ecam['gamma_std']:.1f} (N={res_ecam['n_pairs']})")
    if res_smarts:
        print(f"SMARTS:   Gamma = {res_smarts['gamma_mean']:.1f} +/- {res_smarts['gamma_std']:.1f} (N={res_smarts['n_pairs']})")
        
    if res_ecam and res_smarts:
        diff = abs(res_ecam['gamma_mean'] - res_smarts['gamma_mean'])
        print(f"Difference: {diff:.1f}")
        
        # Save
        out = {"ecam": res_ecam, "smarts": res_smarts, "diff": diff}
        with open(RESULTS_DIR / "step_3_10_wfi2033_consistency.json", 'w') as f:
            # Convert numpy types
            def convert(o):
                if isinstance(o, np.generic): return o.item()
                raise TypeError
            json.dump(out, f, indent=2, default=convert)

if __name__ == "__main__":
    main()
