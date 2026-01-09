#!/usr/bin/env python3
"""
Step 3.8: Telescope Consistency & Multi-Band Analysis

Tests the robustness of the Temporal Shear signal by comparing measurements 
from different telescopes (and potentially bands) for the same lens systems.

TEP Prediction: The signal is physical (gravitational), so Γ should be consistent 
across different instruments and bands (Achromatic).
Systematics/Microlensing Prediction: Signal might vary with instrument (blending/PSF) 
or band (chromaticity).

Data Analyzed:
1. RXJ1131-1231: Split by telescope (Euler, SMARTS, Mercator, etc.)
2. WFI2033-4723: Split by telescope (EulerCAM vs SMARTS)
3. DESJ0408-5354: Split by telescope (if multiple present)

Methodology:
1. Load RDB files.
2. Split light curves by 'telescope' column.
3. Compute Γ for each telescope subset independently.
4. Compare consistency.

Author: TEP Collaboration
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import json
from scipy import stats, interpolate
from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional
from scipy.ndimage import gaussian_filter1d

# Paths
SCRIPT_DIR = Path(__file__).parent
PROJECT_ROOT = SCRIPT_DIR.parent.parent
DATA_DIR = PROJECT_ROOT / "data" / "cosmograil"
RESULTS_DIR = PROJECT_ROOT / "results" / "outputs"
FIGURE_DIR = PROJECT_ROOT / "results" / "figures" / "consistency"

RESULTS_DIR.mkdir(parents=True, exist_ok=True)
FIGURE_DIR.mkdir(parents=True, exist_ok=True)

# -----------------------------------------------------------------------------
# Core Analysis Classes (Adapted from Step 3.0)
# -----------------------------------------------------------------------------

@dataclass
class LightCurve:
    label: str
    t: np.ndarray
    mag: np.ndarray
    magerr: np.ndarray
    telescope: np.ndarray = None

    def __post_init__(self):
        valid = np.isfinite(self.t) & np.isfinite(self.mag) & np.isfinite(self.magerr)
        if self.telescope is not None:
            self.telescope = self.telescope[valid]
        self.t = self.t[valid]
        self.mag = self.mag[valid]
        self.magerr = self.magerr[valid]
        
        # Sort
        idx = np.argsort(self.t)
        self.t = self.t[idx]
        self.mag = self.mag[idx]
        self.magerr = self.magerr[idx]
        if self.telescope is not None:
            self.telescope = self.telescope[idx]

@dataclass
class LensSystem:
    system_id: str
    light_curves: Dict[str, LightCurve]

    def get_telescopes(self) -> List[str]:
        """Get list of unique telescopes present in the data."""
        telescopes = set()
        for lc in self.light_curves.values():
            if lc.telescope is not None:
                telescopes.update(np.unique(lc.telescope))
        return sorted(list(telescopes))

    def filter_by_telescope(self, telescope_name: str) -> 'LensSystem':
        """Return a new LensSystem with data only from the specified telescope."""
        new_lcs = {}
        for label, lc in self.light_curves.items():
            if lc.telescope is None:
                continue
            mask = (lc.telescope == telescope_name)
            if np.sum(mask) < 10:
                continue
            
            new_lcs[label] = LightCurve(
                label=label,
                t=lc.t[mask],
                mag=lc.mag[mask],
                magerr=lc.magerr[mask],
                telescope=lc.telescope[mask]
            )
        
        if len(new_lcs) < 2:
            return None
            
        return LensSystem(f"{self.system_id}_{telescope_name}", new_lcs)

# -----------------------------------------------------------------------------
# File Parsing
# -----------------------------------------------------------------------------

def parse_rdb_file(filepath: Path) -> LensSystem:
    """Parse COSMOGRAIL .rdb file with telescope column."""
    system_id = filepath.stem.split("_")[0]
    
    with open(filepath, "r") as f:
        lines = f.readlines()
    
    # Find header
    header_line = None
    data_start = 0
    for i, line in enumerate(lines):
        line = line.strip()
        if not line: continue
        if "====" in line:
            data_start = i + 1
            continue
        if header_line is None and ("mhjd" in line.lower() or "mjd" in line.lower()):
            header_line = line
            continue
        if header_line is not None and "====" not in line:
            data_start = i
            break
            
    if not header_line:
        print(f"Error: No header found in {filepath}")
        return None

    parts = header_line.split()
    col_map = {}
    image_labels = []
    
    telescope_col = None
    for j, col in enumerate(parts):
        if col.lower() == 'telescope':
            telescope_col = j
        if col.lower().startswith('mag_') and 'err' not in col.lower():
            lbl = col.split('_')[1]
            image_labels.append(lbl)
            col_map[lbl] = {'mag': j}
        elif col.lower().startswith('magerr_'):
            lbl = col.split('_')[1]
            if lbl in col_map:
                col_map[lbl]['err'] = j

    # Parse data
    raw_data = {lbl: {'t': [], 'mag': [], 'err': [], 'tel': []} for lbl in image_labels}
    
    for line in lines[data_start:]:
        line = line.strip()
        if not line or "====" in line: continue
        
        p = line.split()
        if len(p) < len(parts): continue # Skip incomplete lines
        
        try:
            t = float(p[0])
            tel = p[telescope_col] if telescope_col is not None else "Unknown"
            
            for lbl in image_labels:
                m_idx = col_map[lbl]['mag']
                e_idx = col_map[lbl].get('err')
                
                try:
                    m = float(p[m_idx])
                    e = float(p[e_idx]) if e_idx else 0.01
                    
                    if np.isfinite(m):
                        raw_data[lbl]['t'].append(t)
                        raw_data[lbl]['mag'].append(m)
                        raw_data[lbl]['err'].append(e)
                        raw_data[lbl]['tel'].append(tel)
                except:
                    pass
        except:
            continue

    # Build objects
    lcs = {}
    for lbl in image_labels:
        if len(raw_data[lbl]['t']) > 20:
            lcs[lbl] = LightCurve(
                label=lbl,
                t=np.array(raw_data[lbl]['t']),
                mag=np.array(raw_data[lbl]['mag']),
                magerr=np.array(raw_data[lbl]['err']),
                telescope=np.array(raw_data[lbl]['tel'])
            )
            
    return LensSystem(system_id, lcs)

# -----------------------------------------------------------------------------
# Analysis Functions (Simplified/Robust)
# -----------------------------------------------------------------------------

def full_gamma_analysis(system: LensSystem) -> Dict:
    """
    Actually compute Gamma using the multiscale approach.
    Re-implementing minimal filter logic.
    """
    lcs = system.light_curves
    labels = sorted(lcs.keys())
    
    gammas = []
    
    # Gaussian filter
    def bandpass(t, y, tau):
        # Handle gaps properly
        if len(t) < 5: return None, None
        
        dt = 1.0
        t_min, t_max = t.min(), t.max()
        tg = np.arange(t_min, t_max, dt)
        
        # Linear interp
        f_int = interpolate.interp1d(t, y, bounds_error=False, fill_value=np.nan)
        yg = f_int(tg)
        
        # Fill nans
        mask = np.isnan(yg)
        if np.all(mask): return None, None
        yg[mask] = np.interp(np.flatnonzero(mask), np.flatnonzero(~mask), yg[~mask])
        
        s1 = tau * 0.5 / dt
        s2 = tau * 1.5 / dt
        
        f1 = gaussian_filter1d(yg, s1)
        f2 = gaussian_filter1d(yg, s2)
        res = f1 - f2
        return tg, res

    for i, l1 in enumerate(labels):
        for l2 in labels[i+1:]:
            c1 = lcs[l1]
            c2 = lcs[l2]
            
            taus = [20, 40, 80, 160]
            delays = []
            valid_taus = []
            
            for tau in taus:
                res1 = bandpass(c1.t, c1.mag, tau)
                res2 = bandpass(c2.t, c2.mag, tau)
                if not res1[0] is not None or not res2[0] is not None: continue
                
                t1, y1 = res1
                t2, y2 = res2
                
                # Cross corr
                lags = np.arange(-150, 150, 1)
                corrs = []
                for l in lags:
                    shift = int(l)
                    if shift >= 0: 
                        y2s = np.roll(y2, shift)
                        y2s[:shift] = np.nan
                    else: 
                        y2s = np.roll(y2, shift)
                        y2s[shift:] = np.nan
                        
                    # Valid overlap
                    mask = np.isfinite(y1) & np.isfinite(y2s)
                    if np.sum(mask) < 20: 
                        corrs.append(np.nan)
                        continue
                        
                    r = np.corrcoef(y1[mask], y2s[mask])[0,1]
                    corrs.append(r)
                
                corrs = np.array(corrs)
                if np.any(np.isfinite(corrs)) and np.nanmax(corrs) > 0.2:
                    best_lag = lags[np.nanargmax(corrs)]
                    delays.append(best_lag)
                    valid_taus.append(tau)
            
            if len(valid_taus) >= 3:
                # Fit slope
                slope, intercept, r_val, p_val, std_err = stats.linregress(np.log10(valid_taus), delays)
                if np.isfinite(slope):
                    gammas.append(slope)
                
    if not gammas:
        return {'gamma_mean': np.nan, 'gamma_std': np.nan, 'n_pairs': 0}
        
    return {
        'gamma_mean': float(np.mean(gammas)),
        'gamma_std': float(np.std(gammas)),
        'n_pairs': len(gammas)
    }

# -----------------------------------------------------------------------------
# Main Execution
# -----------------------------------------------------------------------------

def main():
    print("Running Telescope Consistency Analysis...")
    
    results = {}
    
    # 1. RXJ1131 & DESJ0408 (Single file with telescope column)
    systems_to_check = {
        'RXJ1131': 'RXJ1131_Tewes2013.rdb',
        'DESJ0408': 'DESJ0408_Courbin2017.rdb', 
    }
    
    for name, fname in systems_to_check.items():
        fpath = DATA_DIR / fname
        if not fpath.exists(): continue
        
        print(f"\nProcessing {name}...")
        sys = parse_rdb_file(fpath)
        if sys:
            telescopes = sys.get_telescopes()
            print(f"  Found telescopes: {telescopes}")
            
            for tel in telescopes:
                if tel == 'Unknown': continue
                
                sub_sys = sys.filter_by_telescope(tel)
                if sub_sys:
                    res = full_gamma_analysis(sub_sys)
                    results[f"{name}_{tel}"] = res
                    print(f"  {tel}: Gamma = {res['gamma_mean']:.1f} +/- {res['gamma_std']:.1f} (N={res['n_pairs']})")

    # 2. WFI2033 (Separate files)
    print("\nProcessing WFI2033...")
    wfi_files = {
        'EulerCAM': DATA_DIR / 'WFI2033_ecam.dat',
        'SMARTS': DATA_DIR / 'WFI2033_smarts.dat'
    }
    
    for tel, fpath in wfi_files.items():
        if not fpath.exists(): continue
        
        # Custom parse for these dat files
        try:
            raw = np.loadtxt(fpath)
            # WFI2033 ecam/smarts format: MJD A errA B errB C errC
            lcs = {}
            if raw.shape[1] >= 7:
                lcs['A'] = LightCurve('A', raw[:,0], raw[:,1], raw[:,2])
                lcs['B'] = LightCurve('B', raw[:,0], raw[:,3], raw[:,4])
                lcs['C'] = LightCurve('C', raw[:,0], raw[:,5], raw[:,6])
                
                sys = LensSystem('WFI2033', lcs)
                res = full_gamma_analysis(sys)
                results[f"WFI2033_{tel}"] = res
                print(f"  {tel}: Gamma = {res['gamma_mean']:.1f} +/- {res['gamma_std']:.1f} (N={res['n_pairs']})")
        except Exception as e:
            print(f"  Failed to parse {tel}: {e}")

    # 3. Compare Results
    print("\nConsistency Summary:")
    print("====================")
    
    # WFI2033
    wfi_res = {k:v for k,v in results.items() if 'WFI2033' in k}
    if len(wfi_res) > 1:
        print("WFI2033:")
        vals = []
        for k, v in wfi_res.items():
            print(f"  {k}: {v['gamma_mean']:.1f}")
            if np.isfinite(v['gamma_mean']): vals.append(v['gamma_mean'])
            
        if len(vals) >= 2:
            diff = abs(vals[0] - vals[1])
            print(f"  Difference: {diff:.1f}")
            if diff < 100: # Heuristic threshold
                print("  -> CONSISTENT")
            else:
                print("  -> INCONSISTENT / NOISY")

    # Save Results
    out_path = RESULTS_DIR / "step_3_8_consistency_analysis.json"
    
    with open(out_path, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\nSaved results to {out_path}")

if __name__ == "__main__":
    main()
