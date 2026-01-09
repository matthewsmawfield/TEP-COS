#!/usr/bin/env python3
"""
Step 3.9: Q2237 Temporal Stability Analysis (OGLE-II vs OGLE-III)

Tests the stability of the Temporal Shear signal over decadal timescales.
TEP Prediction: Γ is a structural property of the potential -> Constant over time.
Microlensing Prediction: Caused by caustic crossings/stars -> Transient (varies between epochs).

Data:
- OGLE-II (1997-2000): V-band
- OGLE-III (2001-2009): V-band

Methodology:
1. Load and align OGLE-II and OGLE-III datasets.
2. Compute Γ independently for both epochs.
3. Check for consistency.
"""

import numpy as np
from pathlib import Path
import json
from scipy import stats, interpolate
from dataclasses import dataclass
from typing import Dict, List, Tuple
from scipy.ndimage import gaussian_filter1d

# Reuse core classes from previous steps
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

def load_ogle2():
    fpath = DATA_DIR / "Q2237_ogle2_phot.dat"
    if not fpath.exists():
        print(f"Missing {fpath}")
        return None
        
    try:
        # Format: DateString HJD Airmass V_A err_A V_B err_B V_C err_C V_D err_D
        # Use pandas to handle the date string column
        import pandas as pd
        df = pd.read_csv(fpath, delim_whitespace=True, header=None, comment='#')
        
        # Col 1 is HJD-2450000
        t = df.iloc[:, 1].values + 2450000
        
        lcs = {}
        # A: 3,4; B: 5,6; C: 7,8; D: 9,10
        lcs['A'] = LightCurve('A', t, df.iloc[:, 3].values, df.iloc[:, 4].values)
        lcs['B'] = LightCurve('B', t, df.iloc[:, 5].values, df.iloc[:, 6].values)
        lcs['C'] = LightCurve('C', t, df.iloc[:, 7].values, df.iloc[:, 8].values)
        lcs['D'] = LightCurve('D', t, df.iloc[:, 9].values, df.iloc[:, 10].values)
        
        return LensSystem("Q2237_OGLE2", lcs)
    except Exception as e:
        print(f"Error loading OGLE2: {e}")
        return None

def load_ogle3():
    # Files: Q2237_A.dat, etc.
    # Format: HJD-2450000, V, err
    lcs = {}
    for img in ['A', 'B', 'C', 'D']:
        fpath = DATA_DIR / f"Q2237_{img}.dat"
        if not fpath.exists():
            print(f"Missing {fpath}")
            continue
            
        try:
            # Use pandas for robustness against bad lines
            # Force reading only the first 3 columns to avoid footer issues
            import pandas as pd
            df = pd.read_csv(fpath, delim_whitespace=True, header=None, comment='#', 
                           on_bad_lines='skip', usecols=[0, 1, 2])
            
            # Ensure we have data
            if df.empty:
                print(f"Empty dataframe for {img}")
                continue

            t = df.iloc[:, 0].values + 2450000
            mag = df.iloc[:, 1].values
            err = df.iloc[:, 2].values
            
            lcs[img] = LightCurve(img, t, mag, err)
        except Exception as e:
            print(f"Error loading OGLE3 {img}: {e}")
            
    if len(lcs) < 2: return None
    return LensSystem("Q2237_OGLE3", lcs)

def load_glendama():
    # File: glendama_J_ApA_616_A118_table16.csv
    # Columns: MJD, mA, e_mA, mB, e_mB, mC, e_mC, mD, e_mD, mS, e_mS
    fpath = DATA_DIR / "glendama_J_ApA_616_A118_table16.csv"
    if not fpath.exists():
        print(f"Missing {fpath}")
        return None
        
    try:
        import pandas as pd
        df = pd.read_csv(fpath)
        
        # MJD in GLENDAMA usually MJD (JD-2400000.5)
        # But sample printed earlier showed ~3900. 
        # If 3900 is MJD, that's year 1869. 
        # If it's JD-2450000, 3900 is year 2006. This matches "1999-2016" overlap with OGLE-III.
        # Let's assume JD-2450000 for consistency with OGLE.
        t = df['MJD'].values + 2450000
        
        lcs = {}
        # Columns mA, e_mA, etc.
        for img in ['A', 'B', 'C', 'D']:
            col_m = f"m{img}"
            col_e = f"e_m{img}"
            if col_m in df.columns:
                # Filter NaNs
                mask = ~df[col_m].isna()
                lcs[img] = LightCurve(img, t[mask], df.loc[mask, col_m].values, df.loc[mask, col_e].values)
                
        return LensSystem("Q2237_GLENDAMA", lcs)
    except Exception as e:
        print(f"Error loading GLENDAMA: {e}")
        return None

def analyze_gamma(system: LensSystem):
    # Simplified Gamma estimation (slope of delay vs log tau)
    # Using 3 timescales: 20, 40, 80 days
    
    taus = [20, 40, 80]
    lcs = system.light_curves
    labels = sorted(lcs.keys())
    
    gammas = []
    
    def get_filtered(lc, tau):
        if len(lc.t) < 5: return None
        
        # Fill gaps
        dt = 1.0
        tg = np.arange(lc.t.min(), lc.t.max(), dt)
        if len(tg) < 50: return None
        
        # Linear interp
        y = interpolate.interp1d(lc.t, lc.mag, bounds_error=False, fill_value=np.nan)(tg)
        mask = np.isnan(y)
        
        # If too many gaps, skip
        if np.sum(mask) / len(mask) > 0.5: return None
        
        if np.all(mask): return None
        # Simple gap filling
        y[mask] = np.interp(np.flatnonzero(mask), np.flatnonzero(~mask), y[~mask])
        
        s = tau / 2.355 / dt
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
                
                # Check overlap length
                n_ov = min(len(y1), len(y2))
                if n_ov < 2 * tau: continue
                
                # Cross corr
                lags = np.arange(-20, 20, 0.5)
                corrs = []
                for l in lags:
                    shift = int(l)
                    # Safe roll/shift
                    if shift > 0:
                        # y2 shifted right means y2(t-shift)
                        # correlation: sum y1(t) * y2(t-tau)
                        # This is tricky with numpy.
                        # Let's just use valid part
                        valid_len = n_ov - abs(shift)
                        v1 = y1[shift:]
                        v2 = y2[:-shift]
                    elif shift < 0:
                        v1 = y1[:shift]
                        v2 = y2[-shift:]
                    else:
                        v1 = y1
                        v2 = y2
                    
                    if len(v1) < 10: 
                        corrs.append(0)
                        continue
                        
                    r = np.corrcoef(v1, v2)[0,1]
                    corrs.append(r)
                
                if np.max(corrs) > 0.3: # Lower threshold
                    best_lag = lags[np.argmax(corrs)]
                    delays.append(best_lag)
                    valid_taus.append(tau)
            
            if len(valid_taus) >= 2:
                # If we have at least 2 points, we can fit a slope
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
    print("Running Q2237 Temporal Stability Analysis...")
    
    ogle2 = load_ogle2()
    ogle3 = load_ogle3()
    glendama = load_glendama()
    
    res2 = analyze_gamma(ogle2) if ogle2 else None
    res3 = analyze_gamma(ogle3) if ogle3 else None
    resG = analyze_gamma(glendama) if glendama else None
    
    print("\nResults:")
    if res2:
        print(f"OGLE-II (1997-2000): Gamma = {res2['gamma_mean']:.1f} +/- {res2['gamma_std']:.1f} (N={res2['n_pairs']})")
    if res3:
        print(f"OGLE-III (2001-2009): Gamma = {res3['gamma_mean']:.1f} +/- {res3['gamma_std']:.1f} (N={res3['n_pairs']})")
    if resG:
        print(f"GLENDAMA (2006-2016): Gamma = {resG['gamma_mean']:.1f} +/- {resG['gamma_std']:.1f} (N={resG['n_pairs']})")
        
    # Interpretation
    results = {}
    if res2: results["ogle2"] = res2
    if res3: results["ogle3"] = res3
    if resG: results["glendama"] = resG
    
    if len(results) >= 2:
        # Calculate max difference
        means = [r['gamma_mean'] for r in results.values()]
        diff = max(means) - min(means)
        print(f"Max Difference: {diff:.1f}")
        results["max_diff"] = diff
        
    with open(RESULTS_DIR / "step_3_9_q2237_stability.json", 'w') as f:
        # Convert numpy types
        def convert(o):
            if isinstance(o, np.generic): return o.item()
            raise TypeError
        json.dump(results, f, indent=2, default=convert)

if __name__ == "__main__":
    main()
