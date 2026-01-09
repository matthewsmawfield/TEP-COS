#!/usr/bin/env python3
"""
Step 6.93: SDSS Test DW - Blue Straggler Fraction (Collisional/Binary Clock)

Hypothesis:
Blue Stragglers (BSS) are formed via binary mass transfer or collisions. They appear 
as main-sequence stars younger than the turn-off. Their formation and lifetime are 
rate-dependent. In deep potentials, if stellar evolution or dynamical encounter rates 
are time-dilated, the steady-state fraction of BSS relative to the underlying 
population (e.g., HB or RGB stars) should vary.

Prediction:
BSS Fraction (N_BSS / N_RGB) decreases in the Inner Galaxy (slower formation/faster evolution?).
(Or varies with potential depth).

Data:
- sppParams: TEFFADOP, LOGGADOP, FEHADOP
- SpecObjAll: glon, glat

Method:
1. Select metal-poor stars (FEHADOP < -1.0) to ensure old population (Halo/Thick Disk).
   This minimizes contamination from young thin disk stars which mimic BSS.
2. Define BSS Region:
   - High Teff (6000 < T < 8000 K)
   - High Logg (> 3.5)
   - Old population (FeH < -1.0) -> These shouldn't exist as single stars.
3. Define Reference Region (RGB):
   - Cool (4000 < T < 5000 K)
   - Low Logg (< 3.0)
4. Calculate Fraction f_BSS = N_BSS / N_RGB.
5. Bin by Galactic Longitude |l| to separate Inner vs Outer Galaxy.
   - Inner: |l| < 45 (Towards center)
   - Outer: |l| > 135 (Anti-center)
   - Intermediate: 45 < |l| < 135
   (Note: SDSS footprint is mostly high latitude, but spans wide longitude).
6. Compare fractions.
"""

import pandas as pd
import numpy as np
from scipy import stats
import os
import json
import matplotlib.pyplot as plt
import requests
import time

# Paths
DATA_DIR = os.path.join(os.path.dirname(__file__), '..', '..', 'data', 'sdss')
RESULTS_DIR = os.path.join(os.path.dirname(__file__), '..', '..', 'results', 'outputs')
FIGURES_DIR = os.path.join(os.path.dirname(__file__), '..', '..', 'results', 'figures')

SDSS_URL = "https://skyserver.sdss.org/dr18/SkyServerWS/SearchTools/SqlSearch"

def query_sdss(sql, max_retries=3):
    for attempt in range(max_retries):
        try:
            response = requests.get(
                SDSS_URL,
                params={"cmd": sql, "format": "json"},
                timeout=300
            )
            if response.status_code == 200:
                data = response.json()
                if data and len(data) > 0 and "Rows" in data[0]:
                    return pd.DataFrame(data[0]["Rows"])
            else:
                print(f"  HTTP {response.status_code}")
        except Exception as e:
            print(f"  Error: {e}")
        if attempt < max_retries - 1:
            time.sleep(2)
    return None

def download_data(limit=10000):
    print(f"Querying SDSS for Test DW (Limit: {limit})...")
    
    # Selecting old, metal-poor stars
    # We join sppParams with SpecObjAll to get coordinates
    
    sql = f"""
    SELECT TOP {limit}
        s.specObjID, 
        s.TEFFADOP as teff, 
        s.LOGGADOP as logg, 
        s.FEHADOP as feh,
        p.cx, p.cy, p.cz -- Cartesian coords on unit sphere (for SkyCoord if needed)
        -- Or just calculate l, b from ra, dec or use if available? 
        -- SpecObjAll usually doesn't have l, b directly columns? 
        -- It has ra, dec.
    FROM sppParams s
    JOIN SpecObjAll p ON s.specObjID = p.specObjID
    WHERE s.FEHADOP < -1.0 
      AND s.TEFFADOP BETWEEN 4000 AND 9000
      AND s.LOGGADOP > 0
    """
    
    # Getting RA/DEC
    sql = f"""
    SELECT TOP {limit}
        s.specObjID, 
        s.TEFFADOP as teff, 
        s.LOGGADOP as logg, 
        s.FEHADOP as feh,
        p.ra, p.dec
    FROM sppParams s
    JOIN SpecObjAll p ON s.specObjID = p.specObjID
    WHERE s.FEHADOP < -1.0 
      AND s.TEFFADOP BETWEEN 4000 AND 9000
      AND s.LOGGADOP > 0
    """
    return query_sdss(sql)

def analyze_bss_fraction(df):
    print("Analyzing Blue Straggler Fraction...")
    
    if df is None or len(df) < 100:
        print("  Insufficient data.")
        return None
        
    df = df.dropna().copy()
    
    # Calculate Galactic Coordinates
    from astropy.coordinates import SkyCoord
    from astropy import units as u
    
    coords = SkyCoord(ra=df['ra'].values*u.degree, dec=df['dec'].values*u.degree, frame='icrs')
    gal = coords.galactic
    df['l'] = gal.l.degree
    df['b'] = gal.b.degree
    
    # Shift l to be 0-360 or -180 to 180? Astropy gives 0-360.
    # Angle from center: min(|l|, 360-l)
    df['l_dist'] = np.minimum(df['l'], 360 - df['l'])
    
    # Define Regions
    # Inner: l_dist < 60
    # Outer: l_dist > 120
    # Mid: 60-120
    
    df['region'] = 'Mid'
    df.loc[df['l_dist'] < 60, 'region'] = 'Inner'
    df.loc[df['l_dist'] > 120, 'region'] = 'Outer'
    
    # Define Populations
    # BSS: Teff 6000-8000, Logg > 3.5
    # RGB: Teff 4000-5000, Logg < 3.0
    
    bss_mask = (df['teff'].between(6000, 8000)) & (df['logg'] > 3.5)
    rgb_mask = (df['teff'].between(4000, 5000)) & (df['logg'] < 3.0)
    
    df['type'] = 'Other'
    df.loc[bss_mask, 'type'] = 'BSS'
    df.loc[rgb_mask, 'type'] = 'RGB'
    
    print("  Population Counts:")
    print(df['type'].value_counts())
    
    # Calculate Fractions per region
    results = []
    regions = ['Inner', 'Mid', 'Outer']
    
    print("\n  BSS Fraction by Region:")
    for reg in regions:
        subset = df[df['region'] == reg]
        n_bss = len(subset[subset['type'] == 'BSS'])
        n_rgb = len(subset[subset['type'] == 'RGB'])
        
        if n_rgb > 0:
            frac = n_bss / n_rgb
            err = frac * np.sqrt(1/n_bss + 1/n_rgb) if n_bss > 0 else 0
            print(f"    {reg}: N_BSS={n_bss}, N_RGB={n_rgb}, f={frac:.4f} +/- {err:.4f}")
            results.append({'region': reg, 'f': frac, 'err': err, 'mean_l': subset['l_dist'].mean()})
        else:
            print(f"    {reg}: N_BSS={n_bss}, N_RGB={n_rgb}, f=NaN")
            
    # Visualize
    if len(results) == 3:
        res_df = pd.DataFrame(results)
        
        plt.figure(figsize=(8, 6))
        plt.errorbar(res_df['mean_l'], res_df['f'], yerr=res_df['err'], fmt='o-', capsize=5, color='blue')
        plt.xlabel('Angular Distance from Galactic Center (deg)')
        plt.ylabel('BSS Fraction (N_BSS / N_RGB)')
        plt.title('Test DW: Blue Straggler Fraction vs Longitude')
        plt.grid(True, alpha=0.3)
        plt.ylim(bottom=0)
        
        out_path = os.path.join(FIGURES_DIR, 'sdss_test_dw_bss_fraction.png')
        plt.savefig(out_path, dpi=150)
        print(f"Figure saved to {out_path}")
        
        # Determine trend
        slope, intercept, r_val, p_val, std_err = stats.linregress(res_df['mean_l'], res_df['f'])
        print(f"  Gradient (f vs l_dist): slope={slope:.5f}, r={r_val:.2f}, p={p_val:.3f}")
        
        return {
            'slope': slope,
            'p_value': p_val,
            'regions': results
        }
    else:
        return None

def main():
    cache_path = os.path.join(DATA_DIR, 'sdss_bss_fraction.csv')
    
    if os.path.exists(cache_path):
        print("Loading cached data...")
        df = pd.read_csv(cache_path)
    else:
        df = download_data()
        if df is not None:
            df.to_csv(cache_path, index=False)
        else:
            print("Download failed.")
            return

    results = analyze_bss_fraction(df)
    
    out_path = os.path.join(RESULTS_DIR, 'sdss_test_dw_results.json')
    if results:
        with open(out_path, 'w') as f:
            json.dump(results, f, indent=2)
            
        print("\nSUMMARY TEST DW:")
        print(f"Slope (Fraction vs Distance from GC): {results['slope']:.5f}")
        
        # Hypothesis: Decreases in Inner Galaxy -> Low f at Low l_dist
        # -> Positive slope (f increases with l_dist)
        
        if results['p_value'] < 0.1: # Relaxed p-value for 3 points
            if results['slope'] > 0:
                print("RESULT: SIGNAL (BSS Fraction lower in Inner Galaxy)")
            else:
                print("RESULT: CONTRADICTED (BSS Fraction higher in Inner Galaxy)")
        else:
            print("RESULT: NULL (No significant gradient)")

if __name__ == "__main__":
    main()
