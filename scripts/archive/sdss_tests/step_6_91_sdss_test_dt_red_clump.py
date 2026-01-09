#!/usr/bin/env python3
"""
Step 6.91: SDSS Test DT - Red Clump Absolute Magnitude

Hypothesis:
The Red Clump (RC) magnitude is a standard candle fixed by He-core burning physics. 
If nuclear rates are time-dilated in the Inner Galaxy, the intrinsic luminosity L 
might differ. Comparing the geometric distance modulus (Gaia) to the apparent 
magnitude should reveal if M_RC varies with potential depth (Galactocentric Radius).

Prediction:
Derived Absolute Magnitude M_G of RC stars fades (or varies) in the Inner Galaxy 
relative to the Outer Galaxy.

Data:
- mos_gaia_dr2_source: phot_g_mean_mag, bp_rp, l, b
- mos_geometric_distances_gaia_dr2: r_est (Distance in pc)

Method:
1. Select RC stars using color (1.1 < BP-RP < 1.3) and absolute magnitude cuts.
2. Calculate Absolute Magnitude M_G = m_G - 5*log10(r) + 5.
   (Ignoring extinction for now? Or select low extinction regions?)
   Using abs(b) > 10 helps, but we want Inner Galaxy which is often high extinction.
   However, we look for a gradient with R_GC.
3. Calculate Galactocentric Radius R_GC.
   R_GC = sqrt(R0^2 + (d cos b)^2 - 2 R0 d cos b cos l)
   Assume R0 = 8000 pc.
4. Correlate M_G with R_GC.
   Control for extinction? A_G varies. 
   If we see a trend, is it extinction or TEP?
   Extinction makes stars fainter -> Derived M_G = m - DM. m is larger (fainter).
   So M_G appears fainter (larger positive) if extinction is underestimated.
   We need to be careful. Ideally use extinction-corrected magnitudes if available.
   Or select low extinction fields.
   
   Let's try to look at high latitudes |b| > 20 where extinction is lower, 
   but varying R_GC (towards vs away from center).
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
    print(f"Querying Gaia for Test DT (Limit: {limit})...")
    
    # Simple RC selection: 1.1 < bp_rp < 1.3
    # Use |b| > 15 to minimize extinction issues
    
    sql = f"""
    SELECT TOP {limit}
        g.source_id,
        g.phot_g_mean_mag as g_mag,
        g.bp_rp,
        g.l, g.b,
        gd.r_est as dist
    FROM mos_gaia_dr2_source g
    JOIN mos_geometric_distances_gaia_dr2 gd ON g.source_id = gd.source_id
    WHERE g.bp_rp BETWEEN 1.1 AND 1.3
      AND abs(g.b) > 15
      AND gd.r_est > 0 AND gd.r_est < 5000 -- Limit distance to reasonably reliable
    """
    return query_sdss(sql)

def analyze_red_clump(df):
    print("Analyzing Red Clump Magnitude...")
    
    if df is None or len(df) < 100:
        print("  Insufficient data.")
        return None
        
    df = df.dropna().copy()
    
    # Calculate Absolute Magnitude (Raw, no extinction corr)
    # M = m - 5 log(d) + 5
    df['M_g'] = df['g_mag'] - 5 * np.log10(df['dist']) + 5
    
    # Calculate Galactocentric Radius
    # R0 = 8000 pc
    R0 = 8000.0
    l_rad = np.radians(df['l'])
    b_rad = np.radians(df['b'])
    d = df['dist']
    
    # Law of cosines
    # Proj dist on plane = d cos b
    d_proj = d * np.cos(b_rad)
    
    # R_GC^2 = R0^2 + d_proj^2 - 2 R0 d_proj cos l
    df['R_gc'] = np.sqrt(R0**2 + d_proj**2 - 2 * R0 * d_proj * np.cos(l_rad))
    
    # Filter for rough RC magnitude range (e.g. M_g ~ 0.5)
    # RC in G is roughly +0.5.
    # Allow some spread.
    rc_candidates = df[(df['M_g'] > -1.0) & (df['M_g'] < 2.0)].copy()
    
    print(f"  RC Candidates: {len(rc_candidates)} (from {len(df)})")
    
    if len(rc_candidates) < 50:
        print("  Not enough RC candidates.")
        return None
        
    # Correlate M_g with R_gc
    # Prediction: In Inner Galaxy (Low R_gc), M_g is fainter (higher value) or brighter?
    # TEP: "Fades in Inner Galaxy" -> M_g increases (fainter) at low R_gc.
    # Correlation (R_gc, M_g): Negative slope (Low R -> High M).
    
    slope, intercept, r_val, p_val, std_err = stats.linregress(rc_candidates['R_gc'], rc_candidates['M_g'])
    
    print(f"  Correlation (R_gc vs M_g): r={r_val:.3f}, p={p_val:.2e}, slope={slope:.5f}")
    
    # Plot
    plt.figure(figsize=(8, 6))
    
    # Scatter
    plt.scatter(rc_candidates['R_gc'], rc_candidates['M_g'], s=1, alpha=0.3, c='red')
    
    # Binned trend
    rc_candidates['bin'] = pd.cut(rc_candidates['R_gc'], bins=10)
    grouped = rc_candidates.groupby('bin').agg({'R_gc': 'mean', 'M_g': 'mean'}).reset_index()
    plt.plot(grouped['R_gc'], grouped['M_g'], 'bo-', label='Binned Mean')
    
    # Fit
    x_range = np.linspace(rc_candidates['R_gc'].min(), rc_candidates['R_gc'].max(), 100)
    plt.plot(x_range, intercept + slope*x_range, 'k--', label=f'Slope={slope*1000:.2f} mag/kpc')
    
    plt.xlabel('Galactocentric Radius (pc)')
    plt.ylabel('Absolute Magnitude M_G')
    plt.title('Test DT: Red Clump Magnitude vs Potential')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.ylim(2.0, -1.0) # Inverted mag axis
    
    out_path = os.path.join(FIGURES_DIR, 'sdss_test_dt_red_clump.png')
    plt.savefig(out_path, dpi=150)
    print(f"Figure saved to {out_path}")
    
    return {
        'correlation_r': r_val,
        'p_value': p_val,
        'slope_mag_pc': slope,
        'slope_mag_kpc': slope * 1000,
        'n_stars': int(len(rc_candidates))
    }

def main():
    cache_path = os.path.join(DATA_DIR, 'sdss_red_clump.csv')
    
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

    results = analyze_red_clump(df)
    
    out_path = os.path.join(RESULTS_DIR, 'sdss_test_dt_results.json')
    if results:
        with open(out_path, 'w') as f:
            json.dump(results, f, indent=2)
            
        print("\nSUMMARY TEST DT:")
        print(f"Slope (mag/kpc): {results['slope_mag_kpc']:.4f}")
        
        # Prediction: Fades in Inner -> High M at Low R -> Negative Slope of M vs R
        # Wait: High M (fainter) at Low R. 
        # Low R -> High M. High R -> Low M.
        # Slope dM/dR should be Negative.
        
        if results['p_value'] < 0.05 and results['slope_mag_kpc'] < -0.01:
             print("RESULT: SIGNAL (RC fades in Inner Galaxy)")
        elif results['p_value'] < 0.05 and results['slope_mag_kpc'] > 0.01:
             print("RESULT: CONTRADICTED (RC brightens in Inner Galaxy)")
        else:
             print("RESULT: NULL (No significant gradient)")

if __name__ == "__main__":
    main()
