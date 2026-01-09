#!/usr/bin/env python3
"""
Step 6.71: SDSS Test CT - Schechter Function High-Mass Cutoff

Hypothesis:
The exponential cutoff in the galaxy mass function is driven by feedback (AGN/Supernovae) 
preventing star formation in massive halos. Feedback is a rate-limited process. 
In the deepest potentials (massive halos, high sigma), these rates are time-dilated. 
Suppression might be less efficient per unit cosmic time, leading to a "softer" cutoff 
(more ultra-massive galaxies) than predicted by standard feedback models tuned to the field.

Prediction:
The high-mass end of the Schechter function is shallower (more massive galaxies) in high-sigma environments.

Data:
- stellarMassFSPSGranWideDust: logMass
- ebossMCPM: mid_dens_1 (Density) - Optional, we focus on Sigma as potential proxy.
- emissionLinesPort: sigma_stars

Method:
1. Select massive galaxies (logMass > 10).
2. Split into High-Sigma (Top 25%) and Low-Sigma (Bottom 25%) bins.
3. Compute the Galaxy Stellar Mass Function (GSMF) for both (normalized).
4. Fit Schechter-like cutoff (or just compare the tail).
   Tail slope: d(log N)/d(log M) at high M.
   Or compare fraction of Ultra-Massive Galaxies (logM > 11.5).
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

def download_data(limit=100):
    print(f"Querying SDSS for Test CT (Limit: {limit})...")
    
    # Use emissionLinesPort for sigma, stellarMass for Mass
    # Need ebossMCPM for density? Maybe not strictly required if Sigma is the TEP variable.
    # But checking if we can join ebossMCPM is good. Previous checks said yes.
    # Let's keep it simple first: Mass + Sigma.
    
    sql = f"""
    SELECT TOP {limit}
        s.specObjID,
        s.logMass,
        p.sigma_stars
    FROM stellarMassFSPSGranWideDust s
    JOIN emissionLinesPort p ON s.specObjID = p.specObjID
    WHERE 
        s.logMass > 10.0
        AND p.sigma_stars > 50 AND p.sigma_stars < 450
    """
    return query_sdss(sql)

def analyze_schechter(df):
    print("Analyzing Schechter Cutoff...")
    
    # Clean
    df = df.dropna().copy()
    
    print(f"  Sample size: {len(df)}")
    
    # Define High/Low Sigma
    q25 = df['sigma_stars'].quantile(0.25)
    q75 = df['sigma_stars'].quantile(0.75)
    
    print(f"  Low Sigma < {q25:.1f} km/s")
    print(f"  High Sigma > {q75:.1f} km/s")
    
    low_sig = df[df['sigma_stars'] <= q25]
    high_sig = df[df['sigma_stars'] >= q75]
    
    # Compute GSMF (Histogram)
    bins = np.linspace(10.0, 12.5, 26) # 0.1 dex bins
    
    hist_low, edges = np.histogram(low_sig['logMass'], bins=bins, density=True)
    hist_high, _ = np.histogram(high_sig['logMass'], bins=bins, density=True)
    centers = (edges[:-1] + edges[1:]) / 2
    
    # Calculate High-Mass Slope (between 11.0 and 12.0)
    # Log counts
    mask_fit = (centers > 11.0) & (centers < 12.0)
    
    def get_slope(hist, name):
        # Filter zero bins for log
        valid = (hist > 0) & mask_fit
        if valid.sum() < 3: return 0
        
        y = np.log10(hist[valid])
        x = centers[valid]
        slope, intercept, r, p, _ = stats.linregress(x, y)
        print(f"  {name} Slope (log N vs log M): {slope:.3f}")
        return slope
        
    slope_low = get_slope(hist_low, "Low Sigma")
    slope_high = get_slope(hist_high, "High Sigma")
    
    delta_slope = slope_high - slope_low
    print(f"  Delta Slope (High - Low): {delta_slope:.3f}")
    
    # Interpretation:
    # Schechter exp cutoff: N ~ exp(-M/M*). log N ~ -M.
    # So slope is negative.
    # "Softer" cutoff means less negative (flatter, closer to 0).
    # So if High Sigma has softer cutoff, Slope_High > Slope_Low (e.g., -2 vs -4).
    # Positive Delta Slope => Signal.
    
    # Plot
    fig, ax = plt.subplots(figsize=(8, 6))
    
    ax.plot(centers, hist_low, 'b-o', label=f'Low Sigma (<{q25:.0f})')
    ax.plot(centers, hist_high, 'r-s', label=f'High Sigma (>{q75:.0f})')
    
    ax.set_yscale('log')
    ax.set_xlabel('log Stellar Mass [M_sun]')
    ax.set_ylabel('Normalized PDF')
    ax.set_title(f'Test CT: GSMF Cutoff (Delta Slope={delta_slope:.2f})')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    out_path = os.path.join(FIGURES_DIR, 'sdss_test_ct_schechter.png')
    plt.savefig(out_path, dpi=150)
    print(f"Figure saved to {out_path}")
    
    return {
        'slope_low': slope_low,
        'slope_high': slope_high,
        'delta_slope': delta_slope,
        'n_sample': int(len(df))
    }

def main():
    cache_path = os.path.join(DATA_DIR, 'sdss_schechter.csv')
    
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

    results = analyze_schechter(df)
    
    out_path = os.path.join(RESULTS_DIR, 'sdss_test_ct_results.json')
    if results:
        with open(out_path, 'w') as f:
            json.dump(results, f, indent=2)
            
        print("\nSUMMARY TEST CT:")
        print("Prediction: High-mass cutoff is shallower (softer) in high-sigma environments.")
        print(f"Observed Delta Slope (High - Low): {results['delta_slope']:.3f}")
        
        if results['delta_slope'] > 0.5:
             print("RESULT: SIGNAL (Cutoff significantly softer in deep potentials)")
        elif results['delta_slope'] < -0.5:
             print("RESULT: CONTRADICTED (Cutoff sharper in deep potentials)")
        else:
             print("RESULT: NULL (Similar cutoff shape)")

if __name__ == "__main__":
    main()
