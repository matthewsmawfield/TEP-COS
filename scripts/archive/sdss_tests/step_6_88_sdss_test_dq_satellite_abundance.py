#!/usr/bin/env python3
"""
Step 6.88: SDSS Test DQ - Satellite Abundance (Halo Mass Proxy)

Hypothesis:
The number of satellite galaxies (N_sat) orbiting a central galaxy scales with the 
Dark Matter halo mass. In TEP, the "halo" is a phantom mass effect. The scaling of 
N_sat with Stellar Mass (baryons) might differ if the phantom mass ratio 
M_phan/M_bar depends on potential depth/concentration differently than LCDM 
abundance matching predicts.

Prediction:
Satellite Abundance N_sat (at fixed M_star) varies with Compactness/Sigma.

Data:
- stellarMassFSPSGranWideDust: logMass
- PhotoObjAll: petroR50_r (Compactness proxy)
- Neighbors: Count satellites
- SpecPhotoAll: z (for isolation)

Method:
1. Select Central Galaxies (spectroscopic sample).
2. Count neighbors within projected radius (e.g. 5 arcmin ~ 250 kpc at z=0.05).
   - Use Photometric neighbors (fainter).
   - Neighbors table stores pairs within 0.5 arcmin? No, usually small search.
   - Standard Neighbors table limit is often 0.5 arcmin. This is too small for satellites.
   - We might need to assume Neighbors table is not sufficient for N_sat > 100kpc.
   - However, query plan suggests using Neighbors table. Let's check max distance in Neighbors.
   - If max distance is small, we can only test "Close Companions" not full satellite halo.
   - Alternative: Use `Neighbors` table if it has larger search, or accept "Close Companion Fraction".
3. Calculate N_sat per central.
4. Bin by Compactness (Size at fixed Mass).
5. Compare Mean N_sat.
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

def download_data(limit=1000):
    print(f"Querying SDSS for Test DQ (Limit: {limit})...")
    
    # We want to count neighbors. Doing this in SQL with GROUP BY is efficient.
    # Neighbors table usually contains all pairs < 0.5 arcmin.
    # 0.5 arcmin at z=0.05 is ~30 kpc. This is "Close Companions".
    # This is still a valid proxy for "Environment" or "Halo Substructure".
    
    sql = f"""
    SELECT TOP {limit} 
        s.specObjID, 
        s.logMass,
        ph.petroR50_r as radius_r,
        ph.modelMag_r,
        count(n.NeighborObjID) as n_neighbors
    FROM stellarMassFSPSGranWideDust s
    JOIN PhotoObjAll ph ON s.specObjID = ph.objID
    JOIN SpecPhotoAll sp ON s.specObjID = sp.specObjID
    LEFT JOIN Neighbors n ON s.specObjID = n.objID -- Note: specObjID usually links to bestObjID?
    -- Actually Neighbors links objID. Need to ensure link.
    -- Join s -> SpecPhoto -> PhotoObj (objID) -> Neighbors
    WHERE sp.z BETWEEN 0.02 AND 0.1
      AND s.logMass > 10
      AND ph.petroR50_r > 0
    GROUP BY s.specObjID, s.logMass, ph.petroR50_r, ph.modelMag_r
    """
    
    # Correct join logic for Neighbors:
    # stellarMass -> SpecObjAll -> PhotoObjAll (via bestObjID) -> Neighbors (via objID)
    
    sql = f"""
    SELECT TOP {limit}
        s.specObjID,
        s.logMass,
        ph.petroR50_r as radius_r,
        count(n.NeighborObjID) as n_neighbors
    FROM stellarMassFSPSGranWideDust s
    JOIN SpecObjAll so ON s.specObjID = so.specObjID
    JOIN PhotoObjAll ph ON so.bestObjID = ph.objID
    LEFT JOIN Neighbors n ON ph.objID = n.objID
    WHERE so.z BETWEEN 0.02 AND 0.1
      AND s.logMass > 10.0
      AND ph.petroR50_r > 0
    GROUP BY s.specObjID, s.logMass, ph.petroR50_r
    """
    
    return query_sdss(sql)

def analyze_satellite_abundance(df):
    print("Analyzing Satellite Abundance...")
    
    if df is None or len(df) < 100:
        print("  Insufficient data.")
        return None
        
    df = df.dropna().copy()
    
    # Control for Mass
    # N_sat increases with Mass.
    # Radius increases with Mass.
    # We want to see if N_sat depends on Radius at fixed Mass.
    
    df['log_r'] = np.log10(df['radius_r'])
    
    from sklearn.linear_model import LinearRegression
    X = df[['logMass']]
    y = df['n_neighbors']
    
    reg = LinearRegression().fit(X, y)
    print(f"  Control Fit (Mass -> N_sat) R2: {reg.score(X, y):.3f}")
    
    df['n_resid'] = y - reg.predict(X)
    
    # Correlate with Radius (Compactness)
    # Compactness ~ 1/Radius.
    # Prediction: Compact (Small R) -> Deep Potential -> More Satellites?
    # Or TEP: Phantom Mass halo is different.
    
    slope, intercept, r_val, p_val, std_err = stats.linregress(df['log_r'], df['n_resid'])
    
    print(f"  Correlation (log Radius vs N_sat Resid): r={r_val:.3f}, p={p_val:.2e}, slope={slope:.4f}")
    
    # Plot
    fig, ax = plt.subplots(1, 2, figsize=(12, 5))
    
    # Raw N vs Mass (color=Radius)
    sc = ax[0].scatter(df['logMass'], df['n_neighbors'], c=df['log_r'], cmap='viridis_r', s=10, alpha=0.6)
    plt.colorbar(sc, ax=ax[0], label='log(Radius) (inverted)')
    ax[0].set_xlabel('log Stellar Mass')
    ax[0].set_ylabel('Number of Close Neighbors')
    ax[0].set_title('Neighbors vs Mass')
    
    # Residuals vs Radius
    ax[1].scatter(df['log_r'], df['n_resid'], alpha=0.5, s=10, c='indigo')
    
    x_range = np.linspace(df['log_r'].min(), df['log_r'].max(), 100)
    ax[1].plot(x_range, intercept + slope*x_range, 'k--', label=f'Slope={slope:.3f}')
    
    ax[1].set_xlabel('log(Petro Radius R50)')
    ax[1].set_ylabel('Neighbor Count Residual')
    ax[1].set_title('Test DQ: Satellites vs Compactness')
    ax[1].legend()
    ax[1].grid(True, alpha=0.3)
    
    out_path = os.path.join(FIGURES_DIR, 'sdss_test_dq_satellites.png')
    plt.savefig(out_path, dpi=150)
    print(f"Figure saved to {out_path}")
    
    return {
        'correlation_r': r_val,
        'p_value': p_val,
        'slope': slope,
        'n_gal': int(len(df))
    }

def main():
    cache_path = os.path.join(DATA_DIR, 'sdss_satellite_abundance.csv')
    
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

    results = analyze_satellite_abundance(df)
    
    out_path = os.path.join(RESULTS_DIR, 'sdss_test_dq_results.json')
    if results:
        with open(out_path, 'w') as f:
            json.dump(results, f, indent=2)
            
        print("\nSUMMARY TEST DQ:")
        print(f"Slope (Radius vs N_sat Resid): {results['slope']:.4f}")
        
        # Interpretation:
        # Negative slope: Larger Radius -> Fewer Satellites (So Compact -> More Satellites).
        # Positive slope: Larger Radius -> More Satellites.
        
        if results['p_value'] < 0.05 and abs(results['slope']) > 0.1:
             print("RESULT: SIGNAL (Significant dependence)")
        else:
             print("RESULT: NULL (No significant dependence)")

if __name__ == "__main__":
    main()
