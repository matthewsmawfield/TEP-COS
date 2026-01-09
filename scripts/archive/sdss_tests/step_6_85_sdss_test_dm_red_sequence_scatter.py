#!/usr/bin/env python3
"""
Step 6.85: SDSS Test DM - Red Sequence Scatter (The Variance of Time)

Hypothesis:
The Red Sequence has a finite width (scatter) in color at fixed mass. Part of this 
scatter is due to age variations. In TEP, variations in potential depth (sigma) 
at fixed mass (due to concentration/profile differences) lead to variations in 
proper time elapsed since formation. This "Time Variance" should contribute to 
the age scatter and thus the color scatter.

Prediction:
Color Residuals (relative to RS mean) correlate with Sigma Residuals (relative to M-Sigma relation).
At fixed Mass, higher Sigma -> Deeper Potential -> Less Proper Time -> Younger/Bluer? 
Or Slower Evolution -> Redder?
TEP Prediction: Time dilation slows down evolution. High Sigma -> "Younger" stellar population 
(less evolved) at fixed cosmic time. BUT, if they formed early, they might be old.
Let's stick to the prediction: Color correlates with Sigma at fixed Mass.

Data:
- stellarMassFSPSGranWideDust: logMass
- SpecPhotoAll: modelMag_u, modelMag_r (Color u-r)
- emissionLinesPort: sigma_stars
- PhotoObjAll: petroR50_r (Size, compactness)

Method:
1. Select Red Sequence Galaxies (Quiescent).
   - u-r vs Mass selection or sSFR selection.
2. Fit the Mean Red Sequence: (u-r) = a * logM + b.
3. Calculate Color Residuals: Delta(u-r).
4. Fit the Mass-Sigma Relation: log(sigma) = c * logM + d.
5. Calculate Sigma Residuals: Delta(log sigma).
6. Correlate Delta(u-r) vs Delta(log sigma).
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
    print(f"Querying SDSS for Test DM (Limit: {limit})...")
    
    sql = f"""
    SELECT TOP {limit}
        s.specObjID, 
        s.logMass, 
        sp.modelMag_u - sp.modelMag_r as color_ur,
        e.sigma_stars
    FROM stellarMassFSPSGranWideDust s
    JOIN SpecPhotoAll sp ON s.specObjID = sp.specObjID
    JOIN emissionLinesPort e ON s.specObjID = e.specObjID
    WHERE sp.class = 'GALAXY' 
      AND s.logMass > 10.0
      AND e.sigma_stars > 50
      AND abs(sp.modelMag_u - sp.modelMag_r) < 5 -- Sanity check
    """
    return query_sdss(sql)

def analyze_rs_scatter(df):
    print("Analyzing Red Sequence Scatter...")
    
    if df is None or len(df) < 100:
        print("  Insufficient data.")
        return None
        
    df = df.dropna().copy()
    
    # 1. Define Red Sequence
    # Simple cut: u-r > 2.0 (roughly)
    # Better: Fit the ridge line.
    # Let's visualize first.
    
    # Isolate RS roughly
    rs_candidates = df[df['color_ur'] > 2.0].copy()
    
    if len(rs_candidates) < 50:
        print("  Not enough red galaxies.")
        return None
        
    # Fit RS: Color = f(Mass)
    from sklearn.linear_model import LinearRegression, RANSACRegressor
    
    X_mass = rs_candidates[['logMass']]
    y_color = rs_candidates['color_ur']
    
    # Use RANSAC to be robust against Green Valley / outliers
    ransac = RANSACRegressor(LinearRegression(), residual_threshold=0.2)
    ransac.fit(X_mass, y_color)
    
    rs_candidates['color_pred'] = ransac.predict(X_mass)
    rs_candidates['color_resid'] = rs_candidates['color_ur'] - rs_candidates['color_pred']
    
    print(f"  RS Fit: slope={ransac.estimator_.coef_[0]:.3f}, intercept={ransac.estimator_.intercept_:.3f}")
    
    # 2. Fit Mass-Sigma Relation
    # log(sigma) = f(logMass)
    rs_candidates['log_sigma'] = np.log10(rs_candidates['sigma_stars'])
    
    X_mass = rs_candidates[['logMass']]
    y_sigma = rs_candidates['log_sigma']
    
    reg_sigma = LinearRegression().fit(X_mass, y_sigma)
    rs_candidates['log_sigma_pred'] = reg_sigma.predict(X_mass)
    rs_candidates['sigma_resid'] = rs_candidates['log_sigma'] - rs_candidates['log_sigma_pred']
    
    print(f"  M-Sigma Fit: slope={reg_sigma.coef_[0]:.3f}")
    
    # 3. Correlate Residuals
    # Delta(Color) vs Delta(log Sigma)
    slope, intercept, r_val, p_val, std_err = stats.linregress(rs_candidates['sigma_resid'], rs_candidates['color_resid'])
    
    print(f"  Correlation (Sigma Resid vs Color Resid): r={r_val:.3f}, p={p_val:.2e}, slope={slope:.4f}")
    
    # Plot
    fig, ax = plt.subplots(1, 3, figsize=(18, 5))
    
    # RS
    ax[0].scatter(df['logMass'], df['color_ur'], s=1, c='gray', alpha=0.3)
    ax[0].scatter(rs_candidates['logMass'], rs_candidates['color_ur'], s=5, c='crimson', alpha=0.5)
    ax[0].plot(rs_candidates['logMass'], rs_candidates['color_pred'], 'k--')
    ax[0].set_xlabel('log Stellar Mass')
    ax[0].set_ylabel('u-r Color')
    ax[0].set_title('Red Sequence')
    
    # M-Sigma
    ax[1].scatter(rs_candidates['logMass'], rs_candidates['log_sigma'], s=5, c='teal', alpha=0.5)
    ax[1].plot(rs_candidates['logMass'], rs_candidates['log_sigma_pred'], 'k--')
    ax[1].set_xlabel('log Stellar Mass')
    ax[1].set_ylabel('log Sigma')
    ax[1].set_title('Mass-Sigma Relation')
    
    # Residuals
    ax[2].scatter(rs_candidates['sigma_resid'], rs_candidates['color_resid'], s=10, c='purple', alpha=0.5)
    
    x_range = np.linspace(rs_candidates['sigma_resid'].min(), rs_candidates['sigma_resid'].max(), 100)
    ax[2].plot(x_range, intercept + slope*x_range, 'k--', label=f'Slope={slope:.3f}')
    
    ax[2].set_xlabel('Residual log(Sigma) (at fixed Mass)')
    ax[2].set_ylabel('Residual u-r Color (at fixed Mass)')
    ax[2].set_title('Test DM: RS Scatter vs Potential')
    ax[2].legend()
    ax[2].grid(True, alpha=0.3)
    
    out_path = os.path.join(FIGURES_DIR, 'sdss_test_dm_rs_scatter.png')
    plt.savefig(out_path, dpi=150)
    print(f"Figure saved to {out_path}")
    
    return {
        'correlation_r': r_val,
        'p_value': p_val,
        'slope': slope,
        'n_gal': int(len(rs_candidates))
    }

def main():
    cache_path = os.path.join(DATA_DIR, 'sdss_rs_scatter.csv')
    
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

    results = analyze_rs_scatter(df)
    
    out_path = os.path.join(RESULTS_DIR, 'sdss_test_dm_results.json')
    if results:
        with open(out_path, 'w') as f:
            json.dump(results, f, indent=2)
            
        print("\nSUMMARY TEST DM:")
        print(f"Slope (Sigma Resid vs Color Resid): {results['slope']:.4f}")
        
        # Interpretation:
        # If High Sigma -> Redder (Slower evolution/Older?): Positive Slope
        # If High Sigma -> Bluer (Younger effective age): Negative Slope
        if results['p_value'] < 0.05 and abs(results['slope']) > 0.05:
             print("RESULT: SIGNAL (Significant Correlation)")
        else:
             print("RESULT: NULL (No significant dependence)")

if __name__ == "__main__":
    main()
