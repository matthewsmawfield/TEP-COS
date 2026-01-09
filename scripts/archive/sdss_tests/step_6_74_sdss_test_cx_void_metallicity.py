#!/usr/bin/env python3
"""
Step 6.74: SDSS Test CX - Void Metallicity (Fast Evolution)

Hypothesis:
Voids are regions of shallow potential (time runs fast). Galaxies there experience more 
proper time than cluster galaxies. Chemical enrichment (SFR integrated over time) should 
proceed further. At fixed stellar mass, void galaxies should be **more metal-rich** 
than field galaxies. (Standard model often predicts the opposite due to retarded infall).

Prediction:
Mass-Metallicity relation in Voids is offset to higher Metallicity.

Data:
- ebossMCPM: mid_dens_1 (Density), PLATE, MJD, FIBERID
- SpecObjAll: specObjID, PLATE, MJD, FIBERID (Linking)
- stellarMassFSPSGranWideDust: logMass
- galSpecExtra: oh_p50 (Metallicity)

Method:
1. Join ebossMCPM to SpecObjAll to get specObjID.
2. Join to stellarMass and galSpecExtra.
3. Classify Void (density < -0.5) vs Field/Cluster (density > 0).
4. Compute Mean Metallicity at fixed Mass bins.
5. Calculate Delta Z (Void - Field).
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
    print(f"Querying SDSS for Test CX (Limit: {limit})...")
    
    # Join chain: ebossMCPM -> SpecObjAll -> stellarMass/galSpecExtra
    
    sql = f"""
    SELECT TOP {limit}
        s.specObjID,
        m.mid_dens_1 as density,
        st.logMass,
        g.oh_p50 as metallicity
        
    FROM ebossMCPM m
    JOIN SpecObjAll s ON m.PLATE = s.PLATE AND m.MJD = s.MJD AND m.FIBERID = s.FIBERID
    JOIN stellarMassFSPSGranWideDust st ON s.specObjID = st.specObjID
    JOIN galSpecExtra g ON s.specObjID = g.specObjID
    
    WHERE 
        st.logMass > 8.5
        AND g.oh_p50 > -9
        AND abs(m.mid_dens_1) < 10 -- Sanity check on density
        AND (m.mid_dens_1 < -0.5 OR m.mid_dens_1 > 0) -- Select Void or Field/Cluster
    """
    return query_sdss(sql)

def analyze_void_metallicity(df):
    print("Analyzing Void Metallicity...")
    
    if df is None or len(df) < 100:
        print("  Insufficient data.")
        return None
        
    # Clean
    df = df.dropna().copy()
    
    # Classify
    df['env'] = 'Field'
    df.loc[df['density'] < -0.5, 'env'] = 'Void'
    
    print("  Environment counts:")
    print(df['env'].value_counts())
    
    voids = df[df['env'] == 'Void']
    field = df[df['env'] == 'Field']
    
    if len(voids) < 20 or len(field) < 20:
        print("  Not enough galaxies in bins.")
        return None
        
    # Fit MZR for Field
    # Z = f(M)
    z_fit = np.polyfit(field['logMass'], field['metallicity'], 2)
    p = np.poly1d(z_fit)
    
    print(f"  Field MZR: {z_fit}")
    
    # Calculate Residuals for Voids
    # Delta Z = Z_void - Z_field_model(M_void)
    # Positive -> Void is Metal Richer
    voids = voids.copy() # Avoid SettingWithCopy
    voids['mzr_resid'] = voids['metallicity'] - p(voids['logMass'])
    
    mean_offset = voids['mzr_resid'].mean()
    sem_offset = voids['mzr_resid'].sem()
    
    print(f"  Mean Void Metallicity Offset: {mean_offset:.4f} +/- {sem_offset:.4f} dex")
    
    # Control comparison: Field residuals (should be 0 by definition of fit)
    field = field.copy()
    field['mzr_resid'] = field['metallicity'] - p(field['logMass'])
    print(f"  Field Mean Residual: {field['mzr_resid'].mean():.4f}")
    
    # T-test
    t_stat, p_val = stats.ttest_ind(voids['mzr_resid'], field['mzr_resid'], equal_var=False)
    print(f"  T-test (Void vs Field): t={t_stat:.2f}, p={p_val:.2e}")
    
    # Plot
    fig, ax = plt.subplots(1, 2, figsize=(12, 5))
    
    # MZR
    ax[0].scatter(field['logMass'], field['metallicity'], alpha=0.1, s=2, c='gray', label='Field')
    ax[0].scatter(voids['logMass'], voids['metallicity'], alpha=0.3, s=5, c='blue', label='Void')
    x_range = np.linspace(8.5, 11.5, 100)
    ax[0].plot(x_range, p(x_range), 'k--', label='Field Fit')
    ax[0].set_xlabel('log Stellar Mass')
    ax[0].set_ylabel('Metallicity (12+log O/H)')
    ax[0].set_title('Mass-Metallicity Relation')
    ax[0].legend()
    ax[0].grid(True, alpha=0.3)
    
    # Residuals Hist
    ax[1].hist(field['mzr_resid'], bins=30, density=True, alpha=0.5, color='gray', label='Field')
    ax[1].hist(voids['mzr_resid'], bins=30, density=True, alpha=0.5, color='blue', label='Void')
    ax[1].axvline(mean_offset, color='blue', linestyle='--', label=f'Void Mean={mean_offset:.3f}')
    ax[1].set_xlabel('MZR Residual (Delta Z)')
    ax[1].set_ylabel('Density')
    ax[1].set_title(f'Test CX: Void Metallicity (Delta={mean_offset:.3f})')
    ax[1].legend()
    ax[1].grid(True, alpha=0.3)
    
    out_path = os.path.join(FIGURES_DIR, 'sdss_test_cx_metallicity.png')
    plt.savefig(out_path, dpi=150)
    print(f"Figure saved to {out_path}")
    
    return {
        'void_offset': mean_offset,
        'void_sem': sem_offset,
        'n_void': int(len(voids)),
        'n_field': int(len(field)),
        'p_value': p_val
    }

def main():
    cache_path = os.path.join(DATA_DIR, 'sdss_void_metals.csv')
    
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

    results = analyze_void_metallicity(df)
    
    out_path = os.path.join(RESULTS_DIR, 'sdss_test_cx_results.json')
    if results:
        with open(out_path, 'w') as f:
            json.dump(results, f, indent=2)
            
        print("\nSUMMARY TEST CX:")
        print("Prediction: Void galaxies are more metal-rich (Positive Offset).")
        print(f"Observed Offset: {results['void_offset']:.4f} dex")
        
        if results['void_offset'] > 0.02 and results['p_value'] < 0.05:
             print("RESULT: SIGNAL (Void galaxies are metal-rich)")
        elif results['void_offset'] < -0.02 and results['p_value'] < 0.05:
             print("RESULT: CONTRADICTED (Void galaxies are metal-poor)")
        else:
             print("RESULT: NULL (No significant difference)")

if __name__ == "__main__":
    main()
