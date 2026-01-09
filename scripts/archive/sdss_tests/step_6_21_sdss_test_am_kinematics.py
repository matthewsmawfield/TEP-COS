#!/usr/bin/env python3
"""
Step 6.21: SDSS Test AM - Gas-Stellar Kinematic Alignment

Hypothesis:
In CDM, dark matter halos can be triaxial and misaligned with the baryonic disk, leading to warps and gas-star misalignments.
TEP predicts that the "phantom mass" is a metric effect tied directly to the baryonic potential (plus the soliton wake).
Prediction: Tighter alignment between gas and stellar kinematic axes (or photometric vs kinematic axes) than CDM predicts, especially in low-mass galaxies where DM usually dominates.

Data:
- mangaDAPall: Photometric Position Angle (nsa_elpetro_phi).
- mangaPipe3D: Kinematic Position Angle (PA).

Observable: Misalignment Angle Delta_PA = |PA_kin - PA_phot|.
Test: r(Delta_PA, Mass) or r(Delta_PA, Sigma).
TEP predicts Delta_PA should be SMALLER (more aligned) in low mass galaxies compared to standard expectations?
Actually, the text says "tighter alignment... especially in low-mass galaxies".
Standard CDM predicts misalignments are common in low mass (dwarfs, irregulars).
So TEP predicts Delta_PA is low everywhere, or specifically low at low mass.

We will just measure the distribution of Delta_PA vs Mass.
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
    print(f"Querying SDSS MaNGA for Test AM (Limit: {limit})...")
    
    # Join DAP and Pipe3D on plateifu or mangaid
    # Pipe3D PA is usually the kinematic PA of the gas/stars (check doc, usually stellar kinematics PA).
    # DAP nsa_elpetro_phi is photometric PA.
    
    sql = f"""
    SELECT TOP {limit}
        d.mangaid,
        d.nsa_elpetro_phi as PA_phot,
        d.nsa_elpetro_ba as ba, -- axis ratio to select disks
        d.stellar_sigma_1re as sigma,
        p.PA as PA_kin,
        p.log_Mass as logMass -- Pipe3D mass
        
    FROM mangaDAPall d
    JOIN mangaPipe3D p ON d.mangaid = p.mangaid
    
    WHERE 
        d.nsa_elpetro_ba < 0.8 -- Defined axis (not round)
        AND d.stellar_sigma_1re > 0
    """
    return query_sdss(sql)

def analyze_alignment(df):
    print("Analyzing Kinematic Alignment...")
    
    # 1. Clean
    df_clean = df.dropna(subset=['PA_phot', 'PA_kin', 'logMass', 'sigma']).copy()
    
    # 2. Calculate Misalignment
    # Angles in degrees. Range 0-180 usually.
    # Misalignment is min(|PA1-PA2|, 180-|PA1-PA2|)
    
    delta = np.abs(df_clean['PA_phot'] - df_clean['PA_kin'])
    df_clean['misalignment'] = np.minimum(delta, 180 - delta) # If range is 0-180
    # Sometimes PA is 0-360.
    # Assuming standard 0-180 for galaxy PA.
    
    print(f"N = {len(df_clean)}")
    print(f"Mean Misalignment: {df_clean['misalignment'].mean():.2f} deg")
    
    # 3. Correlation with Mass
    # TEP: Tighter alignment at low mass? -> Positive correlation (Low mass = Low misalignment).
    # Wait, "tighter alignment... especially in low-mass".
    # Means misalignment should be LOW at LOW Mass.
    # Standard CDM: Low mass galaxies are messy/triaxial/dm dominated -> Higher misalignment?
    # So TEP predicts r(Misalignment, Mass) > 0 ?? (Or just low misalignment everywhere).
    
    # Let's check correlation.
    r_mass, p_mass = stats.pearsonr(df_clean['logMass'], df_clean['misalignment'])
    r_sigma, p_sigma = stats.pearsonr(df_clean['sigma'], df_clean['misalignment'])
    
    print(f"r(Misalignment, Mass): {r_mass:.4f} (p={p_mass:.2e})")
    print(f"r(Misalignment, Sigma): {r_sigma:.4f} (p={p_sigma:.2e})")
    
    # 4. Bin Analysis
    df_clean['mass_bin'] = pd.qcut(df_clean['logMass'], 5)
    binned = df_clean.groupby('mass_bin')['misalignment'].mean()
    print("\nMean Misalignment by Mass Bin:")
    print(binned)
    
    return {
        'r_mass': float(r_mass),
        'p_mass': float(p_mass),
        'mean_misalignment': float(df_clean['misalignment'].mean()),
        'n_sample': int(len(df_clean))
    }, df_clean

def create_figure(df, results):
    print("Generating figure...")
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # Plot Misalignment vs Mass
    ax.scatter(df['logMass'], df['misalignment'], alpha=0.1, s=2, c='black')
    
    # Running mean
    bins = np.linspace(df['logMass'].min(), df['logMass'].max(), 10)
    centers = 0.5 * (bins[1:] + bins[:-1])
    means = []
    for i in range(len(bins)-1):
        mask = (df['logMass'] >= bins[i]) & (df['logMass'] < bins[i+1])
        means.append(df.loc[mask, 'misalignment'].mean())
        
    ax.plot(centers, means, 'r-o', lw=2, label='Mean Misalignment')
    
    ax.set_xlabel(r'$\log(M_*)$')
    ax.set_ylabel(r'Misalignment $|PA_{phot} - PA_{kin}|$ (deg)')
    ax.set_title("Kinematic-Photometric Alignment vs Mass")
    ax.legend()
    ax.set_ylim(0, 90)
    
    plt.tight_layout()
    out_path = os.path.join(FIGURES_DIR, 'sdss_test_am_kinematics.png')
    plt.savefig(out_path, dpi=150)
    print(f"Figure saved to {out_path}")

def main():
    cache_path = os.path.join(DATA_DIR, 'sdss_manga_kinematics.csv')
    
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

    results, df_clean = analyze_alignment(df)
    
    out_path = os.path.join(RESULTS_DIR, 'sdss_test_am_results.json')
    with open(out_path, 'w') as f:
        json.dump(results, f, indent=2)
        
    create_figure(df_clean, results)
    
    print("\nSUMMARY TEST AM:")
    print("TEP Prediction: Tight alignment (Low mean misalignment, especially at low mass)")
    print(f"Mean Misalignment: {results['mean_misalignment']:.2f} deg")
    print(f"Correlation with Mass: {results['r_mass']:.4f}")
    
    # Interpretation is complex, but let's report the stats.

if __name__ == "__main__":
    main()
