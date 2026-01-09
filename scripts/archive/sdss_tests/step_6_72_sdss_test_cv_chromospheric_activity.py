#!/usr/bin/env python3
"""
Step 6.72: SDSS Test CV - Chromospheric Activity (Magnetic Braking Clock)

Hypothesis:
Stellar magnetic activity (e.g., Ca II emission) decays as stars spin down over Gyrs. 
Spin-down is a torque-driven rate process. In the deep potential of the Inner Galaxy, 
spin-down is time-dilated (slower). Stars of a given formation epoch should retain 
higher rotation and higher activity levels than their counterparts in the Outer Galaxy. 
Activity fills in absorption lines, reducing the measured Lick index.

Prediction:
Ca4227 Absorption Index is LOWER (more filled) at small R_gc (at fixed Teff/Met).

Data:
- galSpecIndx: lick_ca4227 (Absorption Index)
- sppParams: TEFFADOP, LOGGADOP, FEHADOP
- sppTargets: DISTV_KPC, RA, DEC

Method:
1. Select FGK Dwarfs (4000 < Teff < 6000, logg > 4.0).
2. Compute R_gc.
3. Regress lick_ca4227 on Teff and Fe/H to remove stellar parameter dependence.
4. Analyze Residuals vs R_gc.
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

def download_data(limit=5000):
    print(f"Querying SDSS for Test CV (Limit: {limit})...")
    
    # Join galSpecIndx (Lick), sppParams (Stellar), sppTargets (Dist)
    # Using specObjID link
    
    sql = f"""
    SELECT TOP {limit}
        g.specObjID,
        g.lick_ca4227,
        g.lick_ca4227_err,
        s.TEFFADOP as teff,
        s.LOGGADOP as logg,
        s.FEHADOP as fe_h,
        t.DISTV_KPC as dist,
        t.RA as ra,
        t.DEC as dec
        
    FROM galSpecIndx g
    JOIN sppParams s ON g.specObjID = s.specObjID
    JOIN sppTargets t ON g.specObjID = t.SPECOBJID -- Verify if SPECOBJID exists in Targets
    
    WHERE 
        s.TEFFADOP BETWEEN 4000 AND 6000 -- FGK
        AND s.LOGGADOP > 4.0 -- Dwarfs
        AND t.DISTV_KPC > 0
        AND g.lick_ca4227 > -10
    """
    # Note: sppTargets might use SPECOBJID or BESTOBJID. 
    # Check script said SPECOBJID is present.
    
    return query_sdss(sql)

def compute_rgc(df):
    R0 = 8.2 # kpc
    df['dist_kpc'] = df['dist']
    
    # Convert RA/Dec to Galactic L/B
    # Approximation or use astropy. 
    # Simple approx for now or assume most SEGUE are high lat?
    # Let's use astropy if available, else approx.
    try:
        from astropy.coordinates import SkyCoord
        import astropy.units as u
        c = SkyCoord(ra=df['ra'].values*u.degree, dec=df['dec'].values*u.degree, frame='icrs')
        g = c.galactic
        l_rad = g.l.radian
        b_rad = g.b.radian
    except ImportError:
        print("  Astropy not found, skipping R_gc calculation (or implementing approx).")
        return df
    
    d_proj = df['dist_kpc'] * np.cos(b_rad)
    
    df['R_plane_sq'] = R0**2 + d_proj**2 - 2 * R0 * d_proj * np.cos(l_rad)
    df['Z'] = df['dist_kpc'] * np.sin(b_rad)
    df['R_gc'] = np.sqrt(df['R_plane_sq'] + df['Z']**2)
    
    return df

def analyze_activity(df):
    print("Analyzing Chromospheric Activity...")
    
    if df is None or len(df) < 50:
        print("  Insufficient data.")
        return None
        
    df = compute_rgc(df)
    if 'R_gc' not in df.columns:
        return None
        
    # Clean
    df = df.dropna().copy()
    # Filter bad Lick measurements
    df = df[df['lick_ca4227_err'] < 0.5]
    
    print(f"  Sample size: {len(df)}")
    
    # 1. Control for Teff and Metallicity
    # Index depends on Temp and Z.
    # Model: Index ~ a*Teff + b*FeH + c
    
    X = np.column_stack([df['teff'], df['fe_h'], np.ones(len(df))])
    y = df['lick_ca4227']
    
    try:
        # Robust fit or simple least squares
        theta, residuals, rank, s = np.linalg.lstsq(X, y, rcond=None)
        model = X @ theta
        df['lick_resid'] = y - model
        print(f"  Parameter dependence removed.")
    except Exception as e:
        print(f"  Fit failed: {e}")
        return None
        
    # 2. Correlate Residuals with R_gc
    r_val, p_val = stats.pearsonr(df['R_gc'], df['lick_resid'])
    print(f"  Correlation r(Residual, R_gc): {r_val:.4f} (p={p_val:.2e})")
    
    slope, intercept, _, _, _ = stats.linregress(df['R_gc'], df['lick_resid'])
    print(f"  Slope (Index Resid vs R): {slope:.5f} Angstrom/kpc")
    
    # Binning
    bins = np.linspace(df['R_gc'].min(), df['R_gc'].max(), 10)
    df['r_bin'] = pd.cut(df['R_gc'], bins=bins)
    binned = df.groupby('r_bin')['lick_resid'].agg(['mean', 'sem', 'count'])
    binned['r_center'] = [i.mid for i in binned.index]
    
    # Plot
    fig, ax = plt.subplots(1, 2, figsize=(12, 5))
    
    # Teff dependence
    ax[0].scatter(df['teff'], df['lick_ca4227'], alpha=0.2, s=5)
    ax[0].set_xlabel('Teff [K]')
    ax[0].set_ylabel('Ca4227 Index [Angstrom]')
    ax[0].set_title('Index vs Temperature')
    ax[0].grid(True, alpha=0.3)
    
    # Residual vs R_gc
    ax[1].errorbar(binned['r_center'], binned['mean'], yerr=binned['sem'], fmt='o-', capsize=5)
    ax[1].set_xlabel('Galactocentric Radius [kpc]')
    ax[1].set_ylabel('Index Residual (Activity Proxy)')
    ax[1].set_title(f'Test CV: Activity vs Position (r={r_val:.2f})')
    ax[1].grid(True, alpha=0.3)
    
    out_path = os.path.join(FIGURES_DIR, 'sdss_test_cv_activity.png')
    plt.savefig(out_path, dpi=150)
    print(f"Figure saved to {out_path}")
    
    return {
        'slope': slope,
        'r_val': r_val,
        'n_sample': int(len(df))
    }

def main():
    cache_path = os.path.join(DATA_DIR, 'sdss_activity.csv')
    
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

    results = analyze_activity(df)
    
    out_path = os.path.join(RESULTS_DIR, 'sdss_test_cv_results.json')
    if results:
        def default(o):
            if isinstance(o, pd.Interval): return str(o)
            raise TypeError
            
        with open(out_path, 'w') as f:
            json.dump(results, f, indent=2, default=default)
            
        print("\nSUMMARY TEST CV:")
        print("Prediction: Index Residuals increase with R_gc (Less filling/Lower activity in outskirts).")
        print("Wait, TEP says Inner Galaxy (Deep Potential) -> Slower spin-down -> Higher Activity -> More Filling -> Lower Index.")
        print("So Index should be LOWER at Small R. Slope should be POSITIVE (Index increases with R).")
        print(f"Observed Slope: {results['slope']:.5f}")
        
        if results['slope'] > 0.005:
             print("RESULT: CONSISTENT (Lower index/Higher activity in inner galaxy)")
        elif results['slope'] < -0.005:
             print("RESULT: CONTRADICTED (Higher index/Lower activity in inner galaxy)")
        else:
             print("RESULT: NULL (No gradient)")

if __name__ == "__main__":
    main()
