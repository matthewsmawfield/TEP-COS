#!/usr/bin/env python3
"""
Step 6.50: SDSS Test BV - BAL Quasar Fraction (Wind Launching Clock)

Hypothesis:
Broad Absorption Line (BAL) outflows are driven by radiation pressure and magnetocentrifugal winds. 
These are rate processes. In deep potential wells (massive BHs), wind launching might be sustained 
for longer periods or launching velocities might be altered by the scalar field. 
We expect the fraction of quasars showing BAL features to correlate with the potential depth (M_BH).

Prediction:
BAL Fraction correlates with M_BH (or Eddington Ratio).
TEP: Higher BAL fraction at high Mass (slower evolution/longer phase).

Data:
- mos_sdss_dr16_qso: bal_prob (BAL probability)
- spiders_quasar: logBHMS_mgII, logBHMA_hb (BH Mass)

Method:
1. Fetch BAL prob from DR16Q and Mass from SPIDERS.
2. Join tables (on SPECOBJID or coordinates).
3. Define BAL flag (bal_prob > 0.5).
4. Bin by BH Mass.
5. Calculate BAL Fraction in each bin.
6. Analyze trend (Slope of Fraction vs Mass).
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
    print(f"Querying SDSS for Test BV (Limit: {limit})...")
    
    # Join DR16Q and SPIDERS on SPECOBJID if possible, or assume DR16Q has mass? 
    # DR16Q usually has mass estimates too, but previous checks showed limited columns.
    # Let's try joining spiders_quasar (Mass) and mos_sdss_dr16_qso (BAL).
    # Common key: specObjID / SPECOBJID? Or coordinate match.
    # DR16Q has 'specObjID' usually? Check script didn't show it explicitly but had 'pk'.
    # Let's try coordinate match or ThingID if available.
    # Both have RA/DEC.
    
    # Try simple join on coordinates (approx) or just assume SPIDERS is a subset.
    # SPIDERS has SPECOBJID. DR16Q might not have it exposed as such?
    # Let's rely on coordinate match in SQL.
    # Actually, simpler: spiders_quasar might have BAL info? No, check said no.
    # DR16Q has 'bal_prob'.
    
    sql = f"""
    SELECT TOP {limit}
        s.logBHMS_mgII,
        s.logBHMA_hb,
        q.bal_prob
        
    FROM spiders_quasar s
    JOIN mos_sdss_dr16_qso q ON s.SPECOBJID = q.specObjID -- Optimistic guess on key
    -- If key fails, we might need spatial join
    
    WHERE 
        (s.logBHMS_mgII > 6 OR s.logBHMA_hb > 6)
    """
    
    # Let's try a safer spatial join if IDs don't match
    # But spatial join is slow.
    # Let's assume DR16Q 'specObjID' exists (standard SDSS).
    # Wait, check_qso_all_cols output didn't show specObjID! It showed 'pk', 'id_number', 'objid' (photometric).
    # But it has 'plate', 'mjd', 'fiberid'.
    # spiders_quasar has 'Plate', 'MJD', 'FiberID'.
    
    sql = f"""
    SELECT TOP {limit}
        s.logBHMS_mgII,
        s.logBHMA_hb,
        q.bal_prob
        
    FROM spiders_quasar s
    JOIN mos_sdss_dr16_qso q ON s.Plate = q.plate AND s.MJD = q.mjd AND s.FiberID = q.fiberid
    
    WHERE 
        (s.logBHMS_mgII > 6 OR s.logBHMA_hb > 6)
        AND q.bal_prob >= 0
    """
    
    return query_sdss(sql)

def analyze_bal_fraction(df):
    print("Analyzing BAL Fraction...")
    
    # Consolidate Mass
    df['logBH'] = df['logBHMS_mgII']
    mask_nan = df['logBH'].isna() | (df['logBH'] == -99)
    df.loc[mask_nan, 'logBH'] = df.loc[mask_nan, 'logBHMA_hb']
    
    df = df.dropna(subset=['logBH', 'bal_prob'])
    df = df[df['logBH'] > 6]
    
    print(f"  Sample size: {len(df)}")
    
    # Define BAL
    # bal_prob is probability. Use threshold > 0.5?
    # Or just use mean probability as fraction.
    
    df['is_bal'] = df['bal_prob'] > 0.5
    
    # 1. Logistic Regression (Raw)
    # logit(p) = a * logBH + b
    # Use scipy or simple binning trend first.
    
    # 2. Binning
    df['mass_bin'] = pd.qcut(df['logBH'], 8)
    binned = df.groupby('mass_bin')['is_bal'].agg(['mean', 'sem', 'count'])
    binned['mass_center'] = [i.mid for i in binned.index]
    
    print("\nBAL Fraction by Mass Bin:")
    print(binned[['mean', 'sem', 'count']])
    
    # Fit trend to binned data
    slope, intercept, r_val, p_val, std_err = stats.linregress(binned['mass_center'], binned['mean'])
    print(f"  Slope (Fraction vs logM): {slope:.4f}")
    print(f"  Correlation r: {r_val:.4f} (p={p_val:.2e})")
    
    # Plot
    fig, ax = plt.subplots(figsize=(8, 6))
    
    ax.errorbar(binned['mass_center'], binned['mean'], yerr=binned['sem'], fmt='o-', capsize=5)
    ax.set_xlabel('log(BH Mass) [M_sun]')
    ax.set_ylabel('BAL Quasar Fraction')
    ax.set_title(f'Test BV: BAL Fraction vs Mass (r={r_val:.2f})')
    ax.grid(True, alpha=0.3)
    
    out_path = os.path.join(FIGURES_DIR, 'sdss_test_bv_bal.png')
    plt.savefig(out_path, dpi=150)
    print(f"Figure saved to {out_path}")
    
    return {
        'slope': slope,
        'r_val': r_val,
        'n_sample': int(len(df))
    }

def main():
    cache_path = os.path.join(DATA_DIR, 'sdss_bal_fraction.csv')
    
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

    results = analyze_bal_fraction(df)
    
    out_path = os.path.join(RESULTS_DIR, 'sdss_test_bv_results.json')
    if results:
        with open(out_path, 'w') as f:
            json.dump(results, f, indent=2)
            
        print("\nSUMMARY TEST BV:")
        print("TEP Prediction: BAL Fraction correlates with Potential Depth (Mass).")
        print(f"Observed Slope: {results['slope']:.4f}")
        
        if abs(results['r_val']) > 0.5:
             print("RESULT: CONSISTENT (Strong correlation)")
        else:
             print("RESULT: NULL (Weak/No correlation)")

if __name__ == "__main__":
    main()
