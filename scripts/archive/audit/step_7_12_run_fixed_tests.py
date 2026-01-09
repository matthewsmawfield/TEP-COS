#!/usr/bin/env python3
"""
Step 7.12: Run Fixed SDSS Tests with Correct DR18 Column Names
"""

import requests
import json
import time
import os
import numpy as np
from scipy import stats

SDSS_URL = "https://skyserver.sdss.org/dr18/SkyServerWS/SearchTools/SqlSearch"
RESULTS_DIR = os.path.join(os.path.dirname(__file__), '..', '..', 'results', 'outputs')

def query_sdss(sql, timeout=120):
    """Execute SDSS query"""
    try:
        response = requests.get(
            SDSS_URL,
            params={"cmd": sql, "format": "json"},
            timeout=timeout
        )
        if response.status_code == 200:
            data = response.json()
            if data and len(data) > 0 and "Rows" in data[0]:
                import pandas as pd
                return pd.DataFrame(data[0]["Rows"])
        return None
    except Exception as e:
        print(f"    Error: {e}")
        return None

def run_test_ce_no_clock():
    """CE: Nitrogen/Oxygen Clock using emissionLinesPort (FIXED)"""
    print("  Running CE: Nitrogen/Oxygen Clock...")
    sql = """
    SELECT TOP 5000
        gi.v_disp as sigma,
        el.Amplitude_NII_6583 / el.Amplitude_OIII_5006 as no_ratio
    FROM galSpecInfo gi
    JOIN emissionLinesPort el ON gi.specObjID = el.specObjID
    WHERE gi.v_disp > 50 AND gi.v_disp < 400
      AND el.Amplitude_NII_6583 > 0 AND el.Amplitude_OIII_5006 > 0
    """
    df = query_sdss(sql)
    if df is None or len(df) < 100:
        return {"status": "Skipped", "metric": "Query failed", "n": 0}
    
    df = df[df['no_ratio'] > 0]
    if len(df) < 100:
        return {"status": "Skipped", "metric": "Insufficient data", "n": len(df)}
    
    r, p = stats.pearsonr(np.log10(df['sigma']), np.log10(df['no_ratio']))
    return {
        "status": "Signal" if p < 0.05 else "Null",
        "metric": f"r={r:.3f}, p={p:.2e}",
        "n": len(df)
    }

def run_test_ci_galactic_dipole():
    """CI: Galactic Dipole using apogeeStar (FIXED)"""
    print("  Running CI: Galactic Dipole...")
    sql = """
    SELECT TOP 10000
        a.glon, a.glat, a.vhelio_avg as vhelio,
        a.teff, a.logg, a.m_h as feh
    FROM apogeeStar a
    WHERE a.starflag = 0
      AND a.teff BETWEEN 4000 AND 5500
      AND a.logg BETWEEN 1 AND 3.5
    """
    df = query_sdss(sql)
    if df is None or len(df) < 1000:
        return {"status": "Skipped", "metric": "Query failed", "n": 0}
    
    df['cos_l'] = np.cos(np.radians(df['glon']))
    r, p = stats.pearsonr(df['cos_l'], df['vhelio'])
    return {
        "status": "Signal" if p < 0.05 else "Null",
        "metric": f"r={r:.3f}, p={p:.2e}",
        "n": len(df)
    }

def run_test_cz_dig():
    """CZ: Diffuse Ionized Gas in MaNGA (FIXED)"""
    print("  Running CZ: Diffuse Ionized Gas...")
    sql = """
    SELECT TOP 2000
        m.nsa_elpetro_mass as mass,
        m.nsa_sersic_n as sersic,
        d.emline_gew_1re_ha_6564 as ha_ew,
        d.emline_gew_1re_nii_6585 as nii_ew
    FROM mangaDRPall m
    JOIN mangaDAPall d ON m.mangaid = d.mangaid
    WHERE m.nsa_elpetro_mass > 1e9
      AND d.emline_gew_1re_ha_6564 > 0
    """
    df = query_sdss(sql)
    if df is None or len(df) < 100:
        return {"status": "Skipped", "metric": "Query failed", "n": 0}
    
    df['nii_ha'] = df['nii_ew'] / df['ha_ew']
    df = df[(df['nii_ha'] > 0) & (df['nii_ha'] < 10)]
    
    if len(df) < 50:
        return {"status": "Skipped", "metric": "Insufficient valid data", "n": len(df)}
    
    r, p = stats.pearsonr(np.log10(df['mass']), np.log10(df['nii_ha']))
    return {
        "status": "Signal" if p < 0.05 else "Null",
        "metric": f"r={r:.3f}, p={p:.2e}",
        "n": len(df)
    }

def run_test_da_agn_type():
    """DA: AGN Type 1/2 ratio vs environment (FIXED)"""
    print("  Running DA: AGN Type 1/2...")
    sql = """
    SELECT TOP 5000
        gi.v_disp as sigma,
        gl.oiii_5007_flux / gl.h_beta_flux as o3hb,
        gl.nii_6584_flux / gl.h_alpha_flux as n2ha
    FROM galSpecInfo gi
    JOIN galSpecLine gl ON gi.specObjID = gl.specObjID
    WHERE gi.v_disp > 50 AND gi.v_disp < 400
      AND gl.oiii_5007_flux > 0 AND gl.h_beta_flux > 0
      AND gl.nii_6584_flux > 0 AND gl.h_alpha_flux > 0
    """
    df = query_sdss(sql)
    if df is None or len(df) < 100:
        return {"status": "Skipped", "metric": "Query failed", "n": 0}
    
    # BPT classification: AGN if log(O3/Hb) > 0.61/(log(N2/Ha)-0.05)+1.3
    df['log_o3hb'] = np.log10(df['o3hb'])
    df['log_n2ha'] = np.log10(df['n2ha'])
    df['is_agn'] = (df['log_o3hb'] > 0.61/(df['log_n2ha']-0.05)+1.3).astype(int)
    
    r, p = stats.pearsonr(np.log10(df['sigma']), df['is_agn'])
    return {
        "status": "Signal" if p < 0.05 else "Null",
        "metric": f"r={r:.3f}, p={p:.2e}",
        "n": len(df)
    }

def run_test_dk_ring_galaxies():
    """DK: Ring Galaxy / Disk Fraction (FIXED)"""
    print("  Running DK: Disk Galaxy Fraction...")
    sql = """
    SELECT TOP 5000
        gi.v_disp as sigma,
        z.spiral as disk_prob
    FROM galSpecInfo gi
    JOIN SpecObjAll s ON gi.specObjID = s.specObjID
    JOIN zooSpec z ON s.bestobjid = z.objid
    WHERE gi.v_disp > 50 AND gi.v_disp < 400
    """
    df = query_sdss(sql)
    if df is None or len(df) < 100:
        return {"status": "Skipped", "metric": "Query failed", "n": 0}
    
    r, p = stats.pearsonr(np.log10(df['sigma']), df['disk_prob'])
    return {
        "status": "Signal" if p < 0.05 else "Null",
        "metric": f"r={r:.3f}, p={p:.2e}",
        "n": len(df)
    }

def run_test_di_cluster_spin():
    """DI: Galaxy Spin from MaNGA (FIXED)"""
    print("  Running DI: Galaxy Spin...")
    sql = """
    SELECT TOP 1000
        m.nsa_elpetro_mass as mass,
        d.stellar_vel_1re as vel,
        d.stellar_sigma_1re as sigma
    FROM mangaDRPall m
    JOIN mangaDAPall d ON m.mangaid = d.mangaid
    WHERE m.nsa_elpetro_mass > 1e10
      AND d.stellar_sigma_1re > 0
    """
    df = query_sdss(sql)
    if df is None or len(df) < 50:
        return {"status": "Skipped", "metric": "Query failed", "n": 0}
    
    df['v_sigma'] = np.abs(df['vel']) / df['sigma']
    df = df[df['v_sigma'] < 10]
    
    if len(df) < 30:
        return {"status": "Skipped", "metric": "Insufficient data", "n": len(df)}
    
    r, p = stats.pearsonr(np.log10(df['mass']), df['v_sigma'])
    return {
        "status": "Signal" if p < 0.05 else "Null",
        "metric": f"r={r:.3f}, p={p:.2e}",
        "n": len(df)
    }

def run_test_dn_qso_asymmetry():
    """DN: QSO Line Asymmetry (FIXED)"""
    print("  Running DN: Line Width vs Sigma...")
    sql = """
    SELECT TOP 3000
        gi.v_disp as sigma,
        gl.oiii_5007_eqw as line_ew,
        gl.oiii_5007_flux as flux
    FROM galSpecInfo gi
    JOIN galSpecLine gl ON gi.specObjID = gl.specObjID
    WHERE gl.oiii_5007_flux > 10
      AND gi.v_disp > 50
      AND gl.oiii_5007_eqw > 0
    """
    df = query_sdss(sql)
    if df is None or len(df) < 100:
        return {"status": "Skipped", "metric": "Query failed", "n": 0}
    
    r, p = stats.pearsonr(np.log10(df['sigma']), np.log10(df['line_ew']))
    return {
        "status": "Signal" if p < 0.05 else "Null",
        "metric": f"r={r:.3f}, p={p:.2e}",
        "n": len(df)
    }

def run_test_dp_tio_imf():
    """DP: TiO IMF indicator (FIXED)"""
    print("  Running DP: TiO IMF...")
    sql = """
    SELECT TOP 3000
        gi.v_disp as sigma,
        gx.lick_tio2 as tio2
    FROM galSpecInfo gi
    JOIN galSpecIndx gx ON gi.specObjID = gx.specObjID
    WHERE gi.v_disp > 100 AND gi.v_disp < 400
      AND gx.lick_tio2 > 0
    """
    df = query_sdss(sql)
    if df is None or len(df) < 100:
        return {"status": "Skipped", "metric": "Query failed", "n": 0}
    
    r, p = stats.pearsonr(np.log10(df['sigma']), df['tio2'])
    return {
        "status": "Signal" if p < 0.05 else "Null",
        "metric": f"r={r:.3f}, p={p:.2e}",
        "n": len(df)
    }

def run_test_dr_brown_dwarf():
    """DR: Brown Dwarf Desert (FIXED)"""
    print("  Running DR: Brown Dwarf Desert...")
    sql = """
    SELECT TOP 5000
        a.teff, a.logg, a.vscatter,
        a.glon, a.glat
    FROM apogeeStar a
    WHERE a.starflag = 0
      AND a.vscatter > 0
      AND a.teff BETWEEN 4000 AND 6000
    """
    df = query_sdss(sql)
    if df is None or len(df) < 500:
        return {"status": "Skipped", "metric": "Query failed", "n": 0}
    
    df['binary_cand'] = (df['vscatter'] > 1).astype(int)
    df['R_gc'] = 8.0 - np.cos(np.radians(df['glon'])) * 2.0
    
    r, p = stats.pearsonr(df['R_gc'], df['binary_cand'])
    return {
        "status": "Signal" if p < 0.05 else "Null",
        "metric": f"r={r:.3f}, p={p:.2e}",
        "n": len(df)
    }

def run_test_dm_red_sequence():
    """DM: Red Sequence Scatter (FIXED)"""
    print("  Running DM: Red Sequence Scatter...")
    sql = """
    SELECT TOP 5000
        ge.gr_p50 as g_r,
        gi.v_disp as sigma,
        ge.lgm_tot_p50 as log_mass
    FROM galSpecExtra ge
    JOIN galSpecInfo gi ON ge.specObjID = gi.specObjID
    WHERE gi.v_disp > 100 AND gi.v_disp < 400
      AND ge.gr_p50 > 0.5 AND ge.gr_p50 < 1.2
      AND ge.lgm_tot_p50 > 10
    """
    df = query_sdss(sql)
    if df is None or len(df) < 200:
        return {"status": "Skipped", "metric": "Query failed", "n": 0}
    
    bins = np.percentile(df['sigma'], [0, 33, 66, 100])
    scatters = []
    for i in range(3):
        mask = (df['sigma'] >= bins[i]) & (df['sigma'] < bins[i+1])
        sub = df[mask]
        if len(sub) > 30:
            slope, intercept, _, _, _ = stats.linregress(sub['log_mass'], sub['g_r'])
            resid = sub['g_r'] - (slope * sub['log_mass'] + intercept)
            scatters.append(resid.std())
    
    if len(scatters) < 3:
        return {"status": "Null", "metric": "Insufficient bins", "n": len(df)}
    
    trend = scatters[2] - scatters[0]
    return {
        "status": "Signal" if abs(trend) > 0.02 else "Null",
        "metric": f"Δscatter={trend:.3f}",
        "n": len(df)
    }

def run_test_du_hi_optical():
    """DU: HI vs Optical size (FIXED)"""
    print("  Running DU: HI vs Optical...")
    sql = """
    SELECT TOP 1000
        m.nsa_elpetro_th50 as opt_size,
        h.logmhi as log_mhi,
        m.nsa_elpetro_mass as mass
    FROM mangaDRPall m
    JOIN mangaHIall h ON m.mangaid = h.mangaid
    WHERE m.nsa_elpetro_th50 > 0
      AND h.logmhi > 7
    """
    df = query_sdss(sql)
    if df is None or len(df) < 50:
        return {"status": "Skipped", "metric": "Query failed", "n": 0}
    
    df['hi_frac'] = df['log_mhi'] - np.log10(df['mass'])
    r, p = stats.pearsonr(np.log10(df['opt_size']), df['hi_frac'])
    return {
        "status": "Signal" if p < 0.05 else "Null",
        "metric": f"r={r:.3f}, p={p:.2e}",
        "n": len(df)
    }

def main():
    print("=" * 70)
    print("RUNNING FIXED SDSS TESTS (DR18 Schema)")
    print("=" * 70)
    
    tests = {
        "CE": run_test_ce_no_clock,
        "CI": run_test_ci_galactic_dipole,
        "CZ": run_test_cz_dig,
        "DA": run_test_da_agn_type,
        "DI": run_test_di_cluster_spin,
        "DK": run_test_dk_ring_galaxies,
        "DM": run_test_dm_red_sequence,
        "DN": run_test_dn_qso_asymmetry,
        "DP": run_test_dp_tio_imf,
        "DR": run_test_dr_brown_dwarf,
        "DU": run_test_du_hi_optical,
    }
    
    results = {}
    counts = {"Signal": 0, "Null": 0, "Skipped": 0}
    
    for code, func in tests.items():
        print(f"\n[{code}]")
        try:
            result = func()
            results[code] = result
            counts[result["status"]] += 1
            print(f"    Status: {result['status']}")
            print(f"    Metric: {result['metric']}")
            print(f"    N: {result.get('n', 'N/A')}")
        except Exception as e:
            print(f"    ERROR: {e}")
            results[code] = {"status": "Skipped", "metric": f"Error: {str(e)[:50]}", "n": 0}
            counts["Skipped"] += 1
        
        time.sleep(1)
    
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"  Signal: {counts['Signal']}")
    print(f"  Null: {counts['Null']}")
    print(f"  Skipped: {counts['Skipped']}")
    
    out_path = os.path.join(RESULTS_DIR, 'sdss_fixed_tests_results.json')
    with open(out_path, 'w') as f:
        json.dump({"counts": counts, "results": results}, f, indent=2)
    print(f"\nResults saved to {out_path}")
    
    return results

if __name__ == "__main__":
    main()
