#!/usr/bin/env python3
"""
Diagnose SDSS test failures and identify fixable queries.
"""

import requests
import json
import time
import os

SDSS_URL = "https://skyserver.sdss.org/dr18/SkyServerWS/SearchTools/SqlSearch"
RESULTS_DIR = os.path.join(os.path.dirname(__file__), '..', '..', 'results', 'outputs')

def test_query(sql, timeout=30):
    """Test if a query works"""
    try:
        response = requests.get(
            SDSS_URL,
            params={"cmd": sql, "format": "json"},
            timeout=timeout
        )
        if response.status_code == 200:
            data = response.json()
            if data and len(data) > 0:
                if "Rows" in data[0]:
                    return True, len(data[0]["Rows"]), None
                elif "error" in str(data).lower():
                    return False, 0, str(data)[:200]
            return True, 0, "Empty result"
        else:
            return False, 0, f"HTTP {response.status_code}"
    except requests.exceptions.Timeout:
        return False, 0, "Timeout"
    except Exception as e:
        return False, 0, str(e)[:100]

def check_table_exists(table_name):
    """Check if a table exists in SDSS"""
    sql = f"SELECT TOP 1 * FROM {table_name}"
    success, count, error = test_query(sql)
    return success

def check_column_exists(table_name, column_name):
    """Check if a column exists in a table"""
    sql = f"SELECT TOP 1 {column_name} FROM {table_name}"
    success, count, error = test_query(sql)
    return success

# Define all skipped tests with their likely issues
SKIPPED_TESTS = {
    "BT": {
        "name": "Galaxy Merger Rate",
        "issue": "zoo2MainSpecz",
        "tables": ["zoo2MainSpecz", "Galaxy"],
        "alt_tables": ["zooSpec", "Zoo1Spec"],
    },
    "BY": {
        "name": "Astrometric Binary Fraction",
        "issue": "Gaia tables",
        "tables": ["gaiadr3_source"],
        "alt_tables": [],  # Gaia not in SDSS
        "skip_reason": "Gaia data not in SDSS - requires external join"
    },
    "CA": {
        "name": "BLR Kinematics",
        "issue": "No overlap",
        "tables": ["SpecObjAll", "galSpecLine"],
    },
    "CD": {
        "name": "Satellite Concentration",
        "issue": "redMaPPer",
        "tables": ["redMaPPer"],
        "alt_tables": [],
        "skip_reason": "redMaPPer cluster catalog not in SDSS DR18"
    },
    "CE": {
        "name": "Nitrogen/Oxygen Clock",
        "issue": "emissionLinesPort",
        "tables": ["emissionLinesPort"],
        "alt_tables": ["galSpecLine", "emissionLinesDR8"],
    },
    "CF": {
        "name": "QSO Gravitational Redshift",
        "issue": "mos_sdss_dr16_qso",
        "tables": ["sdss_dr16_qso"],
        "alt_tables": ["SpecObjAll", "QsoCatalogAll"],
    },
    "CG": {
        "name": "Cluster Sigma Profiles",
        "issue": "redMaPPer",
        "tables": ["redMaPPer"],
        "skip_reason": "redMaPPer not in SDSS DR18"
    },
    "CI": {
        "name": "Galactic Dipole",
        "issue": "aspcapStar",
        "tables": ["aspcapStar"],
        "alt_tables": ["apogeeStar", "allStar"],
    },
    "CM": {
        "name": "TRGB Magnitude",
        "issue": "Gaia",
        "tables": ["gaiadr3_source"],
        "skip_reason": "Requires Gaia parallax"
    },
    "CN": {
        "name": "BCG X-ray Offsets",
        "issue": "redMaPPer",
        "tables": ["redMaPPer"],
        "skip_reason": "redMaPPer not in SDSS DR18"
    },
    "CO": {
        "name": "Exoplanet Yield",
        "issue": "MARVELS",
        "tables": ["marvelsStar"],
        "alt_tables": [],
        "skip_reason": "MARVELS tables may not be in DR18"
    },
    "CQ": {
        "name": "Cluster Cooling Flows",
        "issue": "eFEDS",
        "tables": ["eFEDS"],
        "skip_reason": "eROSITA eFEDS not in SDSS"
    },
    "CS": {
        "name": "Void HI Fraction",
        "issue": "mangaHIall",
        "tables": ["mangaHIall"],
        "alt_tables": ["mangaDRPall", "mangaFirefly"],
    },
    "CT": {
        "name": "Schechter Cutoff",
        "issue": "HTTP 500",
        "tables": ["PhotoObjAll"],
    },
    "CU": {
        "name": "Binary Quasar Fraction",
        "issue": "No Primary Key",
        "tables": ["SpecObjAll"],
    },
    "CW": {
        "name": "Stellar Twins",
        "issue": "HTTP 500",
        "tables": ["apogeeStar"],
    },
    "CX": {
        "name": "Void Metallicity",
        "issue": "HTTP 500",
        "tables": ["galSpecExtra", "galSpecInfo"],
    },
    "CY": {
        "name": "Quasar Clustering",
        "issue": "HTTP 500",
        "tables": ["QsoCatalogAll", "SpecObjAll"],
    },
    "CZ": {
        "name": "Diffuse Ionized Gas",
        "issue": "HTTP 500",
        "tables": ["mangaDAPall", "mangaDRPall"],
    },
    "DA": {
        "name": "AGN Type 1/2",
        "issue": "HTTP 500",
        "tables": ["galSpecLine", "galSpecInfo"],
    },
    "DB": {
        "name": "Void Hubble Drift",
        "issue": "HTTP 500",
        "tables": ["SpecObjAll", "PhotoObjAll"],
    },
    "DC": {
        "name": "Pair Decay Ratio",
        "issue": "Neighbors table",
        "tables": ["Neighbors"],
        "alt_tables": ["PhotoObjAll"],  # Can compute neighbors ourselves
    },
    "DD": {
        "name": "QSO Color-Potential",
        "issue": "HTTP 500",
        "tables": ["QsoCatalogAll", "SpecObjAll"],
    },
    "DE": {
        "name": "Richness-Sigma Tension",
        "issue": "redMaPPer",
        "tables": ["redMaPPer"],
        "skip_reason": "redMaPPer not in SDSS DR18"
    },
    "DF": {
        "name": "Lithium Survival",
        "issue": "No Li columns",
        "tables": ["apogeeStar"],
        "columns": ["li_fe", "LI_FE"],
    },
    "DG": {
        "name": "ICL Growth",
        "issue": "No Cluster Data",
        "tables": ["redMaPPer"],
        "skip_reason": "redMaPPer not in SDSS DR18"
    },
    "DI": {
        "name": "Cluster Stellar Spin",
        "issue": "HTTP 500",
        "tables": ["mangaDAPall", "mangaDRPall"],
    },
    "DK": {
        "name": "Ring Galaxy Fraction",
        "issue": "HTTP 500",
        "tables": ["zooSpec", "Galaxy"],
    },
    "DL": {
        "name": "Tidal Debris",
        "issue": "HTTP 500",
        "tables": ["PhotoObjAll", "SpecObjAll"],
    },
    "DM": {
        "name": "Red Sequence Scatter",
        "issue": "HTTP 500",
        "tables": ["galSpecExtra", "galSpecInfo"],
    },
    "DN": {
        "name": "QSO Line Asymmetry",
        "issue": "HTTP 500",
        "tables": ["SpecObjAll", "galSpecLine"],
    },
    "DO": {
        "name": "Cluster Lx-Sigma",
        "issue": "redMaPPer",
        "tables": ["redMaPPer"],
        "skip_reason": "redMaPPer not in SDSS DR18"
    },
    "DP": {
        "name": "TiO IMF",
        "issue": "HTTP 500",
        "tables": ["galSpecIndx", "galSpecInfo"],
    },
    "DR": {
        "name": "Brown Dwarf Desert",
        "issue": "HTTP 500",
        "tables": ["apogeeStar"],
    },
    "DS": {
        "name": "QSO Variability",
        "issue": "Timeout",
        "tables": ["PhotoObjAll"],  # Time-domain requires Stripe 82
    },
    "DU": {
        "name": "HI vs Optical",
        "issue": "HTTP 500",
        "tables": ["mangaHIall", "mangaDRPall"],
    },
    "DV": {
        "name": "Cannon vs ASPCAP",
        "issue": "cannonStar",
        "tables": ["cannonStar"],
        "alt_tables": [],
        "skip_reason": "cannonStar may not be in DR18"
    },
    "DZ": {
        "name": "Potassium Anomaly",
        "issue": "No K column",
        "tables": ["apogeeStar"],
        "columns": ["k_fe", "K_FE"],
    },
}

def diagnose_all():
    """Run diagnostics on all skipped tests"""
    print("=" * 70)
    print("SDSS TEST FAILURE DIAGNOSTICS")
    print("=" * 70)
    
    # First, check which common tables exist
    common_tables = [
        "galSpecLine", "galSpecInfo", "galSpecExtra", "galSpecIndx",
        "SpecObjAll", "SpecPhotoAll", "PhotoObjAll", "Galaxy",
        "apogeeStar", "aspcapStar", "allStar",
        "mangaDRPall", "mangaDAPall", "mangaFirefly", "mangaHIall",
        "zooSpec", "Zoo1Spec", "zoo2MainSpecz",
        "QsoCatalogAll", "emissionLinesPort", "emissionLinesDR8",
        "redMaPPer", "Neighbors", "marvelsStar", "cannonStar"
    ]
    
    print("\n1. CHECKING TABLE AVAILABILITY")
    print("-" * 50)
    
    table_status = {}
    for table in common_tables:
        exists = check_table_exists(table)
        table_status[table] = exists
        status = "✓" if exists else "✗"
        print(f"  {status} {table}")
        time.sleep(0.3)  # Rate limit
    
    # Categorize tests
    fixable = []
    unfixable = []
    needs_alt = []
    
    print("\n2. CATEGORIZING SKIPPED TESTS")
    print("-" * 50)
    
    for code, info in SKIPPED_TESTS.items():
        if info.get("skip_reason"):
            unfixable.append((code, info["name"], info["skip_reason"]))
            continue
            
        # Check if main tables exist
        main_tables_exist = all(table_status.get(t, False) for t in info.get("tables", []))
        
        if main_tables_exist:
            fixable.append((code, info["name"], "Tables exist - query issue"))
        elif info.get("alt_tables"):
            alt_exist = [t for t in info["alt_tables"] if table_status.get(t, False)]
            if alt_exist:
                needs_alt.append((code, info["name"], f"Use alt tables: {alt_exist}"))
            else:
                unfixable.append((code, info["name"], "No available tables"))
        else:
            unfixable.append((code, info["name"], "Missing tables, no alternatives"))
    
    print("\n  FIXABLE (tables exist, likely query syntax issue):")
    for code, name, reason in fixable:
        print(f"    {code}: {name} - {reason}")
    
    print("\n  NEEDS ALTERNATIVE TABLES:")
    for code, name, reason in needs_alt:
        print(f"    {code}: {name} - {reason}")
    
    print("\n  UNFIXABLE (data not in SDSS DR18):")
    for code, name, reason in unfixable:
        print(f"    {code}: {name} - {reason}")
    
    # Summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"  Fixable (query issues): {len(fixable)}")
    print(f"  Needs alt tables: {len(needs_alt)}")
    print(f"  Unfixable (missing data): {len(unfixable)}")
    print(f"  Total skipped: {len(SKIPPED_TESTS)}")
    
    # Save results
    results = {
        "table_status": table_status,
        "fixable": [{"code": c, "name": n, "reason": r} for c, n, r in fixable],
        "needs_alt": [{"code": c, "name": n, "reason": r} for c, n, r in needs_alt],
        "unfixable": [{"code": c, "name": n, "reason": r} for c, n, r in unfixable],
    }
    
    out_path = os.path.join(RESULTS_DIR, 'sdss_failure_diagnosis.json')
    with open(out_path, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\nResults saved to {out_path}")
    
    return results

if __name__ == "__main__":
    diagnose_all()
