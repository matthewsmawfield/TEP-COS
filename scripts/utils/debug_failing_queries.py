#!/usr/bin/env python3
"""Debug failing SDSS queries"""

import requests
import json

SDSS_URL = "https://skyserver.sdss.org/dr18/SkyServerWS/SearchTools/SqlSearch"

def test_query(name, sql):
    """Test a query and show detailed error"""
    print(f"\n=== {name} ===")
    print(f"SQL: {sql[:200]}...")
    try:
        response = requests.get(
            SDSS_URL,
            params={"cmd": sql, "format": "json"},
            timeout=60
        )
        print(f"Status: {response.status_code}")
        if response.status_code == 200:
            data = response.json()
            if data and len(data) > 0:
                if "Rows" in data[0]:
                    print(f"SUCCESS: {len(data[0]['Rows'])} rows")
                    return True
                else:
                    print(f"Response: {str(data)[:500]}")
        else:
            print(f"Error: {response.text[:500]}")
    except Exception as e:
        print(f"Exception: {e}")
    return False

# Test simplified queries for failing tests
queries = {
    "CE_simple": """
    SELECT TOP 10 el.nii_6584_flux, el.oiii_5007_flux
    FROM emissionLinesPort el
    WHERE el.nii_6584_flux > 0
    """,
    
    "CI_simple": """
    SELECT TOP 10 a.glon, a.glat, a.vhelio_avg
    FROM aspcapStar a
    WHERE a.aspcapflag = 0
    """,
    
    "CZ_simple": """
    SELECT TOP 10 d.emline_gflux_ha_6564
    FROM mangaDAPall d
    WHERE d.emline_gflux_ha_6564 > 0
    """,
    
    "DA_simple": """
    SELECT TOP 10 gi.bptclass, gl.oiii_5007_flux
    FROM galSpecInfo gi
    JOIN galSpecLine gl ON gi.specObjID = gl.specObjID
    WHERE gi.bptclass IN (3, 4, 5)
    """,
    
    "DC_simple": """
    SELECT TOP 10 n.objID, n.distance
    FROM Neighbors n
    WHERE n.distance > 10
    """,
    
    "DK_simple": """
    SELECT TOP 10 z.p_features_or_disk
    FROM zooSpec z
    """,
    
    "DI_simple": """
    SELECT TOP 10 d.stellar_vel, d.stellar_sigma
    FROM mangaDAPall d
    WHERE d.stellar_sigma > 0
    """,
    
    "DP_simple": """
    SELECT TOP 10 gx.tio2sdss
    FROM galSpecIndx gx
    WHERE gx.tio2sdss > 0
    """,
    
    "DR_simple": """
    SELECT TOP 10 a.vscatter
    FROM apogeeStar a
    WHERE a.vscatter > 0
    """,
}

print("Testing simplified queries to identify column issues...")
for name, sql in queries.items():
    test_query(name, sql)
