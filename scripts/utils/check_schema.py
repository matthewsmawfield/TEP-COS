#!/usr/bin/env python3
"""Check actual column names in SDSS DR18 tables"""

import requests

SDSS_URL = "https://skyserver.sdss.org/dr18/SkyServerWS/SearchTools/SqlSearch"

def get_columns(table):
    """Get column names for a table"""
    sql = f"SELECT TOP 1 * FROM {table}"
    try:
        response = requests.get(
            SDSS_URL,
            params={"cmd": sql, "format": "json"},
            timeout=30
        )
        if response.status_code == 200:
            data = response.json()
            if data and len(data) > 0 and "Rows" in data[0] and len(data[0]["Rows"]) > 0:
                return list(data[0]["Rows"][0].keys())
    except:
        pass
    return []

tables = [
    "emissionLinesPort",
    "aspcapStar", 
    "mangaDAPall",
    "galSpecInfo",
    "galSpecLine",
    "galSpecIndx",
    "zooSpec",
]

for table in tables:
    print(f"\n=== {table} ===")
    cols = get_columns(table)
    if cols:
        # Show relevant columns
        relevant = [c for c in cols if any(k in c.lower() for k in 
            ['nii', 'oiii', 'ha', 'hb', 'flux', 'glon', 'glat', 'vhelio', 'bpt', 
             'stellar', 'emline', 'tio', 'disk', 'features', 'sigma', 'disp'])]
        print(f"Relevant columns ({len(relevant)}):")
        for c in sorted(relevant)[:30]:
            print(f"  {c}")
        if not relevant:
            print(f"All columns: {cols[:20]}")
    else:
        print("  Could not fetch columns")
