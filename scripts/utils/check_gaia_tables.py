
import requests
import pandas as pd

SDSS_URL = "https://skyserver.sdss.org/dr18/SkyServerWS/SearchTools/SqlSearch"

def check_tables():
    # In SDSS CasJobs, we can query information_schema or just try to select from likely names
    # But SqlSearch might not let us query info schema easily.
    # Let's try to query a known table or just check system tables if possible.
    # Actually, the best way often is to query a known table name and see if it works.
    
    # Common names in SDSS contexts:
    candidates = [
        "mos_gaia_dr2_source",
        "mos_geometric_distances_gaia_dr2",
        "gaia_dr2_source",
        "gaia_dr3_source",
        "gaiadr2.gaia_source", 
    ]
    
    for table in candidates:
        sql = f"SELECT TOP 1 * FROM {table}"
        print(f"Checking {table}...")
        try:
            response = requests.get(SDSS_URL, params={"cmd": sql, "format": "json"}, timeout=30)
            if response.status_code == 200:
                print(f"  [FOUND] {table}")
            else:
                print(f"  [MISSING] {table} (Status {response.status_code})")
        except Exception as e:
            print(f"  Error checking {table}: {e}")

if __name__ == "__main__":
    check_tables()
