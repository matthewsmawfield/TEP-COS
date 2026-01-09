
import requests
import pandas as pd

SDSS_URL = "https://skyserver.sdss.org/dr18/SkyServerWS/SearchTools/SqlSearch"

def check_tables():
    # Try more candidates including legacy and variations
    candidates = [
        "redMapper",
        "mapper",
        "cluster",
        "members",
        "redmapper_cluster",
        "redmapper_members",
        "run_redmapper_v6.3_lgt20_catalog",
        "run_redmapper_v6.3_members",
        "spiders_cluster_redmapper_v5.2_members",
        "spiders_cluster_redmapper_v5.2_catalog"
    ]
    
    # Also try to query a known schema if possible
    # Note: run_redmapper_v6_3_lgt20_catalog is standard DR17. 
    # Maybe dot vs underscore? v6.3 vs v6_3
    
    candidates_v2 = [
        "run_redmapper_v6_3_lgt20_catalog",
        "run_redmapper_v6_3_members",
        "run_redmapper_v5_10_lgt20_catalog",
        "redmapper_dr8_cluster",
        "redmapper_dr8_member"
    ]
    
    for table in candidates + candidates_v2:
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
