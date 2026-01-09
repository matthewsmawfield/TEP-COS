
import requests
import pandas as pd

SDSS_URL = "https://skyserver.sdss.org/dr18/SkyServerWS/SearchTools/SqlSearch"

def check_table(table_name):
    sql = f"SELECT TOP 1 * FROM {table_name}"
    print(f"Checking {table_name}...")
    try:
        response = requests.get(SDSS_URL, params={"cmd": sql, "format": "json"}, timeout=30)
        if response.status_code == 200:
            data = response.json()
            if "Rows" in data[0]:
                print(f"  [OK] Table {table_name} exists.")
                return True
            else:
                print(f"  [FAIL] Table {table_name} returned no rows (might exist but be empty).")
                return False
        else:
            print(f"  [FAIL] Table {table_name} not found or error (HTTP {response.status_code}).")
            return False
    except Exception as e:
        print(f"  [ERROR] {e}")
        return False

def main():
    tables_to_check = [
        "mangaPipe3D",
        "mangaFirefly",
        "ebossMCPM",
        "MaNGA_GZ2",
        "zoo2MainSpecz", # Alternative for GZ2
        "mos_sagitta",
        "mos_gaia_dr2_source",
        "mos_sdss_dr16_qso",
        "apogeeStar",
        "apogee_starhorse"
    ]
    
    results = {}
    for t in tables_to_check:
        results[t] = check_table(t)
        
    print("\nSummary:")
    for t, exists in results.items():
        print(f"{t}: {'Available' if exists else 'Missing'}")

if __name__ == "__main__":
    main()
