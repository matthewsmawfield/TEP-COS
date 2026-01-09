
import requests
import pandas as pd

SDSS_URL = "https://skyserver.sdss.org/dr18/SkyServerWS/SearchTools/SqlSearch"

def check_tables():
    candidates = [
        "eFEDS_Main_speccomp",
        "eFEDS_Main",
        "eFEDS",
        "spiders_cluster_redmapper_v5.2", # Already checked and missing, but good to double check context
        "spiders_cluster_v5.2",
        "rosat_2rxs",
        "xmm_cluster"
    ]
    
    for table in candidates:
        sql = f"SELECT TOP 1 * FROM {table}"
        print(f"Checking {table}...")
        try:
            response = requests.get(SDSS_URL, params={"cmd": sql, "format": "json"}, timeout=30)
            if response.status_code == 200:
                print(f"  [FOUND] {table}")
                data = response.json()
                if "Rows" in data[0]:
                    df = pd.DataFrame(data[0]["Rows"])
                    cols = sorted(df.columns.tolist())
                    print(f"    Columns: {cols[:5]}...")
            else:
                print(f"  [MISSING] {table} (Status {response.status_code})")
        except Exception as e:
            print(f"  Error checking {table}: {e}")

if __name__ == "__main__":
    check_tables()
