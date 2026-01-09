
import requests
import pandas as pd

SDSS_URL = "https://skyserver.sdss.org/dr18/SkyServerWS/SearchTools/SqlSearch"

def check_tables():
    candidates = [
        "gaia_dr2.gaia_source",
        "gaiadr2.gaia_source",
        "GaiaSource",
        "gaia_source",
        "gaia_dr3.gaia_source",
        "gaiadr3.gaia_source",
        "edr3.gaia_source",
        "gaiadr3", 
        "gaiadr2",
        "apogee_starhorse", # Just to confirm we can still query something
    ]
    
    for table in candidates:
        sql = f"SELECT TOP 1 * FROM {table}"
        print(f"Checking {table}...")
        try:
            response = requests.get(SDSS_URL, params={"cmd": sql, "format": "json"}, timeout=30)
            if response.status_code == 200:
                print(f"  [FOUND] {table}")
                # Print columns if found
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
