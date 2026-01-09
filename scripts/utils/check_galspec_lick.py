
import requests
import pandas as pd

SDSS_URL = "https://skyserver.sdss.org/dr18/SkyServerWS/SearchTools/SqlSearch"

def check_tables():
    candidates = ["galSpecIndx", "sppParams"]
    
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
                    if table == "galSpecIndx":
                        lick_cols = [c for c in cols if 'lick' in c.lower() and 'ca' in c.lower()]
                        print(f"    Lick Ca Cols: {lick_cols}")
                    elif table == "sppParams":
                        param_cols = [c for c in cols if any(x in c.lower() for x in ['teff', 'logg', 'feh', 'elodi'])]
                        print(f"    Params: {param_cols[:10]}...")
            else:
                print(f"  [MISSING] {table} (Status {response.status_code})")
        except Exception as e:
            print(f"  Error checking {table}: {e}")

if __name__ == "__main__":
    check_tables()
