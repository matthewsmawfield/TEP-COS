
import requests
import pandas as pd

SDSS_URL = "https://skyserver.sdss.org/dr18/SkyServerWS/SearchTools/SqlSearch"

def check_tables():
    candidates = [
        "stellarMassFSPSGranWideDust",
        "stellarMassFSPSGranEarlyDust",
        "stellarMassStarformingPort",
        "stellarMassPCAWiscM11",
        "galSpecExtra",
        "galSpecInfo" # For sigma
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
                    # Check specific columns
                    if 'mass' in str(cols).lower():
                        mass_cols = [c for c in cols if 'mass' in c.lower()]
                        print(f"    Mass Cols: {mass_cols}")
                    if 'oh_p50' in cols:
                        print(f"    Has oh_p50: Yes")
            else:
                print(f"  [MISSING] {table} (Status {response.status_code})")
        except Exception as e:
            print(f"  Error checking {table}: {e}")

if __name__ == "__main__":
    check_tables()
