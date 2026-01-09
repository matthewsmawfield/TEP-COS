
import requests
import pandas as pd

SDSS_URL = "https://skyserver.sdss.org/dr18/SkyServerWS/SearchTools/SqlSearch"

def check_tables():
    candidates = [
        "zoo2MainSpecz",
        "zoo2MainPhotoz",
        "mangaGalaxyZoo",
        "zooSpec", # GZ1
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
                    
                    # Check chirality
                    chirality = [c for c in cols if 'clock' in c.lower()]
                    print(f"    Chirality Cols: {chirality}")
            else:
                print(f"  [MISSING] {table} (Status {response.status_code})")
        except Exception as e:
            print(f"  Error checking {table}: {e}")

if __name__ == "__main__":
    check_tables()
