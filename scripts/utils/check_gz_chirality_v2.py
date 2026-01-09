
import requests
import pandas as pd

SDSS_URL = "https://skyserver.sdss.org/dr18/SkyServerWS/SearchTools/SqlSearch"

def check_tables():
    candidates = ["zoo2MainSpecz", "zooSpec"]
    
    for table in candidates:
        sql = f"SELECT TOP 1 * FROM {table}"
        print(f"Checking {table}...")
        try:
            response = requests.get(SDSS_URL, params={"cmd": sql, "format": "json"}, timeout=30)
            if response.status_code == 200:
                data = response.json()
                if "Rows" in data[0]:
                    df = pd.DataFrame(data[0]["Rows"])
                    cols = sorted(df.columns.tolist())
                    
                    # Search for likely chirality columns
                    targets = [c for c in cols if any(x in c.lower() for x in ['cw', 'acw', 'wind', 'dir', 'spiral'])]
                    print(f"  Candidates in {table}: {targets}")
        except Exception as e:
            print(f"  Error checking {table}: {e}")

if __name__ == "__main__":
    check_tables()
