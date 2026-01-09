
import requests
import pandas as pd

SDSS_URL = "https://skyserver.sdss.org/dr18/SkyServerWS/SearchTools/SqlSearch"

def check_table_patterns(table_name, patterns):
    print(f"\nChecking table: {table_name}")
    sql = f"SELECT TOP 1 * FROM {table_name}"
    try:
        response = requests.get(
            SDSS_URL,
            params={"cmd": sql, "format": "json"},
            timeout=30
        )
        if response.status_code == 200:
            data = response.json()
            if "Rows" in data[0]:
                df = pd.DataFrame(data[0]["Rows"])
                cols = df.columns.tolist()
                print(f"  Total columns: {len(cols)}")
                
                for p in patterns:
                    matches = [c for c in cols if p.lower() in c.lower()]
                    print(f"  Matches for '{p}': {matches}")
                
                return
        print(f"  Failed: {response.status_code}")
    except Exception as e:
        print(f"  Error: {e}")

if __name__ == "__main__":
    check_table_patterns("mangaFirefly_miles", ["age", "plate", "ifu", "z", "mass"])
    check_table_patterns("mangaDAPall", ["sigma", "plate", "ifu", "z"])
