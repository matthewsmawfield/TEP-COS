
import requests
import pandas as pd

SDSS_URL = "https://skyserver.sdss.org/dr18/SkyServerWS/SearchTools/SqlSearch"

def check_tables():
    print("Checking APOGEE columns...")
    
    tables = ["apogeeStar", "apogee_starhorse"]
    
    for t in tables:
        sql = f"SELECT TOP 1 * FROM {t}"
        try:
            response = requests.get(SDSS_URL, params={"cmd": sql, "format": "json"}, timeout=30)
            if response.status_code == 200:
                data = response.json()
                if "Rows" in data[0]:
                    print(f"\n{t}: Found.")
                    df = pd.DataFrame(data[0]["Rows"])
                    cols = df.columns.tolist()
                    print(f"  Columns: {cols}")
                    
                    # Check for PM or Gaia ID keywords
                    relevant = [c for c in cols if 'pm' in c.lower() or 'gaia' in c.lower() or 'source_id' in c.lower()]
                    print(f"  Relevant Columns: {relevant}")
                else:
                    print(f"\n{t}: Found but empty.")
            else:
                print(f"\n{t}: Failed {response.status_code}")
        except Exception as e:
            print(f"\n{t}: Error {e}")

if __name__ == "__main__":
    check_tables()
