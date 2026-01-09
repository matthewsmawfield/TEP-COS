
import requests
import pandas as pd

SDSS_URL = "https://skyserver.sdss.org/dr18/SkyServerWS/SearchTools/SqlSearch"

def check_cols():
    # Check galSpecExtra
    sql = "SELECT TOP 1 * FROM galSpecExtra"
    print(f"Checking galSpecExtra...")
    try:
        response = requests.get(SDSS_URL, params={"cmd": sql, "format": "json"}, timeout=30)
        if response.status_code == 200:
            data = response.json()
            if "Rows" in data[0]:
                df = pd.DataFrame(data[0]["Rows"])
                cols = df.columns.tolist()
                print(f"Columns in galSpecExtra: {sorted(cols)}")
                
                # Check for SFR
                sfr = [c for c in cols if "sfr" in c.lower()]
                print(f"SFR columns: {sfr}")
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    check_cols()
