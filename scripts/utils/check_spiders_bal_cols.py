
import requests
import pandas as pd

SDSS_URL = "https://skyserver.sdss.org/dr18/SkyServerWS/SearchTools/SqlSearch"

def check_cols():
    sql = "SELECT TOP 1 * FROM spiders_quasar"
    print(f"Checking spiders_quasar columns for BAL info...")
    try:
        response = requests.get(SDSS_URL, params={"cmd": sql, "format": "json"}, timeout=30)
        if response.status_code == 200:
            data = response.json()
            if "Rows" in data[0]:
                df = pd.DataFrame(data[0]["Rows"])
                cols = sorted(df.columns.tolist())
                
                # Search for BAL related columns
                bal_cols = [c for c in cols if 'bal' in c.lower() or 'civ' in c.lower() or 'absorp' in c.lower()]
                print(f"BAL Columns: {bal_cols}")
                
                # Check for plate/mjd
                link_cols = [c for c in cols if 'plate' in c.lower() or 'mjd' in c.lower()]
                print(f"Link Columns: {link_cols}")
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    check_cols()
