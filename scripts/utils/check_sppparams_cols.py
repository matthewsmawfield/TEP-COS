
import requests
import pandas as pd

SDSS_URL = "https://skyserver.sdss.org/dr18/SkyServerWS/SearchTools/SqlSearch"

def check_cols():
    sql = "SELECT TOP 1 * FROM sppParams"
    print(f"Checking sppParams columns...")
    try:
        response = requests.get(SDSS_URL, params={"cmd": sql, "format": "json"}, timeout=30)
        if response.status_code == 200:
            data = response.json()
            if "Rows" in data[0]:
                df = pd.DataFrame(data[0]["Rows"])
                cols = sorted(df.columns.tolist())
                
                # Search for Teff, Logg, Dist
                teff_cols = [c for c in cols if 'teff' in c.lower()]
                logg_cols = [c for c in cols if 'logg' in c.lower()]
                dist_cols = [c for c in cols if 'dist' in c.lower()]
                
                print(f"Teff Cols: {teff_cols}")
                print(f"Logg Cols: {logg_cols}")
                print(f"Dist Cols: {dist_cols}")
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    check_cols()
