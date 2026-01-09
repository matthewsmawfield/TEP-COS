
import requests
import pandas as pd

SDSS_URL = "https://skyserver.sdss.org/dr18/SkyServerWS/SearchTools/SqlSearch"

def check_cols():
    sql = "SELECT TOP 1 * FROM eFEDS_Main_speccomp"
    print(f"Checking eFEDS_Main_speccomp columns...")
    try:
        response = requests.get(SDSS_URL, params={"cmd": sql, "format": "json"}, timeout=30)
        if response.status_code == 200:
            data = response.json()
            if "Rows" in data[0]:
                df = pd.DataFrame(data[0]["Rows"])
                cols = sorted(df.columns.tolist())
                
                # Search for Flux/Extent
                flux_cols = [c for c in cols if 'flux' in c.lower()]
                ext_cols = [c for c in cols if 'ext' in c.lower()]
                print(f"Flux Cols: {flux_cols}")
                print(f"Extent Cols: {ext_cols}")
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    check_cols()
