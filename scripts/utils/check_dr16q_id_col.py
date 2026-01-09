
import requests
import pandas as pd

SDSS_URL = "https://skyserver.sdss.org/dr18/SkyServerWS/SearchTools/SqlSearch"

def check_cols():
    sql = "SELECT TOP 1 * FROM mos_sdss_dr16_qso"
    print(f"Checking mos_sdss_dr16_qso columns...")
    try:
        response = requests.get(SDSS_URL, params={"cmd": sql, "format": "json"}, timeout=30)
        if response.status_code == 200:
            data = response.json()
            if "Rows" in data[0]:
                df = pd.DataFrame(data[0]["Rows"])
                cols = sorted(df.columns.tolist())
                
                # Search for ID columns
                id_cols = [c for c in cols if 'id' in c.lower() or 'spec' in c.lower()]
                print(f"ID Columns: {id_cols}")
        else:
            print(f"Status: {response.status_code}")
            print(response.text[:200])
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    check_cols()
