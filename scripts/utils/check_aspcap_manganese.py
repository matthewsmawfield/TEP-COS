
import requests
import pandas as pd

SDSS_URL = "https://skyserver.sdss.org/dr18/SkyServerWS/SearchTools/SqlSearch"

def check_cols():
    sql = "SELECT TOP 1 * FROM aspcapStar"
    print(f"Checking aspcapStar columns...")
    try:
        response = requests.get(SDSS_URL, params={"cmd": sql, "format": "json"}, timeout=30)
        if response.status_code == 200:
            data = response.json()
            if "Rows" in data[0]:
                df = pd.DataFrame(data[0]["Rows"])
                cols = sorted(df.columns.tolist())
                
                # Search for Mn columns
                mn_cols = [c for c in cols if 'mn' in c.lower()]
                print(f"Mn Columns: {mn_cols}")
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    check_cols()
