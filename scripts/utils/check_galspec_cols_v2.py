
import requests
import pandas as pd

SDSS_URL = "https://skyserver.sdss.org/dr18/SkyServerWS/SearchTools/SqlSearch"

def check_cols(table):
    sql = f"SELECT TOP 1 * FROM {table}"
    print(f"Checking {table} columns...")
    try:
        response = requests.get(SDSS_URL, params={"cmd": sql, "format": "json"}, timeout=30)
        if response.status_code == 200:
            data = response.json()
            if "Rows" in data[0]:
                df = pd.DataFrame(data[0]["Rows"])
                cols = sorted(df.columns.tolist())
                # Search for sigma/dispersion
                sigma_cols = [c for c in cols if 'sigma' in c.lower() or 'disp' in c.lower() or 'vel' in c.lower()]
                print(f"  Sigma/Vel related: {sigma_cols}")
    except Exception as e:
        print(f"  Error: {e}")

if __name__ == "__main__":
    check_cols("galSpecExtra")
    check_cols("galSpecInfo")
