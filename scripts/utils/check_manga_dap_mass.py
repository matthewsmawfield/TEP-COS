
import requests
import pandas as pd

SDSS_URL = "https://skyserver.sdss.org/dr18/SkyServerWS/SearchTools/SqlSearch"

def check_cols():
    sql = "SELECT TOP 1 * FROM mangaDAPall"
    print(f"Checking mangaDAPall columns for 'mass'...")
    try:
        response = requests.get(SDSS_URL, params={"cmd": sql, "format": "json"}, timeout=30)
        if response.status_code == 200:
            data = response.json()
            if "Rows" in data[0]:
                df = pd.DataFrame(data[0]["Rows"])
                cols = sorted(df.columns.tolist())
                mass_cols = [c for c in cols if 'mass' in c.lower()]
                print(f"Mass Columns in DAPall: {mass_cols}")
                
                nsa_cols = [c for c in cols if 'nsa' in c.lower()]
                print(f"NSA Columns in DAPall: {nsa_cols}")
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    check_cols()
