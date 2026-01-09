
import requests
import pandas as pd

SDSS_URL = "https://skyserver.sdss.org/dr18/SkyServerWS/SearchTools/SqlSearch"

def check_cols():
    sql = "SELECT TOP 1 * FROM mangaDAPall"
    print(f"Checking mangaDAPall columns for properties...")
    try:
        response = requests.get(SDSS_URL, params={"cmd": sql, "format": "json"}, timeout=30)
        if response.status_code == 200:
            data = response.json()
            if "Rows" in data[0]:
                df = pd.DataFrame(data[0]["Rows"])
                cols = sorted(df.columns.tolist())
                
                props = [c for c in cols if 'lambda' in c.lower() or 'eps' in c.lower() or 'ellip' in c.lower()]
                print(f"Property Columns: {props}")
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    check_cols()
