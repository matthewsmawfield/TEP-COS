
import requests
import pandas as pd

SDSS_URL = "https://skyserver.sdss.org/dr18/SkyServerWS/SearchTools/SqlSearch"

def check_cols():
    # Check galSpecLine
    sql = "SELECT TOP 1 * FROM galSpecLine"
    print(f"Checking galSpecLine...")
    try:
        response = requests.get(SDSS_URL, params={"cmd": sql, "format": "json"}, timeout=30)
        if response.status_code == 200:
            data = response.json()
            if "Rows" in data[0]:
                df = pd.DataFrame(data[0]["Rows"])
                cols = df.columns.tolist()
                print(f"Columns in galSpecLine: {sorted(cols)}")
                
                # Check for Paschen
                paschen = [c for c in cols if "Pa" in c or "12818" in c or "18751" in c]
                print(f"Potential Paschen columns: {paschen}")
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    check_cols()
