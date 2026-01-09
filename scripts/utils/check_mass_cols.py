import requests
import pandas as pd

SDSS_URL = "https://skyserver.sdss.org/dr18/SkyServerWS/SearchTools/SqlSearch"

def check_columns(table):
    sql = f"SELECT TOP 1 * FROM {table}"
    try:
        response = requests.get(SDSS_URL, params={"cmd": sql, "format": "json"}, timeout=30)
        if response.status_code == 200:
            data = response.json()
            if len(data) > 0 and "Rows" in data[0]:
                df = pd.DataFrame(data[0]["Rows"])
                print(f"Columns in {table}:")
                for c in df.columns:
                    print(f"  {c}")
            else:
                print(f"No data for {table}")
        else:
            print(f"Error {response.status_code} for {table}")
    except Exception as e:
        print(f"Exception: {e}")

check_columns("stellarMassFSPSGranWideDust")
check_columns("stellarMassStarformingPort")
