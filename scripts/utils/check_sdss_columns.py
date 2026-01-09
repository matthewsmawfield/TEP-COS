import requests
import pandas as pd
import sys

SDSS_URL = "https://skyserver.sdss.org/dr18/SkyServerWS/SearchTools/SqlSearch"

def check_table(table_name):
    print(f"\nChecking table: {table_name}")
    sql = f"SELECT TOP 1 * FROM {table_name}"
    try:
        response = requests.get(
            SDSS_URL,
            params={"cmd": sql, "format": "json"},
            timeout=30
        )
        if response.status_code == 200:
            data = response.json()
            if len(data) > 0 and "Rows" in data[0]:
                df = pd.DataFrame(data[0]["Rows"])
                print(f"Columns: {', '.join(df.columns)}")
                return True
            else:
                print(f"Empty result or error: {data}")
                return False
        else:
            print(f"HTTP {response.status_code}: {response.text}")
            return False
    except Exception as e:
        print(f"Exception: {e}")
        return False

tables = [
    "SpecPhotoAll",
    "emissionLinesPort"
]

for t in tables:
    check_table(t)
