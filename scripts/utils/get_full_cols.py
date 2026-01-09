import requests
import pandas as pd

SDSS_URL = "https://skyserver.sdss.org/dr18/SkyServerWS/SearchTools/SqlSearch"

def get_columns(table):
    print(f"\nGetting columns for {table}...")
    sql = f"SELECT TOP 1 * FROM {table}"
    try:
        response = requests.get(SDSS_URL, params={"cmd": sql, "format": "json"}, timeout=30)
        if response.status_code == 200:
            data = response.json()
            if len(data) > 0 and "Rows" in data[0]:
                df = pd.DataFrame(data[0]["Rows"])
                cols = sorted(df.columns.tolist())
                print(f"Columns in {table}:")
                # Print in a way that doesn't get truncated easily
                for i in range(0, len(cols), 5):
                    print("  " + ", ".join(cols[i:i+5]))
            else:
                print("No rows.")
        else:
            print(f"Error {response.status_code}")
    except Exception as e:
        print(f"Exception: {e}")

get_columns("emissionLinesPort")
get_columns("SpecPhotoAll")
