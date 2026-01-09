
import requests
import pandas as pd
import sys

SDSS_URL = "https://skyserver.sdss.org/dr18/SkyServerWS/SearchTools/SqlSearch"

def list_cols():
    sql = "SELECT TOP 1 * FROM mangaDAPall"
    print(f"Querying mangaDAPall columns...")
    try:
        response = requests.get(SDSS_URL, params={"cmd": sql, "format": "json"}, timeout=30)
        if response.status_code == 200:
            data = response.json()
            if "Rows" in data[0]:
                df = pd.DataFrame(data[0]["Rows"])
                cols = sorted(df.columns.tolist())
                
                with open("manga_cols.txt", "w") as f:
                    for c in cols:
                        f.write(c + "\n")
                print("Columns written to manga_cols.txt")
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    list_cols()
