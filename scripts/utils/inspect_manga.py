
import requests
import pandas as pd
import json

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
            try:
                data = response.json()
                if "Rows" in data[0]:
                    df = pd.DataFrame(data[0]["Rows"])
                    cols = sorted(df.columns.tolist())
                    print(f"  Columns ({len(cols)}): {cols}")
                else:
                    print("  No rows returned.")
            except Exception as e:
                print(f"  JSON decode error or format error: {e}")
                print(response.text[:200])
        else:
            print(f"  Failed: {response.status_code}")
            print(response.text[:200])
    except Exception as e:
        print(f"  Error: {e}")

if __name__ == "__main__":
    check_table("mangaFirefly")
    check_table("mangaDAPall")
