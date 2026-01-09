
import requests
import pandas as pd

SDSS_URL = "https://skyserver.sdss.org/dr18/SkyServerWS/SearchTools/SqlSearch"

def check_table(table_name):
    sql = f"SELECT TOP 1 * FROM {table_name}"
    print(f"Checking {table_name}...")
    try:
        response = requests.get(SDSS_URL, params={"cmd": sql, "format": "json"}, timeout=30)
        if response.status_code == 200:
            data = response.json()
            if "Rows" in data[0]:
                df = pd.DataFrame(data[0]["Rows"])
                cols = sorted(df.columns.tolist())
                print(f"Columns in {table_name} ({len(cols)}): {cols}")
            else:
                print(f"  [Empty] No rows returned for {table_name}.")
        else:
            print(f"  [Error] HTTP {response.status_code}")
    except Exception as e:
        print(f"  [Error] {e}")

def main():
    check_table("mos_sagitta")
    check_table("mangaFirefly")
    check_table("spiders_quasar")

if __name__ == "__main__":
    main()
