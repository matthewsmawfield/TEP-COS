
import requests
import pandas as pd

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
            if "Rows" in data[0]:
                df = pd.DataFrame(data[0]["Rows"])
                cols = df.columns.tolist()
                print(f"  Columns: {cols}")
                
                # Check for specific columns of interest
                if 'Firefly' in table_name:
                    age_cols = [c for c in cols if 'AGE' in c]
                    print(f"  Age Columns: {age_cols}")
                if 'DAP' in table_name:
                    sigma_cols = [c for c in cols if 'sigma' in c.lower()]
                    print(f"  Sigma Columns: {sigma_cols}")
                return
        print(f"  Failed: {response.status_code}")
    except Exception as e:
        print(f"  Error: {e}")

if __name__ == "__main__":
    check_table("mangaFireflyGlobal") # Often the table name is slightly different
    check_table("mangaFirefly")       # checking base name just in case
    check_table("mangaDAPall")
