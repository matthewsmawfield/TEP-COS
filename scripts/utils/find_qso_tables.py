
import requests
import pandas as pd

SDSS_URL = "https://skyserver.sdss.org/dr18/SkyServerWS/SearchTools/SqlSearch"

def check_table(table_name):
    print(f"Checking {table_name}...")
    sql = f"SELECT TOP 1 * FROM {table_name}"
    try:
        response = requests.get(SDSS_URL, params={"cmd": sql, "format": "json"}, timeout=30)
        if response.status_code == 200:
            print(f"  EXISTS: {table_name}")
            data = response.json()
            if "Rows" in data[0]:
                df = pd.DataFrame(data[0]["Rows"])
                print(f"  Columns: {sorted(df.columns.tolist())}")
            return True
        else:
            print(f"  NOT FOUND or Error: {response.status_code}")
            return False
    except Exception as e:
        print(f"  Error: {e}")
        return False

def main():
    tables_to_check = [
        "qsoVarIndexStripe82",
        "MacLeod2010",
        "Stripe82Score",
        "QsoVar",
        "DRW",
        "qso_drw"
    ]
    
    for t in tables_to_check:
        check_table(t)

if __name__ == "__main__":
    main()
