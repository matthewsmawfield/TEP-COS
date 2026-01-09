
import requests
import pandas as pd

SDSS_URL = "https://skyserver.sdss.org/dr18/SkyServerWS/SearchTools/SqlSearch"

def check_table(table_name):
    print(f"\nChecking {table_name}...")
    sql = f"SELECT TOP 1 * FROM {table_name}"
    try:
        response = requests.get(SDSS_URL, params={"cmd": sql, "format": "json"}, timeout=30)
        if response.status_code == 200:
            data = response.json()
            if "Rows" in data[0]:
                df = pd.DataFrame(data[0]["Rows"])
                cols = sorted(df.columns.tolist())
                
                # Filter for PA/Phi/Angle/Kin/Rot
                keywords = ['pa', 'phi', 'angle', 'rot', 'kin', 'axis', 'vel']
                interesting = [c for c in cols if any(k in c.lower() for k in keywords)]
                print(f"Interesting columns in {table_name}: {interesting}")
        else:
            print(f"HTTP {response.status_code}")
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    check_table("mangaDrpAll")
    check_table("mangaPipe3D")
