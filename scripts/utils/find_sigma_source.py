
import requests
import pandas as pd

SDSS_URL = "https://skyserver.sdss.org/dr18/SkyServerWS/SearchTools/SqlSearch"

def check_table_columns(table_name):
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
                
                # Check for sigma/dispersion related columns
                sigma_cols = [c for c in cols if "sigma" in c.lower() or "disp" in c.lower() or "vel" in c.lower()]
                print(f"  Potential velocity dispersion columns: {sigma_cols}")
                
                # Check for objID related columns
                obj_cols = [c for c in cols if "objid" in c.lower()]
                print(f"  ObjID columns: {obj_cols}")
                return
        print(f"  Failed: {response.status_code}")
    except Exception as e:
        print(f"  Error: {e}")

if __name__ == "__main__":
    check_table_columns("galSpecInfo")
    check_table_columns("galSpecExtra")
    check_table_columns("SpecObjAll")
    check_table_columns("emissionLinesPort")
    check_table_columns("SpecPhotoAll")
