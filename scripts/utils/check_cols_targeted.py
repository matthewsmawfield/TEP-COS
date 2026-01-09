import requests
import pandas as pd

SDSS_URL = "https://skyserver.sdss.org/dr18/SkyServerWS/SearchTools/SqlSearch"

def check_columns(table):
    print(f"\nChecking table: {table}")
    sql = f"SELECT TOP 1 * FROM {table}"
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
                cols = df.columns.tolist()
                print(f"Columns found: {len(cols)}")
                # Check for specific columns we care about
                targets = ['sigma_stars', 'sigma_stars_err', 'bestObjID', 'specObjID', 'objID', 'logMass', 'mstellar_median']
                for t in targets:
                    if t in cols:
                        print(f"  FOUND: {t}")
                    else:
                        pass # Don't print missing to keep output clean, unless we want to verify absence
                
                # Print all columns if list is short, or first few if long
                print(f"First 10 columns: {cols[:10]}")
            else:
                print(f"Empty result: {data}")
        else:
            print(f"HTTP {response.status_code}")
    except Exception as e:
        print(f"Error: {e}")

tables = [
    "emissionLinesPort",
    "SpecPhotoAll",
    "stellarMassFSPSGranWideDust",
    "stellarMassPCAWiscBC03",
    "stellarMassStarformingPort"
]

for t in tables:
    check_columns(t)
