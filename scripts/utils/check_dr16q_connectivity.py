
import requests
import pandas as pd

SDSS_URL = "https://skyserver.sdss.org/dr18/SkyServerWS/SearchTools/SqlSearch"

def check_connectivity():
    print(f"Checking mos_sdss_dr16_qso connectivity...")
    sql = "SELECT TOP 10 specObjID, ra, dec, z FROM mos_sdss_dr16_qso"
    
    try:
        response = requests.get(SDSS_URL, params={"cmd": sql, "format": "json"}, timeout=30)
        if response.status_code == 200:
            print(f"  [SUCCESS] Status 200")
            data = response.json()
            if "Rows" in data[0]:
                df = pd.DataFrame(data[0]["Rows"])
                print(f"  Retrieved {len(df)} rows.")
                print(df.head())
            else:
                print("  No rows returned.")
        else:
            print(f"  [FAILURE] Status {response.status_code}")
            print(response.text[:200])
    except Exception as e:
        print(f"  Error: {e}")

if __name__ == "__main__":
    check_connectivity()
