
import requests
import pandas as pd

SDSS_URL = "https://skyserver.sdss.org/dr18/SkyServerWS/SearchTools/SqlSearch"

def check_overlap():
    print("Checking for VALID overlap between galSpecIndx and sppParams...")
    
    # Check if we have valid parameters in the join
    # Use teffadop from sppParams
    sql = """
    SELECT TOP 20
        g.specObjID, s.teffadop, s.loggadop, g.lick_ca4227
    FROM galSpecIndx g
    JOIN sppParams s ON g.specObjID = s.specObjID
    WHERE s.teffadop > 0 AND s.teffadop < 10000
    """
    
    try:
        response = requests.get(SDSS_URL, params={"cmd": sql, "format": "json"}, timeout=30)
        if response.status_code == 200:
            data = response.json()
            if "Rows" in data[0]:
                df = pd.DataFrame(data[0]["Rows"])
                print(f"  Valid overlap found! {len(df)} rows.")
                print(df)
            else:
                print("  No valid overlap found (0 rows with valid Teff).")
        else:
            print(f"  Error: Status {response.status_code}")
            print(response.text[:200])
    except Exception as e:
        print(f"  Exception: {e}")

if __name__ == "__main__":
    check_overlap()
