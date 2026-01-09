
import requests
import pandas as pd

SDSS_URL = "https://skyserver.sdss.org/dr18/SkyServerWS/SearchTools/SqlSearch"

def check_link():
    print("Checking link between ebossMCPM and SpecObjAll...")
    
    sql = """
    SELECT TOP 10
        e.PLATE, e.MJD, e.FIBERID,
        s.specObjID
    FROM ebossMCPM e
    JOIN SpecObjAll s ON e.PLATE = s.PLATE AND e.MJD = s.MJD AND e.FIBERID = s.FIBERID
    """
    
    try:
        response = requests.get(SDSS_URL, params={"cmd": sql, "format": "json"}, timeout=30)
        if response.status_code == 200:
            data = response.json()
            if "Rows" in data[0]:
                df = pd.DataFrame(data[0]["Rows"])
                print(f"  Link successful! {len(df)} rows.")
                print(df)
            else:
                print("  No rows returned (Empty join).")
        else:
            print(f"  Error: Status {response.status_code}")
            print(response.text[:200])
    except Exception as e:
        print(f"  Exception: {e}")

if __name__ == "__main__":
    check_link()
