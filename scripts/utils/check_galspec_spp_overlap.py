
import requests
import pandas as pd

SDSS_URL = "https://skyserver.sdss.org/dr18/SkyServerWS/SearchTools/SqlSearch"

def check_overlap():
    print("Checking overlap between galSpecIndx and sppParams...")
    
    # Try to join on specObjID
    sql = """
    SELECT TOP 10
        g.specObjID, s.teffadop
    FROM galSpecIndx g
    JOIN sppParams s ON g.specObjID = s.specObjID
    """
    
    try:
        response = requests.get(SDSS_URL, params={"cmd": sql, "format": "json"}, timeout=30)
        if response.status_code == 200:
            data = response.json()
            if "Rows" in data[0]:
                df = pd.DataFrame(data[0]["Rows"])
                print(f"  Overlap found! {len(df)} rows.")
                print(df)
            else:
                print("  No overlap found (0 rows).")
        else:
            print(f"  Error: Status {response.status_code}")
            print(response.text[:200])
    except Exception as e:
        print(f"  Exception: {e}")

    # Check if we can find Lick indices elsewhere if no overlap
    print("\nChecking sppLines for Lick indices...")
    sql2 = "SELECT TOP 1 * FROM sppLines"
    try:
        response = requests.get(SDSS_URL, params={"cmd": sql2, "format": "json"}, timeout=30)
        if response.status_code == 200:
            data = response.json()
            if "Rows" in data[0]:
                df = pd.DataFrame(data[0]["Rows"])
                cols = sorted(df.columns.tolist())
                lick = [c for c in cols if 'lick' in c.lower() or 'ca' in c.lower()]
                print(f"  sppLines candidates: {lick[:10]}...")
            else:
                print("  sppLines empty.")
        else:
            print(f"  sppLines Error: Status {response.status_code}")
    except Exception as e:
        print(f"  Exception: {e}")

if __name__ == "__main__":
    check_overlap()
