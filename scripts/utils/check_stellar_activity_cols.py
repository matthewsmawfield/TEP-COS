
import requests
import pandas as pd

SDSS_URL = "https://skyserver.sdss.org/dr18/SkyServerWS/SearchTools/SqlSearch"

def check_cols():
    # Check sppLines for Ca II
    print("Checking sppLines columns...")
    sql1 = "SELECT TOP 1 * FROM sppLines"
    try:
        response = requests.get(SDSS_URL, params={"cmd": sql1, "format": "json"}, timeout=30)
        if response.status_code == 200:
            data = response.json()
            if "Rows" in data[0]:
                df = pd.DataFrame(data[0]["Rows"])
                cols = sorted(df.columns.tolist())
                ca_cols = [c for c in cols if 'ca' in c.lower()]
                k_cols = [c for c in cols if 'k' in c.lower() and len(c) < 5] # Short names like K, CaK
                print(f"  sppLines Ca/K: {ca_cols} {k_cols}")
    except Exception as e:
        print(f"Error 1: {e}")

    # Check galSpecIndx class
    print("\nChecking galSpecIndx object classes...")
    sql2 = """
    SELECT TOP 20 s.class
    FROM galSpecIndx g
    JOIN SpecObjAll s ON g.specObjID = s.specObjID
    """
    try:
        response = requests.get(SDSS_URL, params={"cmd": sql2, "format": "json"}, timeout=30)
        if response.status_code == 200:
            data = response.json()
            if "Rows" in data[0]:
                df = pd.DataFrame(data[0]["Rows"])
                print(f"  Classes found: {df['class'].unique()}")
    except Exception as e:
        print(f"Error 2: {e}")

if __name__ == "__main__":
    check_cols()
