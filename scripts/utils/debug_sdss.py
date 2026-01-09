import requests
import pandas as pd
import time

SDSS_URL = "https://skyserver.sdss.org/dr18/SkyServerWS/SearchTools/SqlSearch"

def query_sdss(sql):
    print(f"Executing: {sql}")
    try:
        response = requests.get(
            SDSS_URL,
            params={"cmd": sql, "format": "json"},
            timeout=30
        )
        print(f"Status Code: {response.status_code}")
        if response.status_code == 200:
            data = response.json()
            if len(data) > 0 and "Rows" in data[0]:
                df = pd.DataFrame(data[0]["Rows"])
                print(f"Rows returned: {len(df)}")
                print(df.head())
            else:
                print("No rows returned or unexpected format.")
                print(data)
        else:
            print(f"Error: {response.text}")
    except Exception as e:
        print(f"Exception: {e}")

if __name__ == "__main__":
    # Test 1: Very simple query
    query_sdss("SELECT TOP 5 p.objid, p.ra, p.dec FROM PhotoObjAll p")
    
    # Test 2: Join
    query_sdss("SELECT TOP 5 s.specObjID, s.z FROM SpecObjAll s JOIN PhotoObjAll p ON s.bestObjID = p.objID")
