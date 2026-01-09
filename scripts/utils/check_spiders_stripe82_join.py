
import requests
import pandas as pd

SDSS_URL = "https://skyserver.sdss.org/dr18/SkyServerWS/SearchTools/SqlSearch"

def check_join():
    # Join spiders_quasar -> SpecObjAll -> qsoVarStripe
    # spiders_quasar.SPECOBJID = SpecObjAll.specObjID
    # SpecObjAll.bestObjID = qsoVarStripe.VAR_OBJID
    
    sql = """
    SELECT count(*) as count
    FROM spiders_quasar s
    JOIN SpecObjAll sp ON s.SPECOBJID = sp.specObjID
    JOIN qsoVarStripe v ON sp.bestObjID = v.VAR_OBJID
    WHERE s.width2_OIII5007 > 0
    """
    print("Checking join: spiders_quasar -> SpecObjAll -> qsoVarStripe...")
    try:
        response = requests.get(SDSS_URL, params={"cmd": sql, "format": "json"}, timeout=30)
        if response.status_code == 200:
            data = response.json()
            if "Rows" in data[0]:
                count = data[0]["Rows"][0]['count']
                print(f"Join count with OIII width: {count}")
        else:
            print(f"Failed: {response.status_code} {response.text}")
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    check_join()
