
import requests
import pandas as pd

SDSS_URL = "https://skyserver.sdss.org/dr18/SkyServerWS/SearchTools/SqlSearch"

def check_join():
    # Check if we can join spiders_quasar and emissionLinesPort
    sql = """
    SELECT count(*) as count
    FROM spiders_quasar q
    JOIN emissionLinesPort e ON q.SPECOBJID = e.specObjID
    """
    print(f"Checking BG Join count...")
    try:
        response = requests.get(SDSS_URL, params={"cmd": sql, "format": "json"}, timeout=60)
        if response.status_code == 200:
            data = response.json()
            if "Rows" in data[0]:
                print(f"Join Count: {data[0]['Rows'][0]['count']}")
            else:
                print("No rows returned.")
        else:
            print(f"HTTP {response.status_code}")
            print(response.text)
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    check_join()
