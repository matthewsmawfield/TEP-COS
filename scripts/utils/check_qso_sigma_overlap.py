
import requests
import pandas as pd

SDSS_URL = "https://skyserver.sdss.org/dr18/SkyServerWS/SearchTools/SqlSearch"

def check_overlap():
    sql = """
    SELECT count(*) as count
    FROM spiders_quasar s
    JOIN emissionLinesPort p ON s.SpecObjID = p.specObjID
    WHERE p.sigmaStars > 0 AND s.logBHMA_hb > 0
    """
    print("Checking overlap between spiders_quasar and emissionLinesPort...")
    try:
        response = requests.get(SDSS_URL, params={"cmd": sql, "format": "json"}, timeout=30)
        if response.status_code == 200:
            data = response.json()
            if "Rows" in data[0]:
                count = data[0]["Rows"][0]['count']
                print(f"Overlap count: {count}")
        else:
            print(f"Failed: {response.status_code}")
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    check_overlap()
