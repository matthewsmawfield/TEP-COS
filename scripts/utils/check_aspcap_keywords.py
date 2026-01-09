
import requests
import pandas as pd

SDSS_URL = "https://skyserver.sdss.org/dr18/SkyServerWS/SearchTools/SqlSearch"

def check_keywords():
    sql = "SELECT TOP 1 * FROM aspcapStar"
    print(f"Checking aspcapStar columns for keywords...")
    try:
        response = requests.get(SDSS_URL, params={"cmd": sql, "format": "json"}, timeout=30)
        if response.status_code == 200:
            data = response.json()
            if "Rows" in data[0]:
                df = pd.DataFrame(data[0]["Rows"])
                cols = df.columns.tolist()
                
                keywords = ['dist', 'plx', 'gaia', 'glon', 'glat', 'ra', 'dec']
                found = []
                for k in keywords:
                    matches = [c for c in cols if k.lower() in c.lower()]
                    found.extend(matches)
                
                print(f"Found keywords: {found}")
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    check_keywords()
