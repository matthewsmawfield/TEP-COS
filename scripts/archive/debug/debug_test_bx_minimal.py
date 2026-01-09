
import requests
import pandas as pd

SDSS_URL = "https://skyserver.sdss.org/dr18/SkyServerWS/SearchTools/SqlSearch"

def query(sql, tag):
    print(f"--- {tag} ---")
    print(f"SQL: {sql}")
    try:
        response = requests.get(SDSS_URL, params={"cmd": sql, "format": "json"}, timeout=30)
        if response.status_code == 200:
            data = response.json()
            if "Rows" in data[0]:
                df = pd.DataFrame(data[0]["Rows"])
                print(f"Success! Got {len(df)} rows.")
                print(df.head(2))
            else:
                print("No rows returned.")
        else:
            print(f"Error: {response.status_code} {response.text[:200]}")
    except Exception as e:
        print(f"Exception: {e}")

def main():
    # 1. Very basic query
    query("SELECT TOP 10 apogee_id, vscatter FROM apogeeStar", "Basic apogeeStar")
    
    # 2. Filter nvisits
    query("SELECT TOP 10 apogee_id, vscatter FROM apogeeStar WHERE nvisits > 2", "Filter nvisits")
    
    # 3. Filter vscatter
    query("SELECT TOP 10 apogee_id, vscatter FROM apogeeStar WHERE vscatter > 0", "Filter vscatter")

if __name__ == "__main__":
    main()
