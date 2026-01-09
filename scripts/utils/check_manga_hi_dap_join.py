
import requests
import pandas as pd

SDSS_URL = "https://skyserver.sdss.org/dr18/SkyServerWS/SearchTools/SqlSearch"

def check_join():
    # Check join and HI mass column
    sql = """
    SELECT TOP 10
        h.mangaid,
        h.logMHI, 
        d.stellar_sigma_1re
    FROM mangaHIall h
    JOIN mangaDAPall d ON h.mangaid = d.mangaid
    """
    print(f"Checking join...")
    try:
        response = requests.get(SDSS_URL, params={"cmd": sql, "format": "json"}, timeout=30)
        if response.status_code == 200:
            data = response.json()
            if "Rows" in data[0]:
                df = pd.DataFrame(data[0]["Rows"])
                print(f"Success! Columns: {df.columns.tolist()}")
                print(df.head())
            else:
                print("No rows returned.")
        else:
            print(f"HTTP {response.status_code}")
            print(response.text)
            
            # Fallback: check columns of mangaHIall again
            print("\nChecking mangaHIall columns...")
            sql2 = "SELECT TOP 1 * FROM mangaHIall"
            r2 = requests.get(SDSS_URL, params={"cmd": sql2, "format": "json"}, timeout=30)
            if r2.status_code == 200:
                print(r2.json()[0]["Rows"][0].keys())

    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    check_join()
