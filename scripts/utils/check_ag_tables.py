
import requests
import pandas as pd

SDSS_URL = "https://skyserver.sdss.org/dr18/SkyServerWS/SearchTools/SqlSearch"

def check_gaia_tables():
    sql = "SELECT table_name FROM information_schema.tables WHERE table_name LIKE '%gaia%'"
    print(f"Searching for Gaia tables...")
    try:
        response = requests.get(SDSS_URL, params={"cmd": sql, "format": "json"}, timeout=30)
        if response.status_code == 200:
            data = response.json()
            if "Rows" in data[0]:
                df = pd.DataFrame(data[0]["Rows"])
                print(df['table_name'].tolist())
            else:
                print("No Gaia tables found.")
        else:
            print(f"HTTP {response.status_code}")
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    check_gaia_tables()
