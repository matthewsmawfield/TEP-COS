
import requests
import pandas as pd

SDSS_URL = "https://skyserver.sdss.org/dr18/SkyServerWS/SearchTools/SqlSearch"

def check_cols():
    # Check mangaFirefly_miles
    sql = "SELECT TOP 1 * FROM mangaFirefly_miles"
    print(f"Checking mangaFirefly_miles...")
    try:
        response = requests.get(SDSS_URL, params={"cmd": sql, "format": "json"}, timeout=30)
        if response.status_code == 200:
            data = response.json()
            if "Rows" in data[0]:
                df = pd.DataFrame(data[0]["Rows"])
                cols = sorted(df.columns.tolist())
                print(f"Columns in mangaFirefly_miles: {cols}")
                
                # Check for Gradient
                grads = [c for c in cols if "grad" in c.lower()]
                print(f"Gradient columns: {grads}")
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    check_cols()
