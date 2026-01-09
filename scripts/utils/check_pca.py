import requests
import pandas as pd

SDSS_URL = "https://skyserver.sdss.org/dr18/SkyServerWS/SearchTools/SqlSearch"

def check_pca_columns():
    sql = "SELECT TOP 1 * FROM stellarMassPCAWiscBC03"
    try:
        response = requests.get(SDSS_URL, params={"cmd": sql, "format": "json"}, timeout=30)
        if response.status_code == 200:
            data = response.json()
            if len(data) > 0 and "Rows" in data[0]:
                df = pd.DataFrame(data[0]["Rows"])
                print("Columns in stellarMassPCAWiscBC03:")
                for c in df.columns:
                    print(f"  {c}")
            else:
                print("No data")
        else:
            print(f"Error: {response.status_code}")
    except Exception as e:
        print(f"Exception: {e}")

if __name__ == "__main__":
    check_pca_columns()
