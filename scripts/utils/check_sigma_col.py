
import requests
import pandas as pd

SDSS_URL = "https://skyserver.sdss.org/dr18/SkyServerWS/SearchTools/SqlSearch"

def check_sigma():
    # Query to get one row to inspect columns
    sql = "SELECT TOP 1 * FROM emissionLinesPort"
    try:
        response = requests.get(
            SDSS_URL,
            params={"cmd": sql, "format": "json"},
            timeout=30
        )
        if response.status_code == 200:
            data = response.json()
            if "Rows" in data[0]:
                df = pd.DataFrame(data[0]["Rows"])
                cols = df.columns.tolist()
                sigma_cols = [c for c in cols if "sigma" in c.lower()]
                print(f"Sigma columns in emissionLinesPort: {sigma_cols}")
                return
        print(f"Failed: {response.status_code} {response.text}")
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    check_sigma()
