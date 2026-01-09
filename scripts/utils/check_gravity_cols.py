
import requests
import pandas as pd

SDSS_URL = "https://skyserver.sdss.org/dr18/SkyServerWS/SearchTools/SqlSearch"

def check_cols():
    # Check aspcapStar
    sql1 = "SELECT TOP 1 * FROM aspcapStar"
    print(f"Checking aspcapStar columns...")
    try:
        response = requests.get(SDSS_URL, params={"cmd": sql1, "format": "json"}, timeout=30)
        if response.status_code == 200:
            data = response.json()
            if "Rows" in data[0]:
                df = pd.DataFrame(data[0]["Rows"])
                cols = sorted(df.columns.tolist())
                logg_cols = [c for c in cols if 'logg' in c.lower()]
                print(f"  aspcapStar logg: {logg_cols}")
    except Exception as e:
        print(f"Error: {e}")

    # Check apogee_starhorse
    sql2 = "SELECT TOP 1 * FROM apogee_starhorse"
    print(f"Checking apogee_starhorse columns...")
    try:
        response = requests.get(SDSS_URL, params={"cmd": sql2, "format": "json"}, timeout=30)
        if response.status_code == 200:
            data = response.json()
            if "Rows" in data[0]:
                df = pd.DataFrame(data[0]["Rows"])
                cols = sorted(df.columns.tolist())
                logg_cols = [c for c in cols if 'logg' in c.lower()]
                print(f"  apogee_starhorse logg: {logg_cols}")
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    check_cols()
