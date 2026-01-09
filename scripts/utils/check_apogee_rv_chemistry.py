
import requests
import pandas as pd

SDSS_URL = "https://skyserver.sdss.org/dr18/SkyServerWS/SearchTools/SqlSearch"

def check_data():
    sql = "SELECT TOP 20 rv_feh, rv_alpha, glon, glat, gaiaedr3_r_med_photogeo FROM apogeeStar WHERE rv_feh > -5"
    print(f"Checking apogeeStar chemistry...")
    try:
        response = requests.get(SDSS_URL, params={"cmd": sql, "format": "json"}, timeout=30)
        if response.status_code == 200:
            data = response.json()
            if "Rows" in data[0]:
                df = pd.DataFrame(data[0]["Rows"])
                print(df)
            else:
                print("No rows returned.")
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    check_data()
