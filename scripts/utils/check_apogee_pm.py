
import requests
import pandas as pd

SDSS_URL = "https://skyserver.sdss.org/dr18/SkyServerWS/SearchTools/SqlSearch"

def check_pm():
    sql = "SELECT TOP 1 * FROM apogeeStar"
    print(f"Checking apogeeStar columns for PM/RV...")
    try:
        response = requests.get(SDSS_URL, params={"cmd": sql, "format": "json"}, timeout=30)
        if response.status_code == 200:
            data = response.json()
            if "Rows" in data[0]:
                df = pd.DataFrame(data[0]["Rows"])
                cols = sorted(df.columns.tolist())
                
                kin_cols = [c for c in cols if 'pm' in c.lower() or 'vhelio' in c.lower()]
                print(f"Kinematic Columns: {kin_cols}")
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    check_pm()
