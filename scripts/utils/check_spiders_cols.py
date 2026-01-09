
import requests
import pandas as pd

SDSS_URL = "https://skyserver.sdss.org/dr18/SkyServerWS/SearchTools/SqlSearch"

def check_cols():
    # Check spiders_quasar
    sql = "SELECT TOP 1 * FROM spiders_quasar"
    print(f"Checking spiders_quasar...")
    try:
        response = requests.get(SDSS_URL, params={"cmd": sql, "format": "json"}, timeout=30)
        if response.status_code == 200:
            data = response.json()
            if "Rows" in data[0]:
                df = pd.DataFrame(data[0]["Rows"])
                cols = df.columns.tolist()
                print(f"Columns in spiders_quasar: {sorted(cols)}")
                
                # Check for BH mass keywords
                keywords = ['bh', 'mass', 'sigma', 'fwhm', 'lum', 'bol']
                interesting = [c for c in cols if any(k in c.lower() for k in keywords)]
                print(f"Interesting columns: {sorted(interesting)}")
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    check_cols()
