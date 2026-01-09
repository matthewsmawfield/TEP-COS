
import requests
import pandas as pd

SDSS_URL = "https://skyserver.sdss.org/dr18/SkyServerWS/SearchTools/SqlSearch"

def check_cols():
    # Check mos_sdss_dr16_qso
    sql = "SELECT TOP 1 * FROM mos_sdss_dr16_qso"
    print(f"Checking mos_sdss_dr16_qso...")
    try:
        response = requests.get(SDSS_URL, params={"cmd": sql, "format": "json"}, timeout=30)
        if response.status_code == 200:
            data = response.json()
            if "Rows" in data[0]:
                df = pd.DataFrame(data[0]["Rows"])
                cols = df.columns.tolist()
                
                # Filter for interesting columns
                keywords = ['mass', 'bh', 'sigma', 'var', 'spec', 'obj', 'z']
                interesting = [c for c in cols if any(k in c.lower() for k in keywords)]
                print(f"Interesting columns in mos_sdss_dr16_qso: {sorted(interesting)}")
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    check_cols()
