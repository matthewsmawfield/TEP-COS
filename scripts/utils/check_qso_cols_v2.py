
import requests
import pandas as pd

SDSS_URL = "https://skyserver.sdss.org/dr18/SkyServerWS/SearchTools/SqlSearch"

def check_qso_cols():
    sql = "SELECT TOP 1 * FROM mos_sdss_dr16_qso"
    print(f"Checking mos_sdss_dr16_qso columns...")
    try:
        response = requests.get(SDSS_URL, params={"cmd": sql, "format": "json"}, timeout=30)
        if response.status_code == 200:
            data = response.json()
            if "Rows" in data[0]:
                df = pd.DataFrame(data[0]["Rows"])
                cols = sorted(df.columns.tolist())
                # Filter to avoid printing too much
                sigma_cols = [c for c in cols if 'sigma' in c.lower()]
                print(f"Sigma related: {sigma_cols}")
                bh_cols = [c for c in cols if 'bh' in c.lower()]
                print(f"BH related: {bh_cols}")
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    check_qso_cols()
