
import requests
import pandas as pd

SDSS_URL = "https://skyserver.sdss.org/dr18/SkyServerWS/SearchTools/SqlSearch"

def check_snr():
    sql = "SELECT TOP 1 * FROM aspcapStar"
    print(f"Checking aspcapStar columns for snr...")
    try:
        response = requests.get(SDSS_URL, params={"cmd": sql, "format": "json"}, timeout=30)
        if response.status_code == 200:
            data = response.json()
            if "Rows" in data[0]:
                df = pd.DataFrame(data[0]["Rows"])
                cols = df.columns.tolist()
                
                snr_cols = [c for c in cols if 'snr' in c.lower()]
                print(f"SNR-related columns: {snr_cols}")
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    check_snr()
