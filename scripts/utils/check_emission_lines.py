
import requests
import pandas as pd

SDSS_URL = "https://skyserver.sdss.org/dr18/SkyServerWS/SearchTools/SqlSearch"

def check_cols():
    sql = "SELECT TOP 1 * FROM emissionLinesPort"
    print(f"Checking emissionLinesPort columns...")
    try:
        response = requests.get(SDSS_URL, params={"cmd": sql, "format": "json"}, timeout=30)
        if response.status_code == 200:
            data = response.json()
            if "Rows" in data[0]:
                df = pd.DataFrame(data[0]["Rows"])
                cols = sorted(df.columns.tolist())
                
                # Search for Flux columns
                flux_cols = [c for c in cols if 'flux' in c.lower()]
                print(f"Flux Columns: {flux_cols[:10]}... ({len(flux_cols)} total)")
                
                # Check specific lines
                targets = ['Flux_NII_6583', 'Flux_OII_3726', 'Flux_OIII_5006', 'Flux_Ha_6562', 'Flux_Hb_4861']
                found = [t for t in targets if any(c.lower() == t.lower() for c in cols)]
                print(f"Found Targets: {found}")
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    check_cols()
