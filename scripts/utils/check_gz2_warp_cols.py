
import requests
import pandas as pd

SDSS_URL = "https://skyserver.sdss.org/dr18/SkyServerWS/SearchTools/SqlSearch"

def check_cols():
    sql = "SELECT TOP 1 * FROM MaNGA_GZ2"
    print(f"Checking MaNGA_GZ2 columns...")
    try:
        response = requests.get(SDSS_URL, params={"cmd": sql, "format": "json"}, timeout=30)
        if response.status_code == 200:
            data = response.json()
            if "Rows" in data[0]:
                df = pd.DataFrame(data[0]["Rows"])
                cols = sorted(df.columns.tolist())
                
                # Check for warp-related keywords
                keywords = ['warp', 'edge', 'odd', 'irr', 'disturb']
                found = []
                for k in keywords:
                    matches = [c for c in cols if k.lower() in c.lower()]
                    found.extend(matches)
                
                print(f"Potential Warp Columns: {found}")
                print(f"All Columns (first 20): {cols[:20]}")
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    check_cols()
