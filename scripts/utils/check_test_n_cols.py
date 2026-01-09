
import requests
import pandas as pd

SDSS_URL = "https://skyserver.sdss.org/dr18/SkyServerWS/SearchTools/SqlSearch"

def check_columns():
    tables = ['spiders_quasar', 'qsoVarStripe', 'mos_sdss_dr16_qso']
    
    for table in tables:
        print(f"\nChecking {table}...")
        sql = f"SELECT TOP 1 * FROM {table}"
        try:
            response = requests.get(SDSS_URL, params={"cmd": sql, "format": "json"}, timeout=30)
            if response.status_code == 200:
                data = response.json()
                if "Rows" in data[0]:
                    df = pd.DataFrame(data[0]["Rows"])
                    cols = sorted(df.columns.tolist())
                    
                    # Filter for interesting columns to avoid massive output
                    interesting = []
                    keywords = ['id', 'obj', 'spec', 'ra', 'dec', 'mass', 'bh', 'log', 'var', 'tau', 'sigma']
                    for c in cols:
                        if any(k in c.lower() for k in keywords):
                            interesting.append(c)
                            
                    print(f"Key columns: {interesting}")
            else:
                print(f"Failed: {response.status_code}")
        except Exception as e:
            print(f"Error: {e}")

if __name__ == "__main__":
    check_columns()
