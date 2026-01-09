
import requests
import pandas as pd

SDSS_URL = "https://skyserver.sdss.org/dr18/SkyServerWS/SearchTools/SqlSearch"

def check_spp_lines():
    print("Checking sppLines table for Lithium...")
    
    # Check table cols
    sql = "SELECT TOP 1 * FROM sppLines"
    
    try:
        response = requests.get(SDSS_URL, params={"cmd": sql, "format": "json"}, timeout=30)
        if response.status_code == 200:
            data = response.json()
            if "Rows" in data[0]:
                df = pd.DataFrame(data[0]["Rows"])
                cols = sorted(df.columns.tolist())
                print("sppLines table found.")
                
                # Check for Li
                # li_cols = [c for c in cols if 'li' in c.lower()]
                # print(f"Lithium columns: {li_cols}")
                print(f"All columns: {cols}")
            else:
                print("sppLines table empty.")
        else:
            print(f"Error: Status {response.status_code}")
            
    except Exception as e:
        print(f"Exception: {e}")

if __name__ == "__main__":
    check_spp_lines()
