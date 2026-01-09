
import requests
import pandas as pd

SDSS_URL = "https://skyserver.sdss.org/dr18/SkyServerWS/SearchTools/SqlSearch"

def check_table(table_name):
    sql = f"SELECT TOP 1 * FROM {table_name}"
    print(f"Checking {table_name} columns...")
    try:
        response = requests.get(SDSS_URL, params={"cmd": sql, "format": "json"}, timeout=30)
        if response.status_code == 200:
            data = response.json()
            if "Rows" in data[0]:
                df = pd.DataFrame(data[0]["Rows"])
                cols = sorted(df.columns.tolist())
                print(f"Columns in {table_name}: {cols}")
                return cols
    except Exception as e:
        print(f"Error checking {table_name}: {e}")
    return []

def main():
    # Check aspcapStar for coordinates
    aspcap_cols = check_table("aspcapStar")
    if 'glon' in aspcap_cols or 'l' in aspcap_cols:
        print("FOUND COORDINATES IN ASPCAPSTAR")
    
    # Check spiders_quasar
    check_table("spiders_quasar")
    
    # Check for lambda in manga tables?
    # We can't search all tables easily via SQL without schema access, 
    # but we can check if a specific kinematic table exists
    check_table("mangaDAPall") # Double check hidden cols?
    
if __name__ == "__main__":
    main()
