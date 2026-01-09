
import requests
import pandas as pd

SDSS_URL = "https://skyserver.sdss.org/dr18/SkyServerWS/SearchTools/SqlSearch"

def check_table(table_name, filter_term=None):
    sql = f"SELECT TOP 1 * FROM {table_name}"
    print(f"Checking {table_name}...")
    try:
        response = requests.get(SDSS_URL, params={"cmd": sql, "format": "json"}, timeout=30)
        if response.status_code == 200:
            data = response.json()
            if "Rows" in data[0]:
                df = pd.DataFrame(data[0]["Rows"])
                cols = sorted(df.columns.tolist())
                if filter_term:
                    filtered = [c for c in cols if filter_term.lower() in c.lower()]
                    print(f"Columns matching '{filter_term}': {filtered}")
                else:
                    print(f"Columns: {cols}")
                return cols
    except Exception as e:
        print(f"Error checking {table_name}: {e}")
    return []

def main():
    # Test AZ: Check aspcapStar for coordinates
    cols = check_table("aspcapStar", "ra")
    check_table("aspcapStar", "glon")
    
    # Test BB: Check spiders_quasar
    check_table("spiders_quasar", "flux")
    
    # Test BA: Check mangaDAPall for lambda
    check_table("mangaDAPall", "lambda")
    check_table("mangaDAPall", "eps")

if __name__ == "__main__":
    main()
