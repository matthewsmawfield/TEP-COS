
import requests
import pandas as pd
import os

SDSS_URL = "https://skyserver.sdss.org/dr18/SkyServerWS/SearchTools/SqlSearch"
OUT_FILE = "lithium_cols.txt"

def check_tables():
    print("Checking tables for Lithium columns...")
    
    with open(OUT_FILE, "w") as f:
        for table in ["sppLines", "sppParams"]:
            sql = f"SELECT TOP 1 * FROM {table}"
            try:
                response = requests.get(SDSS_URL, params={"cmd": sql, "format": "json"}, timeout=30)
                if response.status_code == 200:
                    data = response.json()
                    if "Rows" in data[0]:
                        df = pd.DataFrame(data[0]["Rows"])
                        cols = sorted(df.columns.tolist())
                        f.write(f"--- {table} ---\n")
                        f.write(", ".join(cols) + "\n\n")
                        
                        # Check for Li specific
                        li_cols = [c for c in cols if 'li' in c.lower()]
                        f.write(f"Lithium candidates in {table}: {li_cols}\n\n")
                    else:
                        f.write(f"{table} returned no rows.\n")
                else:
                    f.write(f"{table} error: {response.status_code}\n")
            except Exception as e:
                f.write(f"{table} exception: {e}\n")
    
    print(f"Columns written to {OUT_FILE}")

if __name__ == "__main__":
    check_tables()
