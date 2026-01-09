
import requests
import pandas as pd

SDSS_URL = "https://skyserver.sdss.org/dr18/SkyServerWS/SearchTools/SqlSearch"

def check_tables():
    print("Checking Test DO tables...")
    
    tables = ["redmapperDR12", "eFEDS_Main_speccomp"]
    
    for t in tables:
        sql = f"SELECT TOP 1 * FROM {t}"
        try:
            response = requests.get(SDSS_URL, params={"cmd": sql, "format": "json"}, timeout=30)
            if response.status_code == 200:
                data = response.json()
                if "Rows" in data[0]:
                    print(f"  {t}: Found.")
                    df = pd.DataFrame(data[0]["Rows"])
                    print(f"  Columns: {df.columns.tolist()}")
                else:
                    print(f"  {t}: Found but empty.")
            else:
                print(f"  {t}: Not Found (Status {response.status_code})")
        except Exception as e:
            print(f"  {t}: Error {e}")

if __name__ == "__main__":
    check_tables()
