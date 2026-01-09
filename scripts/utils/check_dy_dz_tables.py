
import requests
import pandas as pd

SDSS_URL = "https://skyserver.sdss.org/dr18/SkyServerWS/SearchTools/SqlSearch"

def check_tables():
    print("Checking Test DY and DZ tables...")
    
    # Test DY Tables
    dy_tables = ["apogeeStar", "apogee_starhorse", "apogeeObject", "mos_gaia_dr2_source", "gaia_dr2_source"]
    
    # Test DZ Tables
    dz_tables = ["aspcapStar"]
    
    print("\n--- Test DY Tables ---")
    for t in dy_tables:
        sql = f"SELECT TOP 1 * FROM {t}"
        try:
            response = requests.get(SDSS_URL, params={"cmd": sql, "format": "json"}, timeout=30)
            if response.status_code == 200:
                data = response.json()
                if "Rows" in data[0]:
                    print(f"  {t}: Found.")
                    df = pd.DataFrame(data[0]["Rows"])
                    # print(f"  Columns: {df.columns.tolist()}") 
                else:
                    print(f"  {t}: Found but empty.")
            else:
                print(f"  {t}: Not Found (Status {response.status_code})")
        except Exception as e:
            print(f"  {t}: Error {e}")

    print("\n--- Test DZ Columns ---")
    # Check for K abundance in aspcapStar
    sql = "SELECT TOP 1 param_k_h, param_fe_h, param_mg_h FROM aspcapStar"
    try:
        response = requests.get(SDSS_URL, params={"cmd": sql, "format": "json"}, timeout=30)
        if response.status_code == 200:
            data = response.json()
            if "Rows" in data[0]:
                print("  aspcapStar columns (K, Fe, Mg): Found.")
                df = pd.DataFrame(data[0]["Rows"])
                print(f"  Sample: {df.iloc[0].to_dict()}")
            else:
                print("  aspcapStar columns: Not Found (Empty result).")
        else:
            print(f"  aspcapStar query failed (Status {response.status_code})")
    except Exception as e:
        print(f"  aspcapStar query Error: {e}")

if __name__ == "__main__":
    check_tables()
