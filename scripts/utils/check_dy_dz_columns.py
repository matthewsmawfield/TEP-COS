
import requests
import pandas as pd

SDSS_URL = "https://skyserver.sdss.org/dr18/SkyServerWS/SearchTools/SqlSearch"

def check_tables():
    print("Checking specific columns...")
    
    # Check mos_gaia_dr2_source columns
    sql_dy = "SELECT TOP 1 * FROM mos_gaia_dr2_source"
    try:
        response = requests.get(SDSS_URL, params={"cmd": sql_dy, "format": "json"}, timeout=30)
        if response.status_code == 200:
            data = response.json()
            if "Rows" in data[0]:
                print("mos_gaia_dr2_source: Found.")
                df = pd.DataFrame(data[0]["Rows"])
                print(f"  Columns: {df.columns.tolist()}")
            else:
                print("mos_gaia_dr2_source: Found but empty.")
        else:
            print(f"mos_gaia_dr2_source: Failed {response.status_code}")
    except Exception as e:
        print(f"mos_gaia_dr2_source Error: {e}")

    # Check aspcapStar columns
    sql_dz = "SELECT TOP 1 * FROM aspcapStar"
    try:
        response = requests.get(SDSS_URL, params={"cmd": sql_dz, "format": "json"}, timeout=30)
        if response.status_code == 200:
            data = response.json()
            if "Rows" in data[0]:
                print("aspcapStar: Found.")
                df = pd.DataFrame(data[0]["Rows"])
                cols = df.columns.tolist()
                # Filter for relevant ones to keep output short
                rel_cols = [c for c in cols if 'param' in c or 'FE_' in c or 'MG_' in c or 'K_' in c]
                print(f"  Relevant Columns: {rel_cols}")
            else:
                print("aspcapStar: Found but empty.")
        else:
            print(f"aspcapStar: Failed {response.status_code}")
    except Exception as e:
        print(f"aspcapStar Error: {e}")

if __name__ == "__main__":
    check_tables()
