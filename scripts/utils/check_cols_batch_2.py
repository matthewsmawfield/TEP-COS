
import requests
import pandas as pd

SDSS_URL = "https://skyserver.sdss.org/dr18/SkyServerWS/SearchTools/SqlSearch"

def check_table(table_name, keywords=None):
    sql = f"SELECT TOP 1 * FROM {table_name}"
    print(f"Checking {table_name}...")
    try:
        response = requests.get(SDSS_URL, params={"cmd": sql, "format": "json"}, timeout=30)
        if response.status_code == 200:
            data = response.json()
            if "Rows" in data[0]:
                df = pd.DataFrame(data[0]["Rows"])
                cols = sorted(df.columns.tolist())
                if keywords:
                    found = []
                    for k in keywords:
                        matches = [c for c in cols if k.lower() in c.lower()]
                        found.extend(matches)
                    print(f"  Keywords {keywords} -> Found: {found}")
                else:
                    print(f"  Cols: {cols[:10]}...") # Print first 10 if no keywords
                return cols
            else:
                print(f"  [Empty] No rows returned.")
                return []
        else:
            print(f"  [Error] HTTP {response.status_code}")
            return []
    except Exception as e:
        print(f"  [Error] {e}")
        return []

def main():
    # Test BC: YSO Contraction
    check_table("mos_sagitta", ["age", "mass", "g_mag", "bp_rp"])
    
    # Test BD: Pipe3D vs Firefly
    check_table("mangaPipe3D", ["age", "mass", "log_mass"])
    check_table("mangaFirefly", ["age", "lw_age"])
    
    # Test BF: AGB Dust
    check_table("mos_allwise", ["w3", "w4", "mpro"]) # Check if table exists too
    
    # Test BG: M-sigma
    check_table("mos_sdss_dr16_qso", ["bh", "sigma", "mass"])
    
    # Test BH: Satellite Quenching
    check_table("ebossMCPM", ["dens", "mid"])
    
    # Test BK: Warp
    check_table("MaNGA_GZ2", ["warp", "disk", "total"])

if __name__ == "__main__":
    main()
