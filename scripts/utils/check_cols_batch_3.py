
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
                    print(f"  Cols ({len(cols)}): {cols[:10]}...")
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
    # Test BF: AGB Dust
    check_table("mos_gaia_dr2_source", ["source_id", "phot_g_mean_mag", "bp_rp", "l", "b"])
    check_table("mos_allwise", ["w3", "w4", "mpro"]) 
    check_table("mos_geometric_distances_gaia_dr2", ["r_est", "dist"])
    
    # Test BG: M-sigma (QSO)
    check_table("mos_sdss_dr16_qso", ["sigma", "bh", "mass"])
    
    # Test BK: Warp
    check_table("MaNGA_GZ2", ["warp", "disk", "inclination"])
    
    # Test BH: Satellite Quenching
    check_table("ebossMCPM", ["dens", "mid_dens"])
    
    # Test BI: FMR
    check_table("galSpecExtra", ["sfr", "oh_p50"])
    
    # Test BJ: Halo Escape
    check_table("apogeeStar", ["vhelio", "rv"])

if __name__ == "__main__":
    main()
