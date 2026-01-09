
import requests
import pandas as pd

SDSS_URL = "https://skyserver.sdss.org/dr18/SkyServerWS/SearchTools/SqlSearch"

def check():
    # 1. Check overlap mos_sdss_dr16_qso and qsoVarStripe
    sql = """
    SELECT count(*) as count
    FROM mos_sdss_dr16_qso q
    JOIN qsoVarStripe v ON q.objid = v.VAR_OBJID
    """
    print("Checking overlap mos_sdss_dr16_qso and qsoVarStripe...")
    try:
        response = requests.get(SDSS_URL, params={"cmd": sql, "format": "json"}, timeout=30)
        if response.status_code == 200:
            data = response.json()
            if "Rows" in data[0]:
                print(f"Overlap count: {data[0]['Rows'][0]['count']}")
        else:
            print(f"Overlap check failed: {response.status_code}")
    except Exception as e:
        print(f"Overlap check error: {e}")

    # 2. Check mos_sdss_dr16_qso for OIII columns
    sql = "SELECT TOP 1 * FROM mos_sdss_dr16_qso"
    print("\nChecking mos_sdss_dr16_qso columns...")
    try:
        response = requests.get(SDSS_URL, params={"cmd": sql, "format": "json"}, timeout=30)
        if response.status_code == 200:
            data = response.json()
            if "Rows" in data[0]:
                df = pd.DataFrame(data[0]["Rows"])
                cols = df.columns.tolist()
                oiii_cols = [c for c in cols if "oiii" in c.lower() or "5007" in c]
                width_cols = [c for c in cols if "fwhm" in c.lower() or "width" in c.lower() or "sigma" in c.lower()]
                print(f"OIII columns: {sorted(oiii_cols)}")
                print(f"Width columns: {sorted(width_cols)}")
    except Exception as e:
        print(f"Column check error: {e}")

if __name__ == "__main__":
    check()
