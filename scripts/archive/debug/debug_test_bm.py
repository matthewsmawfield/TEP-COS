
import requests
import pandas as pd
import time

SDSS_URL = "https://skyserver.sdss.org/dr18/SkyServerWS/SearchTools/SqlSearch"

def query(sql, tag):
    print(f"--- {tag} ---")
    print(f"SQL: {sql}")
    try:
        start = time.time()
        response = requests.get(SDSS_URL, params={"cmd": sql, "format": "json"}, timeout=60)
        duration = time.time() - start
        print(f"Status: {response.status_code}, Time: {duration:.2f}s")
        
        if response.status_code == 200:
            data = response.json()
            if "Rows" in data[0]:
                df = pd.DataFrame(data[0]["Rows"])
                print(f"Success! Got {len(df)} rows.")
                print(df.head(2))
            else:
                print("No rows returned.")
        else:
            print(f"Error: {response.text[:200]}")
            
    except Exception as e:
        print(f"Exception: {e}")

def main():
    # Test 1: fparam columns
    query("SELECT TOP 10 apogee_id, fparam_alpha_m, fparam_m_h FROM aspcapStar", "Test 1: fparam columns")
    
    # Test 2: alpha_m, m_h
    query("SELECT TOP 10 apogee_id, alpha_m, m_h FROM aspcapStar", "Test 2: Short names")
    
    # Test 3: TEFF and LOGG checks (usually valid)
    query("SELECT TOP 10 apogee_id, teff, logg FROM aspcapStar", "Test 3: Teff/Logg")

if __name__ == "__main__":
    main()
