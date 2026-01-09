
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
            else:
                print("No rows returned.")
        else:
            print(f"Error: {response.text[:200]}")
            
    except Exception as e:
        print(f"Exception: {e}")

def main():
    # Test 1: Just Mass
    sql1 = "SELECT TOP 10 specObjID, logMass FROM stellarMassFSPSGranWideDust WHERE logMass > 10.5"
    query(sql1, "Test 1: Mass Table")
    
    # Test 2: Mass + Emission (Sigma)
    sql2 = """
    SELECT TOP 10 s.specObjID, s.logMass, e.sigma_stars
    FROM stellarMassFSPSGranWideDust s
    JOIN emissionLinesPort e ON s.specObjID = e.specObjID
    WHERE s.logMass > 10.5 AND e.sigma_stars > 50
    """
    query(sql2, "Test 2: Mass + Emission")
    
    # Test 3: Mass + Emission + SpecObj (Z, ObjID)
    sql3 = """
    SELECT TOP 10 s.specObjID, s.logMass, e.sigma_stars, sp.z, sp.bestObjID
    FROM stellarMassFSPSGranWideDust s
    JOIN emissionLinesPort e ON s.specObjID = e.specObjID
    JOIN SpecObjAll sp ON s.specObjID = sp.specObjID
    WHERE s.logMass > 10.5 AND e.sigma_stars > 50
    """
    query(sql3, "Test 3: Mass + Emission + SpecObj")
    
    # Test 4: All 4 tables
    sql4 = """
    SELECT TOP 10 s.specObjID, s.logMass, e.sigma_stars, sp.z, ph.petroR50_r
    FROM stellarMassFSPSGranWideDust s
    JOIN emissionLinesPort e ON s.specObjID = e.specObjID
    JOIN SpecObjAll sp ON s.specObjID = sp.specObjID
    JOIN PhotoObjAll ph ON sp.bestObjID = ph.objID
    WHERE s.logMass > 10.5 AND e.sigma_stars > 50
    """
    query(sql4, "Test 4: Full Join")

if __name__ == "__main__":
    main()
