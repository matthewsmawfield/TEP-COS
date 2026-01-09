
import requests
import pandas as pd

SDSS_URL = "https://skyserver.sdss.org/dr18/SkyServerWS/SearchTools/SqlSearch"

def query(sql, tag):
    print(f"--- {tag} ---")
    print(f"SQL: {sql}")
    try:
        response = requests.get(SDSS_URL, params={"cmd": sql, "format": "json"}, timeout=60)
        if response.status_code == 200:
            data = response.json()
            if "Rows" in data[0]:
                df = pd.DataFrame(data[0]["Rows"])
                print(f"Success! Got {len(df)} rows.")
                print(df.head(2))
            else:
                print("No rows returned.")
        else:
            print(f"Error: {response.status_code} {response.text[:200]}")
    except Exception as e:
        print(f"Exception: {e}")

def main():
    # Check counts
    query("SELECT count(*) as count FROM spiders_quasar", "Total Spiders")
    query("SELECT count(*) as count FROM spiders_quasar WHERE logBHMS_mgII > 0", "Spiders with MgII Mass")
    query("SELECT count(*) as count FROM spiders_quasar WHERE logBHMA_hb > 0", "Spiders with Hb Mass")
    
    # Check overlap
    sql_join = """
    SELECT count(*) as count 
    FROM spiders_quasar q
    JOIN emissionLinesPort e ON q.SPECOBJID = e.specObjID
    """
    query(sql_join, "Join Spiders-EmissionLines")
    
    # Check overlap with constraints
    sql_join_constr = """
    SELECT count(*) as count 
    FROM spiders_quasar q
    JOIN emissionLinesPort e ON q.SPECOBJID = e.specObjID
    WHERE e.sigmaStars > 50
    AND (q.logBHMS_mgII > 6 OR q.logBHMA_hb > 6)
    """
    query(sql_join_constr, "Join with Constraints")

if __name__ == "__main__":
    main()
