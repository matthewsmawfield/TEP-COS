
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
                print(df.head(5))
            else:
                print("No rows returned.")
        else:
            print(f"Error: {response.status_code} {response.text[:200]}")
    except Exception as e:
        print(f"Exception: {e}")

def main():
    # Check join size without filters
    sql_join = """
    SELECT count(*) as count 
    FROM spiders_quasar q
    JOIN galSpecInfo g ON q.SPECOBJID = g.specObjID
    """
    query(sql_join, "Join Spiders-GalSpecInfo Total")
    
    # Check v_disp validity in join
    sql_vdisp = """
    SELECT count(*) as count 
    FROM spiders_quasar q
    JOIN galSpecInfo g ON q.SPECOBJID = g.specObjID
    WHERE g.v_disp > 0
    """
    query(sql_vdisp, "Join with v_disp > 0")
    
    # Check v_disp range
    sql_range = """
    SELECT TOP 20 g.v_disp, q.logBHMS_mgII
    FROM spiders_quasar q
    JOIN galSpecInfo g ON q.SPECOBJID = g.specObjID
    WHERE g.v_disp > 0
    """
    query(sql_range, "Sample v_disp values")

if __name__ == "__main__":
    main()
