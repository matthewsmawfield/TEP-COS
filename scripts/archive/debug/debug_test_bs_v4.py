
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
                print(df)
            else:
                print("No rows returned.")
        else:
            print(f"Error: {response.status_code} {response.text[:200]}")
    except Exception as e:
        print(f"Exception: {e}")

def main():
    # Check valid overlap
    sql_valid = """
    SELECT count(*) as count 
    FROM spiders_quasar q
    JOIN galSpecInfo g ON q.SPECOBJID = g.specObjID
    WHERE g.v_disp > 50 
      AND (
          (q.logBHMS_mgII > 0 AND q.logBHMS_mgII < 99) 
          OR 
          (q.logBHMA_hb > 0 AND q.logBHMA_hb < 99)
      )
    """
    query(sql_valid, "Valid BS Sample Size")

if __name__ == "__main__":
    main()
