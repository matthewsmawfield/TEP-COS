
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
    # Check overlap with galSpecExtra (MPA-JHU)
    sql_join_gal = """
    SELECT count(*) as count 
    FROM spiders_quasar q
    JOIN galSpecExtra g ON q.SPECOBJID = g.specObjID
    """
    query(sql_join_gal, "Join Spiders-GalSpecExtra")
    
    # Check spiders_quasar for sigma columns
    sql_cols = "SELECT TOP 1 * FROM spiders_quasar"
    # We already checked cols, but let's check for 'width' or 'fwhm' which might be proxies?
    # Actually, let's just check the join first.
    
    # Also check if there is another table 'dr16q' with mass
    # sql_tables = "SELECT * FROM outputTable WHERE name LIKE '%dr16%'" -- not easy to list tables this way.
    
    pass

if __name__ == "__main__":
    main()
