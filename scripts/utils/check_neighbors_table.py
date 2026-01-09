
import requests
import pandas as pd

SDSS_URL = "https://skyserver.sdss.org/dr18/SkyServerWS/SearchTools/SqlSearch"

def check_neighbors():
    print("Checking Neighbors table...")
    
    # Check if table exists by selecting top 1
    sql = "SELECT TOP 1 * FROM Neighbors"
    
    try:
        response = requests.get(SDSS_URL, params={"cmd": sql, "format": "json"}, timeout=30)
        if response.status_code == 200:
            data = response.json()
            if "Rows" in data[0]:
                df = pd.DataFrame(data[0]["Rows"])
                print("Neighbors table found.")
                print("Columns:", df.columns.tolist())
            else:
                print("Neighbors table found but empty or no rows returned.")
        else:
            print(f"Error: Status {response.status_code}")
            print(response.text[:200])
            
            # Try lowercase just in case (usually case insensitive but good to check)
            print("Trying lowercase 'neighbors'...")
            sql_lower = "SELECT TOP 1 * FROM neighbors"
            response = requests.get(SDSS_URL, params={"cmd": sql_lower, "format": "json"}, timeout=30)
            if response.status_code == 200:
                 print("Lowercase table found.")
    except Exception as e:
        print(f"Exception: {e}")

if __name__ == "__main__":
    check_neighbors()
