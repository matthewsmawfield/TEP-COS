
import requests
import pandas as pd

SDSS_URL = "https://skyserver.sdss.org/dr18/SkyServerWS/SearchTools/SqlSearch"

def check_zoo2():
    print("Checking zoo2MainSpecz access without joins...")
    
    sql = "SELECT TOP 10 * FROM zoo2MainSpecz"
    
    try:
        response = requests.get(SDSS_URL, params={"cmd": sql, "format": "json"}, timeout=30)
        if response.status_code == 200:
            data = response.json()
            if "Rows" in data[0]:
                df = pd.DataFrame(data[0]["Rows"])
                print("zoo2MainSpecz access successful.")
                print("Columns:", df.columns.tolist())
            else:
                print("zoo2MainSpecz returned no rows.")
        else:
            print(f"Error: Status {response.status_code}")
            print(response.text[:200])
            
    except Exception as e:
        print(f"Exception: {e}")

if __name__ == "__main__":
    check_zoo2()
