
import requests
import pandas as pd

SDSS_URL = "https://skyserver.sdss.org/dr18/SkyServerWS/SearchTools/SqlSearch"

def check_manga_hi():
    print("Checking mangaHIall table...")
    
    sql = "SELECT TOP 1 * FROM mangaHIall"
    
    try:
        response = requests.get(SDSS_URL, params={"cmd": sql, "format": "json"}, timeout=30)
        if response.status_code == 200:
            data = response.json()
            if "Rows" in data[0]:
                df = pd.DataFrame(data[0]["Rows"])
                print("mangaHIall table found.")
                print("Columns:", df.columns.tolist())
            else:
                print("mangaHIall table found but empty.")
        else:
            print(f"Error: Status {response.status_code}")
            
    except Exception as e:
        print(f"Exception: {e}")

if __name__ == "__main__":
    check_manga_hi()
