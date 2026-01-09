
import requests
import pandas as pd

SDSS_URL = "https://skyserver.sdss.org/dr18/SkyServerWS/SearchTools/SqlSearch"

def test_post():
    sql = "SELECT TOP 5 ra, dec FROM apogeeStar"
    print(f"Testing POST to SDSS...")
    try:
        # payload for POST
        data = {"cmd": sql, "format": "json"}
        response = requests.post(SDSS_URL, data=data, timeout=30)
        
        if response.status_code == 200:
            print("POST Success!")
            json_data = response.json()
            if "Rows" in json_data[0]:
                df = pd.DataFrame(json_data[0]["Rows"])
                print(df)
            else:
                print("No rows returned.")
        else:
            print(f"POST Failed. HTTP {response.status_code}")
            print(response.text)
            
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    test_post()
