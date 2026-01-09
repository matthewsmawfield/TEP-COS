
import requests
import pandas as pd

SDSS_URL = "https://skyserver.sdss.org/dr18/SkyServerWS/SearchTools/SqlSearch"

def check_join():
    # Test CS requires joining HI data (mangaHIall) with Density data (ebossMCPM)
    # These might not share IDs. Spatial match is robust.
    
    sql = """
    SELECT TOP 10
        h.mangaid,
        h.objra, h.objdec,
        m.MATTERDENS,
        m.MJD, m.PLATE
    FROM mangaHIall h
    JOIN ebossMCPM m ON 
        abs(h.objra - m.RA) < 0.001 AND abs(h.objdec - m.DEC) < 0.001
    """
    print(f"Checking spatial join between mangaHIall and ebossMCPM...")
    try:
        response = requests.get(SDSS_URL, params={"cmd": sql, "format": "json"}, timeout=30)
        if response.status_code == 200:
            data = response.json()
            if "Rows" in data[0]:
                df = pd.DataFrame(data[0]["Rows"])
                print(f"Success! Found {len(df)} matches.")
                print(df.head())
            else:
                print("No rows returned. Spatial match too tight or no overlap.")
        else:
            print(f"HTTP {response.status_code}")
            print(response.text)

    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    check_join()
