
import requests
import pandas as pd

SDSS_URL = "https://skyserver.sdss.org/dr18/SkyServerWS/SearchTools/SqlSearch"

def check_table(table_name, keywords):
    sql = f"SELECT TOP 1 * FROM {table_name}"
    print(f"Checking {table_name}...")
    try:
        response = requests.get(SDSS_URL, params={"cmd": sql, "format": "json"}, timeout=30)
        if response.status_code == 200:
            data = response.json()
            if "Rows" in data[0]:
                df = pd.DataFrame(data[0]["Rows"])
                cols = sorted(df.columns.tolist())
                
                found = []
                for k in keywords:
                    matches = [c for c in cols if k.lower() in c.lower()]
                    found.extend(matches)
                
                print(f"  Keywords {keywords} -> Found: {found}")
                if len(cols) < 20:
                    print(f"  All Cols: {cols}")
    except Exception as e:
        print(f"Error checking {table_name}: {e}")

def main():
    # Test AZ: Coordinates
    check_table("aspcapStar", ["ra", "dec", "glon", "glat", "lon", "lat"])
    
    # Test BB: X-ray Flux
    check_table("spiders_quasar", ["flux", "rate", "cts", "xray", "erosita"])
    
    # Test BA: Rotation
    check_table("mangaDAPall", ["lambda", "ellip", "vel", "sigma", "rot"])

if __name__ == "__main__":
    main()
