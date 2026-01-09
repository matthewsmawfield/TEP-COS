
import requests
import pandas as pd

SDSS_URL = "https://skyserver.sdss.org/dr18/SkyServerWS/SearchTools/SqlSearch"

def check_cols():
    # Check emissionLinesPort for bpt
    sql = "SELECT TOP 1 * FROM emissionLinesPort"
    print(f"Checking emissionLinesPort...")
    try:
        response = requests.get(SDSS_URL, params={"cmd": sql, "format": "json"}, timeout=30)
        if response.status_code == 200:
            data = response.json()
            if "Rows" in data[0]:
                df = pd.DataFrame(data[0]["Rows"])
                cols = df.columns.tolist()
                print(f"Columns in emissionLinesPort: {cols}")
                if 'bpt' in cols:
                    print("FOUND: bpt in emissionLinesPort")
                else:
                    print("NOT FOUND: bpt in emissionLinesPort")
    except Exception as e:
        print(f"Error: {e}")

    # Check galSpecExtra for bptclass
    sql = "SELECT TOP 1 * FROM galSpecExtra"
    print(f"\nChecking galSpecExtra...")
    try:
        response = requests.get(SDSS_URL, params={"cmd": sql, "format": "json"}, timeout=30)
        if response.status_code == 200:
            data = response.json()
            if "Rows" in data[0]:
                df = pd.DataFrame(data[0]["Rows"])
                cols = df.columns.tolist()
                print(f"Columns in galSpecExtra: {cols}")
                if 'bptclass' in cols:
                    print("FOUND: bptclass in galSpecExtra")
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    check_cols()
