from astroquery.vizier import Vizier
from astropy.table import Table
import os

def download_q0957():
    print("Querying VizieR for J/A+A/492/401 (Shalyapin 2008)...")
    Vizier.ROW_LIMIT = -1
    
    try:
        # returns a TableList (list of Table)
        catalogs = Vizier.get_catalogs("J/A+A/492/401")
        
        print(f"Found {len(catalogs)} tables.")
        for table in catalogs:
            # table is an astropy Table
            # meta data is in table.meta (dict)
            name = table.meta.get('name', 'unknown')
            desc = table.meta.get('description', 'no description')
            print(f"Table: {name} - {desc}")
            
            # Identify bands
            # Table 1: g band light curves
            # Table 2: r band light curves
            
            if "table1" in name:
                filename = "data/cosmograil/q0957_g_shalyapin.csv"
                table.write(filename, format="csv", overwrite=True)
                print(f"Saved g-band to {filename}")
            elif "table2" in name:
                filename = "data/cosmograil/q0957_r_shalyapin.csv"
                table.write(filename, format="csv", overwrite=True)
                print(f"Saved r-band to {filename}")
            else:
                safe_name = name.replace("/", "_").replace("+", "p")
                filename = f"data/cosmograil/q0957_{safe_name}.csv"
                table.write(filename, format="csv", overwrite=True)
                print(f"Saved {name} to {filename}")
            
    except Exception as e:
        print(f"Error downloading Q0957: {e}")

def download_glendama():
    print("\nQuerying VizieR for J/A+A/616/A118 (GLENDAMA)...")
    Vizier.ROW_LIMIT = -1
    
    try:
        catalogs = Vizier.get_catalogs("J/A+A/616/A118")
        
        if not catalogs:
            print("No catalogs found for GLENDAMA.")
            return

        print(f"Found {len(catalogs)} tables.")
        for table in catalogs:
            name = table.meta.get('name', 'unknown')
            print(f"Table: {name}")
            
            safe_name = name.replace("/", "_").replace("+", "p")
            filename = f"data/cosmograil/glendama_{safe_name}.csv"
            table.write(filename, format="csv", overwrite=True)
            print(f"Saved to {filename}")
            
    except Exception as e:
        print(f"Error downloading GLENDAMA: {e}")

if __name__ == "__main__":
    # Ensure directory exists
    os.makedirs("data/cosmograil", exist_ok=True)
    
    download_q0957()
    download_glendama()
