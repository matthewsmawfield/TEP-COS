from astroquery.vizier import Vizier
import astropy.units as u
from astropy.table import Table

def download_q0957():
    print("Querying VizieR for J/A+A/492/401...")
    
    # Enable downloading all columns
    Vizier.ROW_LIMIT = -1
    
    try:
        catalogs = Vizier.get_catalogs("J/A+A/492/401")
        
        print(f"Found {len(catalogs)} tables.")
        for name, table in catalogs.items():
            print(f"Table: {name}")
            print(table.info)
            print(table[:5])
            
            # Save tables
            # Table 1: g-band? Table 2: r-band?
            # Usually Vizier names them like J/A+A/492/401/table1
            
            safe_name = name.replace("/", "_").replace("+", "p")
            filename = f"data/cosmograil/q0957_{safe_name}.csv"
            table.write(filename, format="csv", overwrite=True)
            print(f"Saved to {filename}")
            
    except Exception as e:
        print(f"Error downloading catalog: {e}")

if __name__ == "__main__":
    download_q0957()
