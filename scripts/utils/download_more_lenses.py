from astroquery.vizier import Vizier
from astropy.table import Table
import os

def download_glendama():
    print("Querying VizieR for J/A+A/616/A118 (GLENDAMA)...")
    Vizier.ROW_LIMIT = -1
    
    try:
        # GLENDAMA database paper
        catalogs = Vizier.get_catalogs("J/A+A/616/A118")
        
        print(f"Found {len(catalogs)} tables.")
        for name, table in catalogs.items():
            print(f"Table: {name}")
            # print(table.info)
            
            safe_name = name.replace("/", "_").replace("+", "p")
            filename = f"data/cosmograil/glendama_{safe_name}.csv"
            table.write(filename, format="csv", overwrite=True)
            print(f"Saved to {filename}")
            
    except Exception as e:
        print(f"Error downloading GLENDAMA: {e}")

def download_q2237_ogle():
    print("\nQuerying VizieR for OGLE Q2237 data...")
    # OGLE-III Q2237: J/A+A/529/A146 (Udalski+ 2011) ? Or similar.
    # Let's search for Q2237 object
    
    try:
        # Searching by object name might return many catalogs
        # We specifically want multi-band.
        # Let's try to find the OGLE Q2237 monitoring.
        # J/ApJ/659/1040 ? No that's RXJ1131.
        
        # Try to find catalogs for object "Q2237+0305"
        catalogs = Vizier.query_object("Q2237+0305")
        print(f"Found {len(catalogs)} catalogs for Q2237+0305")
        
        for cat in catalogs:
            desc = cat.description
            name = cat.name
            if "light curve" in desc.lower() or "photometry" in desc.lower():
                print(f"Candidate: {name} - {desc}")
                
    except Exception as e:
        print(f"Error searching Q2237: {e}")

if __name__ == "__main__":
    download_glendama()
    download_q2237_ogle()
