import pandas as pd
from pathlib import Path

def identify_tables():
    path = Path("data/cosmograil")
    files = list(path.glob("glendama*.csv"))
    
    for p in files:
        try:
            df = pd.read_csv(p, nrows=5)
            print(f"\n--- {p.name} ---")
            print(df.columns.tolist())
            print(df.head())
            print(f"Mean Mag A: {df['mA'].mean() if 'mA' in df.columns else 'N/A'}")
        except Exception as e:
            print(f"Error reading {p}: {e}")

if __name__ == "__main__":
    identify_tables()
