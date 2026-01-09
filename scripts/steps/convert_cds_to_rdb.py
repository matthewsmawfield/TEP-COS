#!/usr/bin/env python3
"""
Convert CDS-format COSMOGRAIL light curves to RDB format for temporal shear analysis.
"""

import sys
from pathlib import Path


def convert_pg1115(input_path: Path, output_path: Path):
    """Convert PG1115 CDS format (3 columns: A1A2, B, C) to RDB."""
    lines = input_path.read_text().strip().split('\n')
    
    with open(output_path, 'w') as f:
        f.write("mhjd\tmag_A\tmagerr_A\tmag_B\tmagerr_B\tmag_C\tmagerr_C\n")
        f.write("====\t=====\t========\t=====\t========\t=====\t========\n")
        
        for line in lines:
            parts = line.split()
            if len(parts) >= 7:
                mhjd = parts[0]
                mag_a, err_a = parts[1], parts[2]
                mag_b, err_b = parts[3], parts[4]
                mag_c, err_c = parts[5], parts[6]
                f.write(f"{mhjd}\t{mag_a}\t{err_a}\t{mag_b}\t{err_b}\t{mag_c}\t{err_c}\n")
    
    print(f"Converted {len(lines)} epochs to {output_path}")


def convert_wfi2033(input_paths: list, output_path: Path):
    """Convert WFI2033 CDS format files (multiple telescopes) to single RDB."""
    all_data = []
    
    for input_path in input_paths:
        lines = input_path.read_text().strip().split('\n')
        for line in lines:
            parts = line.split()
            if len(parts) >= 6:
                mhjd = float(parts[0])
                mag_a, err_a = float(parts[1]), float(parts[2])
                mag_b, err_b = float(parts[3]), float(parts[4])
                mag_c, err_c = float(parts[5]), float(parts[6]) if len(parts) > 6 else (0.0, 0.01)
                all_data.append((mhjd, mag_a, err_a, mag_b, err_b, mag_c, err_c))
    
    # Sort by time
    all_data.sort(key=lambda x: x[0])
    
    with open(output_path, 'w') as f:
        f.write("mhjd\tmag_A\tmagerr_A\tmag_B\tmagerr_B\tmag_C\tmagerr_C\n")
        f.write("====\t=====\t========\t=====\t========\t=====\t========\n")
        
        for row in all_data:
            f.write(f"{row[0]:.5f}\t{row[1]:.5f}\t{row[2]:.5f}\t{row[3]:.5f}\t{row[4]:.5f}\t{row[5]:.5f}\t{row[6]:.5f}\n")
    
    print(f"Converted {len(all_data)} epochs to {output_path}")


def main():
    data_dir = Path(__file__).parent.parent.parent / 'data' / 'cosmograil'
    
    # Convert PG1115
    pg1115_in = data_dir / 'PG1115_Bonvin2018.dat'
    pg1115_out = data_dir / 'PG1115_Bonvin2018.rdb'
    if pg1115_in.exists():
        convert_pg1115(pg1115_in, pg1115_out)
    
    # Convert WFI2033 (merge all telescope data)
    wfi2033_files = [
        data_dir / 'WFI2033_ecam.dat',
        data_dir / 'WFI2033_smarts.dat',
        data_dir / 'WFI2033_wfi.dat',
    ]
    wfi2033_files = [f for f in wfi2033_files if f.exists()]
    if wfi2033_files:
        wfi2033_out = data_dir / 'WFI2033_Bonvin2019.rdb'
        convert_wfi2033(wfi2033_files, wfi2033_out)
    
    print("\nDone! New RDB files ready for analysis.")


if __name__ == '__main__':
    main()
