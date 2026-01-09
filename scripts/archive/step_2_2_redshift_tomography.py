#!/usr/bin/env python3

import sys
from pathlib import Path

# Ensure repo root is importable when executing this script directly
PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import argparse
import json
import math
import pandas as pd
import numpy as np
from astropy.io import fits
from typing import Dict, List, Tuple, Optional

from scripts.utils.logger import TEPLogger, print_status, set_step_logger

def weighted_linear_fit(x: np.ndarray, y: np.ndarray, w: np.ndarray) -> Tuple[float, float, float, float]:
    wsum = np.sum(w)
    xbar = np.sum(w * x) / wsum
    ybar = np.sum(w * y) / wsum
    xx = np.sum(w * (x - xbar) ** 2)
    xy = np.sum(w * (x - xbar) * (y - ybar))
    if xx == 0:
        return 0.0, 0.0, 0.0, 0.0
    a = xy / xx
    b = ybar - a * xbar
    yhat = a * x + b
    resid = y - yhat
    dof = max(int(x.size) - 2, 1)
    s2 = float(np.sum(w * resid**2) / dof)
    a_se = math.sqrt(s2 / xx)
    b_se = math.sqrt(s2 * (1.0 / wsum + (xbar**2) / xx))
    return float(a), float(b), float(a_se), float(b_se)

def robust_huber_fit(x: np.ndarray, y: np.ndarray, w: np.ndarray) -> Tuple[float, float, float]:
    # Simple Huber loop
    a, b, _, _ = weighted_linear_fit(x, y, w)
    for _ in range(20):
        r = y - (a * x + b)
        mad = np.median(np.abs(r))
        s = max(1.4826 * mad, 1e-12)
        u = r / (s * 1.345)
        wr = np.ones_like(u)
        mask = np.abs(u) > 1.0
        wr[mask] = 1.0 / np.abs(u[mask])
        w_eff = w * wr
        a, b, _, _ = weighted_linear_fit(x, y, w_eff)
    return float(a), float(b), float(s)

def main():
    parser = argparse.ArgumentParser(description="TEP-COS Step 2.2 - Redshift Tomography ('Onion Layer') Analysis")
    parser.add_argument("--stellar-csv", required=True, help="Path to step_2_0_per_galaxy CSV for stellar data")
    parser.add_argument("--gas-csv", help="Path to step_2_0_per_galaxy CSV for gas data (optional)")
    parser.add_argument("--dapall", required=True, help="Path to dapall fits file")
    parser.add_argument("--n-bins", type=int, default=4, help="Number of redshift bins (onion layers)")
    parser.add_argument("--output-dir", default="results/outputs", help="Output directory")
    args = parser.parse_args()

    logger = TEPLogger("step_2_2_redshift_tomography", log_file_path=PROJECT_ROOT / "logs" / "step_2_2_redshift_tomography.log")
    set_step_logger(logger)
    
    # 1. Load Redshift Data
    print_status(f"Loading redshifts from {args.dapall}...", "PROCESS")
    with fits.open(args.dapall) as hdul:
        # Try finding the table. Usually extension 1 or named by DAPTYPE
        # We'll just look for PLATEIFU and Z in the first table extension
        data = None
        for hdu in hdul:
            if isinstance(hdu, fits.BinTableHDU) and "PLATEIFU" in hdu.columns.names:
                data = hdu.data
                break
        
        if data is None:
            raise RuntimeError("Could not find table with PLATEIFU in dapall")
            
        z_map = {}
        for row in data:
            pid = str(row["PLATEIFU"]).strip()
            # Prefer NSA_Z (NASA Sloan Atlas) if available, else Z
            z_val = row["NSA_Z"] if "NSA_Z" in row.array.names else row["Z"]
            z_map[pid] = float(z_val)
            
    print_status(f"Loaded {len(z_map)} redshifts.", "SUCCESS")

    # 2. Analyze Stellar
    analyze_component(args.stellar_csv, z_map, args.n_bins, "stellar", Path(args.output_dir))

    # 3. Analyze Gas (if provided)
    if args.gas_csv:
        analyze_component(args.gas_csv, z_map, args.n_bins, "gas", Path(args.output_dir))

def analyze_component(csv_path: str, z_map: Dict[str, float], n_bins: int, label: str, out_dir: Path):
    print_status(f"Analyzing {label} component from {csv_path}...", "PROCESS")
    df = pd.read_csv(csv_path)
    
    # Add Z column
    df["z"] = df["plateifu"].apply(lambda x: z_map.get(str(x).strip(), np.nan))
    
    # Filter missing Z
    n_total = len(df)
    df = df.dropna(subset=["z"])
    n_valid = len(df)
    if n_valid < n_total:
        print_status(f"Dropped {n_total - n_valid} rows due to missing redshift.", "WARNING")

    # Sort by Z
    df = df.sort_values("z")
    
    # Create Bins
    df["bin"] = pd.qcut(df["z"], n_bins, labels=False)
    
    results = []
    
    for i in range(n_bins):
        bin_data = df[df["bin"] == i]
        z_min = bin_data["z"].min()
        z_max = bin_data["z"].max()
        z_mean = bin_data["z"].mean()
        n = len(bin_data)
        
        # Observable: delta_v_axis (Kinematic Axis Hemisphere Diff)
        # Weights: 1/sigma^2
        y = bin_data["delta_v_axis"].values
        # Handle cases where sigma is missing or 0
        sig = bin_data["delta_v_axis_sigma"].values
        w = 1.0 / (sig**2 + 1e-6)
        x = bin_data["x_cmb"].values
        
        # Robust Fit
        a, b, s = robust_huber_fit(x, y, w)
        
        # Also do Normalized version if available
        a_norm = np.nan
        if "delta_v_axis_norm" in bin_data.columns:
            y_norm = bin_data["delta_v_axis_norm"].values
            # Assuming sigma scales roughly, or just use equal weights for norm check
            w_norm = np.ones_like(y_norm) 
            a_norm, _, _ = robust_huber_fit(x, y_norm, w_norm)

        print_status(f"Bin {i+1}/{n_bins}: z=[{z_min:.4f}, {z_max:.4f}], N={n}, a={a:.4f} (km/s)", "INFO")
        
        results.append({
            "bin": i,
            "z_min": z_min,
            "z_max": z_max,
            "z_mean": z_mean,
            "n_galaxies": n,
            "slope_km_s": a,
            "slope_norm": a_norm
        })

    # Save Results
    out_file = out_dir / f"step_2_2_redshift_tomography_{label}.json"
    with open(out_file, "w") as f:
        json.dump(results, f, indent=2)
    print_status(f"Saved {out_file}", "SUCCESS")
    
    # Quick ASCII Plot
    print(f"\n--- {label.upper()} TOMOGRAPHY ---")
    print(f"{'Bin':<5} {'z_range':<15} {'N':<5} {'Slope (km/s)':<15} {'Slope (Norm)':<15}")
    for r in results:
        z_rng = f"{r['z_min']:.3f}-{r['z_max']:.3f}"
        print(f"{r['bin']:<5} {z_rng:<15} {r['n_galaxies']:<5} {r['slope_km_s']:<15.4f} {r['slope_norm']:<15.4f}")
    print("--------------------------------\n")

if __name__ == "__main__":
    main()
