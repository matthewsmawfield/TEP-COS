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

def vec_from_radec(ra_deg, dec_deg):
    ra = np.radians(ra_deg)
    dec = np.radians(dec_deg)
    return np.array([np.cos(dec)*np.cos(ra), np.cos(dec)*np.sin(ra), np.sin(dec)])

def fit_axis_grid_search(df, n_grid=1000):
    # Quick grid search for best axis in this bin
    # Generate random axes
    rng = np.random.default_rng(42)
    best_a = 0.0
    best_axis = None
    
    # Precompute vectors
    vecs = np.array([vec_from_radec(r, d) for r, d in zip(df["ra_deg"], df["dec_deg"])])
    y = df["delta_v_axis"].values
    sig = df["delta_v_axis_sigma"].values
    w = 1.0 / (sig**2 + 1e-6)
    
    # Try n_grid random directions
    for _ in range(n_grid):
        u = rng.uniform(0, 1)
        v = rng.uniform(0, 1)
        phi = 2 * np.pi * u
        theta = np.arccos(2 * v - 1)
        axis = np.array([np.sin(theta)*np.cos(phi), np.sin(theta)*np.sin(phi), np.cos(theta)])
        
        x_proj = vecs @ axis
        a, _, _ = robust_huber_fit(x_proj, y, w)
        
        if abs(a) > best_a:
            best_a = abs(a)
            best_axis = axis * np.sign(a) # Point in direction of correlation
            
    if best_axis is not None:
        ra = np.degrees(np.arctan2(best_axis[1], best_axis[0])) % 360
        dec = np.degrees(np.arcsin(best_axis[2]))
        return ra, dec, best_a
    return 0.0, 0.0, 0.0

def main():
    parser = argparse.ArgumentParser(description="TEP-COS Step 2.3 - Deep Tomography & Paired Analysis")
    parser.add_argument("--stellar-csv", required=True)
    parser.add_argument("--gas-csv", required=True)
    parser.add_argument("--dapall", required=True)
    parser.add_argument("--n-bins", type=int, default=4)
    parser.add_argument("--output-dir", default="results/outputs")
    args = parser.parse_args()

    logger = TEPLogger("step_2_3_deep_tomography", log_file_path=PROJECT_ROOT / "logs" / "step_2_3_deep_tomography.log")
    set_step_logger(logger)
    
    # 1. Load Data
    print_status("Loading datasets...", "PROCESS")
    z_map = {}
    with fits.open(args.dapall) as hdul:
        data = None
        for hdu in hdul:
            if isinstance(hdu, fits.BinTableHDU) and "PLATEIFU" in hdu.columns.names:
                data = hdu.data
                break
        for row in data:
            pid = str(row["PLATEIFU"]).strip()
            z_val = row["NSA_Z"] if "NSA_Z" in row.array.names else row["Z"]
            z_map[pid] = float(z_val)

    df_star = pd.read_csv(args.stellar_csv)
    df_gas = pd.read_csv(args.gas_csv)
    
    df_star["z"] = df_star["plateifu"].apply(lambda x: z_map.get(str(x).strip(), np.nan))
    df_gas["z"] = df_gas["plateifu"].apply(lambda x: z_map.get(str(x).strip(), np.nan))
    
    df_star = df_star.dropna(subset=["z"])
    df_gas = df_gas.dropna(subset=["z"])
    
    # 2. Per-Bin Axis Stability (Stellar)
    print_status("Analyzing Axis Stability (Stellar)...", "PROCESS")
    df_star["bin"] = pd.qcut(df_star["z"], args.n_bins, labels=False)
    
    axis_results = []
    for i in range(args.n_bins):
        bin_data = df_star[df_star["bin"] == i]
        ra, dec, strength = fit_axis_grid_search(bin_data)
        z_mean = bin_data["z"].mean()
        print_status(f"Bin {i}: z={z_mean:.3f}, Axis=({ra:.1f}, {dec:.1f}), Strength={strength:.2f}", "INFO")
        axis_results.append({"bin": i, "z": z_mean, "ra": ra, "dec": dec, "strength": strength})
        
    # 3. Paired Analysis (Star - Gas)
    print_status("Performing Paired Star-Gas Analysis...", "PROCESS")
    # Merge on plateifu
    df_pair = pd.merge(df_star, df_gas, on="plateifu", suffixes=("_star", "_gas"))
    
    # Compute Difference Observable
    # We expect Star < 0 and Gas > 0 (or less negative)
    # So Difference (Gas - Star) should be POSITIVE and strongly correlated with dipole
    df_pair["diff_v"] = df_pair["delta_v_axis_gas"] - df_pair["delta_v_axis_star"]
    
    # Combine sigmas
    df_pair["diff_sigma"] = np.sqrt(df_pair["delta_v_axis_sigma_star"]**2 + df_pair["delta_v_axis_sigma_gas"]**2)
    
    # Fit Dipole to Difference
    y = df_pair["diff_v"].values
    sig = df_pair["diff_sigma"].values
    w = 1.0 / (sig**2 + 1e-6)
    x = df_pair["x_cmb_star"].values # x_cmb is same for both
    
    a_diff, b_diff, s_diff = robust_huber_fit(x, y, w)
    
    print_status(f"Paired Difference (Gas-Star): N={len(df_pair)}", "INFO")
    print_status(f"Difference Slope: a = {a_diff:.4f} km/s", "SUCCESS")
    
    # Save Results
    out_file = Path(args.output_dir) / "step_2_3_deep_tomography.json"
    with open(out_file, "w") as f:
        json.dump({
            "axis_stability": axis_results,
            "paired_diff": {
                "n": len(df_pair),
                "slope": a_diff,
                "intercept": b_diff,
                "scale": s_diff
            }
        }, f, indent=2)

if __name__ == "__main__":
    main()
