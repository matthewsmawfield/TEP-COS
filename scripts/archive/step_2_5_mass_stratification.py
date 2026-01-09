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

def load_dapall_dispersion(dapall_path: str) -> Dict[str, float]:
    print_status(f"Loading dispersion from {dapall_path}...", "PROCESS")
    meta = {}
    with fits.open(dapall_path) as hdul:
        data = None
        for hdu in hdul:
            if isinstance(hdu, fits.BinTableHDU) and "PLATEIFU" in hdu.columns.names:
                data = hdu.data
                break
        
        if data is None:
            raise RuntimeError("Could not find table with PLATEIFU in dapall")
            
        names = data.columns.names
        # STELLAR_SIGMA_1RE is the flux-weighted stellar velocity dispersion within 1 Re
        col_name = "STELLAR_SIGMA_1RE"
        if col_name not in names:
            print_status(f"{col_name} not found, checking for alternatives...", "WARNING")
            return {}

        for row in data:
            pid = str(row["PLATEIFU"]).strip()
            sigma = float(row[col_name])
            if np.isfinite(sigma) and sigma > 0:
                meta[pid] = sigma
                
    return meta

def main():
    parser = argparse.ArgumentParser(description="TEP-COS Step 2.5 - Mass (Dispersion) Stratification")
    parser.add_argument("--stellar-csv", required=True)
    parser.add_argument("--dapall", required=True)
    parser.add_argument("--output-dir", default="results/outputs")
    args = parser.parse_args()

    logger = TEPLogger("step_2_5_mass_stratification", log_file_path=PROJECT_ROOT / "logs" / "step_2_5_mass_stratification.log")
    set_step_logger(logger)
    
    # 1. Load Data
    sigma_map = load_dapall_dispersion(args.dapall)
    df = pd.read_csv(args.stellar_csv)
    
    df["sigma_star"] = df["plateifu"].apply(lambda x: sigma_map.get(str(x).strip(), np.nan))
    df = df.dropna(subset=["sigma_star"])
    
    print_status(f"Loaded {len(df)} galaxies with dispersion data.", "INFO")
    
    # 2. Stratify by Dispersion (Low Mass vs High Mass)
    # Median split
    median_sig = df["sigma_star"].median()
    low_mass = df[df["sigma_star"] <= median_sig]
    high_mass = df[df["sigma_star"] > median_sig]
    
    results = {}
    
    for subset, label in [(low_mass, "Low Mass (Low Dispersion)"), (high_mass, "High Mass (High Dispersion)")]:
        y = subset["delta_v_axis"].values
        sig = subset["delta_v_axis_sigma"].values
        w = 1.0 / (sig**2 + 1e-6)
        x = subset["x_cmb"].values
        
        a, b, s = robust_huber_fit(x, y, w)
        
        print_status(f"{label}: N={len(subset)}, Median Sigma={subset['sigma_star'].median():.1f} km/s, Slope a={a:.3f} km/s", "INFO")
        results[label] = {
            "n": len(subset),
            "median_sigma": subset['sigma_star'].median(),
            "slope": a
        }

    # 3. Quartile check for trend
    print_status("Checking Quartiles...", "PROCESS")
    df["q"] = pd.qcut(df["sigma_star"], 4, labels=False)
    quartiles = []
    for i in range(4):
        subset = df[df["q"] == i]
        y = subset["delta_v_axis"].values
        sig = subset["delta_v_axis_sigma"].values
        w = 1.0 / (sig**2 + 1e-6)
        x = subset["x_cmb"].values
        a, _, _ = robust_huber_fit(x, y, w)
        quartiles.append({
            "quartile": i,
            "range_sigma": f"{subset['sigma_star'].min():.1f}-{subset['sigma_star'].max():.1f}",
            "slope": a
        })
        print_status(f"Q{i+1}: Sigma=[{subset['sigma_star'].min():.1f}, {subset['sigma_star'].max():.1f}], Slope={a:.3f}", "INFO")

    # Save
    out_data = {
        "median_split": results,
        "quartiles": quartiles
    }
    out_file = Path(args.output_dir) / "step_2_5_mass_stratification.json"
    with open(out_file, "w") as f:
        json.dump(out_data, f, indent=2)
    print_status(f"Saved {out_file}", "SUCCESS")

if __name__ == "__main__":
    main()
