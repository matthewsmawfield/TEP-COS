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
            return {}
        names = data.columns.names
        col_name = "STELLAR_SIGMA_1RE"
        if col_name not in names:
            return {}
        for row in data:
            pid = str(row["PLATEIFU"]).strip()
            sigma = float(row[col_name])
            if np.isfinite(sigma) and sigma > 0:
                meta[pid] = sigma
    return meta

def bootstrap_paired_diff(df: pd.DataFrame, n_boot: int = 2000) -> Dict:
    diffs = []
    n = len(df)
    
    y = df["diff_v"].values
    sig = df["diff_sigma"].values
    w = 1.0 / (sig**2 + 1e-6)
    x = df["x_cmb_star"].values
    a_orig, _, _ = robust_huber_fit(x, y, w)
    
    rng = np.random.default_rng(42)
    
    for _ in range(n_boot):
        idx = rng.choice(n, n, replace=True)
        x_b = x[idx]
        y_b = y[idx]
        w_b = w[idx]
        try:
            a, _, _ = robust_huber_fit(x_b, y_b, w_b)
            diffs.append(a)
        except:
            continue
            
    diffs = np.array(diffs)
    return {
        "mean": float(np.mean(diffs)),
        "std": float(np.std(diffs)),
        "ci_low": float(np.percentile(diffs, 2.5)),
        "ci_high": float(np.percentile(diffs, 97.5)),
        "orig": a_orig,
        "n_samples": n
    }

def main():
    parser = argparse.ArgumentParser(description="TEP-COS Step 2.6 - Refined Paired Analysis (Screening Control)")
    parser.add_argument("--stellar-csv", required=True)
    parser.add_argument("--gas-csv", required=True)
    parser.add_argument("--dapall", required=True)
    parser.add_argument("--output-dir", default="results/outputs")
    args = parser.parse_args()

    logger = TEPLogger("step_2_6_refined_paired", log_file_path=PROJECT_ROOT / "logs" / "step_2_6_refined_paired.log")
    set_step_logger(logger)
    
    # 1. Load Data
    sigma_map = load_dapall_dispersion(args.dapall)
    df_star = pd.read_csv(args.stellar_csv)
    df_gas = pd.read_csv(args.gas_csv)
    
    df_pair = pd.merge(df_star, df_gas, on="plateifu", suffixes=("_star", "_gas"))
    df_pair["sigma_star"] = df_pair["plateifu"].apply(lambda x: sigma_map.get(str(x).strip(), np.nan))
    df_pair = df_pair.dropna(subset=["sigma_star"])
    
    df_pair["diff_v"] = df_pair["delta_v_axis_gas"] - df_pair["delta_v_axis_star"]
    df_pair["diff_sigma"] = np.sqrt(df_pair["delta_v_axis_sigma_star"]**2 + df_pair["delta_v_axis_sigma_gas"]**2)
    
    print_status(f"Total Paired Galaxies: {len(df_pair)}", "INFO")
    
    # 2. Define Subsets based on Screening (Sigma)
    # We saw earlier that High Sigma (> 165 km/s) is screened (slope ~ 0).
    # So we should look for the Two-Fluid effect in the Unscreened (Low/Med Sigma) population.
    
    median_sig = df_pair["sigma_star"].median()
    cutoff_sig = 160.0 # From previous analysis Q3/Q4 boundary approx
    
    subsets = {
        "All": df_pair,
        "Unscreened (Sigma < 160 km/s)": df_pair[df_pair["sigma_star"] < cutoff_sig],
        "Screened (Sigma >= 160 km/s)": df_pair[df_pair["sigma_star"] >= cutoff_sig],
        "Deep Unscreened (Sigma < 100 km/s)": df_pair[df_pair["sigma_star"] < 100.0]
    }
    
    results = {}
    
    print("\n--- REFINED PAIRED ANALYSIS ---")
    for name, subset in subsets.items():
        if len(subset) < 10:
            print_status(f"Skipping {name}: N={len(subset)} too small", "WARNING")
            continue
            
        boot = bootstrap_paired_diff(subset, n_boot=2000)
        results[name] = boot
        
        signif = ""
        if boot['ci_low'] > 0: signif = "** SIGNIFICANT **"
        if boot['ci_high'] < 0: signif = "** SIGNIFICANT (NEGATIVE) **"
        
        print(f"{name:<35} | N={boot['n_samples']:<4} | Slope={boot['orig']:>7.3f} | CI=[{boot['ci_low']:>6.2f}, {boot['ci_high']:>6.2f}] {signif}")

    # 3. Save
    out_file = Path(args.output_dir) / "step_2_6_refined_paired.json"
    with open(out_file, "w") as f:
        json.dump(results, f, indent=2)
    print_status(f"Saved {out_file}", "SUCCESS")

if __name__ == "__main__":
    main()
