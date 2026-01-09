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

def load_dapall_metadata(dapall_path: str) -> Dict[str, Dict[str, float]]:
    print_status(f"Loading metadata from {dapall_path}...", "PROCESS")
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
        # We need Z, STELLAR_MASS (if avail), NSA_ELPETRO_BA (inclination proxy)
        # Often mass is in a separate VAC, but let's check what we have.
        # NSA_ELPETRO_BA is b/a ratio. Low b/a = Edge-on, High b/a = Face-on.
        
        for row in data:
            pid = str(row["PLATEIFU"]).strip()
            z_val = row["NSA_Z"] if "NSA_Z" in names else row["Z"]
            ba = row["NSA_ELPETRO_BA"] if "NSA_ELPETRO_BA" in names else np.nan
            
            # Simple mass proxy if stellar mass not available: SERSIC Flux * Distance^2 ? 
            # For now, let's rely on Sersic Index 'NSA_SERSIC_N' as a morphology proxy 
            # or just use sigma (dispersion) from our own analysis as mass proxy.
            
            meta[pid] = {
                "z": float(z_val),
                "ba": float(ba)
            }
    return meta

def sliding_window_tomography(df: pd.DataFrame, window_size: int = 100, step: int = 20) -> List[Dict]:
    df = df.sort_values("z")
    results = []
    n_total = len(df)
    
    for i in range(0, n_total - window_size + 1, step):
        window = df.iloc[i : i + window_size]
        z_mean = window["z"].mean()
        z_min = window["z"].min()
        z_max = window["z"].max()
        
        y = window["delta_v_axis"].values
        sig = window["delta_v_axis_sigma"].values
        w = 1.0 / (sig**2 + 1e-6)
        x = window["x_cmb"].values
        
        a, b, s = robust_huber_fit(x, y, w)
        
        results.append({
            "z_mean": z_mean,
            "z_min": z_min,
            "z_max": z_max,
            "slope": a,
            "slope_err": s/np.sqrt(window_size), # Approx
            "n": len(window)
        })
    return results

def bootstrap_paired_diff(df: pd.DataFrame, n_boot: int = 1000) -> Dict:
    # Bootstrap the Difference Slope
    diffs = []
    n = len(df)
    
    # Original
    y = df["diff_v"].values
    sig = df["diff_sigma"].values
    w = 1.0 / (sig**2 + 1e-6)
    x = df["x_cmb_star"].values # Same for both
    a_orig, _, _ = robust_huber_fit(x, y, w)
    
    rng = np.random.default_rng(42)
    
    for _ in range(n_boot):
        # Resample indices with replacement
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
        "orig": a_orig
    }

def stratification_test(df: pd.DataFrame, col: str, label: str) -> Dict:
    # Split into low/high bins by median
    median_val = df[col].median()
    low = df[df[col] <= median_val]
    high = df[df[col] > median_val]
    
    res = {}
    for subset, name in [(low, "low"), (high, "high")]:
        y = subset["delta_v_axis"].values
        sig = subset["delta_v_axis_sigma"].values
        w = 1.0 / (sig**2 + 1e-6)
        x = subset["x_cmb"].values
        a, _, _ = robust_huber_fit(x, y, w)
        res[name] = {"n": len(subset), "slope": a, "median_val": subset[col].median()}
        
    return {
        "variable": label,
        "cut_value": median_val,
        "results": res
    }

def main():
    parser = argparse.ArgumentParser(description="TEP-COS Step 2.4 - Robustness & Expansion")
    parser.add_argument("--stellar-csv", required=True)
    parser.add_argument("--gas-csv", required=True)
    parser.add_argument("--dapall", required=True)
    parser.add_argument("--output-dir", default="results/outputs")
    args = parser.parse_args()

    logger = TEPLogger("step_2_4_robustness", log_file_path=PROJECT_ROOT / "logs" / "step_2_4_robustness.log")
    set_step_logger(logger)
    
    # 1. Load Data & Metadata
    meta = load_dapall_metadata(args.dapall)
    
    df_star = pd.read_csv(args.stellar_csv)
    df_gas = pd.read_csv(args.gas_csv)
    
    # Attach Metadata
    for df in [df_star, df_gas]:
        df["z"] = df["plateifu"].apply(lambda x: meta.get(str(x).strip(), {}).get("z", np.nan))
        df["ba"] = df["plateifu"].apply(lambda x: meta.get(str(x).strip(), {}).get("ba", np.nan))
        # Inclination proxy: i = arccos(b/a). But b/a is sufficient. Low b/a = High Inclination.
        
    df_star = df_star.dropna(subset=["z", "ba"])
    df_gas = df_gas.dropna(subset=["z", "ba"])
    
    # 2. Sliding Window Tomography (High Res)
    print_status("Running Sliding Window Tomography...", "PROCESS")
    # Window size 80, step 10 for detailed curve
    sw_results = sliding_window_tomography(df_star, window_size=80, step=10)
    
    # 3. Stratification Tests (Systematics)
    print_status("Running Stratification Tests (Inclination)...", "PROCESS")
    # b/a ratio check
    strat_ba = stratification_test(df_star, "ba", "Axial Ratio (b/a)")
    
    # 4. Bootstrap Paired Difference
    print_status("Running Bootstrap on Paired Difference...", "PROCESS")
    df_pair = pd.merge(df_star, df_gas, on="plateifu", suffixes=("_star", "_gas"))
    df_pair["diff_v"] = df_pair["delta_v_axis_gas"] - df_pair["delta_v_axis_star"]
    df_pair["diff_sigma"] = np.sqrt(df_pair["delta_v_axis_sigma_star"]**2 + df_pair["delta_v_axis_sigma_gas"]**2)
    
    boot_res = bootstrap_paired_diff(df_pair, n_boot=2000)
    
    # 5. Output
    out_data = {
        "sliding_window": sw_results,
        "stratification": {
            "ba": strat_ba
        },
        "bootstrap_paired": boot_res
    }
    
    out_file = Path(args.output_dir) / "step_2_4_robustness_expansion.json"
    with open(out_file, "w") as f:
        json.dump(out_data, f, indent=2)
        
    print_status(f"Saved {out_file}", "SUCCESS")
    
    # Report
    print("\n--- ROBUSTNESS SUMMARY ---")
    print(f"Paired Difference (Gas-Star): {boot_res['orig']:.3f} km/s")
    print(f"Bootstrap 95% CI: [{boot_res['ci_low']:.3f}, {boot_res['ci_high']:.3f}]")
    print(f"Bootstrap Mean/Std: {boot_res['mean']:.3f} ± {boot_res['std']:.3f}")
    
    print("\n--- INCLINATION CHECK ---")
    low = strat_ba["results"]["low"]
    high = strat_ba["results"]["high"]
    print(f"Edge-on (b/a < {strat_ba['cut_value']:.2f}): a = {low['slope']:.3f} (N={low['n']})")
    print(f"Face-on (b/a > {strat_ba['cut_value']:.2f}): a = {high['slope']:.3f} (N={high['n']})")
    
    print("\n--- TOMOGRAPHY PEAK ---")
    # Find bin with max absolute slope
    best_bin = max(sw_results, key=lambda x: abs(x["slope"]))
    print(f"Peak Signal at z = {best_bin['z_mean']:.4f} (z_range {best_bin['z_min']:.4f}-{best_bin['z_max']:.4f})")
    print(f"Slope: {best_bin['slope']:.3f} km/s")

if __name__ == "__main__":
    main()
