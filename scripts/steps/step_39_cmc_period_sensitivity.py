#!/usr/bin/env python3
"""
Step 39: CMC Period Sensitivity Sweep
=========================================

Test whether the fixed 5 ms MSP period assumption suppresses or inflates
the CMC amplitude prediction. Sweeps across period models from 2 ms to
observed cluster-period distributions.

Outputs:
- step_39_period_sensitivity.json

Author: M. Smawfield
Date: July 2026
"""

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from step_01_cmc_parser import CMCParser

REPO_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = REPO_ROOT / "data" / "cmc"
RESULTS_DIR = REPO_ROOT / "results" / "outputs"
OUT_JSON = RESULTS_DIR / "step_39_period_sensitivity.json"

CLUSTER_NAMES = ["M15", "47_Tuc", "Terzan_5"]
C = 3e8
FIELD_LOG_PDOT = -19.7


def load_cluster(name: str) -> pd.DataFrame:
    d = DATA_DIR / name
    if not d.exists():
        return pd.DataFrame()
    try:
        return CMCParser(d).parse_morepulsars()
    except Exception:
        return pd.DataFrame()


def compute_excess(df: pd.DataFrame, period_s: float) -> dict:
    """Compute predicted excess for a given fixed MSP period."""
    if df.empty or "a_grav_ms2" not in df.columns:
        return {"error": "No data"}

    # Apply default spatial filter (3 rc) to match main analysis
    if "r_pc" in df.columns and "r_core_pc" in df.columns:
        rc = df["r_core_pc"].iloc[0]
        if rc > 0:
            df = df[df["r_pc"] <= 3.0 * rc].copy()

    field_pdot = 10**FIELD_LOG_PDOT
    pdot_accel = np.abs(df["a_grav_ms2"].values) * period_s / C
    pdot_total = field_pdot + pdot_accel
    log_pdot = np.log10(pdot_total)
    cmc_mean = float(np.mean(log_pdot))
    return {
        "cmc_mean_logpdot": cmc_mean,
        "predicted_excess_dex": float(cmc_mean - FIELD_LOG_PDOT),
        "n_pulsars": len(df),
    }


def observed_period_distribution() -> list[float]:
    """Load observed MSP period distribution from Freire catalog if available."""
    freire_path = RESULTS_DIR / "step_02_pulsar_population_controls.json"
    periods = [0.002, 0.003, 0.005, 0.008, 0.010]
    if freire_path.exists():
        try:
            with open(freire_path) as f:
                data = json.load(f)
            # Try to extract period distribution from catalog data if present
            if "catalog" in data:
                cat = data["catalog"]
                if isinstance(cat, list):
                    p_ms = [r.get("P_ms", 5.0) for r in cat if r.get("P_ms")]
                    if p_ms:
                        periods = [p / 1000.0 for p in p_ms]
        except Exception:
            pass
    return periods


def main():
    print("Step 39: CMC Period Sensitivity Sweep")
    results = {"fixed_periods": {}, "period_distribution": {}, "per_cluster": {}}

    # Fixed period sweep
    for period_ms in [2, 3, 5, 8, 10]:
        period_s = period_ms / 1000.0
        per_cluster = {}
        for name in CLUSTER_NAMES:
            df = load_cluster(name)
            res = compute_excess(df, period_s)
            per_cluster[name] = res
        results["fixed_periods"][f"{period_ms}ms"] = {
            "period_s": period_s,
            "per_cluster": per_cluster,
            "combined_excess_dex": float(
                np.mean([v["predicted_excess_dex"] for v in per_cluster.values() if "predicted_excess_dex" in v])
            ) if per_cluster else None,
        }

    # Observed period distribution (Monte Carlo draw)
    obs_periods = observed_period_distribution()
    print(f"  Observed period distribution: N={len(obs_periods)}, median={np.median(obs_periods)*1000:.2f} ms")

    per_cluster_dist = {}
    for name in CLUSTER_NAMES:
        df = load_cluster(name)
        if df.empty or "a_grav_ms2" not in df.columns:
            per_cluster_dist[name] = {"error": "No data"}
            continue
        # Monte Carlo: draw periods for each NS from observed distribution
        rng = np.random.default_rng(42)
        drawn_periods = rng.choice(obs_periods, size=len(df), replace=True)
        field_pdot = 10**FIELD_LOG_PDOT
        pdot_accel = np.abs(df["a_grav_ms2"].values) * drawn_periods / C
        log_pdot = np.log10(field_pdot + pdot_accel)
        per_cluster_dist[name] = {
            "predicted_excess_dex": float(np.mean(log_pdot) - FIELD_LOG_PDOT),
            "std_excess_dex": float(np.std(log_pdot)),
            "n_pulsars": len(df),
        }

    results["period_distribution"] = {
        "source": "observed_MSP_distribution_or_default",
        "median_ms": float(np.median(obs_periods) * 1000),
        "per_cluster": per_cluster_dist,
        "combined_excess_dex": float(
            np.mean([v["predicted_excess_dex"] for v in per_cluster_dist.values() if isinstance(v, dict) and "predicted_excess_dex" in v])
        ) if per_cluster_dist else None,
    }

    # Per-cluster fixed 5 ms baseline (for direct comparison)
    for name in CLUSTER_NAMES:
        df = load_cluster(name)
        results["per_cluster"][name] = compute_excess(df, 0.005)

    with open(OUT_JSON, "w") as f:
        json.dump(results, f, indent=2)

    print(f"\nBaseline (5 ms combined excess): {results['fixed_periods']['5ms']['combined_excess_dex']:.3f} dex")
    print(f"Observed-distribution combined excess: {results['period_distribution']['combined_excess_dex']:.3f} dex")
    print(f"Output written to {OUT_JSON}")


if __name__ == "__main__":
    main()
