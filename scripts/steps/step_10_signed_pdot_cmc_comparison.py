#!/usr/bin/env python3
"""
Step 10: Signed Pdot CMC Comparison
=========================================

Turn the signed Ṗ/P analysis into a stronger discriminator by comparing
observed signed distributions to CMC predictions.

Tests:
1. Observed signed Ṗ/P distribution by cluster
2. CMC signed Ṗ/P distribution (from acceleration direction)
3. Predicted vs observed negative fraction
4. Predicted vs observed extreme-tail fraction

Outputs:
- step_10_signed_pdot_cmc.json

Author: M. Smawfield
Date: July 2026
"""

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).parent))
from step_01_cmc_parser import CMCParser

REPO_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = REPO_ROOT / "data" / "cmc"
RESULTS_DIR = REPO_ROOT / "results" / "outputs"
OUT_JSON = RESULTS_DIR / "step_10_signed_pdot_cmc.json"

CLUSTER_NAMES = ["M15", "47_Tuc", "Terzan_5"]
MSP_PERIOD_S = 0.005
C = 3e8


def load_cluster(name: str) -> pd.DataFrame:
    d = DATA_DIR / name
    if not d.exists():
        return pd.DataFrame()
    try:
        return CMCParser(d).parse_morepulsars()
    except Exception:
        return pd.DataFrame()


def compute_signed_pdot(df: pd.DataFrame, period_s: float = MSP_PERIOD_S) -> np.ndarray:
    """Compute signed log(Pdot) from CMC line-of-sight acceleration."""
    if df.empty or "a_grav_ms2" not in df.columns:
        return np.array([])
    # Intrinsic field reference
    field_pdot = 10**(-19.7)
    # Signed acceleration contribution (sign = direction of a_los)
    a_los = df["a_grav_ms2"].values
    pdot_accel = a_los * period_s / C
    pdot_total = field_pdot + pdot_accel
    # Return signed log10(Pdot) — negative when pdot_total < 0
    with np.errstate(divide="ignore", invalid="ignore"):
        signed_log = np.log10(np.abs(pdot_total)) * np.sign(pdot_total)
    return signed_log


def cmc_signed_distributions() -> dict:
    """Compute signed Pdot distributions for each CMC cluster."""
    out = {}
    for name in CLUSTER_NAMES:
        df = load_cluster(name)
        sp = compute_signed_pdot(df)
        if len(sp) == 0:
            out[name] = {"error": "No data"}
            continue
        # Apply spatial filter
        if "r_pc" in df.columns and "r_core_pc" in df.columns:
            rc = df["r_core_pc"].iloc[0]
            if rc > 0:
                mask = df["r_pc"].values <= 3.0 * rc
                sp = sp[mask]
        out[name] = {
            "n": len(sp),
            "mean_signed_logpdot": float(np.mean(sp)),
            "std_signed_logpdot": float(np.std(sp)),
            "frac_negative": float(np.mean(sp < 0)),
            "frac_extreme_tail": float(np.mean(np.abs(sp - np.mean(sp)) > 2 * np.std(sp))) if np.std(sp) > 0 else 0.0,
        }
    return out


def load_observed_signed() -> dict:
    """Load observed signed Pdot results from step_08."""
    path = RESULTS_DIR / "step_08_signed_pdot_analysis.json"
    if path.exists():
        with open(path) as f:
            return json.load(f)
    return {"error": "step_08_signed_pdot_analysis.json not found"}


def compare_fractions(cmc: dict, observed: dict) -> dict:
    """Compare predicted vs observed negative fractions."""
    # Extract observed negative fraction
    obs_neg = observed.get("overall", {}).get("frac_negative", None)
    if obs_neg is None:
        return {"status": "observed_data_missing"}

    cmc_fracs = [v["frac_negative"] for v in cmc.values() if isinstance(v, dict) and "frac_negative" in v]
    if not cmc_fracs:
        return {"status": "cmc_data_missing"}

    cmc_mean = float(np.mean(cmc_fracs))
    cmc_std = float(np.std(cmc_fracs))

    return {
        "observed_frac_negative": obs_neg,
        "cmc_mean_frac_negative": cmc_mean,
        "cmc_std_frac_negative": cmc_std,
        "difference": float(obs_neg - cmc_mean),
        "interpretation": (
            "If CMC predicts a much higher negative fraction than observed, "
            "the signed distribution is more compressed than Newtonian dynamics allows, "
            "consistent with TEP saturation."
        ),
    }


def main():
    print("Step 10: Signed Pdot CMC Comparison")

    cmc_signed = cmc_signed_distributions()
    observed = load_observed_signed()
    comparison = compare_fractions(cmc_signed, observed)

    result = {
        "cmc_signed_distributions": cmc_signed,
        "observed_signed_summary": observed,
        "negative_fraction_comparison": comparison,
        "notes": [
            "CMC signed Pdot uses line-of-sight acceleration sign from initial positions.",
            "Present-day snapshot epochs may differ; full validation requires output.window.snapshot.h5.",
            "If observed negative fraction is lower than CMC prediction, the distribution is compressed (TEP signature).",
            "If observed negative fraction is higher, CMC may underpredict acceleration variance.",
        ],
    }

    with open(OUT_JSON, "w") as f:
        json.dump(result, f, indent=2)

    print(f"Comparison complete. Output: {OUT_JSON}")
    if "cmc_mean_frac_negative" in comparison:
        print(f"  Observed negative fraction: {comparison['observed_frac_negative']:.3f}")
        print(f"  CMC mean negative fraction: {comparison['cmc_mean_frac_negative']:.3f} ± {comparison['cmc_std_frac_negative']:.3f}")


if __name__ == "__main__":
    main()
