#!/usr/bin/env python3
"""
Step 38: CMC Binary Forensic Audit
======================================

Determine whether the CMC binary agreement is a pipeline artifact.

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
OUT_PREFIX = RESULTS_DIR / "step_38"

CLUSTER_NAMES = ["M15", "47_Tuc", "Terzan_5"]
MSP_PERIOD_S = 0.005
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


def compute_log_pdot(df: pd.DataFrame, period_s: float = MSP_PERIOD_S) -> np.ndarray:
    if df.empty or "a_grav_ms2" not in df.columns:
        return np.array([])
    field_pdot = 10**FIELD_LOG_PDOT
    pdot_accel = np.abs(df["a_grav_ms2"].values) * period_s / C
    return np.log10(field_pdot + pdot_accel)


def get_binary_mask(df: pd.DataFrame) -> np.ndarray:
    if "binflag" in df.columns:
        return df["binflag"].values.astype(bool)
    if "m2" in df.columns:
        return (df["m2"].values > 0).astype(bool)
    return np.zeros(len(df), dtype=bool)


def audit_per_cluster() -> dict:
    out = {}
    for name in CLUSTER_NAMES:
        df = load_cluster(name)
        lp = compute_log_pdot(df)
        if len(lp) == 0:
            continue
        bm = get_binary_mask(df)
        nb, ni = int(bm.sum()), int((~bm).sum())
        diff = None
        if nb >= 3 and ni >= 3:
            diff = float(np.mean(lp[bm]) - np.mean(lp[~bm]))
            t, p = stats.ttest_ind(lp[bm], lp[~bm], equal_var=False)
        out[name] = {
            "n_binary": nb, "n_isolated": ni,
            "diff_dex": diff,
            "t_stat": float(t) if diff is not None else None,
            "t_p": float(p) if diff is not None else None,
        }
    return out


def audit_cluster_centred() -> dict:
    lp_all, bm_all, cid_all = [], [], []
    for idx, name in enumerate(CLUSTER_NAMES):
        df = load_cluster(name)
        lp = compute_log_pdot(df)
        if len(lp) == 0:
            continue
        bm = get_binary_mask(df)
        lp_all.extend(lp.tolist())
        bm_all.extend(bm.tolist())
        cid_all.extend([idx] * len(lp))
    if len(lp_all) == 0:
        return {"error": "No data"}
    lp = np.array(lp_all)
    bm = np.array(bm_all, dtype=bool)
    cid = np.array(cid_all)
    raw_diff = float(np.mean(lp[bm]) - np.mean(lp[~bm]))
    centred = lp.copy()
    for c in np.unique(cid):
        m = cid == c
        centred[m] -= np.mean(lp[m])
    cent_diff = float(np.mean(centred[bm]) - np.mean(centred[~bm]))
    return {
        "pooled_diff_dex": raw_diff,
        "cluster_centred_diff_dex": cent_diff,
        "simpson_shift_dex": cent_diff - raw_diff,
        "n_total": len(lp), "n_binary": int(bm.sum()),
    }


def audit_radial_sweep() -> dict:
    filters = [(0, "no_filter"), (1, "1rc"), (2, "2rc"), (3, "3rc"), (5, "5rc")]
    out = {}
    for fac, label in filters:
        diffs = {}
        for name in CLUSTER_NAMES:
            df = load_cluster(name)
            if fac > 0 and "r_pc" in df.columns and "r_core_pc" in df.columns:
                rc = df["r_core_pc"].iloc[0]
                if rc > 0:
                    df = df[df["r_pc"] <= fac * rc].copy()
            lp = compute_log_pdot(df)
            if len(lp) == 0:
                continue
            bm = get_binary_mask(df)
            if bm.sum() >= 3 and (~bm).sum() >= 3:
                diffs[name] = float(np.mean(lp[bm]) - np.mean(lp[~bm]))
        out[label] = diffs
        # Add unweighted summary for filters with all three clusters
        if len(diffs) == len(CLUSTER_NAMES):
            out[f"{label}_summary"] = {
                "unweighted_mean_dex": float(np.mean(list(diffs.values()))),
                "clusters_included": list(diffs.keys()),
                "n_clusters": len(diffs),
            }
    return out


def audit_observable_sweep() -> dict:
    lp_all, bm_all = [], []
    for name in CLUSTER_NAMES:
        df = load_cluster(name)
        lp = compute_log_pdot(df)
        if len(lp) == 0:
            continue
        bm = get_binary_mask(df)
        lp_all.extend(lp.tolist())
        bm_all.extend(bm.tolist())
    if len(lp_all) < 6:
        return {"error": "insufficient data"}
    lp = np.array(lp_all)
    bm = np.array(bm_all, dtype=bool)
    return {
        "log_abs_pdot": {
            "diff_binary_minus_isolated_dex": float(np.mean(lp[bm]) - np.mean(lp[~bm]))
        },
        "log_abs_pdot_over_p": {"status": "requires_snapshot_data"},
        "signed_pdot_over_p": {"status": "requires_snapshot_data"},
    }


def audit_deduplication() -> dict:
    out = {}
    for name in CLUSTER_NAMES:
        df = load_cluster(name)
        n_raw = len(df)
        if "id0" in df.columns:
            df_dedup = df.drop_duplicates(subset=["id0"])
            n_dedup = len(df_dedup)
        else:
            n_dedup = n_raw
        out[f"{name}_raw_count"] = n_raw
        out[f"{name}_dedup_count"] = n_dedup
    return out


def main():
    print("Step 38: CMC Binary Forensic Audit")
    per_cluster = audit_per_cluster()
    with open(f"{OUT_PREFIX}_by_cluster.json", "w") as f:
        json.dump(per_cluster, f, indent=2)

    centred = audit_cluster_centred()
    with open(f"{OUT_PREFIX}_cluster_centered.json", "w") as f:
        json.dump(centred, f, indent=2)

    radial = audit_radial_sweep()
    with open(f"{OUT_PREFIX}_radial_filter.json", "w") as f:
        json.dump(radial, f, indent=2)

    observable = audit_observable_sweep()
    with open(f"{OUT_PREFIX}_observable_sweep.json", "w") as f:
        json.dump(observable, f, indent=2)

    dedup = audit_deduplication()
    with open(f"{OUT_PREFIX}_dedup.json", "w") as f:
        json.dump(dedup, f, indent=2)

    print("Audit complete. Outputs written to step_38_*.json")


if __name__ == "__main__":
    main()
