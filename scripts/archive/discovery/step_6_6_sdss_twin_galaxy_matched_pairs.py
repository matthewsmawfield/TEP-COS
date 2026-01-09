#!/usr/bin/env python3

import json
import os
import time
from dataclasses import dataclass

import numpy as np
import pandas as pd
import requests
from scipy import stats
from sklearn.neighbors import NearestNeighbors
import matplotlib.pyplot as plt


SDSS_URL = "https://skyserver.sdss.org/dr18/SkyServerWS/SearchTools/SqlSearch"

_C_KM_S = 299792.458
_H0_KM_S_MPC = 70.0


def query_sdss(sql: str, max_retries: int = 3) -> pd.DataFrame | None:
    for attempt in range(max_retries):
        try:
            response = requests.get(
                SDSS_URL,
                params={"cmd": sql, "format": "json"},
                timeout=300,
            )
            if response.status_code == 200:
                data = response.json()
                if data and len(data) > 0 and "Rows" in data[0]:
                    return pd.DataFrame(data[0]["Rows"])
            else:
                print(f"  HTTP {response.status_code}")
        except requests.exceptions.Timeout:
            print(f"  Timeout (attempt {attempt + 1})")
        except Exception as e:
            print(f"  Error: {e}")
        if attempt < max_retries - 1:
            time.sleep(2)
    return None


def download_base_sample(z_min: float, z_max: float, limit: int = 120000) -> pd.DataFrame | None:
    sql = f"""
    SELECT TOP {limit}
        i.specobjid,
        g.ra, g.dec,
        g.z as redshift,
        g.v_disp as veldisp,
        g.v_disp_err as veldisp_err,
        i.lick_mgb as mgb,
        i.lick_fe5270 as fe5270,
        i.lick_fe5335 as fe5335,
        i.d4000_n as d4000,
        i.lick_hb as hbeta,
        e.lgm_tot_p50 as log_mass,
        p.petroR50_r as petroR50_r_arcsec
    FROM galSpecIndx i
    JOIN galSpecInfo g ON i.specobjid = g.specobjid
    JOIN galSpecExtra e ON i.specobjid = e.specobjid
    JOIN SpecObj s ON i.specobjid = s.specobjid
    JOIN PhotoObjAll p ON s.bestObjID = p.objID
    WHERE g.reliable = 1
        AND g.z BETWEEN {z_min} AND {z_max}
        AND g.z_err < 0.001
        AND g.v_disp > 60 AND g.v_disp < 350
        AND g.v_disp_err > 0 AND g.v_disp_err < 50
        AND i.lick_mgb > 0.5 AND i.lick_mgb < 8
        AND i.lick_fe5270 > 0.5 AND i.lick_fe5270 < 5
        AND i.lick_fe5335 > 0.5 AND i.lick_fe5335 < 5
        AND i.d4000_n > 1.0 AND i.d4000_n < 2.5
        AND i.lick_hb > 0 AND i.lick_hb < 6
        AND e.lgm_tot_p50 > 8.5 AND e.lgm_tot_p50 < 12.5
        AND p.petroR50_r > 0.2 AND p.petroR50_r < 50
    ORDER BY g.z
    """
    return query_sdss(sql)


def prepare(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    for c in [
        "ra",
        "dec",
        "redshift",
        "veldisp",
        "veldisp_err",
        "mgb",
        "fe5270",
        "fe5335",
        "d4000",
        "hbeta",
        "log_mass",
        "petroR50_r_arcsec",
    ]:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    valid = (
        np.isfinite(df["redshift"])
        & np.isfinite(df["veldisp"])
        & np.isfinite(df["d4000"])
        & np.isfinite(df["hbeta"])
        & np.isfinite(df["mgb"])
        & np.isfinite(df["fe5270"])
        & np.isfinite(df["fe5335"])
        & np.isfinite(df["log_mass"])
        & np.isfinite(df["petroR50_r_arcsec"])
    )
    df = df[valid].copy()

    # Geometric size proxy (kpc) from angular size and redshift distance.
    # Use low-z approximation: Dc ≈ (c/H0) z; d_A = Dc / (1+z)
    Dc_mpc = (_C_KM_S / _H0_KM_S_MPC) * df["redshift"].values
    dA_kpc = (Dc_mpc / (1.0 + df["redshift"].values)) * 1000.0
    theta_rad = df["petroR50_r_arcsec"].values * (np.pi / 180.0 / 3600.0)
    df["R50_kpc"] = dA_kpc * theta_rad
    df["log_R50_kpc"] = np.log10(np.clip(df["R50_kpc"].values, 1e-4, None))

    df["fe_avg"] = (df["fe5270"] + df["fe5335"]) / 2.0
    df["mg_fe_ratio"] = df["mgb"] / df["fe_avg"]
    df["log_mg_fe"] = np.log10(df["mg_fe_ratio"])

    df["spec_age_proxy"] = df["d4000"] / (df["hbeta"] + 0.5)
    df["log_spec_age"] = np.log10(df["spec_age_proxy"])

    df["log_sigma"] = np.log10(df["veldisp"])

    finite2 = np.isfinite(df["log_mg_fe"]) & np.isfinite(df["log_spec_age"]) & np.isfinite(df["log_sigma"])
    df = df[finite2].copy()

    df = df[(df["log_mg_fe"] > -0.6) & (df["log_mg_fe"] < 0.6)]
    df = df[(df["log_spec_age"] > -0.8) & (df["log_spec_age"] < 0.8)]

    return df


def comoving_xyz(df: pd.DataFrame) -> np.ndarray:
    z = df["redshift"].values
    ra = np.deg2rad(df["ra"].values)
    dec = np.deg2rad(df["dec"].values)

    Dc = (_C_KM_S / _H0_KM_S_MPC) * z

    x = Dc * np.cos(dec) * np.cos(ra)
    y = Dc * np.cos(dec) * np.sin(ra)
    zz = Dc * np.sin(dec)
    return np.vstack([x, y, zz]).T


def estimate_local_density(df: pd.DataFrame, k: int = 10) -> pd.Series:
    xyz = comoving_xyz(df)
    nn = NearestNeighbors(n_neighbors=k + 1, algorithm="ball_tree")
    nn.fit(xyz)
    dists, _ = nn.kneighbors(xyz)
    rk = dists[:, -1]
    rk = np.clip(rk, 1e-6, None)
    density = k / ((4.0 / 3.0) * np.pi * rk**3)
    return pd.Series(density, index=df.index, name="rho_k")


@dataclass
class MatchConfig:
    min_delta_log_sigma: float = 0.12
    max_feature_distance: float = 0.35
    include_log_mass: bool = False
    include_log_size: bool = True


def match_twins(df: pd.DataFrame, cfg: MatchConfig) -> pd.DataFrame:
    feature_cols = ["redshift", "log_mg_fe", "rho_k"]
    if cfg.include_log_mass:
        feature_cols.insert(1, "log_mass")
    if cfg.include_log_size:
        # Geometric size proxy, avoids relying on stellar-population M/L assumptions.
        feature_cols.append("log_R50_kpc")

    features = df[feature_cols].copy()
    features = (features - features.mean()) / features.std()

    nn = NearestNeighbors(n_neighbors=25, algorithm="ball_tree")
    nn.fit(features.values)
    dists, inds = nn.kneighbors(features.values)

    used = set()
    pairs = []

    for i in range(len(df)):
        idx_i = df.index[i]
        if idx_i in used:
            continue

        cand_inds = inds[i, 1:]
        cand_dists = dists[i, 1:]

        best = None
        best_d = None

        for j_pos, j in enumerate(cand_inds):
            idx_j = df.index[j]
            if idx_j in used:
                continue

            delta_log_sigma = float(df.loc[idx_j, "log_sigma"] - df.loc[idx_i, "log_sigma"])
            if abs(delta_log_sigma) < cfg.min_delta_log_sigma:
                continue

            if cand_dists[j_pos] > cfg.max_feature_distance:
                continue

            best = idx_j
            best_d = float(cand_dists[j_pos])
            break

        if best is None:
            continue

        used.add(idx_i)
        used.add(best)

        if df.loc[idx_i, "log_sigma"] >= df.loc[best, "log_sigma"]:
            hi, lo = idx_i, best
        else:
            hi, lo = best, idx_i

        pairs.append(
            {
                "specobjid_hi": int(df.loc[hi, "specobjid"]),
                "specobjid_lo": int(df.loc[lo, "specobjid"]),
                "match_dist": best_d,
                "z": float(df.loc[hi, "redshift"]),
                "log_mass": float(df.loc[hi, "log_mass"]),
                "log_mg_fe": float(df.loc[hi, "log_mg_fe"]),
                "rho_k": float(df.loc[hi, "rho_k"]),
                "log_sigma_hi": float(df.loc[hi, "log_sigma"]),
                "log_sigma_lo": float(df.loc[lo, "log_sigma"]),
                "log_age_hi": float(df.loc[hi, "log_spec_age"]),
                "log_age_lo": float(df.loc[lo, "log_spec_age"]),
            }
        )

    out = pd.DataFrame(pairs)
    if len(out) == 0:
        return out

    out["delta_log_sigma"] = out["log_sigma_hi"] - out["log_sigma_lo"]
    out["delta_log_age"] = out["log_age_hi"] - out["log_age_lo"]
    return out


def summarize_pairs(pairs: pd.DataFrame) -> dict:
    if len(pairs) == 0:
        return {"n_pairs": 0}

    delta = pairs["delta_log_age"].values

    n = len(delta)
    mean = float(np.mean(delta))
    med = float(np.median(delta))
    se = float(np.std(delta, ddof=1) / np.sqrt(n)) if n > 1 else float("nan")

    n_neg = int(np.sum(delta < 0))
    n_pos = int(np.sum(delta > 0))

    p_sign = stats.binomtest(n_neg, n_neg + n_pos, p=0.5, alternative="greater").pvalue if (n_neg + n_pos) > 0 else float("nan")

    t_stat, p_t = stats.ttest_1samp(delta, popmean=0.0, alternative="less")

    return {
        "n_pairs": int(n),
        "mean_delta_log_age": mean,
        "median_delta_log_age": med,
        "se_delta_log_age": se,
        "n_delta_neg": n_neg,
        "n_delta_pos": n_pos,
        "p_sign_test_neg": float(p_sign),
        "t_stat": float(t_stat),
        "p_ttest_less": float(p_t),
    }


def plot_pairs(pairs: pd.DataFrame, outdir: str) -> None:
    if len(pairs) == 0:
        return

    os.makedirs(outdir, exist_ok=True)

    plt.figure(figsize=(7, 4.5))
    plt.hist(pairs["delta_log_age"], bins=40, color="#1f77b4", alpha=0.85)
    plt.axvline(0, color="black", linewidth=1)
    plt.xlabel(r"$\Delta\log(\mathrm{age})$ (high-$\sigma$ minus low-$\sigma$)")
    plt.ylabel("Pairs")
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "sdss_twin_pairs_delta_log_age.png"), dpi=160)
    plt.close()

    plt.figure(figsize=(6, 5))
    plt.scatter(pairs["delta_log_sigma"], pairs["delta_log_age"], s=10, alpha=0.25)
    plt.axhline(0, color="black", linewidth=1)
    plt.xlabel(r"$\Delta\log\sigma$ (high-low)")
    plt.ylabel(r"$\Delta\log(\mathrm{age})$ (high-low)")
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "sdss_twin_pairs_scatter.png"), dpi=160)
    plt.close()


def main() -> None:
    out_csv = os.path.join("results", "outputs", "sdss_twin_pairs.csv")
    out_json = os.path.join("results", "outputs", "sdss_twin_pairs_summary.json")
    out_fig = os.path.join("results", "figures", "consistency")
    os.makedirs(os.path.dirname(out_csv), exist_ok=True)
    os.makedirs(os.path.dirname(out_json), exist_ok=True)

    cache_path = os.path.join("data", "sdss", "sdss_twin_base_sample_with_size.csv")
    os.makedirs(os.path.dirname(cache_path), exist_ok=True)

    if os.path.exists(cache_path):
        df_raw = pd.read_csv(cache_path)
        print(f"Loaded cached base sample: {len(df_raw)}")
    else:
        z_ranges = [(0.02, 0.05), (0.05, 0.08), (0.08, 0.10), (0.10, 0.12)]
        chunks = []
        for z_min, z_max in z_ranges:
            print(f"Querying z={z_min:.2f}-{z_max:.2f}...")
            df = download_base_sample(z_min, z_max)
            if df is not None and len(df) > 0:
                print(f"  Retrieved {len(df)}")
                chunks.append(df)
            else:
                print("  No data")
        if len(chunks) == 0:
            raise RuntimeError("No data retrieved from SDSS")
        df_raw = pd.concat(chunks, ignore_index=True).drop_duplicates(subset=["specobjid"])
        df_raw.to_csv(cache_path, index=False)
        print(f"Cached base sample: {len(df_raw)}")

    df = prepare(df_raw)
    print(f"Prepared sample: {len(df)}")

    df["rho_k"] = estimate_local_density(df, k=10)

    cfg = MatchConfig(
        min_delta_log_sigma=0.12,
        max_feature_distance=0.35,
        include_log_mass=False,
        include_log_size=True,
    )
    pairs = match_twins(df, cfg)
    print(f"Matched pairs: {len(pairs)}")

    pairs.to_csv(out_csv, index=False)

    summary = summarize_pairs(pairs)
    summary["match_config"] = cfg.__dict__
    with open(out_json, "w") as f:
        json.dump(summary, f, indent=2)

    plot_pairs(pairs, out_fig)

    print("Summary:")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
