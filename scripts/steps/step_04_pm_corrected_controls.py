#!/usr/bin/env python3
"""
Step 04: PM-Corrected Field Control Analysis
==================================================

Recompute the GC–field residual using proper-motion-corrected (intrinsic)
Pdot for field MSPs. The Shklovskii term is subtracted:

    Pdot_shk = P * (mu^2 * D / c)

where mu is the total proper motion (mas/yr) and D is the distance (kpc).
In convenient units:

    Pdot_shk ≈ 2.43 × 10⁻²¹ × P(s) × mu(mas/yr)² × D(kpc)

Pulsars without PMRA, PMDEC, and DIST are excluded from the field pool.

Author: M. Smawfield
Date: July 2026
"""

import json
import math
import re
import sys
from pathlib import Path

import numpy as np
from scipy import stats

sys.path.insert(0, str(Path(__file__).parent))
from step_02_pulsar_population_controls import (
    _parse_freire_gcpsr,
    _parse_atnf_psrcat_db,
    _period_matched_bootstrap,
    _two_dim_match_bootstrap,
    _ttest_logpdot,
    BOOTSTRAP_ITERATIONS,
    RANDOM_SEED,
    MSP_P0_CUT_SECONDS,
    REPO_ROOT,
    RESULTS_DIR,
    DATA_DIR,
)

OUT_JSON = RESULTS_DIR / "step_04_pm_corrected_controls.json"

# Shklovskii constant in convenient units
# Pdot_shk = K * P(s) * mu(mas/yr)^2 * D(kpc)
# K = 2.43e-21 s⁻¹
SHK_K = 2.43e-21
C_MS = 3e8


def _parse_atnf_with_pm(atnf_text: str) -> list[dict]:
    """Parse ATNF db extracting PMRA, PMDEC, DIST, PX in addition to standard fields."""
    rows = []
    current = {}
    _NUM_RE = re.compile(r"^[+-]?(?:\d+\.?\d*|\d*\.?\d+)(?:[eE][+-]?\d+)?$")

    def flush():
        nonlocal current
        if not current:
            return

        name = current.get("PSRJ") or current.get("PSRB")
        p0 = current.get("P0")
        p1 = current.get("P1")
        f0 = current.get("F0")
        f1 = current.get("F1")
        if p0 is None and f0 is not None and f0 != 0:
            p0 = 1.0 / f0
        if p1 is None and f0 is not None and f1 is not None and f0 != 0:
            p1 = -f1 / (f0 * f0)

        if name is None or p0 is None or p1 is None:
            current = {}
            return

        # Parse PM and distance fields
        pmra = current.get("PMRA")
        pmdec = current.get("PMDEC")
        px = current.get("PX")
        dist_amn = current.get("DIST_AMN")
        dist_dm = current.get("DIST_DM")
        dist_dm1 = current.get("DIST_DM1")

        # Compute intrinsic Pdot if PM data available
        p1_intrinsic = None
        shk = None
        dist_kpc = None
        has_pm = pmra is not None and pmdec is not None and all(v is not None for v in [pmra, pmdec])
        has_dist = any(v is not None for v in [dist_amn, dist_dm, dist_dm1, px])

        if has_pm and has_dist:
            # Total proper motion in mas/yr
            mu_total = math.sqrt(pmra**2 + pmdec**2)

            # Distance in kpc (prefer PX -> DIST_AMN -> DIST_DM -> DIST_DM1)
            if px is not None and px > 0:
                dist_kpc = 1.0 / px  # kpc from parallax in mas
            elif dist_amn is not None:
                dist_kpc = float(dist_amn)
            elif dist_dm is not None:
                dist_kpc = float(dist_dm)
            elif dist_dm1 is not None:
                dist_kpc = float(dist_dm1)

            if dist_kpc is not None and dist_kpc > 0 and mu_total > 0:
                shk = SHK_K * p0 * (mu_total ** 2) * dist_kpc
                p1_intrinsic = p1 - shk

        rows.append({
            "source": "atnf",
            "environment": None,
            "cluster": None,
            "name": name,
            "P0_s": p0,
            "P_ms": p0 * 1000.0,
            "P1_sps": p1,
            "P1_intrinsic_sps": p1_intrinsic,
            "P1_e20": p1 / 1e-20,
            "P1_intrinsic_e20": p1_intrinsic / 1e-20 if p1_intrinsic is not None else None,
            "assoc": current.get("ASSOC", ""),
            "pmra": pmra,
            "pmdec": pmdec,
            "px": px,
            "dist_kpc": dist_kpc if dist_kpc is not None else None,
            "shk_e20": shk / 1e-20 if shk is not None else None,
        })
        current = {}

    for raw_line in atnf_text.splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith("#"):
            continue
        if line.startswith("@"):
            flush()
            continue

        parts = line.split()
        if len(parts) < 2:
            continue

        key = parts[0].strip()
        val0 = parts[1].strip()
        val_rest = " ".join(parts[1:]).strip()

        if key in {"P0", "P1", "F0", "F1", "PMRA", "PMDEC", "PX", "DIST_AMN", "DIST_DM", "DIST_DM1"}:
            if _NUM_RE.match(val0):
                current[key] = float(val0)
            continue

        if key in {"PSRJ", "PSRB", "ASSOC"}:
            if key in {"PSRJ", "PSRB"}:
                current[key] = val0
            else:
                current[key] = val_rest  # ASSOC may contain spaces
            continue

    flush()
    return rows


def main():
    print("Step 04: PM-Corrected Field Control Analysis")

    # Load Freire GC catalog (GC MSPs come from Freire)
    freire_path = DATA_DIR / "freire_GCpsr.txt"
    if not freire_path.exists():
        print("Error: Freire catalog not found.")
        return

    with open(freire_path) as f:
        freire_rows = _parse_freire_gcpsr(f.read())

    # Cluster tokens for GC exclusion from ATNF field sample
    import re
    freire_cluster_tokens = set()
    for r in freire_rows:
        if r.get("cluster"):
            for tok in re.split(r"[^A-Za-z0-9]+", r["cluster"]):
                tok = tok.strip()
                if tok:
                    freire_cluster_tokens.add(tok)

    # Load ATNF with PM data for field MSPs
    atnf_path = DATA_DIR / "atnf_psrcat.db"
    if not atnf_path.exists():
        print("Error: ATNF database not found.")
        return

    with open(atnf_path) as f:
        atnf_rows = _parse_atnf_with_pm(f.read())

    # Apply MSP cut and compute proxies
    def is_msp_cut(r):
        return r.get("P0_s") is not None and r["P0_s"] < MSP_P0_CUT_SECONDS

    def compute_proxy(r, use_intrinsic=False):
        p0 = r.get("P0_s")
        if p0 is None:
            return None
        p_ms = p0 * 1000.0
        if p_ms >= MSP_P0_CUT_SECONDS * 1000:
            return None

        if use_intrinsic and r.get("P1_intrinsic_sps") is not None:
            p1 = r["P1_intrinsic_sps"]
        else:
            p1 = r.get("P1_sps")

        if p1 is None:
            return None

        logpdot = math.log10(abs(p1)) if abs(p1) > 0 else None
        if logpdot is None:
            return None

        b_proxy = math.sqrt(p0 * abs(p1)) if p0 > 0 and p1 != 0 else None
        log_b_proxy = math.log10(b_proxy) if b_proxy and b_proxy > 0 else None
        return {
            **r,
            "logPdot_abs": logpdot,
            "logP": math.log10(p0),
            "b_proxy": float(b_proxy) if b_proxy else None,
            "log_b_proxy": float(log_b_proxy) if log_b_proxy else None,
        }

    # GC rows from Freire
    gc_rows = [r for r in freire_rows if is_msp_cut(r)]
    gc_rows_p = [r2 for r in gc_rows if (r2 := compute_proxy(r)) is not None]

    # Field rows from ATNF (exclude GC associations and name matches)
    field_rows_all = []
    field_rows_pm = []
    for r in atnf_rows:
        if not is_msp_cut(r):
            continue
        # Exclude GC associations
        assoc = r.get("assoc", "").lower()
        if "gc" in assoc:
            continue
        if any(tok and tok.lower() in assoc for tok in freire_cluster_tokens):
            continue
        # Exclude name matches with Freire GC pulsars
        if any(fr["name"] == r.get("name") for fr in freire_rows):
            continue

        r2 = compute_proxy(r, use_intrinsic=False)
        if r2 is not None:
            field_rows_all.append(r2)
        r2_pm = compute_proxy(r, use_intrinsic=True)
        if r2_pm is not None:
            field_rows_pm.append(r2_pm)

    print(f"  GC MSPs (Freire): {len(gc_rows_p)}")
    print(f"  Field MSPs (all ATNF): {len(field_rows_all)}")
    print(f"  Field MSPs (PM-corrected ATNF): {len(field_rows_pm)}")

    if len(gc_rows_p) == 0 or len(field_rows_all) == 0:
        print("  Insufficient data for comparison")
        return

    gc_logpdot = np.array([r["logPdot_abs"] for r in gc_rows_p])
    field_logpdot_all = np.array([r["logPdot_abs"] for r in field_rows_all])

    base = _ttest_logpdot(gc_logpdot, field_logpdot_all)
    print(f"  Base residual (all field): {base['diff_dex']:.3f} dex")

    if field_rows_pm:
        field_logpdot_pm = np.array([r["logPdot_abs"] for r in field_rows_pm])
        base_pm = _ttest_logpdot(gc_logpdot, field_logpdot_pm)
        print(f"  Base residual (PM-corrected field): {base_pm['diff_dex']:.3f} dex")
    else:
        base_pm = None
        print("  No PM-corrected field MSPs available")

    # Matching analyses
    period_match_all = _period_matched_bootstrap(gc_rows_p, field_rows_all, n_boot=BOOTSTRAP_ITERATIONS, seed=RANDOM_SEED)
    two_dim_match_all = _two_dim_match_bootstrap(gc_rows_p, field_rows_all, n_boot=BOOTSTRAP_ITERATIONS, seed=RANDOM_SEED)

    print(f"  Period-matched (all field): {period_match_all['diff_mean']:.3f} dex")
    print(f"  2D-matched (all field):     {two_dim_match_all['diff_mean']:.3f} dex")

    if field_rows_pm:
        period_match_pm = _period_matched_bootstrap(gc_rows_p, field_rows_pm, n_boot=BOOTSTRAP_ITERATIONS, seed=RANDOM_SEED)
        two_dim_match_pm = _two_dim_match_bootstrap(gc_rows_p, field_rows_pm, n_boot=BOOTSTRAP_ITERATIONS, seed=RANDOM_SEED)
        print(f"  Period-matched (PM field):  {period_match_pm['diff_mean']:.3f} dex")
        print(f"  2D-matched (PM field):      {two_dim_match_pm['diff_mean']:.3f} dex")
    else:
        period_match_pm = None
        two_dim_match_pm = None

    # Output
    result = {
        "meta": {
            "description": "PM-corrected field control analysis",
            "shklovskii_constant": SHK_K,
            "msp_cut_ms": MSP_P0_CUT_SECONDS * 1000,
        },
        "counts": {
            "gc_msps": len(gc_rows_p),
            "field_msps_all": len(field_rows_all),
            "field_msps_pm_corrected": len(field_rows_pm),
        },
        "base_residual_all_field": base,
        "base_residual_pm_field": base_pm,
        "period_matched_all": period_match_all,
        "two_dim_matched_all": two_dim_match_all,
        "period_matched_pm": period_match_pm,
        "two_dim_matched_pm": two_dim_match_pm,
        "interpretation": (
            "If PM-corrected field MSPs show lower |Pdot| (more negative P1_intrinsic), "
            "the GC–field residual should INCREASE, strengthening the TEP signal. "
            "If the residual decreases, Shklovskii contamination was partially responsible for the observed excess."
        ),
    }

    with open(OUT_JSON, "w") as f:
        json.dump(result, f, indent=2)

    print(f"\nOutput written to {OUT_JSON}")


if __name__ == "__main__":
    main()
