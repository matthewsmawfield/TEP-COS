#!/usr/bin/env python3
"""
Step 18: Binary Model Validation - Test Period/Mass Predictions

This script tests whether the binary MSP sample shows the predicted
correlations between orbital parameters and spin-down rates.

Model Predictions:
1. Longer orbital period → less Temporal Shear suppression → higher log Ṗ (positive correlation)
2. More massive companion → more suppressed Shear → lower log Ṗ (negative correlation)

Potential Confounders:
- Selection effects: different binary types have different typical masses/periods
- Cluster environment: position within cluster affects baseline enhancement
- Evolutionary effects: MSPs with different histories have different intrinsic Ṗ

Author: M. Smawfield
Date: March 2026
"""

import json
import sys
from pathlib import Path

import numpy as np

# Add project root to path
PROJECT_ROOT = Path(__file__).parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

# Physical constants
G = 6.674e-11  # m³ kg⁻¹ s⁻²
M_sun = 1.989e30  # kg

# Paths
SCRIPT_DIR = Path(__file__).parent
RESULTS_DIR = PROJECT_ROOT / "results" / "outputs"


def compute_shielding_fraction(
    Pb_days: float, Mc_msun: float, Mp_msun: float = 1.4
) -> float:
    """
    Compute predicted shielding fraction from density-based Temporal Shear competition.

    Uses the continuous screening formula from step_16:
    f_shield = (M_c/d) / (M_clust/R_core + M_c/d)

    where:
    - M_c = companion mass
    - d = distance from companion to pulsar saturation radius
    - M_clust = cluster mass
    - R_core = cluster core radius

    This derives from the non-linear superposition of Temporal Shear, where competing
    density wells contribute to the topology transition in proportion to their
    gravitational influence.
    """
    # Physical constants
    G = 6.674e-11  # m³ kg⁻¹ s⁻²
    M_sun = 1.989e30  # kg
    PC = 3.086e16  # m

    # TEP parameters
    rho_c = 20.0  # g/cm³
    rho_c_kg_m3 = rho_c * 1000.0

    # Pulsar saturation radius (Temporal Topology transition scale)
    Mp_kg = Mp_msun * M_sun
    R_sol = (3.0 * Mp_kg / (4.0 * np.pi * rho_c_kg_m3)) ** (1 / 3)

    # Orbital separation from Kepler's third law
    Mc_kg = Mc_msun * M_sun
    M_total = Mp_kg + Mc_kg
    Pb_s = Pb_days * 86400
    a_m = (G * M_total * Pb_s**2 / (4 * np.pi**2)) ** (1 / 3)

    # Distance from companion to pulsar saturation radius
    d_edge = a_m - R_sol
    if d_edge <= 0:
        return 1.0  # Complete overlap

    # Gravitational influence weights (mass/distance)
    M_clust = 1e6 * M_sun
    R_core = 0.5 * PC
    weight_cluster = M_clust / R_core

    Mc_kg = Mc_msun * M_sun
    weight_comp = Mc_kg / d_edge

    # Shielding fraction from Temporal Shear competition
    f_shield = weight_comp / (weight_cluster + weight_comp)

    return float(f_shield)


def load_binary_data():
    """Load binary MSP data from the Freire catalog."""
    # Load from Freire catalog
    data_path = PROJECT_ROOT / "data" / "freire_GCpsr.txt"

    binaries = []
    current_cluster = None

    with open(data_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("# Total"):
                continue

            # Check if this is a cluster header (no tab at start, contains parentheses)
            if not line.startswith("J") and "(" in line:
                current_cluster = line.split("(")[0].strip()
                continue

            # Skip comment lines
            if line.startswith("#"):
                continue

            parts = line.split()
            if len(parts) < 6:
                continue

            try:
                # Parse relevant fields based on actual format:
                # J0024-7205E  0.65  3.53633  +9.8510(6)  24.236(2)  2.25684  1.98184  0.00031  0.18
                # Name  Offset  P(ms)  Pdot  DM  Pb(days)  x(s)  e  m2
                psr_name = parts[0]
                P_ms = float(parts[2])  # Period in milliseconds
                P_s = P_ms / 1000.0  # Convert to seconds

                # Parse Pdot (may have uncertainty in parentheses)
                pdot_str = parts[3]
                if "(" in pdot_str:
                    pdot_val = float(pdot_str.split("(")[0])
                else:
                    pdot_val = float(pdot_str)
                log_Pdot = np.log10(abs(pdot_val) * 1e-20)  # Convert from 10^-20 units

                # Parse binary parameters (columns 5-8)
                Pb_days = None
                m2_msun = None

                # Check if Pb column exists and is not 'i' (isolated)
                if len(parts) > 5 and parts[5] not in ["i", "*"]:
                    try:
                        Pb_days = float(parts[5])
                    except ValueError:
                        Pb_days = None

                # Check if m2 column exists
                if len(parts) > 8 and parts[8] not in ["i", "*"]:
                    try:
                        m2_msun = float(parts[8])
                    except ValueError:
                        m2_msun = None

                # Filter for MSPs (P < 30 ms) and binaries
                if P_s < 0.030 and Pb_days is not None and Pb_days > 0:
                    binaries.append(
                        {
                            "name": psr_name,
                            "cluster": current_cluster or "Unknown",
                            "P_s": P_s,
                            "log_Pdot": log_Pdot,
                            "Pb_days": Pb_days,
                            "m2_msun": m2_msun
                            if m2_msun
                            else 0.2,  # Default companion mass
                        }
                    )
            except (ValueError, IndexError):
                continue

    return binaries


def main():
    print("=" * 80)
    print("BINARY MODEL VALIDATION: PERIOD/MASS PREDICTIONS")
    print("=" * 80)

    # Load data
    binaries = load_binary_data()

    print(f"\n1. SAMPLE STATISTICS")
    print(f"-" * 40)
    print(f"Total binary MSPs: {len(binaries)}")

    # Extract arrays
    Pb = np.array([b["Pb_days"] for b in binaries])
    log_Pdot = np.array([b["log_Pdot"] for b in binaries])
    m2 = np.array([b["m2_msun"] for b in binaries])

    print(f"Orbital period range: {Pb.min():.3f} - {Pb.max():.1f} days")
    print(f"Median period: {np.median(Pb):.2f} days")
    print(f"Companion mass range: {m2.min():.3f} - {m2.max():.2f} M_sun")
    print(f"Median companion mass: {np.median(m2):.3f} M_sun")

    # Compute predicted shielding for each binary
    shielding = np.array(
        [compute_shielding_fraction(b["Pb_days"], b["m2_msun"]) for b in binaries]
    )

    print(f"\n2. PREDICTED SHIELDING")
    print(f"-" * 40)
    print(f"Shielding range: {shielding.min():.2f} - {shielding.max():.2f}")
    print(f"Median shielding: {np.median(shielding):.2f}")

    # Test correlations
    from scipy import stats

    print(f"\n3. CORRELATION TESTS")
    print(f"-" * 40)

    # Period vs log Pdot
    r_pb, p_pb = stats.pearsonr(Pb, log_Pdot)
    print(f"Period vs log Ṗ:")
    print(f"  r = {r_pb:.4f}, p = {p_pb:.4f}")
    print(f"  Prediction: positive (longer Pb → less shielding → higher log Ṗ)")
    print(f"  Observed: {'consistent' if r_pb > 0 else 'OPPOSITE'}")

    # Companion mass vs log Pdot
    r_m2, p_m2 = stats.pearsonr(m2, log_Pdot)
    print(f"\nCompanion mass vs log Ṗ:")
    print(f"  r = {r_m2:.4f}, p = {p_m2:.4f}")
    print(f"  Prediction: negative (higher m2 → more shielding → lower log Ṗ)")
    print(f"  Observed: {'consistent' if r_m2 < 0 else 'OPPOSITE'}")

    # Shielding vs log Pdot (the key test)
    r_shield, p_shield = stats.pearsonr(shielding, log_Pdot)
    print(f"\nPredicted shielding vs log Ṗ:")
    print(f"  r = {r_shield:.4f}, p = {p_shield:.4f}")
    print(f"  Prediction: negative (more shielding → lower log Ṗ)")
    print(f"  Observed: {'consistent' if r_shield < 0 else 'OPPOSITE'}")

    # Partial correlations controlling for period
    print(f"\n4. PARTIAL CORRELATIONS")
    print(f"-" * 40)

    # Shielding vs log Pdot, controlling for period
    # Method: regress out period from both, then correlate residuals
    from sklearn.linear_model import LinearRegression

    # Regress period out of shielding
    reg1 = LinearRegression().fit(Pb.reshape(-1, 1), shielding)
    shield_resid = shielding - reg1.predict(Pb.reshape(-1, 1))

    # Regress period out of log Pdot
    reg2 = LinearRegression().fit(Pb.reshape(-1, 1), log_Pdot)
    pdot_resid = log_Pdot - reg2.predict(Pb.reshape(-1, 1))

    r_partial, p_partial = stats.pearsonr(shield_resid, pdot_resid)
    print(f"Shielding vs log Ṗ (controlling for period):")
    print(f"  r = {r_partial:.4f}, p = {p_partial:.4f}")

    # Binned analysis by period
    print(f"\n5. BINNED ANALYSIS BY PERIOD")
    print(f"-" * 40)

    period_bins = [(0, 0.5), (0.5, 2), (2, 10), (10, 400)]
    print(f"{'Period range':<15} | {'N':<5} | {'Mean log Ṗ':<12} | {'Mean shield':<12}")
    print("-" * 50)

    for pmin, pmax in period_bins:
        mask = (Pb >= pmin) & (Pb < pmax)
        if mask.sum() > 0:
            mean_pdot = log_Pdot[mask].mean()
            mean_shield = shielding[mask].mean()
            print(
                f"{pmin:.1f}-{pmax:.1f} d     | {mask.sum():<5} | {mean_pdot:<12.3f} | {mean_shield:<12.2f}"
            )

    # Binned analysis by companion mass
    print(f"\n6. BINNED ANALYSIS BY COMPANION MASS")
    print(f"-" * 40)

    mass_bins = [(0, 0.05), (0.05, 0.3), (0.3, 0.6), (0.6, 2.0)]
    print(f"{'Mass range':<15} | {'N':<5} | {'Mean log Ṗ':<12} | {'Mean shield':<12}")
    print("-" * 50)

    for mmin, mmax in mass_bins:
        mask = (m2 >= mmin) & (m2 < mmax)
        if mask.sum() > 0:
            mean_pdot = log_Pdot[mask].mean()
            mean_shield = shielding[mask].mean()
            print(
                f"{mmin:.2f}-{mmax:.2f} M☉   | {mask.sum():<5} | {mean_pdot:<12.3f} | {mean_shield:<12.2f}"
            )

    # Summary
    print(f"\n" + "=" * 80)
    print("VALIDATION SUMMARY")
    print("=" * 80)

    issues = []

    if r_pb < 0:
        issues.append("Period correlation opposite to prediction")
    if r_m2 > 0:
        issues.append("Companion mass correlation opposite to prediction")
    if r_shield > 0:
        issues.append("Shielding correlation opposite to prediction")

    if issues:
        print("MODEL ISSUES DETECTED:")
        for issue in issues:
            print(f"  - {issue}")
        print("\nPossible explanations:")
        print("  1. Confounding by cluster environment (position within cluster)")
        print("  2. Selection effects (binary type correlates with MSP properties)")
        print("  3. Model oversimplification (ignores binary evolution)")
        print("  4. Statistical noise (sample size limitations)")
    else:
        print("All correlations consistent with model predictions")

    # Save results
    results = {
        "sample": {
            "n_binaries": len(binaries),
            "period_range_days": [float(Pb.min()), float(Pb.max())],
            "mass_range_msun": [float(m2.min()), float(m2.max())],
        },
        "correlations": {
            "period_vs_logPdot": {"r": float(r_pb), "p": float(p_pb)},
            "mass_vs_logPdot": {"r": float(r_m2), "p": float(p_m2)},
            "shielding_vs_logPdot": {"r": float(r_shield), "p": float(p_shield)},
            "partial_shielding_vs_logPdot": {
                "r": float(r_partial),
                "p": float(p_partial),
            },
        },
        "predictions": {
            "period": "positive correlation expected",
            "mass": "negative correlation expected",
            "shielding": "negative correlation expected",
        },
        "validation": {
            "period_consistent": bool(r_pb > 0),
            "mass_consistent": bool(r_m2 < 0),
            "shielding_consistent": bool(r_shield < 0),
        },
    }

    out_path = RESULTS_DIR / "step_18_model_validation.json"
    out_path.parent.mkdir(exist_ok=True, parents=True)
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)

    print(f"\nResults saved to {out_path}")


if __name__ == "__main__":
    main()
