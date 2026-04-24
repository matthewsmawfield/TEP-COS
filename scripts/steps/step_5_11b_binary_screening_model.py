#!/usr/bin/env python3
"""
Step 5.11b: Quantitative Binary Screening Model (Temporal Topology Competition)

This script develops a quantitative model for the binary inversion mechanism,
based on the continuous geometric screening framework established in TEP v0.7.

The physical insight:
1. A neutron star's extreme self-gravity causes the TEP scalar field to saturate,
   flattening the Temporal Topology (φ-profile) and suppressing Temporal Shear (∇φ)
   within the saturation radius.
2. The saturation radius scales as R_sol ∝ M^(1/3) ρ_c^(-1/3).
3. The absolute clock rate inside the pulsar is determined by the field value at the
   topology transition region, where the flattened profile recovers toward the
   ambient cluster field.
4. For isolated pulsars, the topology anchors directly to the cluster's unscreened
   field with active Temporal Shear.
5. For binary pulsars, the companion's deep potential well introduces a competing
   region of suppressed Shear, partially decoupling the pulsar clock from the
   cluster enhancement through continuous interpolation rather than step-function
   boundary crossing.

Author: M. Smawfield
Date: March 2026
"""

import json
from pathlib import Path

import numpy as np

# Physical constants
G = 6.674e-11  # m³ kg⁻¹ s⁻²
c = 2.998e8  # m/s
M_sun = 1.989e30  # kg
R_sun = 6.957e8  # m

# Paths
SCRIPT_DIR = Path(__file__).parent
PROJECT_ROOT = SCRIPT_DIR.parent.parent
RESULTS_DIR = PROJECT_ROOT / "results" / "outputs"


def compute_screening_radius(M_obj_msun: float, rho_c_g_cm3: float = 20.0) -> float:
    """
    Compute the TEP Temporal Topology saturation radius for a compact object.

    In TEP theory, the scalar field saturates at density ρ_c, flattening the
    Temporal Topology and suppressing Temporal Shear within this radius.
    R_sol = (3M / 4πρ_c)^(1/3)

    Parameters:
    -----------
    M_obj_msun : float
        Mass of the object in solar masses
    rho_c_g_cm3 : float
        Critical density in g/cm³ (Universal constant from Paper 7 ~ 20.0)

    Returns:
    --------
    R_sol_m : float
        Saturation radius in meters (Temporal Topology transition scale)
    """
    M_kg = M_obj_msun * M_sun
    rho_c_kg_m3 = rho_c_g_cm3 * 1000.0  # Convert g/cm³ to kg/m³

    R_screen_m = (3.0 * M_kg / (4.0 * np.pi * rho_c_kg_m3)) ** (1 / 3)
    return R_screen_m


def compute_binary_orbital_parameters(
    Pb_days: float, Mc_msun: float, Mp_msun: float = 1.4
) -> dict:
    """
    Compute the physical scale of the binary system.

    Parameters:
    -----------
    Pb_days : float
        Orbital period in days
    Mc_msun : float
        Companion mass in solar masses
    Mp_msun : float
        Pulsar mass in solar masses

    Returns:
    --------
    dict with physical scales
    """
    Pb_s = Pb_days * 86400
    Mc_kg = Mc_msun * M_sun
    Mp_kg = Mp_msun * M_sun
    M_total = Mp_kg + Mc_kg

    # Semi-major axis from Kepler's third law
    a_m = (G * M_total * Pb_s**2 / (4 * np.pi**2)) ** (1 / 3)

    # Distance from pulsar to companion
    a_p = a_m * Mc_kg / M_total  # pulsar's distance from CoM
    a_c = a_m * Mp_kg / M_total  # companion's distance from CoM
    separation_m = a_p + a_c  # = a_m

    return {
        "Pb_days": Pb_days,
        "Mc_msun": Mc_msun,
        "Mp_msun": Mp_msun,
        "separation_m": separation_m,
        "separation_km": separation_m / 1000.0,
        "separation_Rsun": separation_m / R_sun,
    }


def compute_shielding_fraction(
    Mc_msun: float,
    separation_m: float,
    M_clust_msun: float = 1e6,
    R_core_pc: float = 0.5,
) -> float:
    """
    Compute the shielding fraction for a binary pulsar via Temporal Shear competition.

    In the TEP continuous screening picture, the scalar field value at the
    pulsar is determined by competing Temporal Shear contributions from the
    cluster (active Shear, φ ≈ φ_cluster) and the companion (suppressed Shear,
    φ ≈ φ_min(ρ_c)). The companion's region of flattened Temporal Topology
    partially shields the pulsar from the cluster's enhancement.

    The shielding fraction equals the companion's gravitational influence
    relative to the total:

    f_shield = (M_c/d) / (M_clust/R_core + M_c/d)

    This derives from the non-linear superposition of Temporal Shear, where
    competing density wells contribute in proportion to their mass-to-distance
    ratio. The interpolation is continuous across the spatial profile; there
    is no step-function boundary.

    Parameters:
    -----------
    Mc_msun : float
        Companion mass in solar masses
    separation_m : float
        Pulsar-companion separation in meters
    M_clust_msun : float
        Cluster mass in solar masses (default: 10^6)
    R_core_pc : float
        Cluster core radius in parsecs (default: 0.5)

    Returns:
    --------
    f_shield : float
        Shielding fraction (0 to 1), where 1 means fully shielded
        (no cluster enhancement) and 0 means unshielded (full enhancement)
    """
    # Cluster influence: M_clust / R_core
    R_core_m = R_core_pc * 3.086e16  # Convert pc to meters
    M_clust_kg = M_clust_msun * M_sun
    cluster_influence = M_clust_kg / R_core_m

    # Companion influence: M_c / d (distance to pulsar saturation radius)
    # For typical separations >> R_sol, d ≈ separation
    Mc_kg = Mc_msun * M_sun
    companion_influence = Mc_kg / separation_m

    # Shielding fraction is companion's relative influence
    f_shield = companion_influence / (cluster_influence + companion_influence)

    return float(f_shield)


def main():
    print("=" * 80)
    print("TEP BINARY SCREENING: TEMPORAL TOPOLOGY COMPETITION MODEL")
    print("=" * 80)

    # 1. Compute Pulsar Screening Radius
    Mp = 1.4  # typical MSP
    rho_c = 20.0  # g/cm³ from TEP-GNSS unification

    R_sol_p = compute_screening_radius(Mp, rho_c)

    print(f"\n1. PULSAR TEMPORAL TOPOLOGY")
    print(f"--------------------------")
    print(f"Pulsar Mass: {Mp} M_sun")
    print(f"Critical Density (TEP): {rho_c} g/cm³")
    print(f"Saturation Radius (R_sol): {R_sol_p / 1000.0:.1f} km")

    # 2. Typical Binary Scales
    Pb_typical = 1.0  # days
    Mc_typical = 0.2  # M_sun

    binary_scales = compute_binary_orbital_parameters(Pb_typical, Mc_typical, Mp)
    sep_km = binary_scales["separation_km"]

    print(f"\n2. TYPICAL MSP BINARY")
    print(f"---------------------")
    print(f"Companion Mass: {Mc_typical} M_sun")
    print(f"Orbital Period: {Pb_typical} days")
    print(f"Orbital Separation (a): {sep_km:.1f} km")
    print(f"Ratio (a / R_sol): {sep_km / (R_sol_p / 1000.0):.1f}")

    # Check intersection
    R_sol_c = compute_screening_radius(Mc_typical, rho_c)
    print(f"Companion Saturation Radius: {R_sol_c / 1000.0:.1f} km")

    if sep_km < (R_sol_p + R_sol_c) / 1000.0:
        print("Saturation regions overlap (highly nonlinear regime)")
    else:
        print("Saturation regions are distinct, with interacting topology transitions")

    # 3. Disruption Model
    # Target observation: 0.32 dex suppression out of 0.58 dex total enhancement
    # This implies a shielding fraction of ~ 0.32 / 0.58 = 0.55
    target_shielding = 0.32 / 0.58

    shielding_f = compute_shielding_fraction(Mc_typical, binary_scales["separation_m"])

    print(f"\n3. TEMPORAL TOPOLOGY COMPETITION EFFECT")
    print(f"-----------------------------")
    print(f"Observed total cluster enhancement: +0.58 dex")
    print(f"Observed binary suppression:        -0.32 dex")
    print(f"Required shielding fraction:        {target_shielding:.2f} (55%)")

    print(f"\nModel Prediction:")
    print(f"Predicted shielding fraction:       {shielding_f:.2f}")

    # 4. Range Sweep
    print(f"\n4. POPULATION SWEEP")
    print(f"-------------------")
    print(
        f"{'Pb (d)':<8} | {'Mc (M_sun)':<10} | {'Sep (km)':<12} | {'Shielding %':<12}"
    )
    print("-" * 48)

    sweep_results = []

    periods = [0.1, 0.5, 1.0, 5.0, 10.0, 100.0]
    masses = [0.05, 0.2, 0.5]

    for m in masses:
        for p in periods:
            b_scales = compute_binary_orbital_parameters(p, m, Mp)
            f_d = compute_shielding_fraction(m, b_scales["separation_m"])
            print(
                f"{p:<8.1f} | {m:<10.2f} | {b_scales['separation_km']:<12.1e} | {f_d * 100:>8.1f}%"
            )

            sweep_results.append(
                {
                    "period_days": p,
                    "companion_mass": m,
                    "separation_km": b_scales["separation_km"],
                    "shielding_fraction": f_d,
                }
            )

    # Compile output
    results = {
        "constants": {"rho_c_g_cm3": rho_c, "M_pulsar": Mp},
        "pulsar_topology": {"saturation_radius_km": R_sol_p / 1000.0},
        "typical_binary": {
            "period_days": Pb_typical,
            "companion_mass": Mc_typical,
            "separation_km": sep_km,
            "shielding_predicted": shielding_f,
            "shielding_observed": target_shielding,
        },
        "sweep": sweep_results,
    }

    out_path = RESULTS_DIR / "step_5_11b_screening_boundary_model.json"
    out_path.parent.mkdir(exist_ok=True, parents=True)
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)

    print(f"\nResults saved to {out_path}")


if __name__ == "__main__":
    main()
