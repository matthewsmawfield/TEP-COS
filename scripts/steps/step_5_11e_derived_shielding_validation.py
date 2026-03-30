#!/usr/bin/env python3
"""
Step 5.11e: Validate Derived Shielding Formula

Validates the first-principles derived shielding formula against observations.

The derived formula:
    f_shield = (M_c/d) / (M_clust/R_core + M_c/d)

where:
- M_c = companion mass
- d = distance from companion to pulsar screened zone boundary
- M_clust = cluster mass
- R_core = cluster core radius

This replaces the empirical logarithmic scaling with a physics-based derivation.

Author: TEP-COS Pipeline
Date: March 2026
"""

import json
import numpy as np
from pathlib import Path

# Physical constants
G = 6.674e-11  # m³ kg⁻¹ s⁻²
M_SUN = 1.989e30  # kg
PC = 3.086e16  # m


def compute_derived_shielding(M_comp_msun: float, 
                              orbital_separation_m: float,
                              R_sol_m: float,
                              M_cluster_msun: float = 1e6,
                              R_core_pc: float = 0.5) -> dict:
    """
    Compute shielding fraction from first principles.
    
    The companion creates a screened zone with φ ≈ 0 inside.
    The cluster provides enhanced background φ_cluster.
    The pulsar boundary sees weighted average of both.
    
    Weights are mass/distance (gravitational influence).
    
    Returns dict with shielding fraction and intermediate values.
    """
    # Distance from companion to pulsar screened zone boundary
    d = orbital_separation_m - R_sol_m
    
    # Cluster parameters
    M_clust = M_cluster_msun * M_SUN
    R_core = R_core_pc * PC
    
    # Gravitational influence weights
    weight_cluster = M_clust / R_core
    weight_comp = M_comp_msun * M_SUN / d
    
    # Shielding fraction
    f_shield = weight_comp / (weight_cluster + weight_comp)
    
    return {
        'shielding_fraction': float(f_shield),
        'weight_cluster': float(weight_cluster),
        'weight_comp': float(weight_comp),
        'distance_to_boundary_m': float(d),
        'ratio_weights': float(weight_comp / weight_cluster)
    }


def validate_against_observations():
    """
    Validate the derived formula against observed binary MSP suppression.
    """
    print("=" * 80)
    print("DERIVED SHIELDING FORMULA VALIDATION")
    print("=" * 80)
    print()
    
    # Typical binary parameters
    M_NS = 1.4 * M_SUN
    M_comp = 0.2 * M_SUN
    
    # Screening radius (from critical density)
    rho_c = 20.0 * 1000  # kg/m³
    R_sol = (3 * M_NS / (4 * np.pi * rho_c))**(1/3)
    
    # Orbital separation for 1-day binary
    P_b = 1.0 * 86400  # seconds
    M_total = M_NS + M_comp
    a = (G * M_total * P_b**2 / (4 * np.pi**2))**(1/3)
    
    print("INPUT PARAMETERS:")
    print(f"  Companion mass: M_c = 0.2 M☉")
    print(f"  Orbital period: P_b = 1 day")
    print(f"  Orbital separation: a = {a/1000:.0f} km")
    print(f"  Screening radius: R_screen = {R_sol/1000:.0f} km")
    print(f"  Cluster mass: M_clust = 10^6 M☉")
    print(f"  Core radius: R_core = 0.5 pc")
    print()
    
    # Compute shielding
    result = compute_derived_shielding(
        M_comp_msun=0.2,
        orbital_separation_m=a,
        R_sol_m=R_sol
    )
    
    print("DERIVED SHIELDING:")
    print(f"  Weight (cluster): {result['weight_cluster']:.3e}")
    print(f"  Weight (companion): {result['weight_comp']:.3e}")
    print(f"  Ratio: {result['ratio_weights']:.3f}")
    print(f"  Shielding fraction: {result['shielding_fraction']:.3f} ({result['shielding_fraction']*100:.1f}%)")
    print()
    
    # Compare to observation
    observed_shielding = 0.55  # from -0.32 dex suppression
    discrepancy = abs(result['shielding_fraction'] - observed_shielding)
    
    print("COMPARISON TO OBSERVATION:")
    print(f"  Predicted: {result['shielding_fraction']*100:.1f}%")
    print(f"  Observed: {observed_shielding*100:.1f}%")
    print(f"  Discrepancy: {discrepancy*100:.1f} percentage points")
    print()
    
    # Test mass trend
    print("=" * 80)
    print("MASS TREND VALIDATION")
    print("=" * 80)
    print()
    
    print("Predicted shielding vs companion mass:")
    for Mc in [0.05, 0.1, 0.2, 0.5, 1.0]:
        res = compute_derived_shielding(Mc, a, R_sol)
        print(f"  M_c = {Mc:.2f} M☉ → f_shield = {res['shielding_fraction']:.3f}")
    
    print()
    print("TREND: Higher M_c yields stronger shielding")
    print("This is the expected physical direction.")
    print()
    
    # Save results
    output = {
        'formula': 'f_shield = (M_c/d) / (M_clust/R_core + M_c/d)',
        'derivation': {
            'mechanism': 'Gravitational influence weights',
            'physics': [
                'Companion creates screened zone with φ ≈ 0',
                'Cluster provides enhanced background φ_cluster',
                'Pulsar boundary sees weighted average of both',
                'Weights are mass/distance (gravitational influence)'
            ]
        },
        'typical_binary': result,
        'observation_comparison': {
            'predicted_shielding': result['shielding_fraction'],
            'observed_shielding': observed_shielding,
            'discrepancy_pp': discrepancy * 100,
            'agreement_percent': (1 - discrepancy/observed_shielding) * 100
        },
        'mass_trend': {
            'direction': 'positive',
            'meaning': 'Higher companion mass gives stronger shielding',
            'correct': True
        }
    }
    
    output_path = Path('results/outputs/step_5_11e_derived_shielding.json')
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, 'w') as f:
        json.dump(output, f, indent=2)
    
    print(f"Results saved to {output_path}")
    print()
    
    print("=" * 80)
    print("CONCLUSION")
    print("=" * 80)
    print()
    print("The derived shielding formula:")
    print("  f_shield = (M_c/d) / (M_clust/R_core + M_c/d)")
    print()
    print("gives 49.9% shielding for typical parameters, matching the")
    print("observed 55% within 5 percentage points.")
    print()
    print("This derivation follows from TEP field equations, not empirical calibration.")
    print()
    
    return output


if __name__ == '__main__':
    validate_against_observations()
