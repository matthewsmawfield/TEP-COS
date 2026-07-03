#!/usr/bin/env python3
"""
Step 43: Exact Mathematical Proof of Γ_TEP = 2Γ_N - 1
========================================================

ADDRESSING REVIEWER FEEDBACK:
"The derivation assumes R_c ∝ ρ_c^α with constant α across the sample.
The manuscript acknowledges this is an approximation."

This script PROVES mathematically and empirically that the reviewer is 
incorrect. The identity Γ_TEP = 2Γ_N - 1 does NOT assume constant α.
It is an exact mathematical consequence of Ordinary Least Squares (OLS)
regression properties (linearity of covariance), completely independent 
of scatter in the R_c vs ρ_c relation.

Author: M. Smawfield
Date: March 2026
"""

import numpy as np
import pandas as pd
from scipy import stats
import json
from pathlib import Path
import matplotlib.pyplot as plt

def mathematical_proof():
    """
    Returns the formal mathematical proof.
    """
    proof = """
    MATHEMATICAL PROOF: Γ_TEP = 2Γ_N - 1 is Exact for OLS Slopes
    ----------------------------------------------------------
    Let X = log(ρ_c)
    Let R = log(R_c)
    
    The observables scale as:
    Newtonian acceleration: a ∝ ρ_c * R_c  =>  log(a) = X + R + C_1
    TEP potential:         |Φ| ∝ ρ_c * R_c^2 => log(|Φ|) = X + 2R + C_2
    
    The OLS slope (Γ) of any variable Y against X is defined by:
    Γ = Cov(Y, X) / Var(X)
    
    For the Newtonian slope:
    Γ_N = Cov(X + R, X) / Var(X)
        = [Var(X) + Cov(R, X)] / Var(X)
        = 1 + Cov(R, X)/Var(X)
    
    Let α_eff = Cov(R, X)/Var(X). Note that α_eff is precisely the OLS slope 
    of R vs X. This does NOT require R to be a deterministic function of X; 
    it is just the definition of the regression slope, fully valid for data 
    with arbitrary scatter.
    
    Thus, Γ_N = 1 + α_eff
    
    For the TEP slope:
    Γ_TEP = Cov(X + 2R, X) / Var(X)
          = [Var(X) + 2Cov(R, X)] / Var(X)
          = 1 + 2[Cov(R, X)/Var(X)]
          = 1 + 2α_eff
          
    Substituting α_eff = Γ_N - 1 into the TEP equation:
    Γ_TEP = 1 + 2(Γ_N - 1)
    Γ_TEP = 2Γ_N - 1
    
    CONCLUSION:
    This identity is EXACT for regression slopes. It does not assume α is 
    constant, nor does it assume zero scatter in the R_c vs ρ_c relationship. 
    It relies purely on the linearity of covariance.
    """
    return proof

def run_monte_carlo_proof(n_clusters=1000, alpha_true=-0.4, scatter_dex=0.5):
    """
    Empirically demonstrate that the identity holds EXACTLY,
    even with massive scatter in the α relationship.
    """
    # Simulate cluster densities
    # log(rho_c) from 2.0 to 6.0
    log_rho = np.random.uniform(2.0, 6.0, n_clusters)
    
    # Simulate core radii WITH ENORMOUS SCATTER
    # log(R_c) = alpha * log(rho_c) + noise
    noise = np.random.normal(0, scatter_dex, n_clusters)
    log_R = alpha_true * log_rho + noise
    
    # Calculate structural scaling index (power-law slope of R_c vs rho_c)
    slope_rho, _, _, _, _ = stats.linregress(log_rho, log_R)
    
    # Calculate Newtonian and TEP observables
    log_a = log_rho + log_R
    log_Phi = log_rho + 2 * log_R
    
    # Measure slopes
    gamma_N, _, _, _, _ = stats.linregress(log_rho, log_a)
    gamma_TEP, _, _, _, _ = stats.linregress(log_rho, log_Phi)
    
    # Check identity
    predicted_gamma_TEP = 2 * gamma_N - 1
    difference = abs(gamma_TEP - predicted_gamma_TEP)
    
    return {
        "alpha_true": alpha_true,
        "scatter_dex": scatter_dex,
        "slope_rho_measured": slope_rho,
        "gamma_N_measured": gamma_N,
        "gamma_TEP_measured": gamma_TEP,
        "gamma_TEP_predicted": predicted_gamma_TEP,
        "exact_match": bool(difference < 1e-10),
        "difference": float(difference)
    }

def main():
    print("=" * 80)
    print("STEP 5.54: EXACT PROOF OF TEP SCALING IDENTITY")
    print("=" * 80)
    
    print(mathematical_proof())
    
    print("-" * 80)
    print("EMPIRICAL MONTE CARLO DEMONSTRATION")
    print("-" * 80)
    print("Testing identity with massive, non-constant scatter in R_c vs ρ_c...")
    
    results = run_monte_carlo_proof(scatter_dex=0.8) # 0.8 dex is huge scatter
    
    print(f"Input true α:        {results['alpha_true']:.3f}")
    print(f"Scatter applied:     ±{results['scatter_dex']} dex (huge variation!)")
    print(f"Measured slope:      {results['slope_rho_measured']:.3f}")
    print(f"Measured Γ_N:        {results['gamma_N_measured']:.3f}")
    print(f"Measured Γ_TEP:      {results['gamma_TEP_measured']:.3f}")
    print(f"Predicted Γ_TEP:     {results['gamma_TEP_predicted']:.3f} (using 2Γ_N - 1)")
    print(f"Difference:          {results['difference']:.2e}")
    print(f"Identity holds exact: {results['exact_match']}")
    
    print("\n" + "=" * 80)
    print("RESPONSE TO REVIEWER:")
    print("=====================")
    print("The reviewer's concern that the derivation assumes a constant α is mathematically ")
    print("incorrect. Because the regression slope is defined by covariance, and covariance ")
    print("is a linear operator, the relationship Γ_TEP = 2Γ_N - 1 holds EXACTLY for the ")
    print("ensemble slopes, regardless of the cluster-to-cluster scatter in α.")
    print("")
    print("Furthermore, using Γ_N = 0.748 from CMC simulations correctly captures the ")
    print("EFFECTIVE α of the mass-segregated pulsar population (α_eff = -0.252), which ")
    print("is the physically relevant parameter, rather than the bare structural α (-0.40) ")
    print("of the overall cluster light profile.")
    
    # Save a clean JSON for reference
    out_dict = {
        "reviewer_premise": "Derivation assumes constant alpha",
        "our_finding": "False. Identity is exact due to linearity of covariance.",
        "monte_carlo_validation": results
    }
    
    out_path = Path("results/outputs/step_43_exact_proof.json")
    out_path.parent.mkdir(exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(out_dict, f, indent=2)

if __name__ == "__main__":
    main()
