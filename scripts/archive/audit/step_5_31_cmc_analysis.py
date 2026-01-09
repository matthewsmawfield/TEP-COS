#!/usr/bin/env python3
"""
Step 5.31: CMC "Gold Standard" Analysis

This script implements the "Gold Standard" test by comparing observed pulsar 
residuals directly against synthetic pulsars from Cluster Monte Carlo (CMC) models.

It supports two modes:
1. REAL: Loads actual CMC output files (e.g., initial.morepulsars.dat) if available.
2. SYNTHETIC: Generates "CMC-like" data including Mass Segregation and Binary Hardening
   if real files are missing, to validate the analysis pipeline.

Key Physics captured:
- Mass Segregation: Heavier MSPs/Binaries sink to core (r < rc).
- Binary Hardening: Super-thermal velocity kicks from 3-body interactions.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
import json
from pathlib import Path
from datetime import datetime, timezone

# Constants
G = 6.674e-8  # cgs
c = 2.998e10  # cm/s
Msun = 1.989e33  # g
pc_to_cm = 3.086e18

REPO_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = REPO_ROOT / "data" / "cmc"
RESULTS_DIR = REPO_ROOT / "results" / "outputs"
OUT_JSON = RESULTS_DIR / "step_5_31_cmc_analysis.json"
OUT_MD = RESULTS_DIR / "step_5_31_cmc_analysis.md"

def load_cmc_model(filename):
    """
    Load CMC output file.
    Expected columns (approximate): id, mass, r, v_r, v_t, type, ...
    """
    path = DATA_DIR / filename
    if not path.exists():
        return None
    
    print(f"Loading CMC model: {path}")
    # Placeholder for actual column names based on CMC docs
    # Usually: time, id, m, R, r, z, vx, vy, vz ...
    try:
        data = pd.read_csv(path, sep='\s+', comment='#', names=[
            "id", "m", "r", "vr", "vt", "type", "bin_type"
        ])
        return data
    except Exception as e:
        print(f"Error loading {path}: {e}")
        return None

def generate_synthetic_cmc(n_pulsars=1000, cluster_name="Terzan_5"):
    """
    Generate synthetic data that mimics CMC results:
    1. Strong Mass Segregation (r ~ r_core / 3)
    2. Non-Gaussian Velocity Tails (Binary Hardening)
    """
    rng = np.random.default_rng(42)
    
    # Cluster params
    if cluster_name == "Terzan_5":
        M = 2e6 * Msun
        rc = 0.16 * pc_to_cm
    elif cluster_name == "47_Tuc":
        M = 1e6 * Msun
        rc = 0.36 * pc_to_cm
    else:
        M = 5e5 * Msun
        rc = 0.5 * pc_to_cm

    # 1. Mass Segregation: MSPs are centrally concentrated
    # Standard King/Plummer: n(r) ~ r^-2
    # Mass Segregated: n(r) ~ r^-3 or steeper inside core
    
    # Sample radii: Concentrated in core
    # u = (r/rc)^3  => r = rc * u^(1/3) for constant density? 
    # Let's use a specialized segregated profile: 
    # standard stars: r ~ rc
    # MSPs (1.4 Msun): r ~ rc * (m_avg / m_msp)^0.5 ~ rc * (0.4/1.4)^0.5 ~ 0.5 rc
    
    sigma_r = 0.5 * rc # Segregated scale length
    r = np.abs(rng.normal(0, sigma_r, n_pulsars))
    
    # 2. Acceleration (Mean Field)
    # a ~ GM(<r)/r^2. In core, M(<r) ~ r^3, so a ~ r (harmonic)
    # Outside core, a ~ 1/r^2
    
    acc_mean = np.zeros_like(r)
    mask_core = r < rc
    
    # Core: Linear acceleration (harmonic)
    # g_max at rc ~ GM/rc^2
    g_max = G * M / rc**2
    acc_mean[mask_core] = g_max * (r[mask_core] / rc)
    
    # Envelope: Keplerian
    acc_mean[~mask_core] = G * M / r[~mask_core]**2
    
    # Project to Line of Sight
    cos_theta = rng.uniform(-1, 1, n_pulsars)
    a_los = acc_mean * cos_theta
    
    # 3. Binary Hardening (Super-thermal kicks)
    # Add a Lorentzian (Cauchy) component to represent close encounters
    # Standard velocity dispersion
    sigma_v = np.sqrt(G * M / rc)
    
    # Gaussian thermal component
    v_thermal = rng.normal(0, sigma_v, n_pulsars)
    
    # Hardening component (10% of population)
    n_kicked = int(0.1 * n_pulsars)
    # Velocity kick from exchange interaction ~ orbital velocity of binary
    # v_kick ~ 10 * sigma_v
    v_kick = stats.cauchy.rvs(loc=0, scale=2*sigma_v, size=n_pulsars)
    
    # Combine
    v_tot = v_thermal + 0.2 * v_kick # Mixed population
    
    # Total acceleration/jerk equivalent (Pdot/P ~ a/c + v^2/cd ...)
    # Here we just look at a_los/c effect
    
    # Ideally, CMC gives us Pdot directly or we compute it from full state.
    # Here we approximate Pdot_obs = Pdot_int + P * (a_los/c + Shklovskii)
    
    P_s = 0.005 # 5 ms
    
    # Distance to cluster (cm)
    if cluster_name == "Terzan_5":
        D_kpc = 5.9
    elif cluster_name == "47_Tuc":
        D_kpc = 4.5
    else:
        D_kpc = 5.0
        
    D_cm = D_kpc * 1000 * pc_to_cm
    
    pdot_acc = P_s * a_los / c
    pdot_shk = P_s * (v_tot**2) / (c * D_cm) # Centrifugal/Shklovskii approx
    
    # Intrinsic Pdot
    log_pdot_int = rng.normal(-19.76, 0.6, n_pulsars)
    pdot_int = 10**log_pdot_int
    
    pdot_obs = pdot_int + pdot_acc + pdot_shk
    
    return {
        "pdot_obs": pdot_obs,
        "a_los": a_los,
        "v_tot": v_tot,
        "r": r
    }

def analyze_cmc(data, cluster_name):
    """
    Compare CMC/Synthetic data to Observed Residual.
    """
    pdot_obs = data["pdot_obs"]
    
    # Filter negative Pdots (acceleration dominated)
    frac_neg = np.mean(pdot_obs < 0)
    
    # Log|Pdot|
    log_abs_pdot = np.log10(np.abs(pdot_obs))
    
    mu_obs = np.mean(log_abs_pdot)
    mu_field = -19.76
    
    shift = mu_obs - mu_field
    
    return {
        "cluster": cluster_name,
        "shift_dex": float(shift),
        "frac_negative": float(frac_neg),
        "n_pulsars": len(pdot_obs)
    }

def main():
    print("="*70)
    print("CMC GOLD STANDARD ANALYSIS")
    print("="*70)
    
    results = {
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "clusters": {}
    }
    
    # 1. Terzan 5
    print("\nAnalyzing Terzan 5...")
    # Try load real
    df = load_cmc_model("terzan5_msp.dat")
    if df is not None:
        print("✓ Loaded REAL CMC data")
        # Process df...
        # For now, fallback to synthetic for the prototype
        sim_data = generate_synthetic_cmc(2000, "Terzan_5")
        mode = "REAL (Mocked Logic)"
    else:
        print("⚠ CMC data not found. Using SYNTHETIC 'Best-Guess' Model.")
        print("  (Includes Mass Segregation + Binary Hardening)")
        sim_data = generate_synthetic_cmc(2000, "Terzan_5")
        mode = "SYNTHETIC"
        
    res_t5 = analyze_cmc(sim_data, "Terzan_5")
    res_t5["mode"] = mode
    results["clusters"]["Terzan_5"] = res_t5
    print(f"  Shift: {res_t5['shift_dex']:+.3f} dex (Observed: +0.13 dex)")
    
    # 2. 47 Tuc
    print("\nAnalyzing 47 Tuc...")
    sim_data = generate_synthetic_cmc(2000, "47_Tuc")
    res_47 = analyze_cmc(sim_data, "47_Tuc")
    res_47["mode"] = "SYNTHETIC" # Force for now
    results["clusters"]["47_Tuc"] = res_47
    print(f"  Shift: {res_47['shift_dex']:+.3f} dex (Observed: +0.13 dex)")
    
    # Conclusion
    print("\n" + "="*70)
    print("CONCLUSION")
    print("="*70)
    print(f"The 'Gold Standard' CMC/N-body comparison requires:")
    print(f"1. Access to CMC catalog files (e.g. initial.morepulsars.dat)")
    print(f"2. Parsing of full state vectors (r, v, a)")
    print(f"\nCurrent Status:")
    print(f"- Terzan 5: {res_t5['shift_dex']:.3f} dex (Synthetic)")
    print(f"- 47 Tuc:   {res_47['shift_dex']:.3f} dex (Synthetic)")
    print(f"- Observed: +0.13 dex")
    
    if abs(res_t5['shift_dex'] - 0.13) < 0.1:
        print("\nSUCCESS: Synthetic CMC model matches observed residual!")
    else:
        print("\nDISCREPANCY: Synthetic CMC model predicts larger/smaller shift.")
    
    # Save
    with open(OUT_JSON, 'w') as f:
        json.dump(results, f, indent=2)
        
    # Markdown
    md = f"""# CMC 'Gold Standard' Analysis
    
**Generated:** {results['timestamp']}
**Method:** {mode} (Comparison with CMC models)

## Results

| Cluster | Mode | Shift (dex) | Observed |
|---------|------|-------------|----------|
| Terzan 5 | {res_t5['mode']} | {res_t5['shift_dex']:+.3f} | +0.13 |
| 47 Tuc | {res_47['mode']} | {res_47['shift_dex']:+.3f} | +0.13 |

## Interpretation
The inclusion of mass segregation and binary hardening (simulated) results in a shift of {res_t5['shift_dex']:+.3f} dex.
This compares to the simple Mean-Field prediction of ~ +0.65 dex.
The "messy" dynamics [reduce/increase] the discrepancy... (Requires real data for final verdict).

## Action Items
- [ ] Download `initial.morepulsars.dat` for Terzan 5 from CMC Catalog
- [ ] Download `initial.morepulsars.dat` for 47 Tuc from CMC Catalog
"""
    with open(OUT_MD, 'w') as f:
        f.write(md)
    
    print(f"\nResults saved to {OUT_MD}")

if __name__ == "__main__":
    main()
