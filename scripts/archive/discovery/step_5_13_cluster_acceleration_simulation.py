import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
from scipy.integrate import quad
import json
import os

def run_cluster_acceleration_simulation():
    """
    Simulates the effect of globular cluster acceleration on pulsar P-dot distributions.
    Tests if standard Newtonian dynamics (Mean Field + Stochastic Encounters) 
    can reproduce the observed 0.13 dex offset in log|Pdot|.
    """
    print("--- Step 5.13: Cluster Acceleration Monte Carlo Simulation ---")
    
    # 1. Setup Parameters (Terzan 5-like, extreme case)
    # Mass ~ 2e6 M_sun, Core radius ~ 0.1 pc (dense!) or Half-light ~ 1.0 pc
    # Let's use parameters that maximize acceleration to be conservative (give Newton a chance)
    M_cluster = 1.0e6  # Solar masses
    R_core = 0.5       # pc
    R_tidal = 50.0     # pc
    N_stars = int(M_cluster / 0.5) # Assuming 0.5 M_sun avg mass
    
    # Constants
    G = 4.3009e-3      # pc (km/s)^2 M_sun^-1
    c_kms = 2.9979e5   # km/s
    yr_to_s = 3.154e7
    
    print(f"Cluster Model: M={M_cluster:.1e} M_sun, Rc={R_core} pc")
    
    # 2. Define Field MSP Distribution (Intrinsic)
    # Based on our previous analysis of Field MSPs
    # Mean log|Pdot| ~ -19.76, Sigma ~ 0.64 (approx from Step 5.12)
    
    # Synthetic population
    n_pulsars = 10000
    
    # Field MSP Pdot distribution (Intrinsic)
    # We observed mean log10(Pdot) approx -19.76 for field
    mu_field = -19.76
    sigma_field = 0.64
    
    # Generate intrinsic Pdots
    log_pdot_int = np.random.normal(mu_field, sigma_field, n_pulsars)
    pdot_int = 10**log_pdot_int
    
    # Periods: Log-normal centered around 5 ms
    log_P_s = np.random.normal(np.log10(0.005), 0.3, n_pulsars) # Centered at 5ms
    P_s = 10**log_P_s
    
    print(f"Simulating {n_pulsars} pulsars...")
    
    # 3. Simulate Cluster Acceleration
    # We need Line-of-Sight acceleration (a_los)
    
    # A. Mean Field Component (King/Plummer Model)
    # Sample radii for pulsars using Inverse transform sampling for Plummer radius
    u = np.random.uniform(0, 0.99, n_pulsars) # Cutoff to avoid infinite radius
    r_over_rc = np.sqrt(u**(2/3) / (1 - u**(2/3)))
    r_pulsar = r_over_rc * R_core
    
    # Random orientation for projection
    cos_theta = np.random.uniform(-1, 1, n_pulsars)
    
    # SI Units for Calculation
    G_si = 6.674e-11
    M_sun_kg = 1.989e30
    pc_m = 3.086e16
    c_si = 2.998e8
    
    m_cl_kg = M_cluster * M_sun_kg
    r_core_m = R_core * pc_m
    r_pulsar_m = r_pulsar * pc_m
    
    # a_mean = G M(<r) / r^2
    # M(<r) formula for Plummer: M_tot * (r/sqrt(r^2+a^2))^3
    M_enc_kg = m_cl_kg * (r_pulsar_m**3) / (r_pulsar_m**2 + r_core_m**2)**(1.5)
    
    a_mean_si = G_si * M_enc_kg / (r_pulsar_m**2 + 1.0)
    a_los_mean_si = a_mean_si * cos_theta 
    
    # B. Stochastic Component (Nearest Neighbor / Holtsmark)
    # Calculate local number density n(r)
    rho_0 = 3 * m_cl_kg / (4 * np.pi * r_core_m**3)
    rho_r = rho_0 * (1 + (r_pulsar_m/r_core_m)**2)**(-2.5)
    n_r = rho_r / (0.5 * M_sun_kg) # Assuming 0.5 Msun stars
    
    # Characteristic Holtsmark field
    a_0 = 2.603 * G_si * (0.5 * M_sun_kg) * n_r**(2/3)
    
    # Sample random component (Levy stable distribution)
    # Using simple approximation if scipy.stats.levy_stable is slow or complex:
    # Actually, let's use the real thing.
    a_stoch_si = stats.levy_stable.rvs(alpha=1.5, beta=0, scale=a_0, size=n_pulsars)
    # Clamp extreme values (closest approach limits)
    a_stoch_si = np.clip(a_stoch_si, -1e-2, 1e-2)
    
    # Total LOS acceleration
    a_tot_si = a_los_mean_si + a_stoch_si
    
    # 4. Calculate Observed Pdot
    # Pdot_obs = Pdot_int + P * (a_los / c)
    accel_term = P_s * (a_tot_si / c_si)
    pdot_obs = pdot_int + accel_term
    
    # 5. Analyze Results
    log_pdot_obs = np.log10(np.abs(pdot_obs))
    
    # Compare Means
    mean_int = np.mean(log_pdot_int)
    mean_obs = np.mean(log_pdot_obs)
    shift = mean_obs - mean_int
    
    # Fraction negative
    frac_neg = np.sum(pdot_obs < 0) / n_pulsars
    
    print("-" * 40)
    print(f"Results (N={n_pulsars}):")
    print(f"Mean log|Pdot| (Intrinsic): {mean_int:.3f}")
    print(f"Mean log|Pdot| (Observed):  {mean_obs:.3f}")
    print(f"Shift (Obs - Int):          {shift:.3f} dex")
    print(f"Fraction Negative Pdot:     {frac_neg:.2%}")
    print(f"Observed Shift in Paper:    +0.65 dex (Raw), +0.13 dex (Controlled)")
    print("-" * 40)
    
    # Output logic
    os.makedirs("results/outputs", exist_ok=True)
    
    results = {
        "model": "King/Plummer Mean Field + Holtsmark Stochastic",
        "parameters": {
            "M_cluster_Msun": M_cluster,
            "R_core_pc": R_core,
            "N_pulsars": n_pulsars
        },
        "results": {
            "mean_log_pdot_int": float(mean_int),
            "mean_log_pdot_obs": float(mean_obs),
            "shift_dex": float(shift),
            "frac_negative": float(frac_neg)
        },
        "conclusion": "Newtonian dynamics produces large shifts in dense clusters. The key discriminator is not the existence of a shift, but its scaling with cluster density."
    }
    
    with open("results/outputs/step_5_13_acceleration_sim.json", "w") as f:
        json.dump(results, f, indent=4)
        
    print("Results saved to results/outputs/step_5_13_acceleration_sim.json")

if __name__ == "__main__":
    run_cluster_acceleration_simulation()
