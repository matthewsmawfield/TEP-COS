import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
import json
import os

def run_density_scaling_simulation():
    """
    Simulates Newtonian acceleration shift for clusters of varying density.
    Demonstrates that dynamical noise scales with density, whereas the TEP signal should track potential.
    """
    print("--- Step 5.14: Density Scaling Simulation ---")
    
    # Cluster Models (Updated from Baumgardt 2018 / step_5_30 N-body)
    # M in M_sun, Rc in pc
    # Added Mass Segregation Factor (alpha): Rc_pulsar = Rc / alpha
    # alpha=1.0 (No segregation), alpha=1.5 (Moderate), alpha=2.0 (Strong)
    
    clusters = [
        {"name": "Terzan 5", "M": 2.0e6, "Rc": 0.16}, # High Density
        {"name": "47 Tuc",   "M": 1.0e6, "Rc": 0.36}, # High Density
        {"name": "M62",      "M": 1.0e6, "Rc": 0.18}, # High Density
        {"name": "M28",      "M": 5.0e5, "Rc": 0.24}, # High Density
        {"name": "M15",      "M": 5.0e5, "Rc": 0.14}, # Core Collapsed
        {"name": "M2",       "M": 6.0e5, "Rc": 0.32}, # Moderate
        {"name": "M5",       "M": 5.0e5, "Rc": 0.42}, # Moderate
        {"name": "M53",      "M": 3.0e5, "Rc": 0.65}, # Low Density
        {"name": "M13",      "M": 6.0e5, "Rc": 0.62}, # Low Density
        {"name": "M71",      "M": 2.0e4, "Rc": 0.63}  # Very Low Density
    ]
    
    # Simulation Parameters
    n_pulsars = 5000
    G_si = 6.674e-11
    M_sun_kg = 1.989e30
    pc_m = 3.086e16
    c_si = 2.998e8
    
    # Intrinsic Population (Field)
    mu_field = -19.76
    sigma_field = 0.64
    log_pdot_int = np.random.normal(mu_field, sigma_field, n_pulsars)
    pdot_int = 10**log_pdot_int
    
    # Periods (Mean 5ms)
    log_P_s = np.random.normal(np.log10(0.005), 0.3, n_pulsars)
    P_s = 10**log_P_s
    
    # Sensitivity Analysis on Mass Segregation
    alpha_values = [1.0, 1.5, 2.0, 2.5]
    colors = ['#A0C4FF', '#729BFF', '#2E86AB', '#004E89'] # Light to Dark Blue
    
    plt.figure(figsize=(9, 6))
    
    # Plot Observed Data (Red Points) - Fixed
    obs_data = [
        {"name": "Terzan 5", "rho": 3.6e7, "res": 0.28, "err": 0.03}, 
        {"name": "47 Tuc",   "rho": 5.0e4, "res": 0.12, "err": 0.03}, 
        {"name": "M62",      "rho": 4.0e5, "res": 0.33, "err": 0.08},
        {"name": "M28",      "rho": 8.0e4, "res": 0.28, "err": 0.09},
        {"name": "M15",      "rho": 4.0e5, "res": 0.28, "err": 0.10},
        {"name": "M5",       "rho": 1.6e4, "res": 0.02, "err": 0.04},
        {"name": "M53",      "rho": 2.5e2, "res": 0.02, "err": 0.01},
        {"name": "M13",      "rho": 6.0e2, "res": 0.02, "err": 0.14}
    ]
    obs_rhos = [d['rho'] for d in obs_data]
    obs_res = [d['res'] for d in obs_data]
    obs_err = [d['err'] for d in obs_data]
    
    plt.errorbar(obs_rhos, obs_res, yerr=obs_err, fmt='o', color='#E94F37', ecolor='#E94F37', 
                 capsize=3, markersize=5, elinewidth=0.8, markeredgewidth=0.8, label='Observed (Slope 0.33)')

    # Run Simulation for each alpha
    for idx, alpha in enumerate(alpha_values):
        print(f"\n--- Simulating Mass Segregation alpha={alpha} ---")
        sim_rhos = []
        sim_shifts = []
        
        for cl in clusters:
            M_cluster = cl['M']
            R_core = cl['Rc']
            
            # Core Density (approx)
            rho_0_si = 3 * (M_cluster * M_sun_kg) / (4 * np.pi * (R_core * pc_m)**3)
            rho_0_solar = rho_0_si * ((pc_m)**3) / M_sun_kg
            sim_rhos.append(rho_0_solar)

            # Positions with Segregation
            R_core_pulsar = R_core / alpha
            
            u = np.random.uniform(0, 0.99, n_pulsars)
            r_over_rc = np.sqrt(u**(2/3) / (1 - u**(2/3)))
            r_pulsar_m = r_over_rc * R_core_pulsar * pc_m
            cos_theta = np.random.uniform(-1, 1, n_pulsars)
            
            # Mean Field (Cluster Potential)
            M_enc_kg = (M_cluster * M_sun_kg) * (r_pulsar_m**3) / (r_pulsar_m**2 + (R_core * pc_m)**2)**(1.5)
            a_mean_si = G_si * M_enc_kg / (r_pulsar_m**2)
            a_los_mean_si = a_mean_si * cos_theta
            
            # Stochastic Field
            n_r = rho_0_si * (1 + (r_pulsar_m/(R_core*pc_m))**2)**(-2.5) / (0.5 * M_sun_kg)
            a_0 = 2.603 * G_si * (0.5 * M_sun_kg) * n_r**(2/3)
            a_stoch_si = stats.levy_stable.rvs(alpha=1.5, beta=0, scale=a_0, size=n_pulsars)
            a_stoch_si = np.clip(a_stoch_si, -1e-2, 1e-2)
            
            # Total
            a_tot_si = a_los_mean_si + a_stoch_si
            pdot_obs = pdot_int + P_s * (a_tot_si / c_si)
            
            shift = np.mean(np.log10(np.abs(pdot_obs))) - np.mean(log_pdot_int)
            sim_shifts.append(shift)

        # Fit Slope
        slope, intercept, _, _, _ = stats.linregress(np.log10(sim_rhos), sim_shifts)
        print(f"Alpha={alpha}: Slope={slope:.3f}")
        
        # Plot Line
        # Sort for plotting
        sort_idx = np.argsort(sim_rhos)
        plt.plot(np.array(sim_rhos)[sort_idx], np.array(sim_shifts)[sort_idx], 'o-', 
                 color=colors[idx], linewidth=0.8, markersize=4, 
                 label=f'Newtonian $\\alpha={alpha}$ (Slope {slope:.2f})')

    plt.xscale('log')
    plt.xlabel(r'Core Density $\rho_{core}$ ($M_\odot/\mathrm{pc}^3$)')
    plt.ylabel(r'Shift in $\log_{10}|\dot{P}|$ (dex)')
    plt.title('Impact of Mass Segregation on Density Scaling')
    plt.legend(frameon=True, framealpha=0.95, edgecolor='none', loc='upper left', fontsize=10)
    plt.grid(True, alpha=0.3)
    
    os.makedirs('site/figures', exist_ok=True)
    plt.savefig('site/figures/density_scaling_sensitivity.png', dpi=300)
    print("Saved site/figures/density_scaling_sensitivity.png")
    
    return

if __name__ == "__main__":
    run_density_scaling_simulation()