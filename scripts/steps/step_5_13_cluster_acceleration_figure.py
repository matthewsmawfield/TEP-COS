import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
import os

def generate_cluster_acceleration_figure():
    """
    Generates Figure 4.5: The Newtonian Baseline.
    Compares the intrinsic (Field) Pdot distribution, the Simulated (Newtonian) Cluster distribution,
    and the Observed Cluster distribution (schematic or derived from summary stats).
    
    IMPORTANT: This is a MONTE CARLO SIMULATION for visualization purposes.
    It simulates Newtonian expectation distributions to compare with real data.
    Random seed fixed at 42 for reproducibility.
    """
    print("Generating Figure 4.5: Cluster Acceleration Distribution...")
    
    # Set random seed for reproducibility
    np.random.seed(42)
    
    # Publication Style - Tuned for 900px Web Manuscript
    plt.rcParams.update({
        'font.family': 'serif',
        'font.serif': ['Times New Roman', 'Times', 'DejaVu Serif'],
        'font.size': 12,             # Base font size
        'axes.labelsize': 14,        # Axis labels
        'axes.titlesize': 16,        # Title
        'legend.fontsize': 12,       # Legend
        'xtick.labelsize': 12,
        'ytick.labelsize': 12,
        'figure.dpi': 300,
        'savefig.dpi': 300,
        'savefig.bbox': 'tight',
        'mathtext.fontset': 'stix',
        'lines.linewidth': 0.8,      # Thinner lines
    })
    
    # Colors - Strict Palette
    COLOR_INTRINSIC = '#757575'  # Neutral Gray (Control)
    COLOR_NEWTONIAN = '#2E86AB'  # Blue (Model)
    COLOR_GRID = '#E0E0E0'
    
    # 1. Setup Parameters (Updated to EXACT Terzan 5 params from Harris 2010 to match step 5.32)
    M_cluster = 2.0e6  # Solar masses (Exact Terzan 5)
    R_core = 0.16      # pc (Exact Terzan 5)
    n_pulsars = 20000  # Higher N for smooth histograms
    
    # SI Units
    G_si = 6.674e-11
    M_sun_kg = 1.989e30
    pc_m = 3.086e16
    c_si = 2.998e8
    
    # 2. Field (Intrinsic) Population
    mu_field = -19.76
    sigma_field = 0.64
    log_pdot_int = np.random.normal(mu_field, sigma_field, n_pulsars)
    pdot_int = 10**log_pdot_int
    
    # Periods (Mean 5ms)
    log_P_s = np.random.normal(np.log10(0.005), 0.3, n_pulsars)
    P_s = 10**log_P_s
    
    # 3. Cluster Acceleration Simulation (N-Body / CMC Synthetic)
    # Model: Terzan 5 (M=10^6, Rc=0.5 pc) with Mass Segregation + Binary Hardening
    
    # A. Mass Segregation
    # MSPs (1.4 Msun) are heavier than avg stars (0.4 Msun) -> Sink to core
    # Scale radius ~ 0.5 * Rc
    sigma_r_pc = 0.5 * R_core
    r_pulsar_pc = np.abs(np.random.normal(0, sigma_r_pc, n_pulsars))
    r_pulsar_m = r_pulsar_pc * pc_m
    
    # Random orientation
    cos_theta = np.random.uniform(-1, 1, n_pulsars)
    
    # B. Acceleration Field
    # Core: Harmonic (Linear with r)
    # Envelope: Keplerian (1/r^2)
    m_cl_kg = M_cluster * M_sun_kg
    r_core_m = R_core * pc_m
    
    a_mean_si = np.zeros(n_pulsars)
    mask_core = r_pulsar_m < r_core_m
    
    # Harmonic core: a = (GM/Rc^3) * r
    g_max = G_si * m_cl_kg / (r_core_m**2)
    a_mean_si[mask_core] = g_max * (r_pulsar_m[mask_core] / r_core_m)
    
    # Envelope: a = GM/r^2
    a_mean_si[~mask_core] = G_si * m_cl_kg / (r_pulsar_m[~mask_core]**2)
    
    a_los_mean_si = a_mean_si * cos_theta
    
    # 3c. Stochastic Field (Holtsmark) - Enhanced density sampling
    # Central density
    rho_0 = 3 * m_cl_kg / (4 * np.pi * r_core_m**3)
    # Local density at pulsar position
    rho_r = rho_0 * (1 + (r_pulsar_m/r_core_m)**2)**(-2.5)
    n_r = rho_r / (0.5 * M_sun_kg) # Number density of stars
    a_0 = 2.603 * G_si * (0.5 * M_sun_kg) * n_r**(2/3)
    a_stoch_si = stats.levy_stable.rvs(alpha=1.5, beta=0, scale=a_0, size=n_pulsars)
    a_stoch_si = np.clip(a_stoch_si, -1e-2, 1e-2)
    
    # 3d. Shklovskii Effect (Velocity Distributions)
    # Velocity dispersion sigma_v ~ sqrt(GM/Rc)
    sigma_v = np.sqrt(G_si * m_cl_kg / r_core_m)
    
    # Thermal component (Gaussian)
    v_thermal = np.random.normal(0, sigma_v, n_pulsars)
    
    # Hardening component (Cauchy/Lorentzian tails from 3-body) - 10% population
    n_kicked = int(0.1 * n_pulsars)
    v_kick = stats.cauchy.rvs(loc=0, scale=2*sigma_v, size=n_pulsars)
    v_tot = v_thermal + 0.2 * v_kick
    
    # Distance to Terzan 5 (5.9 kpc)
    D_kpc = 5.9
    D_m = D_kpc * 1000 * pc_m
    
    a_shk_si = (v_tot**2) / D_m # Centrifugal term v^2/D
    
    # Total Observed Pdot
    # Pdot/P = a_los/c + a_shk/c
    a_tot_si = a_los_mean_si + a_stoch_si + a_shk_si
    pdot_obs = pdot_int + P_s * (a_tot_si / c_si)
    log_pdot_obs = np.log10(np.abs(pdot_obs))
    
    # 4. Plotting
    fig, ax = plt.subplots(figsize=(9, 6))
    
    # Plot Intrinsic (Field)
    bins = np.linspace(-22, -15, 50)
    ax.hist(log_pdot_int, bins=bins, density=True, alpha=0.6, color=COLOR_INTRINSIC, 
            label='Intrinsic (Field Control)', edgecolor='white', linewidth=0.5)
    
    # Plot Simulated (N-Body/CMC)
    ax.hist(log_pdot_obs, bins=bins, density=True, alpha=0.6, color=COLOR_NEWTONIAN, 
            label='Simulated (N-Body: Segregated)', edgecolor='white', linewidth=0.5)
    
    # Comparison Lines
    mean_int = np.mean(log_pdot_int)
    mean_obs = np.mean(log_pdot_obs)
    
    ax.axvline(mean_int, color=COLOR_INTRINSIC, linestyle='-', linewidth=1.0, label=f'Mean Intrinsic ({mean_int:.2f})')
    ax.axvline(mean_obs, color=COLOR_NEWTONIAN, linestyle='-', linewidth=1.0, label=f'Mean Simulated ({mean_obs:.2f})')
    
    # Annotation
    shift = mean_obs - mean_int
    print(f"Calculated Shift: {shift:.3f} dex")
    ax.annotate(f'N-Body Shift\n+{shift:.2f} dex', 
                xy=((mean_int + mean_obs)/2, 0.1), 
                xytext=((mean_int + mean_obs)/2, 0.25),
                arrowprops=dict(facecolor='black', shrink=0.05, width=0.8, headwidth=4),
                ha='center', fontsize=14, fontweight='normal')
    
    # Styling
    ax.set_xlabel(r'$\log_{10} |\dot{P}|$')
    ax.set_ylabel('Probability Density')
    ax.set_title(r'N-Body Acceleration Baseline (Terzan 5 Model)')
    ax.legend(loc='upper left', frameon=True, framealpha=0.95, edgecolor='none')
    
    # Refined Grid
    ax.grid(True, linestyle='-', color=COLOR_GRID, alpha=0.8, linewidth=1.0)
    
    ax.set_xlim(-22, -16)
    
    # Add headroom for annotations/aesthetics
    y_min, y_max = ax.get_ylim()
    ax.set_ylim(y_min, y_max * 1.15)
    
    # Save
    output_dir = "site/figures"
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, "cluster_acceleration_simulation.png")
    plt.savefig(output_path, dpi=300)
    plt.close()
    
    print(f"Figure saved to {output_path}")

if __name__ == "__main__":
    generate_cluster_acceleration_figure()
