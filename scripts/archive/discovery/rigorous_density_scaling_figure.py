import numpy as np
import matplotlib.pyplot as plt
import json
import os
from scipy import stats

def generate_rigorous_density_scaling():
    print("Generating Rigorous Density Scaling Figure...")
    
    # 1. Load Data
    # Newtonian (Exact Parameters)
    with open('results/outputs/step_5_32_full_density_scaling.json', 'r') as f:
        newton_data = json.load(f)
        
    # Observed (Controlled Residuals)
    with open('results/outputs/step_5_31_per_cluster_controlled_residuals.json', 'r') as f:
        obs_data = json.load(f)

    # 2. Process Data
    # Newtonian
    n_rho = [c['rho_c_log'] for c in newton_data['clusters']]
    n_shift = [c['shift'] for c in newton_data['clusters']]
    n_slope = newton_data['slope']
    n_intercept = newton_data['intercept']
    
    # Observed
    # Map rho from Newtonian data to observed clusters
    rho_map = {c['name']: c['rho_c_log'] for c in newton_data['clusters']}
    
    o_rho = []
    o_shift = []
    o_err = []
    o_names = []
    
    for name, data in obs_data['clusters'].items():
        if name in rho_map:
            o_rho.append(rho_map[name])
            o_shift.append(data['controlled_residual'])
            o_err.append(data['residual_sem'])
            o_names.append(name)
            
    # Re-calculate observed regression for consistency
    o_slope, o_intercept, _, _, _ = stats.linregress(o_rho, o_shift)
    
    # 3. Plotting Setup
    plt.rcParams.update({
        'font.family': 'serif',
        'font.serif': ['Times New Roman', 'Times', 'DejaVu Serif'],
        'font.size': 14,
        'axes.labelsize': 16,
        'axes.titlesize': 18,
        'legend.fontsize': 12,
        'xtick.labelsize': 14,
        'ytick.labelsize': 14,
        'figure.dpi': 300,
        'lines.linewidth': 1.5,
    })
    
    fig, ax = plt.subplots(figsize=(9, 6))
    
    # Colors
    COLOR_NEWTON = '#2E86AB'  # Blue
    COLOR_OBS    = '#D62828'  # Red
    
    # 4. Plot Newtonian Data (Blue)
    ax.scatter(n_rho, n_shift, color=COLOR_NEWTON, alpha=0.5, s=60, edgecolors='none', label='Newtonian Simulation (Exact Params)')
    
    # Regression Line
    x_range = np.linspace(2.0, 6.0, 100)
    y_newton = n_slope * x_range + n_intercept
    ax.plot(x_range, y_newton, color=COLOR_NEWTON, linestyle='--', linewidth=2, label=f'Newtonian Trend (Slope = {n_slope:.2f})')
    
    # 5. Plot Observed Data (Red)
    ax.errorbar(o_rho, o_shift, yerr=o_err, fmt='o', color=COLOR_OBS, ecolor=COLOR_OBS, 
                elinewidth=1.0, capsize=3, markersize=6, label='Observed Residuals')
    
    # Regression Line
    y_obs = o_slope * x_range + o_intercept
    ax.plot(x_range, y_obs, color=COLOR_OBS, linestyle='-', linewidth=2, label=f'Observed Scaling (Slope = {o_slope:.2f})')
    
    # 6. Styling
    ax.set_xlabel(r'Log Central Density $\log_{10}(\rho_c) [L_\odot/pc^3]$')
    ax.set_ylabel(r'Acceleration Shift (dex)')
    ax.set_title('N-Body/CMC Density Scaling (Mass Segregated)')
    
    ax.grid(True, linestyle='-', alpha=0.2)
    ax.legend(loc='upper left', frameon=True, framealpha=0.95)
    
    # Annotations for key clusters
    for i, name in enumerate(o_names):
        if name in ["Terzan 5", "47 Tuc (NGC 104)", "M15 (NGC 7078)", "Omega Centauri (NGC 5139)"]:
            short_name = name.split('(')[0].strip()
            ax.annotate(short_name, (o_rho[i], o_shift[i]), 
                       xytext=(0, 10), textcoords='offset points',
                       ha='center', fontsize=10, color=COLOR_OBS)

    # 7. Save
    output_dir = "site/figures"
    os.makedirs(output_dir, exist_ok=True)
    save_path = os.path.join(output_dir, "density_scaling.png")
    plt.savefig(save_path, bbox_inches='tight')
    print(f"Figure saved to {save_path}")

if __name__ == "__main__":
    generate_rigorous_density_scaling()
