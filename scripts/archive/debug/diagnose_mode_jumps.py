import json
import numpy as np
import matplotlib.pyplot as plt

def plot_mode_jumping():
    # Load Main Data
    with open('results/outputs/step_3_0_cosmograil_temporal_shear_v3_expanded.json', 'r') as f:
        data = json.load(f)
        
    # Load Q0957 Data
    try:
        with open('results/outputs/step_3_11_q0957_glendama_temporal_shear.json', 'r') as f:
            q0957_data = json.load(f)
    except:
        q0957_data = None

    systems = [
        ('DESJ0408', 'A-D', data['systems']['DESJ0408']['pairs']['A-D']),
        ('PG1115', 'B-C', data['systems']['PG1115']['pairs']['B-C']),
        ('J1206', 'A-B', data['systems']['J1206']['pairs']['A-B'])
    ]
    
    if q0957_data:
        systems.insert(0, ('Q0957', 'A-B', q0957_data['pairs']['A-B']))

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    axes = axes.flatten()

    for ax, (sys_name, pair_name, pair_data) in zip(axes, systems):
        ms = pair_data['multiscale']
        taus = sorted([int(k) for k in ms.keys()])
        log_taus = np.log10(taus)
        delays = []
        errs = []
        
        for t in taus:
            d = ms[str(t)].get('delay_days', np.nan)
            e = ms[str(t)].get('uncertainty_days', np.nan)
            delays.append(d)
            errs.append(e)
            
        ax.errorbar(log_taus, delays, yerr=errs, fmt='o-', capsize=5, label='Measured Delay')
        
        # Plot reported Gamma fit if available
        gamma = pair_data.get('gamma', {}).get('value')
        intercept = pair_data.get('gamma', {}).get('intercept')
        
        if gamma is not None:
            x_fit = np.linspace(min(log_taus), max(log_taus), 100)
            y_fit = gamma * x_fit + intercept
            ax.plot(x_fit, y_fit, 'r--', label=f'Fit: $\Gamma={gamma:.0f}$')
            
        ax.set_title(f"{sys_name} {pair_name}")
        ax.set_xlabel("log10(Tau)")
        ax.set_ylabel("Delay (days)")
        ax.legend()
        ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('results/figures/mode_jumping_diagnosis.png')
    print("Saved results/figures/mode_jumping_diagnosis.png")

if __name__ == "__main__":
    plot_mode_jumping()
