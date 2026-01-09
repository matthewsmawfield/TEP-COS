#!/usr/bin/env python3
"""
Generate publication-quality figures for COSMOGRAIL temporal shear analysis.
"""

import json
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# Set publication style
plt.rcParams.update({
    'font.family': 'serif',
    'font.size': 10,
    'axes.labelsize': 11,
    'axes.titlesize': 11,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'legend.fontsize': 9,
    'figure.dpi': 150,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
})

# Colors
COLORS = {
    'significant': '#2980B9',  # Blue
    'marginal': '#F39C12',     # Orange
    'null': '#95A5A6',         # Grey
    'fit': '#C0392B',          # Red
}


def load_results():
    """Load temporal shear results."""
    results_path = Path(__file__).parent.parent.parent / 'results' / 'outputs' / 'step_3_0_cosmograil_temporal_shear_opB_modejump.json'
    with open(results_path) as f:
        return json.load(f)


def plot_gamma_distribution(data, output_dir):
    """Plot distribution of gamma values across all pairs."""
    fig, ax = plt.subplots(figsize=(7, 4))
    
    gammas = []
    sigmas = []
    labels = []
    colors = []
    
    for sys_id, sys_data in data['systems'].items():
        for pair_id, pair_data in sys_data['pairs'].items():
            gamma = pair_data['gamma']
            if gamma['value'] is None or not np.isfinite(gamma['value']):
                continue
            
            g = gamma['value']
            u = gamma['uncertainty']
            s = abs(g / u) if u and u > 0 else 0
            
            gammas.append(g)
            sigmas.append(s)
            labels.append(f"{sys_id}\n{pair_id}")
            
            if s >= 3:
                colors.append(COLORS['significant'])
            elif s >= 2:
                colors.append(COLORS['marginal'])
            else:
                colors.append(COLORS['null'])
    
    # Sort by sigma
    order = np.argsort(sigmas)[::-1]
    gammas = [gammas[i] for i in order]
    sigmas = [sigmas[i] for i in order]
    labels = [labels[i] for i in order]
    colors = [colors[i] for i in order]
    
    x = np.arange(len(gammas))
    bars = ax.bar(x, sigmas, color=colors, edgecolor='black', linewidth=0.5)
    
    # Add significance thresholds
    ax.axhline(3, color='red', linestyle='--', linewidth=1, label='3σ threshold')
    ax.axhline(2, color='orange', linestyle='--', linewidth=1, label='2σ threshold')
    
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45, ha='right', fontsize=8)
    ax.set_ylabel('Significance (|Γ|/σ)')
    ax.set_xlabel('Lens System / Image Pair')
    ax.set_title('COSMOGRAIL Temporal Shear Detection Significance')
    ax.legend(loc='upper right')
    ax.set_ylim(0, max(sigmas) * 1.1)
    
    # Annotate top 3
    for i in range(min(3, len(sigmas))):
        if sigmas[i] >= 3:
            ax.annotate(f'{sigmas[i]:.1f}σ', (x[i], sigmas[i] + 0.2), 
                       ha='center', fontsize=8, fontweight='bold')
    
    plt.tight_layout()
    fig.savefig(output_dir / 'cosmograil_gamma_significance.png')
    fig.savefig(output_dir / 'cosmograil_gamma_significance.pdf')
    plt.close(fig)
    print(f"Saved: cosmograil_gamma_significance.png/pdf")


def plot_multiscale_delays(data, output_dir):
    """Plot delay vs timescale for significant detections."""
    fig, axes = plt.subplots(2, 3, figsize=(10, 6), sharey=False)
    axes = axes.flatten()
    
    significant_pairs = [
        ('DESJ0408', 'A-D'),
        ('DESJ0408', 'B-D'),
        ('PG1115', 'B-C'),
        ('PG1115', 'A-B'),
        ('J1206', 'A-B'),
    ]
    
    tau_values = [5, 10, 20, 40, 80, 160]
    log_tau = np.log10(tau_values)
    
    for idx, (sys_id, pair_id) in enumerate(significant_pairs):
        ax = axes[idx]
        pair_data = data['systems'][sys_id]['pairs'][pair_id]
        
        delays = []
        valid_log_tau = []
        
        for tau in tau_values:
            ms = pair_data['multiscale'].get(str(tau), {})
            d = ms.get('delay_days')
            if d is not None and np.isfinite(d):
                delays.append(d)
                valid_log_tau.append(np.log10(tau))
        
        if len(delays) < 2:
            continue
        
        # Plot data points
        ax.scatter(valid_log_tau, delays, s=60, c=COLORS['significant'], 
                  edgecolor='black', linewidth=0.5, zorder=5)
        
        # Fit line
        gamma = pair_data['gamma']['value']
        intercept = pair_data['gamma']['intercept']
        x_fit = np.linspace(min(valid_log_tau) - 0.1, max(valid_log_tau) + 0.1, 100)
        y_fit = gamma * x_fit + intercept
        ax.plot(x_fit, y_fit, color=COLORS['fit'], linewidth=2, 
               label=f'Γ = {gamma:.0f} days/decade')
        
        # Broadband delay reference
        bb_delay = pair_data['broadband']['delay_days']
        ax.axhline(bb_delay, color='gray', linestyle=':', linewidth=1, alpha=0.7)
        
        ax.set_xlabel('log₁₀(τ) [days]')
        if idx % 3 == 0:
            ax.set_ylabel('Time Delay [days]')
        ax.set_title(f'{sys_id} {pair_id}')
        ax.legend(loc='best', fontsize=8)
        ax.grid(True, alpha=0.3)
    
    # Hide unused subplot
    if len(significant_pairs) < 6:
        axes[5].axis('off')
    
    plt.tight_layout()
    fig.savefig(output_dir / 'cosmograil_multiscale_delays.png')
    fig.savefig(output_dir / 'cosmograil_multiscale_delays.pdf')
    plt.close(fig)
    print(f"Saved: cosmograil_multiscale_delays.png/pdf")


def plot_summary_panel(data, output_dir):
    """Create a summary panel figure."""
    fig = plt.figure(figsize=(10, 8))
    
    # Panel A: Gamma vs broadband delay
    ax1 = fig.add_subplot(2, 2, 1)
    
    gammas = []
    bb_delays = []
    sigmas = []
    
    for sys_id, sys_data in data['systems'].items():
        for pair_id, pair_data in sys_data['pairs'].items():
            gamma = pair_data['gamma']
            if gamma['value'] is None or not np.isfinite(gamma['value']):
                continue
            
            g = gamma['value']
            u = gamma['uncertainty']
            s = abs(g / u) if u and u > 0 else 0
            bb = pair_data['broadband']['delay_days']
            
            gammas.append(g)
            bb_delays.append(bb)
            sigmas.append(s)
    
    colors = [COLORS['significant'] if s >= 3 else COLORS['marginal'] if s >= 2 else COLORS['null'] 
              for s in sigmas]
    
    ax1.scatter(bb_delays, gammas, c=colors, s=60, edgecolor='black', linewidth=0.5)
    ax1.axhline(0, color='gray', linestyle='--', linewidth=1)
    ax1.set_xlabel('Broadband Delay [days]')
    ax1.set_ylabel('Γ [days/decade]')
    ax1.set_title('A. Temporal Shear vs Broadband Delay')
    ax1.grid(True, alpha=0.3)
    
    # Panel B: Histogram of gamma
    ax2 = fig.add_subplot(2, 2, 2)
    valid_gammas = [g for g in gammas if np.isfinite(g)]
    ax2.hist(valid_gammas, bins=15, color=COLORS['null'], edgecolor='black', alpha=0.7)
    ax2.axvline(0, color='red', linestyle='--', linewidth=2, label='GR prediction (Γ=0)')
    ax2.set_xlabel('Γ [days/decade]')
    ax2.set_ylabel('Count')
    ax2.set_title('B. Distribution of Temporal Shear')
    ax2.legend()
    
    # Panel C: Significance histogram
    ax3 = fig.add_subplot(2, 2, 3)
    ax3.hist(sigmas, bins=np.arange(0, 8, 0.5), color=COLORS['significant'], 
            edgecolor='black', alpha=0.7)
    ax3.axvline(3, color='red', linestyle='--', linewidth=2, label='3σ threshold')
    ax3.axvline(2, color='orange', linestyle='--', linewidth=2, label='2σ threshold')
    ax3.set_xlabel('Significance (|Γ|/σ)')
    ax3.set_ylabel('Count')
    ax3.set_title('C. Detection Significance Distribution')
    ax3.legend()
    
    # Panel D: Summary text
    ax4 = fig.add_subplot(2, 2, 4)
    ax4.axis('off')
    
    n_sig_3 = sum(1 for s in sigmas if s >= 3)
    n_sig_2 = sum(1 for s in sigmas if s >= 2)
    n_total = len(sigmas)
    
    summary_text = f"""
COSMOGRAIL Temporal Shear Analysis
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

Systems analyzed: {data['summary']['n_systems']}
Image pairs: {n_total}

Significant detections:
  • >3σ: {n_sig_3} pairs
  • >2σ: {n_sig_2} pairs

Top detections:
  • DESJ0408 A-D: Γ = -333 ± 53 (6.3σ)
  • DESJ0408 B-D: Γ = -129 ± 21 (6.1σ)
  • J1206 A-B:    Γ = -103 ± 30 (3.4σ)

Interpretation:
Under GR, Γ = 0 (constant delay).
Under TEP-GL, Γ ≠ 0 (scale-dependent).

Status: CANDIDATE DETECTION
"""
    
    ax4.text(0.05, 0.95, summary_text, transform=ax4.transAxes, 
            fontsize=10, verticalalignment='top', fontfamily='monospace',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    ax4.set_title('D. Summary')
    
    plt.tight_layout()
    fig.savefig(output_dir / 'cosmograil_summary_panel.png')
    fig.savefig(output_dir / 'cosmograil_summary_panel.pdf')
    plt.close(fig)
    print(f"Saved: cosmograil_summary_panel.png/pdf")


def main():
    output_dir = Path(__file__).parent.parent.parent / 'results' / 'figures'
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print("Loading COSMOGRAIL temporal shear results...")
    data = load_results()
    
    print("Generating figures...")
    plot_gamma_distribution(data, output_dir)
    plot_multiscale_delays(data, output_dir)
    plot_summary_panel(data, output_dir)
    
    print("\nDone!")


if __name__ == '__main__':
    main()
