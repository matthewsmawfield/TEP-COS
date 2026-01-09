#!/usr/bin/env python3
"""
Generate figure comparing observed temporal shear vs microlensing simulation.
"""

import json
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# Set publication style
plt.rcParams.update({
    'font.family': 'serif',
    'font.size': 12,
    'axes.labelsize': 14,
    'axes.titlesize': 14,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
    'legend.fontsize': 12,
    'figure.dpi': 150,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
})

def load_results():
    """Load validation results."""
    results_path = Path(__file__).parent.parent.parent / 'results' / 'outputs' / 'step_3_2_validation_results.json'
    with open(results_path) as f:
        return json.load(f)

def plot_microlensing_comparison(data, output_dir):
    """Plot Observed vs Microlensing Gamma."""
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # Data extraction
    obs_gamma = data['scrambled_residuals']['gamma_orig']
    obs_sigma = data['scrambled_residuals']['sigma_orig']
    
    ml_gamma = data['microlensing_injection']['gamma_ml']
    ml_sigma = data['microlensing_injection']['sigma_ml']
    
    # We want to compare magnitudes effectively, but also signs.
    # The observed signal is large and negative. 
    # The microlensing signal is small and positive (in this specific injection).
    # Standard microlensing is stochastic, but usually creates small fluctuations.
    # We will plot them as points with error bars.
    
    # X positions
    x = [0, 1]
    labels = ['Observed Signal\n(DESJ0408 A-D)', 'Simulated Microlensing\n(0.3 mag, 2000d)']
    
    values = [obs_gamma, ml_gamma]
    errors = [obs_sigma, ml_sigma] # For ML, we might want to show the spread if we had multiple trials, but we have sigma from fit.
    
    # Colors
    colors = ['#C0392B', '#2980B9'] # Red for signal, Blue for simulation
    
    # Plot bars? Or points? Points with large error bars are better for "values".
    # But here the difference is huge.
    
    # Let's use bar chart to show magnitude comparison
    bars = ax.bar(x, values, yerr=errors, capsize=10, color=colors, alpha=0.7, edgecolor='black', width=0.5)
    
    # Add value labels
    for i, v in enumerate(values):
        offset = 20 if v > 0 else -30
        ax.text(i, v + offset, f"{v:.1f} days/dec", ha='center', fontweight='bold', color=colors[i])

    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel('Temporal Shear Γ [days/decade]')
    ax.set_title('Observed Signal vs. Microlensing Prediction')
    
    # Add a horizontal line at 0
    ax.axhline(0, color='black', linewidth=1)
    
    # Add a shaded region for what might be considered "reasonable" microlensing noise floor if we had it
    # But we rely on the point that it is 2 orders of magnitude different.
    
    # Annotate the factor difference
    diff_factor = abs(obs_gamma / ml_gamma)
    ax.text(0.5, -100, f"Discrepancy Factor: ~{diff_factor:.0f}x", ha='center', 
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="black", alpha=0.8))
    
    plt.tight_layout()
    fig.savefig(output_dir / 'cosmograil_microlensing_comparison.png')
    fig.savefig(output_dir / 'cosmograil_microlensing_comparison.pdf')
    plt.close(fig)
    print(f"Saved: cosmograil_microlensing_comparison.png/pdf")

def main():
    output_dir = Path(__file__).parent.parent.parent / 'results' / 'figures'
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print("Loading validation results...")
    data = load_results()
    
    print("Generating figures...")
    plot_microlensing_comparison(data, output_dir)
    
    print("\nDone!")

if __name__ == '__main__':
    main()
