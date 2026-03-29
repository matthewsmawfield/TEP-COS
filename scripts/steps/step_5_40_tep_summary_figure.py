#!/usr/bin/env python3
"""
TEP Cosmology Summary Figure

Creates a publication-ready figure summarizing the TEP cosmology test results
across APOGEE, SDSS, and MaNGA surveys.

Author: M. Smawfield
Date: January 2026
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# Set style
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.size'] = 11
plt.rcParams['axes.linewidth'] = 1.2

# Results data
surveys = ['APOGEE\n(41k stars)', 'SDSS indices\n(361k galaxies)', 'SDSS ages\n(100k galaxies)', 'MaNGA\n(8.6k galaxies)']
correlations = [-0.123, -0.038, -0.270, -0.092]
errors = [0.01, 0.005, 0.01, 0.015]  # Approximate standard errors

# Create figure
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Panel A: Bar chart of correlations
ax = axes[0]
colors = ['#2E86AB', '#A23B72', '#F18F01', '#C73E1D']
x = np.arange(len(surveys))
bars = ax.bar(x, correlations, yerr=errors, capsize=5, color=colors, edgecolor='black', linewidth=1.2, alpha=0.85)

ax.axhline(0, color='gray', linestyle='--', linewidth=1)
ax.set_ylabel('Correlation r(age, Φ) at fixed formation epoch', fontsize=12)
ax.set_xticks(x)
ax.set_xticklabels(surveys, fontsize=10)
ax.set_ylim(-0.35, 0.05)
ax.set_title('A. TEP Cosmology Test Results', fontsize=13, fontweight='bold')

# Add significance stars
for i, (corr, err) in enumerate(zip(correlations, errors)):
    if abs(corr) > 3 * err:
        ax.text(i, corr - 0.03, '***', ha='center', fontsize=14, fontweight='bold')

# Add TEP prediction arrow
ax.annotate('TEP prediction:\nyounger at deeper Φ', xy=(3.5, -0.15), xytext=(3.5, 0.02),
            fontsize=9, ha='center', va='bottom',
            arrowprops=dict(arrowstyle='->', color='darkgreen', lw=2))

# Panel B: Schematic of the test
ax = axes[1]
ax.set_xlim(0, 10)
ax.set_ylim(0, 10)
ax.axis('off')
ax.set_title('B. The TEP Age-Nucleosynthesis Test', fontsize=13, fontweight='bold')

# Draw schematic
# Two galaxies/stars
circle1 = plt.Circle((2.5, 7), 1.2, color='#2E86AB', alpha=0.7)
circle2 = plt.Circle((7.5, 7), 0.8, color='#A23B72', alpha=0.7)
ax.add_patch(circle1)
ax.add_patch(circle2)

ax.text(2.5, 7, 'Deep Φ\n(high σ)', ha='center', va='center', fontsize=9, color='white', fontweight='bold')
ax.text(7.5, 7, 'Shallow Φ\n(low σ)', ha='center', va='center', fontsize=8, color='white', fontweight='bold')

# Arrows and labels
ax.annotate('', xy=(2.5, 5.5), xytext=(2.5, 5.8), arrowprops=dict(arrowstyle='->', lw=2, color='darkblue'))
ax.annotate('', xy=(7.5, 5.5), xytext=(7.5, 5.8), arrowprops=dict(arrowstyle='->', lw=2, color='darkred'))

# Time flow indicators
ax.text(2.5, 5.2, 'Slower τ', ha='center', fontsize=10, color='darkblue', fontweight='bold')
ax.text(7.5, 5.2, 'Faster τ', ha='center', fontsize=10, color='darkred', fontweight='bold')

# Observables box
box_props = dict(boxstyle='round,pad=0.3', facecolor='lightyellow', edgecolor='black', linewidth=1.5)

ax.text(5, 3.5, 'Same [Mg/Fe] or [α/M]\n(same formation epoch)', ha='center', va='center', 
        fontsize=10, bbox=box_props)

ax.text(2.5, 1.5, 'Appears\nYOUNGER', ha='center', va='center', fontsize=11, 
        color='darkblue', fontweight='bold')
ax.text(7.5, 1.5, 'Appears\nOLDER', ha='center', va='center', fontsize=11, 
        color='darkred', fontweight='bold')

# Connecting lines
ax.plot([2.5, 2.5], [4.0, 2.2], 'b--', lw=1.5, alpha=0.7)
ax.plot([7.5, 7.5], [4.0, 2.2], 'r--', lw=1.5, alpha=0.7)

# Legend/explanation
explanation = """TEP Prediction:
At fixed nucleosynthesis ratio
(same coordinate-time formation),
objects in deeper potentials
experience slower proper time
→ less stellar evolution
→ younger spectroscopic age"""

ax.text(5, 0.3, explanation, ha='center', va='bottom', fontsize=9,
        bbox=dict(boxstyle='round', facecolor='lightgray', alpha=0.5))

plt.tight_layout()

# Save
fig_path = 'results/figures/tep_cosmology_summary.png'
plt.savefig(fig_path, dpi=150, bbox_inches='tight', facecolor='white')
plt.savefig(fig_path.replace('.png', '.pdf'), bbox_inches='tight')
print(f'Figure saved: {fig_path}')

plt.close()
