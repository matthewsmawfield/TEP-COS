#!/usr/bin/env python3
"""
Step 54: Environment Axis Scan

Replace central density (rho_c) with alternative environmental variables and
recompute the density-scaling slope Gamma for each.  This tests whether TEP
tracks potential depth, acceleration scale, escape velocity, or relaxation time
rather than raw density.

Author: M. Smawfield
Date: 2026-06-23
"""

import csv
import json
import math
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from scipy import stats

REPO_ROOT = Path(__file__).resolve().parents[2]
RESULTS_DIR = REPO_ROOT / "results" / "outputs"
STEP_5_10_CSV = RESULTS_DIR / "step_02_pulsar_population_controls.csv"

OUT_JSON = RESULTS_DIR / "step_54_environment_axis_scan.json"
OUT_MD = RESULTS_DIR / "step_54_environment_axis_scan.md"

G = 4.302e-3  # pc (km/s)^2 / M_sun

CLUSTER_PARAMS = {
    "Terzan 5":         {"M": 2.0e6, "Rc": 0.16, "Rt": 5.0,  "rho_c": 5.50},
    "47 Tuc (NGC 104)": {"M": 1.0e6, "Rc": 0.36, "Rt": 42.0, "rho_c": 4.88},
    "NGC 6517":         {"M": 2.0e5, "Rc": 0.06, "Rt": 5.0,  "rho_c": 5.80},
    "M28 (NGC 6626)":   {"M": 5.0e5, "Rc": 0.24, "Rt": 12.0, "rho_c": 4.52},
    "M62 (NGC 6266)":   {"M": 1.0e6, "Rc": 0.18, "Rt": 8.0,  "rho_c": 5.16},
    "M13 (NGC 6205)":   {"M": 6.0e5, "Rc": 0.62, "Rt": 25.0, "rho_c": 3.79},
    "M15 (NGC 7078)":   {"M": 5.0e5, "Rc": 0.14, "Rt": 21.0, "rho_c": 5.05},
    "M5 (NGC 5904)":    {"M": 5.0e5, "Rc": 0.42, "Rt": 28.0, "rho_c": 3.53},
    "Terzan 1":         {"M": 1.5e5, "Rc": 0.10, "Rt": 4.0,  "rho_c": 5.00},
    "NGC 6752":         {"M": 3.0e5, "Rc": 0.17, "Rt": 25.0, "rho_c": 4.30},
    "M2 (NGC 7089)":    {"M": 6.0e5, "Rc": 0.32, "Rt": 21.0, "rho_c": 4.15},
    "Omega Centauri (NGC 5139)": {"M": 4.0e6, "Rc": 2.37, "Rt": 57.0, "rho_c": 3.12},
    "M53 (NGC 5024)":   {"M": 3.0e5, "Rc": 0.65, "Rt": 22.0, "rho_c": 2.96},
    "M3 (NGC 5272)":    {"M": 5.0e5, "Rc": 0.37, "Rt": 38.0, "rho_c": 3.68},
    "M71 (NGC 6838)":   {"M": 2.0e4, "Rc": 0.63, "Rt": 8.0,  "rho_c": 2.29},
    "NGC 6397":         {"M": 1.0e5, "Rc": 0.05, "Rt": 15.0, "rho_c": 5.68},
    "NGC 1851":         {"M": 3.0e5, "Rc": 0.09, "Rt": 11.0, "rho_c": 5.09},
    "NGC 6522":         {"M": 2.0e5, "Rc": 0.05, "Rt": 5.0,  "rho_c": 5.50},
    "NGC 6544":         {"M": 5.0e4, "Rc": 0.05, "Rt": 3.0,  "rho_c": 5.20},
    "NGC 6624":         {"M": 2.0e5, "Rc": 0.06, "Rt": 6.0,  "rho_c": 5.60},
    "NGC 6760":         {"M": 2.0e5, "Rc": 0.34, "Rt": 8.0,  "rho_c": 3.80},
    "M22 (NGC 6656)":   {"M": 5.0e5, "Rc": 1.33, "Rt": 32.0, "rho_c": 2.97},
    "M80 (NGC 6093)":   {"M": 4.0e5, "Rc": 0.15, "Rt": 13.0, "rho_c": 4.79},
    "M92 (NGC 6341)":   {"M": 3.0e5, "Rc": 0.26, "Rt": 15.0, "rho_c": 4.30},
    "NGC 6712":         {"M": 1.5e5, "Rc": 0.33, "Rt": 7.0,  "rho_c": 3.70},
    "NGC 6652":         {"M": 1.0e5, "Rc": 0.10, "Rt": 5.0,  "rho_c": 4.50},
    "M14 (NGC 6402)":   {"M": 1.0e6, "Rc": 0.78, "Rt": 18.0, "rho_c": 3.44},
    "NGC 6539":         {"M": 3.0e5, "Rc": 0.60, "Rt": 10.0, "rho_c": 3.30},
    "M4 (NGC 6121)":    {"M": 1.0e5, "Rc": 0.83, "Rt": 33.0, "rho_c": 2.85},
    "NGC 6440":         {"M": 1.5e5, "Rc": 0.12, "Rt": 5.5,  "rho_c": 5.10},
    "NGC 6441":         {"M": 8.0e5, "Rc": 0.20, "Rt": 12.0, "rho_c": 5.00},
    "NGC 6316":         {"M": 1.2e5, "Rc": 0.15, "Rt": 5.0,  "rho_c": 4.80},
    "M30 (NGC 7099)":   {"M": 2.5e5, "Rc": 0.25, "Rt": 18.0, "rho_c": 4.20},
}


def load_data():
    gc_rows = []
    field_rows = []
    with open(STEP_5_10_CSV, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                r = {
                    'environment': row['environment'],
                    'cluster': row['cluster'],
                    'logPdot_abs': float(row['logPdot_abs']),
                }
                if r['environment'] == 'globular_cluster':
                    gc_rows.append(r)
                elif r['environment'] == 'field':
                    field_rows.append(r)
            except (ValueError, KeyError):
                continue
    return gc_rows, field_rows


def compute_cluster_means(gc_rows):
    clusters = defaultdict(list)
    for r in gc_rows:
        clusters[r['cluster']].append(r)
    means = {}
    for name, rows in clusters.items():
        means[name] = float(np.mean([r['logPdot_abs'] for r in rows]))
    return means


def compute_env_variables():
    """Compute alternative environmental variables for each cluster."""
    env = {}
    for name, p in CLUSTER_PARAMS.items():
        M = p['M']
        Rc = p['Rc']
        Rt = p['Rt']
        rho_c = p['rho_c']

        # Mass density proxy (M / Rc^3 in Msun / pc^3)
        mass_density = M / (Rc ** 3) if Rc > 0 else 0

        # Potential depth (M / Rc)
        potential_depth = M / Rc if Rc > 0 else 0

        # Acceleration scale (M / Rc^2)
        accel_scale = M / (Rc ** 2) if Rc > 0 else 0

        # Escape velocity (km/s)
        v_esc = math.sqrt(2 * G * M / Rt) if Rt > 0 else 0

        # Velocity dispersion (km/s)
        sigma_v = math.sqrt(G * M / Rc) if Rc > 0 else 0

        # Relaxation time proxy (R_c / sigma_v in Myr-like units, arbitrary scale)
        # Using a simple proxy: t_relax ~ N^(1/2) * R_c^(3/2) / (M^(1/2) * ln(N))
        # For simplicity, use Rc / sigma_v
        relax_proxy = Rc / sigma_v if sigma_v > 0 else 0

        env[name] = {
            'rho_c': rho_c,
            'log_rho_c': rho_c,  # already log10 scale in CLUSTER_PARAMS
            'log_mass_density': math.log10(mass_density) if mass_density > 0 else None,
            'log_potential_depth': math.log10(potential_depth) if potential_depth > 0 else None,
            'log_accel_scale': math.log10(accel_scale) if accel_scale > 0 else None,
            'log_v_esc': math.log10(v_esc) if v_esc > 0 else None,
            'log_sigma_v': math.log10(sigma_v) if sigma_v > 0 else None,
            'log_relax_proxy': math.log10(relax_proxy) if relax_proxy > 0 else None,
        }
    return env


def run_regression(cluster_means, env, x_key):
    x = []
    y = []
    for name in cluster_means:
        if name not in env:
            continue
        val = env[name].get(x_key)
        if val is None:
            continue
        x.append(val)
        y.append(cluster_means[name])

    if len(x) < 3:
        return None

    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
    return {
        'x_key': x_key,
        'n': len(x),
        'slope': float(slope),
        'intercept': float(intercept),
        'r_value': float(r_value),
        'r_squared': float(r_value ** 2),
        'p_value': float(p_value),
        'std_err': float(std_err),
    }


def main():
    print("=" * 78)
    print("STEP 5.74: ENVIRONMENT AXIS SCAN")
    print("=" * 78)
    print()

    gc_rows, field_rows = load_data()
    cluster_means = compute_cluster_means(gc_rows)
    env = compute_env_variables()

    print(f"Loaded {len(gc_rows)} GC pulsars across {len(cluster_means)} clusters")
    print()

    x_vars = [
        ('log_rho_c', 'Central luminosity density (current)'),
        ('log_mass_density', 'Mass density M/Rc^3'),
        ('log_potential_depth', 'Potential depth M/Rc'),
        ('log_accel_scale', 'Acceleration scale M/Rc^2'),
        ('log_v_esc', 'Escape velocity v_esc'),
        ('log_sigma_v', 'Velocity dispersion sigma_v'),
        ('log_relax_proxy', 'Relaxation time proxy Rc/sigma_v'),
    ]

    print(f"{'Variable':<35} {'Gamma':>8} {'StdErr':>8} {'R^2':>8} {'p-value':>10}")
    print("-" * 78)

    results = {
        'timestamp_utc': datetime.now(timezone.utc).isoformat(),
        'regressions': [],
    }

    for x_key, label in x_vars:
        reg = run_regression(cluster_means, env, x_key)
        if reg:
            results['regressions'].append({**reg, 'label': label})
            print(f"{label:<35} {reg['slope']:>+7.3f} {reg['std_err']:>8.3f} {reg['r_squared']:>8.3f} {reg['p_value']:>10.4f}")
        else:
            print(f"{label:<35} {'N/A':>8} {'N/A':>8} {'N/A':>8} {'N/A':>10}")

    # Sort by R^2 to find the best predictor
    sorted_regs = sorted([r for r in results['regressions'] if r is not None], key=lambda r: -r['r_squared'])

    print()
    print("Ranked by R^2 (best predictor first):")
    for i, reg in enumerate(sorted_regs, 1):
        print(f"  {i}. {reg['label']}: Gamma = {reg['slope']:+.3f}, R^2 = {reg['r_squared']:.3f}")

    # Save
    with open(OUT_JSON, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\nSaved JSON: {OUT_JSON}")

    md = f"""# Step 54: Environment Axis Scan

**Generated:** {results['timestamp_utc']}

## Results

| Variable | Gamma | Std Err | R^2 | p-value |
|----------|-------|---------|-----|---------|
"""
    for reg in results['regressions']:
        md += f"| {reg['label']} | {reg['slope']:+.3f} | {reg['std_err']:.3f} | {reg['r_squared']:.3f} | {reg['p_value']:.4f} |\n"

    md += """
## Interpretation

- If TEP is **density-driven**, the strongest correlation should be with
  **central luminosity density** or **mass density**.
- If TEP is **potential-driven** (coherence / potential depth), the signal
  should strengthen against **M/Rc** or **escape velocity**.
- If TEP is **dynamical-state-driven**, look for correlation with
  **velocity dispersion** or **relaxation time**.
- The axis with the highest R^2 and steepest positive slope is the preferred
  environmental variable for TEP.
"""

    with open(OUT_MD, 'w') as f:
        f.write(md)
    print(f"Saved Markdown: {OUT_MD}")
    print("Done.")


if __name__ == '__main__':
    main()
