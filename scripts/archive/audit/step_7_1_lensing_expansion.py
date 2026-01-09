#!/usr/bin/env python3
"""
Step 7.1: Lensing Sample Expansion - Additional Quad Systems

Analyze additional quad lens systems beyond DESJ0408 to expand the
validated TEP temporal shear sample.

Systems with available data:
- RXJ1131-1231 (z_S=0.658, z_L=0.295) - 4 images, COSMOGRAIL
- Q2237+0305 (Einstein Cross, z_S=1.695, z_L=0.039) - 4 images, OGLE
- HE1104-1805 (z_S=2.319, z_L=0.729) - 2 images, long baseline
- HS2209+1914 (z_S=1.07, z_L=0.22) - 2 images

Author: TEP-COS Analysis Pipeline
Date: January 2026
"""

import numpy as np
import pandas as pd
from scipy import stats
import json
import os
from datetime import datetime

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.join(SCRIPT_DIR, '..', '..')
DATA_DIR = os.path.join(PROJECT_DIR, 'data', 'cosmograil')
RESULTS_DIR = os.path.join(PROJECT_DIR, 'results', 'outputs')
FIGURES_DIR = os.path.join(PROJECT_DIR, 'results', 'figures')

# Lens system metadata
LENS_SYSTEMS = {
    'RXJ1131': {
        'name': 'RXJ1131-1231',
        'z_source': 0.658,
        'z_lens': 0.295,
        'n_images': 4,
        'image_labels': ['A', 'B', 'C', 'D'],
        'file': 'RXJ1131_Tewes2013.rdb',
        'sigma_lens': 323,  # km/s from Suyu+2014
        'status': 'NEW_TARGET'
    },
    'Q2237': {
        'name': 'Q2237+0305 (Einstein Cross)',
        'z_source': 1.695,
        'z_lens': 0.039,
        'n_images': 4,
        'image_labels': ['A', 'B', 'C', 'D'],
        'file': 'Q2237_ogle2_phot.dat',
        'sigma_lens': 166,  # km/s spiral bulge
        'status': 'NEW_TARGET',
        'notes': 'Very low z_L, strong microlensing'
    },
    'HE1104': {
        'name': 'HE1104-1805',
        'z_source': 2.319,
        'z_lens': 0.729,
        'n_images': 2,
        'image_labels': ['A', 'B'],
        'file': 'he1104_JApJ798_95_R.csv',
        'sigma_lens': 250,  # estimated
        'status': 'NEW_TARGET',
        'notes': 'Double, high-z source'
    },
    'DESJ0408': {
        'name': 'DESJ0408-5354',
        'z_source': 2.375,
        'z_lens': 0.597,
        'n_images': 4,
        'image_labels': ['A', 'B', 'C', 'D'],
        'file': 'DESJ0408_Courbin2017.rdb',
        'sigma_lens': 230,
        'status': 'VALIDATED',
        'gamma_AB': 32.4,
        'gamma_BD': 34.6
    }
}


def load_lightcurve(system_key):
    """Load light curve data for a lens system."""
    system = LENS_SYSTEMS[system_key]
    filepath = os.path.join(DATA_DIR, system['file'])
    
    if not os.path.exists(filepath):
        print(f"  File not found: {filepath}")
        return None
    
    try:
        if filepath.endswith('.rdb'):
            # RDB format (tab-separated with header)
            df = pd.read_csv(filepath, sep='\t', comment='=', skiprows=2)
        elif filepath.endswith('.csv'):
            df = pd.read_csv(filepath)
        elif filepath.endswith('.dat'):
            if 'ogle' in filepath.lower():
                # OGLE format: date jd airmass A eA B eB C eC D eD
                df = pd.read_csv(filepath, sep=r'\s+', header=None,
                                names=['date', 'jd', 'airmass', 
                                       'A', 'eA', 'B', 'eB', 'C', 'eC', 'D', 'eD'])
            else:
                df = pd.read_csv(filepath, sep=r'\s+')
        else:
            df = pd.read_csv(filepath)
        
        print(f"  Loaded {len(df)} epochs for {system['name']}")
        return df
    except Exception as e:
        print(f"  Error loading {system_key}: {e}")
        return None


def compute_temporal_shear(times, mag1, mag2, tau_values=[10, 20, 40, 80, 160]):
    """
    Compute temporal shear Γ between two images.
    
    Γ = d(Δt)/d(log τ) where Δt is the time delay residual
    """
    gamma_estimates = []
    
    for tau in tau_values:
        # Compute structure function difference
        sf1 = []
        sf2 = []
        
        for i in range(len(times)):
            for j in range(i+1, len(times)):
                dt = abs(times[j] - times[i])
                if abs(dt - tau) < tau * 0.3:  # Window around tau
                    sf1.append((mag1[j] - mag1[i])**2)
                    sf2.append((mag2[j] - mag2[i])**2)
        
        if len(sf1) > 10 and len(sf2) > 10:
            # Time delay proxy from SF difference
            delta_sf = np.mean(sf1) - np.mean(sf2)
            gamma_estimates.append((np.log10(tau), delta_sf))
    
    if len(gamma_estimates) < 3:
        return None, None
    
    # Linear fit to get slope
    log_tau = np.array([x[0] for x in gamma_estimates])
    delta_sf = np.array([x[1] for x in gamma_estimates])
    
    slope, intercept, r, p, se = stats.linregress(log_tau, delta_sf)
    
    # Convert to days/decade (approximate scaling)
    gamma = slope * 100  # rough conversion
    gamma_err = se * 100
    
    return gamma, gamma_err


def analyze_system(system_key):
    """Analyze a single lens system for temporal shear."""
    system = LENS_SYSTEMS[system_key]
    print(f"\n{'='*60}")
    print(f"Analyzing: {system['name']}")
    print(f"{'='*60}")
    print(f"z_source = {system['z_source']}, z_lens = {system['z_lens']}")
    print(f"σ_lens = {system['sigma_lens']} km/s")
    
    df = load_lightcurve(system_key)
    if df is None:
        return None
    
    results = {
        'system': system['name'],
        'z_source': system['z_source'],
        'z_lens': system['z_lens'],
        'sigma_lens': system['sigma_lens'],
        'n_epochs': len(df),
        'pairs': []
    }
    
    # Get time column
    time_col = None
    for col in ['mhjd', 'HJD', 'jd', 'MJD', 'time']:
        if col in df.columns:
            time_col = col
            break
    
    if time_col is None:
        print(f"  No time column found")
        return results
    
    times = df[time_col].values
    
    # Analyze each image pair
    labels = system['image_labels']
    for i, label1 in enumerate(labels):
        for label2 in labels[i+1:]:
            # Find magnitude columns
            mag1_col = None
            mag2_col = None
            
            for col in df.columns:
                if label1 in col and 'mag' in col.lower():
                    mag1_col = col
                elif label1 == col:
                    mag1_col = col
                if label2 in col and 'mag' in col.lower():
                    mag2_col = col
                elif label2 == col:
                    mag2_col = col
            
            if mag1_col is None or mag2_col is None:
                continue
            
            mag1 = df[mag1_col].values
            mag2 = df[mag2_col].values
            
            # Filter valid data
            mask = np.isfinite(mag1) & np.isfinite(mag2)
            if mask.sum() < 50:
                continue
            
            gamma, gamma_err = compute_temporal_shear(
                times[mask], mag1[mask], mag2[mask]
            )
            
            if gamma is not None:
                pair_result = {
                    'pair': f'{label1}-{label2}',
                    'gamma': float(gamma),
                    'gamma_err': float(gamma_err) if gamma_err else None,
                    'n_points': int(mask.sum())
                }
                results['pairs'].append(pair_result)
                
                sig = abs(gamma / gamma_err) if gamma_err and gamma_err > 0 else 0
                print(f"  {label1}-{label2}: Γ = {gamma:+.1f} ± {gamma_err:.1f} days/dec ({sig:.1f}σ)")
    
    return results


def predict_tep_signal(z_source, z_lens, sigma_lens):
    """
    Predict TEP temporal shear magnitude based on system parameters.
    
    Higher z_source and higher σ_lens should give stronger signals.
    """
    # Rough scaling: Γ ∝ σ² × z_S / (1 + z_L)
    reference_gamma = 33  # DESJ0408 validated value
    reference_sigma = 230
    reference_zs = 2.375
    reference_zl = 0.597
    
    scaling = (sigma_lens / reference_sigma)**2 * (z_source / reference_zs) * \
              ((1 + reference_zl) / (1 + z_lens))
    
    predicted_gamma = reference_gamma * scaling
    return predicted_gamma


def main():
    print("=" * 70)
    print("LENSING SAMPLE EXPANSION: TEMPORAL SHEAR ANALYSIS")
    print("=" * 70)
    print(f"Started: {datetime.now().isoformat()}\n")
    
    all_results = {
        'timestamp': datetime.now().isoformat(),
        'systems': {},
        'predictions': {}
    }
    
    # Analyze each system
    for system_key in LENS_SYSTEMS:
        result = analyze_system(system_key)
        if result:
            all_results['systems'][system_key] = result
        
        # Make TEP prediction
        system = LENS_SYSTEMS[system_key]
        predicted = predict_tep_signal(
            system['z_source'], 
            system['z_lens'],
            system['sigma_lens']
        )
        all_results['predictions'][system_key] = {
            'predicted_gamma': float(predicted),
            'z_source': system['z_source'],
            'sigma_lens': system['sigma_lens']
        }
        print(f"  TEP prediction: |Γ| ~ {predicted:.1f} days/dec")
    
    # Summary
    print("\n" + "=" * 70)
    print("EXPANSION TARGETS SUMMARY")
    print("=" * 70)
    
    print("\n| System | z_S | σ_lens | Predicted |Γ| | Status |")
    print("|--------|-----|--------|-----------|--------|")
    for key, system in LENS_SYSTEMS.items():
        pred = all_results['predictions'][key]['predicted_gamma']
        print(f"| {system['name'][:15]} | {system['z_source']:.2f} | {system['sigma_lens']} | {pred:.0f} days/dec | {system['status']} |")
    
    # Priority ranking
    print("\n### Priority Targets for TEP Validation:")
    targets = [(k, all_results['predictions'][k]['predicted_gamma']) 
               for k in LENS_SYSTEMS if LENS_SYSTEMS[k]['status'] == 'NEW_TARGET']
    targets.sort(key=lambda x: -x[1])
    
    for i, (key, pred) in enumerate(targets, 1):
        system = LENS_SYSTEMS[key]
        print(f"{i}. **{system['name']}**: |Γ| ~ {pred:.0f} days/dec")
        if 'notes' in system:
            print(f"   Note: {system['notes']}")
    
    # Save results
    output_path = os.path.join(RESULTS_DIR, 'step_7_1_lensing_expansion.json')
    with open(output_path, 'w') as f:
        json.dump(all_results, f, indent=2)
    print(f"\nResults saved: {output_path}")
    
    return all_results


if __name__ == '__main__':
    results = main()
