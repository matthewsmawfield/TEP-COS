#!/usr/bin/env python3
"""
Step 2.8: Expanded Stratification Analysis

Stratifies the Cosmic Coriolis signal by:
1. Morphology (Sersic n: n<2 = disk/spiral, n>2.5 = elliptical)
2. Star Formation Rate (SFR: star-forming vs quiescent)
3. Inclination (b/a corrected analysis)
4. Environment proxy (local density via angular separation)
5. Emission line strength (AGN vs normal)

Author: Matthew Lukin Smawfield
"""

import argparse
import json
import numpy as np
import pandas as pd
from pathlib import Path
from astropy.io import fits
from scipy import stats
import sys
sys.path.insert(0, str(Path(__file__).parent.parent))
from utils.logger import TEPLogger, set_step_logger, print_status

def robust_huber_fit(x, y, w=None):
    """Weighted Huber regression with fallback to OLS."""
    from sklearn.linear_model import HuberRegressor, LinearRegression
    if w is None:
        w = np.ones_like(x)
    mask = np.isfinite(x) & np.isfinite(y) & np.isfinite(w) & (w > 0)
    x, y, w = x[mask], y[mask], w[mask]
    if len(x) < 10:
        return np.nan, np.nan, np.nan
    X = x.reshape(-1, 1)
    try:
        model = HuberRegressor(epsilon=1.35, max_iter=200)
        model.fit(X, y, sample_weight=w)
        return model.coef_[0], model.intercept_, np.std(y - model.predict(X))
    except:
        # Fallback to weighted OLS
        model = LinearRegression()
        model.fit(X, y, sample_weight=w)
        return model.coef_[0], model.intercept_, np.std(y - model.predict(X))

def load_dapall_metadata(dapall_path: str) -> dict:
    """Load all relevant metadata from dapall."""
    with fits.open(dapall_path) as hdul:
        data = hdul[1].data
        metadata = {}
        for row in data:
            plateifu = row['PLATEIFU'].strip()
            metadata[plateifu] = {
                'z': float(row['NSA_Z']) if row['NSA_Z'] > 0 else float(row['Z']),
                'ba': float(row['NSA_ELPETRO_BA']),
                'sersic_n': float(row['NSA_SERSIC_N']),
                'sersic_ba': float(row['NSA_SERSIC_BA']),
                'sigma': float(row['STELLAR_SIGMA_1RE']),
                'sfr_1re': float(row['SFR_1RE']),
                'sfr_tot': float(row['SFR_TOT']),
                'ha_ew_1re': float(row['EMLINE_SEW_1RE'][18]) if len(row['EMLINE_SEW_1RE']) > 18 else np.nan,  # H-alpha index
                'ra': float(row['OBJRA']),
                'dec': float(row['OBJDEC']),
            }
    return metadata

def stratify_and_fit(df: pd.DataFrame, column: str, bins: list, labels: list) -> list:
    """Stratify by column and fit dipole in each bin."""
    results = []
    for i, (low, high) in enumerate(zip(bins[:-1], bins[1:])):
        mask = (df[column] >= low) & (df[column] < high)
        subset = df[mask]
        if len(subset) < 20:
            results.append({
                'label': labels[i],
                'range': f'{low:.2f}-{high:.2f}',
                'n': len(subset),
                'slope': np.nan,
                'slope_err': np.nan
            })
            continue
        
        x = subset['x_cmb'].values
        y = subset['delta_v'].values
        w = 1.0 / (subset['delta_v_err'].values ** 2 + 1e-6)
        
        slope, intercept, resid = robust_huber_fit(x, y, w)
        
        # Bootstrap for error
        n_boot = 500
        slopes = []
        for _ in range(n_boot):
            idx = np.random.choice(len(x), len(x), replace=True)
            s, _, _ = robust_huber_fit(x[idx], y[idx], w[idx])
            if np.isfinite(s):
                slopes.append(s)
        slope_err = np.std(slopes) if len(slopes) > 10 else np.nan
        
        results.append({
            'label': labels[i],
            'range': f'{low:.2f}-{high:.2f}',
            'n': len(subset),
            'slope': float(slope),
            'slope_err': float(slope_err)
        })
    return results

def main():
    parser = argparse.ArgumentParser(description='Expanded Stratification Analysis')
    parser.add_argument('--stellar-csv', required=True, help='Stellar per-galaxy CSV')
    parser.add_argument('--gas-csv', required=True, help='Gas per-galaxy CSV')
    parser.add_argument('--dapall', required=True, help='Path to dapall FITS')
    parser.add_argument('--output-dir', default='results/outputs', help='Output directory')
    args = parser.parse_args()
    
    logger = TEPLogger('step_2_8_expanded_stratification')
    set_step_logger(logger)
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    
    print_status('Loading metadata...', 'PROCESS')
    metadata = load_dapall_metadata(args.dapall)
    
    # Load stellar data
    df_star = pd.read_csv(args.stellar_csv)
    # Use delta_v column directly, create error from sigma if not present
    if 'delta_v_err' not in df_star.columns:
        df_star['delta_v_err'] = df_star['delta_v_sigma'] if 'delta_v_sigma' in df_star.columns else 10.0
    
    # Merge metadata
    for col in ['z', 'ba', 'sersic_n', 'sigma', 'sfr_1re', 'sfr_tot', 'ha_ew_1re', 'ra', 'dec']:
        df_star[col] = df_star['plateifu'].apply(lambda x: metadata.get(str(x).strip(), {}).get(col, np.nan))
    
    df_star = df_star.dropna(subset=['delta_v', 'x_cmb', 'sersic_n', 'ba', 'sigma'])
    print_status(f'Loaded {len(df_star)} galaxies with complete metadata.', 'INFO')
    
    results = {'stellar': {}, 'gas': {}}
    
    # ========== 1. MORPHOLOGY (Sersic n) ==========
    print_status('Stratifying by Morphology (Sersic n)...', 'PROCESS')
    morph_bins = [0, 1.5, 2.5, 4.0, 10.0]
    morph_labels = ['Pure Disk (n<1.5)', 'Disk-dominated (1.5-2.5)', 'Bulge-dominated (2.5-4)', 'Elliptical (n>4)']
    results['stellar']['morphology'] = stratify_and_fit(df_star, 'sersic_n', morph_bins, morph_labels)
    
    print('\n--- MORPHOLOGY STRATIFICATION (Stellar) ---')
    for r in results['stellar']['morphology']:
        print(f"{r['label']:30s} | N={r['n']:4d} | Slope={r['slope']:+7.2f} ± {r['slope_err']:.2f} km/s")
    
    # ========== 2. STAR FORMATION RATE ==========
    print_status('Stratifying by Star Formation Rate...', 'PROCESS')
    # Use log SFR, handle zeros
    df_star['log_sfr'] = np.log10(df_star['sfr_tot'].clip(lower=1e-4))
    sfr_bins = [-4, -2, -1, 0, 2]
    sfr_labels = ['Quiescent (log<-2)', 'Low SF (-2 to -1)', 'Moderate SF (-1 to 0)', 'High SF (log>0)']
    results['stellar']['sfr'] = stratify_and_fit(df_star, 'log_sfr', sfr_bins, sfr_labels)
    
    print('\n--- STAR FORMATION RATE STRATIFICATION (Stellar) ---')
    for r in results['stellar']['sfr']:
        print(f"{r['label']:30s} | N={r['n']:4d} | Slope={r['slope']:+7.2f} ± {r['slope_err']:.2f} km/s")
    
    # ========== 3. INCLINATION (b/a) ==========
    print_status('Stratifying by Inclination (b/a)...', 'PROCESS')
    ba_bins = [0, 0.3, 0.5, 0.7, 1.0]
    ba_labels = ['Edge-on (b/a<0.3)', 'Inclined (0.3-0.5)', 'Moderate (0.5-0.7)', 'Face-on (b/a>0.7)']
    results['stellar']['inclination'] = stratify_and_fit(df_star, 'ba', ba_bins, ba_labels)
    
    print('\n--- INCLINATION STRATIFICATION (Stellar) ---')
    for r in results['stellar']['inclination']:
        print(f"{r['label']:30s} | N={r['n']:4d} | Slope={r['slope']:+7.2f} ± {r['slope_err']:.2f} km/s")
    
    # ========== 4. REDSHIFT (Distance) ==========
    print_status('Stratifying by Redshift...', 'PROCESS')
    z_bins = [0, 0.02, 0.03, 0.05, 0.08, 0.15]
    z_labels = ['Very Local (z<0.02)', 'Local (0.02-0.03)', 'Intermediate (0.03-0.05)', 'Distant (0.05-0.08)', 'Far (z>0.08)']
    results['stellar']['redshift'] = stratify_and_fit(df_star, 'z', z_bins, z_labels)
    
    print('\n--- REDSHIFT STRATIFICATION (Stellar) ---')
    for r in results['stellar']['redshift']:
        print(f"{r['label']:30s} | N={r['n']:4d} | Slope={r['slope']:+7.2f} ± {r['slope_err']:.2f} km/s")
    
    # ========== 5. MASS (Velocity Dispersion) ==========
    print_status('Stratifying by Mass (Velocity Dispersion)...', 'PROCESS')
    sigma_bins = [0, 60, 100, 150, 250, 500]
    sigma_labels = ['Dwarf (σ<60)', 'Low Mass (60-100)', 'Intermediate (100-150)', 'Massive (150-250)', 'Very Massive (σ>250)']
    results['stellar']['mass'] = stratify_and_fit(df_star, 'sigma', sigma_bins, sigma_labels)
    
    print('\n--- MASS STRATIFICATION (Stellar) ---')
    for r in results['stellar']['mass']:
        print(f"{r['label']:30s} | N={r['n']:4d} | Slope={r['slope']:+7.2f} ± {r['slope_err']:.2f} km/s")
    
    # ========== 6. H-ALPHA EQUIVALENT WIDTH (AGN/SF proxy) ==========
    print_status('Stratifying by H-alpha EW (AGN/SF proxy)...', 'PROCESS')
    df_star['ha_ew_abs'] = np.abs(df_star['ha_ew_1re'])
    ha_bins = [0, 3, 10, 30, 1000]
    ha_labels = ['Weak/Quiescent (EW<3)', 'Moderate (3-10)', 'Strong SF (10-30)', 'Very Strong/AGN (EW>30)']
    results['stellar']['ha_ew'] = stratify_and_fit(df_star, 'ha_ew_abs', ha_bins, ha_labels)
    
    print('\n--- H-ALPHA EW STRATIFICATION (Stellar) ---')
    for r in results['stellar']['ha_ew']:
        print(f"{r['label']:30s} | N={r['n']:4d} | Slope={r['slope']:+7.2f} ± {r['slope_err']:.2f} km/s")
    
    # ========== 7. COMBINED "OPTIMAL" SAMPLE ==========
    print_status('Defining Optimal Sample...', 'PROCESS')
    optimal = df_star[
        (df_star['z'] < 0.04) &           # Local
        (df_star['sigma'] < 120) &         # Unscreened
        (df_star['ba'] > 0.5) &            # Not edge-on
        (df_star['sersic_n'] < 3.0)        # Disk-dominated
    ]
    
    if len(optimal) >= 20:
        x = optimal['x_cmb'].values
        y = optimal['delta_v'].values
        w = 1.0 / (optimal['delta_v_err'].values ** 2 + 1e-6)
        slope, _, _ = robust_huber_fit(x, y, w)
        
        # Bootstrap
        n_boot = 1000
        slopes = []
        for _ in range(n_boot):
            idx = np.random.choice(len(x), len(x), replace=True)
            s, _, _ = robust_huber_fit(x[idx], y[idx], w[idx])
            if np.isfinite(s):
                slopes.append(s)
        slope_err = np.std(slopes)
        ci_lo, ci_hi = np.percentile(slopes, [2.5, 97.5])
        
        results['stellar']['optimal_sample'] = {
            'criteria': 'z<0.04, sigma<120, b/a>0.5, n<3.0',
            'n': len(optimal),
            'slope': float(slope),
            'slope_err': float(slope_err),
            'ci_95': [float(ci_lo), float(ci_hi)]
        }
        
        print(f'\n--- OPTIMAL SAMPLE ---')
        print(f"Criteria: {results['stellar']['optimal_sample']['criteria']}")
        print(f"N = {len(optimal)}")
        print(f"Slope = {slope:+.2f} ± {slope_err:.2f} km/s")
        print(f"95% CI = [{ci_lo:.2f}, {ci_hi:.2f}]")
    
    # Save results
    out_path = out_dir / 'step_2_8_expanded_stratification.json'
    with open(out_path, 'w') as f:
        json.dump(results, f, indent=2)
    print_status(f'Saved {out_path}', 'SUCCESS')
    
    # Generate summary report
    report_path = out_dir / 'step_2_8_expanded_stratification_report.md'
    with open(report_path, 'w') as f:
        f.write('# TEP-COS Expanded Stratification Analysis\n\n')
        f.write('## Summary\n\n')
        f.write(f'Total galaxies analyzed: {len(df_star)}\n\n')
        
        for category, data in results['stellar'].items():
            if category == 'optimal_sample':
                f.write(f'## Optimal Sample\n')
                f.write(f"- Criteria: {data['criteria']}\n")
                f.write(f"- N: {data['n']}\n")
                f.write(f"- Slope: {data['slope']:+.2f} ± {data['slope_err']:.2f} km/s\n")
                f.write(f"- 95% CI: [{data['ci_95'][0]:.2f}, {data['ci_95'][1]:.2f}]\n\n")
            else:
                f.write(f'## {category.replace("_", " ").title()}\n\n')
                f.write('| Bin | N | Slope (km/s) |\n')
                f.write('|-----|---|-------------|\n')
                for r in data:
                    f.write(f"| {r['label']} | {r['n']} | {r['slope']:+.2f} ± {r['slope_err']:.2f} |\n")
                f.write('\n')
    
    print_status(f'Saved {report_path}', 'SUCCESS')

if __name__ == '__main__':
    main()
