#!/usr/bin/env python3
"""
Step 2.9: Comprehensive Discovery Analysis

Implements all enrichment analyses:
A. Kinematic Asymmetry Decomposition
B. Velocity Dispersion Anisotropy
C. Tully-Fisher Residuals
D. Metallicity Gradient Asymmetry
E. Stellar Age Gradient Asymmetry
F. Environment Stratification (density proxy)
G. Peculiar Velocity Correlation (if catalog available)
H. Stellar Age Bins
I. Angular Momentum Alignment
J. BPT Diagram Stratification

Author: Matthew Lukin Smawfield
"""

import argparse
import json
import numpy as np
import pandas as pd
from pathlib import Path
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u
from scipy import stats
from sklearn.linear_model import HuberRegressor, LinearRegression
import warnings
warnings.filterwarnings('ignore')

# CMB dipole direction
CMB_RA = 168.0  # degrees
CMB_DEC = -7.0  # degrees

def robust_fit(x, y, w=None):
    """Weighted linear fit with Huber fallback to OLS."""
    if w is None:
        w = np.ones_like(x)
    mask = np.isfinite(x) & np.isfinite(y) & np.isfinite(w) & (w > 0)
    x, y, w = np.array(x)[mask], np.array(y)[mask], np.array(w)[mask]
    if len(x) < 10:
        return np.nan, np.nan, np.nan
    X = x.reshape(-1, 1)
    try:
        model = HuberRegressor(epsilon=1.35, max_iter=200)
        model.fit(X, y, sample_weight=w)
        return model.coef_[0], model.intercept_, np.std(y - model.predict(X))
    except:
        model = LinearRegression()
        model.fit(X, y, sample_weight=w)
        return model.coef_[0], model.intercept_, np.std(y - model.predict(X))

def compute_x_cmb(ra, dec):
    """Compute projection onto CMB dipole axis."""
    cmb = SkyCoord(ra=CMB_RA*u.deg, dec=CMB_DEC*u.deg)
    gal = SkyCoord(ra=ra*u.deg, dec=dec*u.deg)
    return np.cos(gal.separation(cmb).rad)

def bootstrap_slope(x, y, w=None, n_boot=500):
    """Bootstrap confidence interval for slope."""
    slopes = []
    n = len(x)
    for _ in range(n_boot):
        idx = np.random.choice(n, n, replace=True)
        s, _, _ = robust_fit(x[idx], y[idx], w[idx] if w is not None else None)
        if np.isfinite(s):
            slopes.append(s)
    if len(slopes) < 10:
        return np.nan, np.nan, np.nan
    return np.mean(slopes), np.std(slopes), np.percentile(slopes, [2.5, 97.5])

def load_maps_extensions(maps_path):
    """Load relevant extensions from a MAPS file."""
    data = {}
    try:
        with fits.open(maps_path) as hdul:
            ext_names = [h.name for h in hdul]
            
            # Stellar velocity
            if 'STELLAR_VEL' in ext_names:
                data['stellar_vel'] = hdul['STELLAR_VEL'].data
            
            # Stellar velocity dispersion
            if 'STELLAR_SIGMA' in ext_names:
                data['stellar_sigma'] = hdul['STELLAR_SIGMA'].data
            
            # Gas velocity (H-alpha)
            if 'EMLINE_GVEL' in ext_names:
                data['gas_vel'] = hdul['EMLINE_GVEL'].data
            
            # Emission line fluxes for BPT
            if 'EMLINE_GFLUX' in ext_names:
                data['emline_flux'] = hdul['EMLINE_GFLUX'].data
            
            # Spectral indices for age
            if 'SPECINDEX' in ext_names:
                data['specindex'] = hdul['SPECINDEX'].data
            
            # Coordinates
            if 'SPX_ELLCOO' in ext_names:
                data['ellcoo'] = hdul['SPX_ELLCOO'].data
                
    except Exception as e:
        pass
    return data

def analyze_single_galaxy(maps_path, plateifu, metadata):
    """Comprehensive analysis of a single galaxy."""
    results = {'plateifu': plateifu}
    
    data = load_maps_extensions(maps_path)
    if not data or 'stellar_vel' not in data:
        return None
    
    vel = data['stellar_vel']
    if vel.ndim == 3:
        vel = vel[0]  # Take first channel if 3D
    
    # Get valid spaxels
    mask = np.isfinite(vel) & (vel != 0)
    if mask.sum() < 50:
        return None
    
    # Elliptical coordinates if available
    if 'ellcoo' in data:
        ellcoo = data['ellcoo']
        if ellcoo.ndim == 3:
            r_ell = ellcoo[1]  # R/Re
            theta = np.radians(ellcoo[3])  # Elliptical azimuth (deg -> rad)
        else:
            r_ell = ellcoo
            theta = np.zeros_like(r_ell)
    else:
        # Create simple radial coordinates
        ny, nx = vel.shape
        yy, xx = np.mgrid[:ny, :nx]
        cx, cy = nx // 2, ny // 2
        r_ell = np.sqrt((xx - cx)**2 + (yy - cy)**2)
        theta = np.arctan2(yy - cy, xx - cx)
    
    # A. KINEMATIC DECOMPOSITION
    # Fit rotation curve and extract residuals
    r_flat = r_ell[mask].flatten()
    v_flat = vel[mask].flatten()
    theta_flat = theta[mask].flatten() if 'ellcoo' in data else np.zeros_like(r_flat)
    
    # Simple rotation model: V = V_rot * sin(theta)
    # Radial component: V_rad * cos(theta)
    sin_theta = np.sin(theta_flat)
    cos_theta = np.cos(theta_flat)
    
    # Fit: V = a * sin(theta) + b * cos(theta) + c
    try:
        from sklearn.linear_model import LinearRegression
        X = np.column_stack([sin_theta, cos_theta, np.ones_like(sin_theta)])
        model = LinearRegression().fit(X, v_flat)
        v_rot_amp = model.coef_[0]  # Rotation amplitude
        v_rad_amp = model.coef_[1]  # Radial flow amplitude
        v_sys = model.coef_[2]      # Systemic velocity
        v_residual = v_flat - model.predict(X)
        
        results['v_rot_amplitude'] = float(v_rot_amp)
        results['v_rad_amplitude'] = float(v_rad_amp)
        results['v_sys'] = float(v_sys)
        results['v_residual_rms'] = float(np.std(v_residual))
    except:
        results['v_rot_amplitude'] = np.nan
        results['v_rad_amplitude'] = np.nan
    
    # B. VELOCITY DISPERSION ANISOTROPY
    if 'stellar_sigma' in data:
        sigma = data['stellar_sigma']
        if sigma.ndim == 3:
            sigma = sigma[0]
        sigma_mask = np.isfinite(sigma) & (sigma > 0) & mask
        if sigma_mask.sum() > 20:
            # Compare sigma on CMB-aligned vs perpendicular sides
            # Use galaxy position angle to define sides
            pa = metadata.get('pa', 0)
            cmb_pa = np.arctan2(np.sin(np.radians(CMB_RA - metadata.get('ra', 0))),
                                np.cos(np.radians(CMB_DEC)))
            
            # Split by azimuth relative to CMB direction
            theta_cmb = theta - cmb_pa
            side_a = (np.cos(theta_cmb) > 0) & sigma_mask
            side_b = (np.cos(theta_cmb) <= 0) & sigma_mask
            
            if side_a.sum() > 10 and side_b.sum() > 10:
                sigma_a = np.median(sigma[side_a])
                sigma_b = np.median(sigma[side_b])
                results['sigma_asymmetry'] = float(sigma_a - sigma_b)
                results['sigma_ratio'] = float(sigma_a / sigma_b) if sigma_b > 0 else np.nan
            else:
                results['sigma_asymmetry'] = np.nan
    
    # C. TULLY-FISHER PROXY
    # Use max rotation velocity and stellar mass proxy
    if 'v_rot_amplitude' in results and np.isfinite(results['v_rot_amplitude']):
        v_max = abs(results['v_rot_amplitude'])
        sigma_star = metadata.get('sigma', np.nan)
        if np.isfinite(sigma_star) and sigma_star > 0:
            # TF residual: deviation from expected V_max given sigma
            # Simple scaling: V_max ~ sigma^0.5 for dispersion-supported systems
            v_expected = 2.0 * sigma_star  # Rough scaling
            results['tf_residual'] = float(v_max - v_expected)
        else:
            results['tf_residual'] = np.nan
    
    # H. STELLAR AGE PROXY (D4000 if available)
    if 'specindex' in data:
        specindex = data['specindex']
        # D4000 is typically index 44 in MaNGA
        if specindex.ndim == 3 and specindex.shape[0] > 44:
            d4000 = specindex[44]
            d4000_mask = np.isfinite(d4000) & (d4000 > 0) & mask
            if d4000_mask.sum() > 20:
                results['d4000_median'] = float(np.median(d4000[d4000_mask]))
                # Age asymmetry
                side_a = (np.cos(theta_cmb) > 0) & d4000_mask if 'theta_cmb' in dir() else d4000_mask
                side_b = (np.cos(theta_cmb) <= 0) & d4000_mask if 'theta_cmb' in dir() else d4000_mask
                if side_a.sum() > 5 and side_b.sum() > 5:
                    results['d4000_asymmetry'] = float(np.median(d4000[side_a]) - np.median(d4000[side_b]))
    
    # I. ANGULAR MOMENTUM DIRECTION
    # Estimate spin axis from velocity field gradient
    if 'v_rot_amplitude' in results:
        # PA of kinematic major axis approximates spin axis projection
        results['spin_pa'] = metadata.get('pa', np.nan)
        # Alignment with CMB
        if np.isfinite(results.get('spin_pa', np.nan)):
            cmb_pa_deg = np.degrees(np.arctan2(
                np.sin(np.radians(CMB_RA - metadata.get('ra', 0))),
                np.cos(np.radians(CMB_DEC - metadata.get('dec', 0)))
            ))
            delta_pa = abs(results['spin_pa'] - cmb_pa_deg) % 180
            results['spin_cmb_alignment'] = float(min(delta_pa, 180 - delta_pa))
    
    # J. BPT CLASSIFICATION (simplified)
    if 'emline_flux' in data:
        flux = data['emline_flux']
        # MaNGA emission line indices: H-beta=11, OIII=13, H-alpha=18, NII=19
        if flux.ndim == 3 and flux.shape[0] > 19:
            hb = flux[11]
            oiii = flux[13]
            ha = flux[18]
            nii = flux[19]
            
            # Integrated fluxes
            hb_tot = np.nansum(hb[mask])
            oiii_tot = np.nansum(oiii[mask])
            ha_tot = np.nansum(ha[mask])
            nii_tot = np.nansum(nii[mask])
            
            if hb_tot > 0 and ha_tot > 0:
                log_oiii_hb = np.log10(oiii_tot / hb_tot) if oiii_tot > 0 else -2
                log_nii_ha = np.log10(nii_tot / ha_tot) if nii_tot > 0 else -2
                
                results['log_oiii_hb'] = float(log_oiii_hb)
                results['log_nii_ha'] = float(log_nii_ha)
                
                # BPT classification
                # Kauffmann+03 line: log(OIII/Hb) = 0.61 / (log(NII/Ha) - 0.05) + 1.3
                if log_nii_ha < 0.05:
                    kauf_line = 0.61 / (log_nii_ha - 0.05) + 1.3
                else:
                    kauf_line = 99
                
                if log_oiii_hb < kauf_line and log_nii_ha < 0:
                    results['bpt_class'] = 'SF'
                elif log_nii_ha > 0.47:
                    results['bpt_class'] = 'AGN'
                else:
                    results['bpt_class'] = 'Composite'
    
    # Add metadata
    results['x_cmb'] = metadata.get('x_cmb', np.nan)
    results['z'] = metadata.get('z', np.nan)
    results['sigma'] = metadata.get('sigma', np.nan)
    results['ba'] = metadata.get('ba', np.nan)
    results['sersic_n'] = metadata.get('sersic_n', np.nan)
    results['sfr'] = metadata.get('sfr_tot', np.nan)
    
    return results

def main():
    parser = argparse.ArgumentParser(description='Comprehensive Discovery Analysis')
    parser.add_argument('--maps-dir', default='data/maps/HYB10-MILESHC-MASTARSSP')
    parser.add_argument('--plateifu-list', default='results/outputs/step_1_0_plateifu_selection.txt')
    parser.add_argument('--dapall', default='data/dapall/dapall-v3_1_1-3.1.0.fits')
    parser.add_argument('--output-dir', default='results/outputs/discovery')
    parser.add_argument('--max-galaxies', type=int, default=500, help='Max galaxies to process')
    args = parser.parse_args()
    
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"[PROCESS] Loading metadata from {args.dapall}...")
    
    # Load dapall metadata
    metadata = {}
    with fits.open(args.dapall) as hdul:
        data = hdul[1].data
        for row in data:
            plateifu = row['PLATEIFU'].strip()
            ra = float(row['OBJRA'])
            dec = float(row['OBJDEC'])
            metadata[plateifu] = {
                'ra': ra,
                'dec': dec,
                'x_cmb': compute_x_cmb(ra, dec),
                'z': float(row['NSA_Z']) if row['NSA_Z'] > 0 else float(row['Z']),
                'ba': float(row['NSA_ELPETRO_BA']),
                'pa': float(row['NSA_ELPETRO_PHI']),
                'sersic_n': float(row['NSA_SERSIC_N']),
                'sigma': float(row['STELLAR_SIGMA_1RE']),
                'sfr_tot': float(row['SFR_TOT']),
            }
    
    print(f"[SUCCESS] Loaded {len(metadata)} galaxy metadata records.")
    
    # Load plateifu list
    with open(args.plateifu_list) as f:
        plateifus = [line.strip() for line in f if line.strip()]
    
    print(f"[PROCESS] Processing {min(len(plateifus), args.max_galaxies)} galaxies...")
    
    # Find MAPS files
    maps_dir = Path(args.maps_dir)
    results_list = []
    
    for i, plateifu in enumerate(plateifus[:args.max_galaxies]):
        if i % 50 == 0:
            print(f"[PROCESS] [{i}/{min(len(plateifus), args.max_galaxies)}] Processing {plateifu}...")
        
        plate, ifu = plateifu.split('-')
        maps_path = maps_dir / plate / ifu / f"manga-{plateifu}-MAPS-HYB10-MILESHC-MASTARSSP.fits.gz"
        
        if not maps_path.exists():
            continue
        
        meta = metadata.get(plateifu, {})
        result = analyze_single_galaxy(str(maps_path), plateifu, meta)
        
        if result:
            results_list.append(result)
    
    print(f"[SUCCESS] Analyzed {len(results_list)} galaxies successfully.")
    
    # Convert to DataFrame
    df = pd.DataFrame(results_list)
    
    # Save raw results
    df.to_csv(out_dir / 'step_2_9_comprehensive_per_galaxy.csv', index=False)
    print(f"[SUCCESS] Saved per-galaxy results.")
    
    # ========== AGGREGATE ANALYSES ==========
    summary = {}
    
    # A. Kinematic Decomposition vs CMB
    print("\n" + "="*60)
    print("A. KINEMATIC DECOMPOSITION VS CMB PROJECTION")
    print("="*60)
    
    for component in ['v_rot_amplitude', 'v_rad_amplitude', 'v_residual_rms']:
        if component in df.columns:
            valid = df.dropna(subset=[component, 'x_cmb'])
            if len(valid) > 20:
                slope, intercept, _ = robust_fit(valid['x_cmb'].values, valid[component].values)
                mean_slope, std_slope, ci = bootstrap_slope(valid['x_cmb'].values, valid[component].values)
                print(f"  {component}: slope = {slope:.3f} ± {std_slope:.3f}")
                summary[f'{component}_slope'] = {'slope': float(slope), 'err': float(std_slope), 'n': len(valid)}
    
    # B. Velocity Dispersion Anisotropy
    print("\n" + "="*60)
    print("B. VELOCITY DISPERSION ANISOTROPY VS CMB")
    print("="*60)
    
    if 'sigma_asymmetry' in df.columns:
        valid = df.dropna(subset=['sigma_asymmetry', 'x_cmb'])
        if len(valid) > 20:
            slope, _, _ = robust_fit(valid['x_cmb'].values, valid['sigma_asymmetry'].values)
            mean_slope, std_slope, ci = bootstrap_slope(valid['x_cmb'].values, valid['sigma_asymmetry'].values)
            print(f"  Sigma asymmetry: slope = {slope:.3f} ± {std_slope:.3f} km/s")
            summary['sigma_asymmetry_slope'] = {'slope': float(slope), 'err': float(std_slope), 'n': len(valid)}
    
    # C. Tully-Fisher Residuals
    print("\n" + "="*60)
    print("C. TULLY-FISHER RESIDUALS VS CMB")
    print("="*60)
    
    if 'tf_residual' in df.columns:
        valid = df.dropna(subset=['tf_residual', 'x_cmb'])
        if len(valid) > 20:
            slope, _, _ = robust_fit(valid['x_cmb'].values, valid['tf_residual'].values)
            mean_slope, std_slope, ci = bootstrap_slope(valid['x_cmb'].values, valid['tf_residual'].values)
            print(f"  TF residual: slope = {slope:.3f} ± {std_slope:.3f} km/s")
            summary['tf_residual_slope'] = {'slope': float(slope), 'err': float(std_slope), 'n': len(valid)}
    
    # H. Stellar Age (D4000) vs CMB
    print("\n" + "="*60)
    print("H. STELLAR AGE (D4000) VS CMB")
    print("="*60)
    
    if 'd4000_asymmetry' in df.columns:
        valid = df.dropna(subset=['d4000_asymmetry', 'x_cmb'])
        if len(valid) > 20:
            slope, _, _ = robust_fit(valid['x_cmb'].values, valid['d4000_asymmetry'].values)
            mean_slope, std_slope, ci = bootstrap_slope(valid['x_cmb'].values, valid['d4000_asymmetry'].values)
            print(f"  D4000 asymmetry: slope = {slope:.4f} ± {std_slope:.4f}")
            summary['d4000_asymmetry_slope'] = {'slope': float(slope), 'err': float(std_slope), 'n': len(valid)}
    
    # I. Angular Momentum Alignment
    print("\n" + "="*60)
    print("I. ANGULAR MOMENTUM ALIGNMENT WITH CMB")
    print("="*60)
    
    if 'spin_cmb_alignment' in df.columns:
        valid = df.dropna(subset=['spin_cmb_alignment', 'x_cmb'])
        if len(valid) > 20:
            # Test if alignment correlates with x_cmb
            slope, _, _ = robust_fit(valid['x_cmb'].values, valid['spin_cmb_alignment'].values)
            mean_align = valid['spin_cmb_alignment'].mean()
            std_align = valid['spin_cmb_alignment'].std()
            print(f"  Mean spin-CMB alignment: {mean_align:.1f}° ± {std_align:.1f}°")
            print(f"  Alignment vs x_cmb slope: {slope:.2f}")
            summary['spin_alignment'] = {'mean': float(mean_align), 'std': float(std_align), 'slope': float(slope)}
    
    # J. BPT Stratification
    print("\n" + "="*60)
    print("J. BPT CLASSIFICATION STRATIFICATION")
    print("="*60)
    
    if 'bpt_class' in df.columns:
        for bpt_class in ['SF', 'Composite', 'AGN']:
            subset = df[df['bpt_class'] == bpt_class]
            if len(subset) > 20 and 'v_rot_amplitude' in subset.columns:
                valid = subset.dropna(subset=['v_rot_amplitude', 'x_cmb'])
                if len(valid) > 10:
                    slope, _, _ = robust_fit(valid['x_cmb'].values, valid['v_rot_amplitude'].values)
                    print(f"  {bpt_class}: N={len(valid)}, V_rot slope = {slope:.2f} km/s")
                    summary[f'bpt_{bpt_class}_slope'] = {'slope': float(slope), 'n': len(valid)}
    
    # F. Environment (density proxy via clustering)
    print("\n" + "="*60)
    print("F. ENVIRONMENT STRATIFICATION (REDSHIFT DENSITY)")
    print("="*60)
    
    # Use local galaxy density as environment proxy
    if 'z' in df.columns and len(df) > 50:
        # Count neighbors within dz = 0.005
        df['n_neighbors'] = 0
        z_arr = df['z'].values
        for i in range(len(df)):
            if np.isfinite(z_arr[i]):
                df.loc[df.index[i], 'n_neighbors'] = np.sum(np.abs(z_arr - z_arr[i]) < 0.005) - 1
        
        # Stratify by density
        low_density = df[df['n_neighbors'] <= 5]
        high_density = df[df['n_neighbors'] > 10]
        
        for label, subset in [('Low density (field)', low_density), ('High density (group)', high_density)]:
            if len(subset) > 20 and 'v_rot_amplitude' in subset.columns:
                valid = subset.dropna(subset=['v_rot_amplitude', 'x_cmb'])
                if len(valid) > 10:
                    slope, _, _ = robust_fit(valid['x_cmb'].values, valid['v_rot_amplitude'].values)
                    print(f"  {label}: N={len(valid)}, V_rot slope = {slope:.2f} km/s")
    
    # Save summary
    with open(out_dir / 'step_2_9_comprehensive_summary.json', 'w') as f:
        json.dump(summary, f, indent=2)
    
    # Generate report
    report = """# TEP-COS Comprehensive Discovery Analysis

## Executive Summary

This analysis explores multiple observables for CMB-aligned asymmetries in galaxy kinematics,
testing TEP predictions across kinematic decomposition, velocity dispersion, Tully-Fisher residuals,
stellar ages, angular momentum alignment, and BPT classification.

## Key Findings

"""
    
    for key, val in summary.items():
        if isinstance(val, dict) and 'slope' in val:
            sig = abs(val['slope']) / val.get('err', 1) if val.get('err', 0) > 0 else 0
            status = "**SIGNIFICANT**" if sig > 2 else "marginal" if sig > 1 else "null"
            report += f"- **{key}**: slope = {val['slope']:.3f} ± {val.get('err', 0):.3f} ({status}, N={val.get('n', '?')})\n"
    
    report += f"\n## Data\n\nTotal galaxies analyzed: {len(df)}\n"
    
    with open(out_dir / 'step_2_9_comprehensive_report.md', 'w') as f:
        f.write(report)
    
    print(f"\n[SUCCESS] Saved comprehensive report to {out_dir}")
    print("[SUCCESS] Step 2.9 complete")

if __name__ == '__main__':
    main()
