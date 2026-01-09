#!/usr/bin/env python3
"""
Step 2.10: Metallicity Gradient & Peculiar Velocity Analysis

Additional discovery analyses:
D. Metallicity Gradient Asymmetry (O3N2 method)
G. Peculiar Velocity Correlation (from redshift-distance relation)
K. Lopsidedness Parameter (m=1 Fourier mode)
L. Kinematic Position Angle Offset from Photometric PA

Author: Matthew Lukin Smawfield
"""

import argparse
import json
import numpy as np
import pandas as pd
from pathlib import Path
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.cosmology import Planck18 as cosmo
import astropy.units as u
from scipy import stats
from sklearn.linear_model import HuberRegressor, LinearRegression
import warnings
warnings.filterwarnings('ignore')

CMB_RA = 168.0
CMB_DEC = -7.0

def robust_fit(x, y, w=None):
    """Weighted linear fit."""
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

def bootstrap_slope(x, y, w=None, n_boot=500):
    """Bootstrap confidence interval."""
    slopes = []
    n = len(x)
    for _ in range(n_boot):
        idx = np.random.choice(n, n, replace=True)
        s, _, _ = robust_fit(x[idx], y[idx], w[idx] if w is not None else None)
        if np.isfinite(s):
            slopes.append(s)
    if len(slopes) < 10:
        return np.nan, np.nan, [np.nan, np.nan]
    return np.mean(slopes), np.std(slopes), np.percentile(slopes, [2.5, 97.5])

def compute_x_cmb(ra, dec):
    """Compute projection onto CMB dipole axis."""
    cmb = SkyCoord(ra=CMB_RA*u.deg, dec=CMB_DEC*u.deg)
    gal = SkyCoord(ra=ra*u.deg, dec=dec*u.deg)
    return np.cos(gal.separation(cmb).rad)

def compute_metallicity_o3n2(oiii_flux, nii_flux, ha_flux, hb_flux):
    """
    Compute gas-phase metallicity using O3N2 method (Pettini & Pagel 2004).
    12 + log(O/H) = 8.73 - 0.32 * O3N2
    where O3N2 = log10((OIII/Hb) / (NII/Ha))
    Works with arrays.
    """
    oiii_flux = np.asarray(oiii_flux)
    nii_flux = np.asarray(nii_flux)
    ha_flux = np.asarray(ha_flux)
    hb_flux = np.asarray(hb_flux)
    
    result = np.full_like(ha_flux, np.nan, dtype=float)
    valid = (ha_flux > 0) & (hb_flux > 0) & (oiii_flux > 0) & (nii_flux > 0)
    
    if valid.any():
        o3n2 = np.log10((oiii_flux[valid] / hb_flux[valid]) / (nii_flux[valid] / ha_flux[valid]))
        result[valid] = 8.73 - 0.32 * o3n2
    
    return result

def compute_peculiar_velocity(z, ra, dec):
    """
    Estimate peculiar velocity from Hubble residual.
    V_pec = c * z - H0 * D_L (simplified)
    """
    if z <= 0 or z > 0.2:
        return np.nan
    
    # Luminosity distance from cosmology
    d_L = cosmo.luminosity_distance(z).to(u.Mpc).value
    
    # Expected recession velocity from Hubble flow
    v_hubble = cosmo.H0.value * d_L  # km/s
    
    # Observed recession velocity
    v_obs = 299792.458 * z  # km/s
    
    # Peculiar velocity (simplified - ignores relativistic corrections)
    v_pec = v_obs - v_hubble
    
    return v_pec

def analyze_galaxy_extended(maps_path, plateifu, metadata):
    """Extended analysis for metallicity and lopsidedness."""
    results = {'plateifu': plateifu}
    
    try:
        with fits.open(maps_path) as hdul:
            ext_names = [h.name for h in hdul]
            
            # Get emission line fluxes
            if 'EMLINE_GFLUX' not in ext_names:
                return None
            
            flux = hdul['EMLINE_GFLUX'].data
            # MaNGA indices: H-beta=11, OIII-5007=13, H-alpha=18, NII-6584=19
            
            if flux.ndim != 3 or flux.shape[0] < 20:
                return None
            
            hb = flux[11]
            oiii = flux[13]
            ha = flux[18]
            nii = flux[19]
            
            # Stellar velocity for lopsidedness
            if 'STELLAR_VEL' in ext_names:
                vel = hdul['STELLAR_VEL'].data
                if vel.ndim == 3:
                    vel = vel[0]
            else:
                vel = None
            
            # Elliptical coordinates
            if 'SPX_ELLCOO' in ext_names:
                ellcoo = hdul['SPX_ELLCOO'].data
                if ellcoo.ndim == 3:
                    r_ell = ellcoo[1]
                    theta = np.radians(ellcoo[3])
                else:
                    r_ell = ellcoo
                    theta = None
            else:
                r_ell = None
                theta = None
    except:
        return None
    
    # D. METALLICITY GRADIENT ASYMMETRY
    # Compute metallicity in each spaxel
    mask = (ha > 0) & (hb > 0) & (oiii > 0) & (nii > 0) & np.isfinite(ha)
    
    if mask.sum() > 50 and r_ell is not None:
        # Compute metallicity map
        z_metal = np.full_like(ha, np.nan)
        z_metal[mask] = compute_metallicity_o3n2(oiii[mask], nii[mask], ha[mask], hb[mask]).flatten()
        
        # Compute radial gradient
        r_flat = r_ell[mask & np.isfinite(z_metal)]
        z_flat = z_metal[mask & np.isfinite(z_metal)]
        
        if len(r_flat) > 30:
            # Fit gradient
            slope, intercept, _ = robust_fit(r_flat, z_flat)
            results['metallicity_gradient'] = float(slope)  # dex/Re
            results['metallicity_central'] = float(intercept)
            
            # Asymmetry: compare metallicity on CMB-aligned sides
            if theta is not None:
                ra_gal = metadata.get('ra', 0)
                dec_gal = metadata.get('dec', 0)
                cmb_pa = np.arctan2(
                    np.sin(np.radians(CMB_RA - ra_gal)),
                    np.cos(np.radians(CMB_DEC - dec_gal))
                )
                theta_cmb = theta - cmb_pa
                
                side_a = (np.cos(theta_cmb) > 0) & mask & np.isfinite(z_metal)
                side_b = (np.cos(theta_cmb) <= 0) & mask & np.isfinite(z_metal)
                
                if side_a.sum() > 20 and side_b.sum() > 20:
                    z_a = np.median(z_metal[side_a])
                    z_b = np.median(z_metal[side_b])
                    results['metallicity_asymmetry'] = float(z_a - z_b)  # dex
    
    # K. LOPSIDEDNESS (m=1 Fourier mode)
    if vel is not None and theta is not None:
        vel_mask = np.isfinite(vel) & np.isfinite(r_ell)
        sys_sel = vel_mask & (r_ell <= 0.1)
        ann_sel = vel_mask & (r_ell >= 0.8) & (r_ell <= 1.2) & np.isfinite(theta)
        if sys_sel.sum() > 5 and ann_sel.sum() > 50:
            v0 = float(np.nanmedian(vel[sys_sel]))
            v_ann = vel[ann_sel] - v0
            theta_ann = theta[ann_sel]

            speed = np.abs(v_ann)
            speed_centered = speed - float(np.mean(speed))

            cos_theta = np.cos(theta_ann)
            sin_theta = np.sin(theta_ann)

            A1 = 2 * float(np.mean(speed_centered * cos_theta))
            B1 = 2 * float(np.mean(speed_centered * sin_theta))

            lopsidedness = float(np.sqrt(A1**2 + B1**2))
            results['lopsidedness'] = lopsidedness

            lop_phase = float(np.arctan2(B1, A1))
            results['lopsidedness_phase'] = float(np.degrees(lop_phase))
    
    # G. PECULIAR VELOCITY
    z = metadata.get('z', np.nan)
    ra = metadata.get('ra', np.nan)
    dec = metadata.get('dec', np.nan)
    
    if np.isfinite(z) and z > 0:
        v_pec = compute_peculiar_velocity(z, ra, dec)
        results['v_peculiar'] = float(v_pec) if np.isfinite(v_pec) else np.nan
    
    # Add metadata
    results['x_cmb'] = metadata.get('x_cmb', np.nan)
    results['z'] = z
    results['ra'] = ra
    results['dec'] = dec
    results['sigma'] = metadata.get('sigma', np.nan)
    results['ba'] = metadata.get('ba', np.nan)
    
    return results

def main():
    parser = argparse.ArgumentParser(description='Metallicity & Peculiar Velocity Analysis')
    parser.add_argument('--maps-dir', default='data/maps/HYB10-MILESHC-MASTARSSP')
    parser.add_argument('--plateifu-list', default='results/outputs/step_1_0_plateifu_selection.txt')
    parser.add_argument('--dapall', default='data/dapall/dapall-v3_1_1-3.1.0.fits')
    parser.add_argument('--output-dir', default='results/outputs/discovery')
    parser.add_argument('--max-galaxies', type=int, default=2000)
    args = parser.parse_args()
    
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"[PROCESS] Loading metadata...")
    
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
                'sigma': float(row['STELLAR_SIGMA_1RE']),
            }
    
    print(f"[SUCCESS] Loaded {len(metadata)} metadata records.")
    
    # Load plateifu list
    with open(args.plateifu_list) as f:
        plateifus = [line.strip() for line in f if line.strip()]
    
    print(f"[PROCESS] Processing {min(len(plateifus), args.max_galaxies)} galaxies...")
    
    maps_dir = Path(args.maps_dir)
    results_list = []
    
    for i, plateifu in enumerate(plateifus[:args.max_galaxies]):
        if i % 100 == 0:
            print(f"[PROCESS] [{i}/{min(len(plateifus), args.max_galaxies)}]...")
        
        plate, ifu = plateifu.split('-')
        maps_path = maps_dir / plate / ifu / f"manga-{plateifu}-MAPS-HYB10-MILESHC-MASTARSSP.fits.gz"
        
        if not maps_path.exists():
            continue
        
        meta = metadata.get(plateifu, {})
        result = analyze_galaxy_extended(str(maps_path), plateifu, meta)
        
        if result:
            results_list.append(result)
    
    print(f"[SUCCESS] Analyzed {len(results_list)} galaxies.")
    
    df = pd.DataFrame(results_list)
    df.to_csv(out_dir / 'step_2_10_metallicity_peculiar.csv', index=False)
    
    # ========== AGGREGATE ANALYSES ==========
    summary = {}
    
    print("\n" + "="*60)
    print("D. METALLICITY GRADIENT ASYMMETRY VS CMB")
    print("="*60)
    
    if 'metallicity_asymmetry' in df.columns:
        valid = df.dropna(subset=['metallicity_asymmetry', 'x_cmb'])
        if len(valid) > 20:
            slope, _, _ = robust_fit(valid['x_cmb'].values, valid['metallicity_asymmetry'].values)
            mean_slope, std_slope, ci = bootstrap_slope(valid['x_cmb'].values, valid['metallicity_asymmetry'].values)
            print(f"  Metallicity asymmetry: slope = {slope:.4f} ± {std_slope:.4f} dex")
            print(f"  95% CI: [{ci[0]:.4f}, {ci[1]:.4f}]")
            summary['metallicity_asymmetry_slope'] = {
                'slope': float(slope), 'err': float(std_slope), 
                'ci_95': [float(ci[0]), float(ci[1])], 'n': len(valid)
            }
    
    if 'metallicity_gradient' in df.columns:
        valid = df.dropna(subset=['metallicity_gradient', 'x_cmb'])
        if len(valid) > 20:
            slope, _, _ = robust_fit(valid['x_cmb'].values, valid['metallicity_gradient'].values)
            mean_slope, std_slope, ci = bootstrap_slope(valid['x_cmb'].values, valid['metallicity_gradient'].values)
            print(f"  Metallicity gradient: slope = {slope:.4f} ± {std_slope:.4f} dex/Re")
            summary['metallicity_gradient_slope'] = {
                'slope': float(slope), 'err': float(std_slope), 'n': len(valid)
            }
    
    print("\n" + "="*60)
    print("G. PECULIAR VELOCITY CORRELATION WITH CMB")
    print("="*60)
    
    if 'v_peculiar' in df.columns:
        valid = df.dropna(subset=['v_peculiar', 'x_cmb'])
        # Filter extreme values
        valid = valid[np.abs(valid['v_peculiar']) < 2000]
        if len(valid) > 20:
            slope, _, _ = robust_fit(valid['x_cmb'].values, valid['v_peculiar'].values)
            mean_slope, std_slope, ci = bootstrap_slope(valid['x_cmb'].values, valid['v_peculiar'].values)
            print(f"  Peculiar velocity: slope = {slope:.1f} ± {std_slope:.1f} km/s")
            print(f"  95% CI: [{ci[0]:.1f}, {ci[1]:.1f}]")
            print(f"  Mean V_pec: {valid['v_peculiar'].mean():.1f} km/s")
            summary['v_peculiar_slope'] = {
                'slope': float(slope), 'err': float(std_slope),
                'ci_95': [float(ci[0]), float(ci[1])], 'n': len(valid),
                'mean_v_pec': float(valid['v_peculiar'].mean())
            }
    
    print("\n" + "="*60)
    print("K. LOPSIDEDNESS VS CMB")
    print("="*60)
    
    if 'lopsidedness' in df.columns:
        valid = df.dropna(subset=['lopsidedness', 'x_cmb'])
        if len(valid) > 20:
            slope, _, _ = robust_fit(valid['x_cmb'].values, valid['lopsidedness'].values)
            mean_slope, std_slope, ci = bootstrap_slope(valid['x_cmb'].values, valid['lopsidedness'].values)
            print(f"  Lopsidedness amplitude: slope = {slope:.2f} ± {std_slope:.2f} km/s")
            print(f"  Mean lopsidedness: {valid['lopsidedness'].mean():.2f} km/s")
            summary['lopsidedness_slope'] = {
                'slope': float(slope), 'err': float(std_slope), 'n': len(valid)
            }
    
    if 'lopsidedness_phase' in df.columns:
        valid = df.dropna(subset=['lopsidedness_phase', 'x_cmb'])
        if len(valid) > 20:
            # Test if lopsidedness phase correlates with CMB direction
            slope, _, _ = robust_fit(valid['x_cmb'].values, valid['lopsidedness_phase'].values)
            print(f"  Lopsidedness phase: slope = {slope:.1f}°")
            summary['lopsidedness_phase_slope'] = {'slope': float(slope), 'n': len(valid)}
    
    # Save summary
    with open(out_dir / 'step_2_10_metallicity_peculiar_summary.json', 'w') as f:
        json.dump(summary, f, indent=2)
    
    print(f"\n[SUCCESS] Saved results to {out_dir}")

if __name__ == '__main__':
    main()
