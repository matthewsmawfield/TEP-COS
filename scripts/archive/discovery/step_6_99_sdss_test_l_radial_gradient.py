#!/usr/bin/env python3
"""
Test L: Radial Age Gradient within Galaxies (MaNGA)

Hypothesis:
Under TEP, galaxy centers (deeper potential) experience slower time flow.
This predicts centers should appear YOUNGER at fixed local [Z/H].
Standard inside-out formation predicts centers are OLDER.

TEP Prediction:
  dAge/dR > 0 (center appears younger)
  
Standard Prediction:
  dAge/dR < 0 (center formed first)
  
The age GRADIENT SIGN is the key discriminator.
"""

import os
import sys
import json
import requests
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt

# Paths
BASE_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATA_DIR = os.path.join(BASE_DIR, 'data', 'sdss')
OUTPUT_DIR = os.path.join(BASE_DIR, 'results', 'outputs')
FIGURE_DIR = os.path.join(BASE_DIR, 'results', 'figures')

os.makedirs(DATA_DIR, exist_ok=True)
os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(FIGURE_DIR, exist_ok=True)

SDSS_URL = "https://skyserver.sdss.org/dr18/SkyServerWS/SearchTools/SqlSearch"

def query_sdss(sql, max_rows=500000):
    """Query SDSS SkyServer."""
    params = {'cmd': sql, 'format': 'csv'}
    try:
        response = requests.get(SDSS_URL, params=params, timeout=300)
        response.raise_for_status()
        from io import StringIO
        df = pd.read_csv(StringIO(response.text))
        if len(df) == 0 or 'error' in df.columns[0].lower():
            return None
        return df
    except Exception as e:
        print(f"Query failed: {e}")
        return None


def download_gradient_data():
    """Download MaNGA radial gradient data."""
    cache_file = os.path.join(DATA_DIR, 'manga_radial_gradients.csv')
    
    if os.path.exists(cache_file):
        print(f"Loading cached gradient data from {cache_file}")
        return pd.read_csv(cache_file)
    
    # Query MaNGA Firefly for age gradients
    sql = """
    SELECT TOP 10000
        f.PLATEIFU,
        f.REDSHIFT AS redshift,
        
        -- Central vs outer ages
        f.LW_AGE_05RE AS age_center,
        f.LW_AGE_1RE AS age_1re,
        f.LW_AGE_15RE AS age_outer,
        
        -- Age gradient
        f.LW_AGE_GRADIENT AS age_gradient,
        f.LW_AGE_GRADIENT_ERROR AS age_gradient_err,
        
        -- Metallicity gradient (control)
        f.LW_Z_05RE AS z_center,
        f.LW_Z_1RE AS z_1re,
        f.LW_Z_15RE AS z_outer,
        f.LW_Z_GRADIENT AS z_gradient,
        f.LW_Z_GRADIENT_ERROR AS z_gradient_err,
        
        -- Mass-weighted ages
        f.MW_AGE_05RE AS mw_age_center,
        f.MW_AGE_1RE AS mw_age_1re,
        f.MW_AGE_15RE AS mw_age_outer,
        f.MW_AGE_GRADIENT AS mw_age_gradient,
        
        -- Alpha enhancement if available
        f.ALPHA_FE_1RE AS alpha_fe,
        
        -- Stellar mass
        f.PHOTOMETRIC_MASS AS log_mass
        
    FROM mangaFirefly f
    
    WHERE 
        f.REDSHIFT BETWEEN 0.01 AND 0.15
        AND f.LW_AGE_GRADIENT_ERROR > 0 AND f.LW_AGE_GRADIENT_ERROR < 3.0
        AND f.LW_AGE_05RE > 0 AND f.LW_AGE_15RE > 0
        AND f.PHOTOMETRIC_MASS > 9.0
    """
    
    print("Querying MaNGA Firefly for gradient data...")
    df = query_sdss(sql)
    
    if df is not None and len(df) > 50:
        df.to_csv(cache_file, index=False)
        print(f"Saved {len(df)} galaxies to {cache_file}")
        return df
    
    print("Firefly query failed. Trying alternative...")
    
    # Alternative: Use mangaDAPall for velocity dispersion gradients
    sql_alt = """
    SELECT TOP 10000
        d.plateifu,
        d.z AS redshift,
        d.stellar_sigma_1re AS sigma_center,
        d.stellar_sigma_elo AS sigma_outer,
        d.nsa_elpetro_mass AS log_mass,
        d.nsa_sersic_n AS sersic_n,
        d.nsa_elpetro_ba AS ba_ratio
        
    FROM mangaDAPall d
    
    WHERE 
        d.z BETWEEN 0.01 AND 0.15
        AND d.stellar_sigma_1re > 50 AND d.stellar_sigma_1re < 400
        AND d.drp3qual = 0
        AND d.nsa_elpetro_mass > 9.0
    """
    
    print("Querying mangaDAPall...")
    df_alt = query_sdss(sql_alt)
    
    if df_alt is not None:
        return df_alt
    
    return None


def load_existing_manga_data():
    """Load existing MaNGA age data as fallback."""
    manga_file = os.path.join(DATA_DIR, 'manga_age_data.csv')
    
    if os.path.exists(manga_file):
        print(f"Loading existing MaNGA data from {manga_file}")
        df = pd.read_csv(manga_file)
        # Remove duplicates
        df = df.drop_duplicates(subset=['PLATEIFU'])
        return df
    return None


def analyze_radial_gradients(df):
    """Analyze radial age gradients."""
    results = {}
    
    # Check what columns we have
    has_gradients = 'age_gradient' in df.columns
    has_lw_mw = 'LW_AGE_1RE' in df.columns and 'MW_AGE_1RE' in df.columns
    
    print(f"\nData columns: {list(df.columns)}")
    print(f"Has explicit gradients: {has_gradients}")
    print(f"Has LW/MW ages: {has_lw_mw}")
    
    if has_gradients:
        # Direct gradient analysis
        clean = df[df['age_gradient'].notna() & df['age_gradient_err'].notna()].copy()
        print(f"\nAnalyzing {len(clean)} galaxies with gradient data")
        results['n_galaxies'] = len(clean)
        
        # 1. Mean gradient sign
        mean_gradient = clean['age_gradient'].mean()
        std_gradient = clean['age_gradient'].std()
        sem_gradient = std_gradient / np.sqrt(len(clean))
        
        # Test if gradient is significantly different from zero
        t_stat = mean_gradient / sem_gradient
        p_value = 2 * stats.t.sf(abs(t_stat), len(clean) - 1)
        
        print(f"\n1. Age Gradient Distribution:")
        print(f"   Mean gradient: {mean_gradient:.4f} ± {sem_gradient:.4f} Gyr/Re")
        print(f"   t = {t_stat:.2f}, p = {p_value:.2e}")
        print(f"   Standard expects: < 0 (center older)")
        print(f"   TEP expects: > 0 (center younger)")
        
        results['mean_gradient'] = {
            'value': float(mean_gradient),
            'stderr': float(sem_gradient),
            't_stat': float(t_stat),
            'p_value': float(p_value),
            'sign': 'positive' if mean_gradient > 0 else 'negative'
        }
        
        # 2. Fraction with positive vs negative gradients
        n_positive = (clean['age_gradient'] > 0).sum()
        n_negative = (clean['age_gradient'] < 0).sum()
        frac_positive = n_positive / len(clean)
        
        # Binomial test against 50%
        binom_p = stats.binom_test(n_positive, len(clean), 0.5, alternative='two-sided')
        
        print(f"\n2. Gradient Sign Distribution:")
        print(f"   Positive (center younger): {n_positive} ({frac_positive*100:.1f}%)")
        print(f"   Negative (center older): {n_negative} ({(1-frac_positive)*100:.1f}%)")
        print(f"   Binomial test p = {binom_p:.2e}")
        
        results['sign_distribution'] = {
            'n_positive': int(n_positive),
            'n_negative': int(n_negative),
            'frac_positive': float(frac_positive),
            'binom_p': float(binom_p)
        }
        
        # 3. Gradient vs central potential (mass proxy)
        if 'log_mass' in clean.columns:
            mask = clean['log_mass'].notna()
            r_mass, p_mass = stats.pearsonr(clean.loc[mask, 'log_mass'], 
                                            clean.loc[mask, 'age_gradient'])
            print(f"\n3. Gradient vs Stellar Mass:")
            print(f"   r = {r_mass:.4f}, p = {p_mass:.2e}")
            
            results['gradient_vs_mass'] = {
                'r': float(r_mass),
                'p': float(p_mass)
            }
        
        # 4. Control for metallicity gradient
        if 'z_gradient' in clean.columns:
            mask = clean['z_gradient'].notna()
            if mask.sum() > 30:
                # Partial correlation
                from sklearn.linear_model import LinearRegression
                
                X = clean.loc[mask, ['z_gradient']].values
                y_age = clean.loc[mask, 'age_gradient'].values
                
                reg = LinearRegression().fit(X, y_age)
                age_resid = y_age - reg.predict(X)
                
                # Mean of residuals
                mean_resid = age_resid.mean()
                sem_resid = age_resid.std() / np.sqrt(len(age_resid))
                t_resid = mean_resid / sem_resid
                p_resid = 2 * stats.t.sf(abs(t_resid), len(age_resid) - 1)
                
                print(f"\n4. Age Gradient (controlling for Z gradient):")
                print(f"   Mean residual: {mean_resid:.4f} ± {sem_resid:.4f}")
                print(f"   t = {t_resid:.2f}, p = {p_resid:.2e}")
                
                results['controlled_for_Z'] = {
                    'mean_residual': float(mean_resid),
                    'stderr': float(sem_resid),
                    't_stat': float(t_resid),
                    'p_value': float(p_resid)
                }
        
        # Verdict
        if mean_gradient > 0 and p_value < 0.05:
            results['verdict'] = 'Signal'
            results['interpretation'] = ('Age gradients are POSITIVE (centers younger). '
                                        'TEP-consistent: deeper potential = younger appearance.')
        elif mean_gradient < 0 and p_value < 0.05:
            results['verdict'] = 'Contradicted'
            results['interpretation'] = ('Age gradients are NEGATIVE (centers older). '
                                        'Standard inside-out formation, contradicting TEP.')
        else:
            results['verdict'] = 'Null'
            results['interpretation'] = 'No significant age gradient detected.'
            
    elif has_lw_mw:
        # Use LW vs MW age discrepancy as proxy
        clean = df.dropna(subset=['LW_AGE_1RE', 'MW_AGE_1RE']).copy()
        print(f"\nAnalyzing LW vs MW age discrepancy for {len(clean)} galaxies")
        results['n_galaxies'] = len(clean)
        
        # LW age emphasizes young stars, MW emphasizes old stars
        # Difference: LW - MW indicates presence of young component
        clean['age_diff'] = clean['LW_AGE_1RE'] - clean['MW_AGE_1RE']
        
        mean_diff = clean['age_diff'].mean()
        sem_diff = clean['age_diff'].std() / np.sqrt(len(clean))
        
        print(f"\nLW - MW Age Difference:")
        print(f"   Mean: {mean_diff:.4f} ± {sem_diff:.4f} log(Gyr)")
        
        # Correlation with sigma
        if 'stellar_sigma_1re' in clean.columns:
            mask = clean['stellar_sigma_1re'].notna()
            r_sigma, p_sigma = stats.pearsonr(clean.loc[mask, 'stellar_sigma_1re'],
                                              clean.loc[mask, 'age_diff'])
            print(f"\nAge diff vs σ: r = {r_sigma:.4f}, p = {p_sigma:.2e}")
            
            results['lw_mw_diff_vs_sigma'] = {
                'r': float(r_sigma),
                'p': float(p_sigma),
                'interpretation': ('Positive r means high-σ galaxies have younger '
                                  'light-weighted ages relative to mass-weighted')
            }
        
        results['lw_mw_analysis'] = {
            'mean_diff': float(mean_diff),
            'stderr': float(sem_diff),
            'note': 'Analyzed LW vs MW age discrepancy as gradient proxy'
        }
        
        results['verdict'] = 'Indirect'
        results['interpretation'] = 'Used LW-MW age comparison (no explicit gradient data available).'
    
    else:
        results['verdict'] = 'Skipped'
        results['interpretation'] = 'Insufficient data for gradient analysis.'
    
    print(f"\n=== VERDICT: {results['verdict']} ===")
    print(results['interpretation'])
    
    return results, df


def create_figure(df, results):
    """Create visualization."""
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    has_gradients = 'age_gradient' in df.columns
    
    if has_gradients:
        clean = df[df['age_gradient'].notna()].copy()
        
        # 1. Histogram of age gradients
        ax1 = axes[0, 0]
        ax1.hist(clean['age_gradient'], bins=50, edgecolor='black', alpha=0.7)
        ax1.axvline(0, color='red', linestyle='--', linewidth=2, label='Zero')
        ax1.axvline(clean['age_gradient'].mean(), color='blue', linestyle='-', 
                   linewidth=2, label=f"Mean = {clean['age_gradient'].mean():.3f}")
        ax1.set_xlabel('Age Gradient (Gyr/Re)')
        ax1.set_ylabel('Count')
        ax1.set_title('Distribution of Radial Age Gradients')
        ax1.legend()
        
        # 2. Age gradient vs mass
        ax2 = axes[0, 1]
        if 'log_mass' in clean.columns:
            mask = clean['log_mass'].notna()
            ax2.scatter(clean.loc[mask, 'log_mass'], clean.loc[mask, 'age_gradient'],
                       alpha=0.3, s=10)
            ax2.axhline(0, color='red', linestyle='--')
            ax2.set_xlabel('log(M*/M☉)')
            ax2.set_ylabel('Age Gradient (Gyr/Re)')
            ax2.set_title('Age Gradient vs Stellar Mass')
        
        # 3. Age gradient vs Z gradient
        ax3 = axes[1, 0]
        if 'z_gradient' in clean.columns:
            mask = clean['z_gradient'].notna()
            ax3.scatter(clean.loc[mask, 'z_gradient'], clean.loc[mask, 'age_gradient'],
                       alpha=0.3, s=10)
            ax3.axhline(0, color='red', linestyle='--')
            ax3.axvline(0, color='red', linestyle='--')
            ax3.set_xlabel('Metallicity Gradient (dex/Re)')
            ax3.set_ylabel('Age Gradient (Gyr/Re)')
            ax3.set_title('Age vs Metallicity Gradient')
        
        # 4. Summary
        ax4 = axes[1, 1]
        ax4.axis('off')
        
        summary = f"""
TEST L: Radial Age Gradient (MaNGA)

HYPOTHESIS:
TEP: Centers appear YOUNGER (gradient > 0)
Standard: Centers are OLDER (gradient < 0)

RESULTS:
• N = {results['n_galaxies']:,} galaxies
• Mean gradient: {results['mean_gradient']['value']:.4f} ± {results['mean_gradient']['stderr']:.4f}
• t = {results['mean_gradient']['t_stat']:.2f}, p = {results['mean_gradient']['p_value']:.2e}
• {results['sign_distribution']['frac_positive']*100:.1f}% have positive gradients

VERDICT: {results['verdict']}
"""
        ax4.text(0.05, 0.95, summary, transform=ax4.transAxes,
                fontsize=10, verticalalignment='top', fontfamily='monospace',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    else:
        # LW vs MW analysis
        if 'LW_AGE_1RE' in df.columns:
            clean = df.dropna(subset=['LW_AGE_1RE', 'MW_AGE_1RE'])
            
            ax1 = axes[0, 0]
            ax1.scatter(clean['MW_AGE_1RE'], clean['LW_AGE_1RE'], alpha=0.3, s=10)
            lims = [clean[['LW_AGE_1RE', 'MW_AGE_1RE']].min().min(),
                   clean[['LW_AGE_1RE', 'MW_AGE_1RE']].max().max()]
            ax1.plot(lims, lims, 'r--', label='1:1')
            ax1.set_xlabel('MW Age (log Gyr)')
            ax1.set_ylabel('LW Age (log Gyr)')
            ax1.set_title('Light-Weighted vs Mass-Weighted Ages')
            ax1.legend()
            
            if 'stellar_sigma_1re' in clean.columns:
                ax2 = axes[0, 1]
                clean['age_diff'] = clean['LW_AGE_1RE'] - clean['MW_AGE_1RE']
                mask = clean['stellar_sigma_1re'].notna()
                ax2.scatter(clean.loc[mask, 'stellar_sigma_1re'], 
                           clean.loc[mask, 'age_diff'], alpha=0.3, s=10)
                ax2.axhline(0, color='red', linestyle='--')
                ax2.set_xlabel('Central σ (km/s)')
                ax2.set_ylabel('LW - MW Age')
                ax2.set_title('Age Difference vs Velocity Dispersion')
        
        ax4 = axes[1, 1]
        ax4.axis('off')
        ax4.text(0.5, 0.5, f"VERDICT: {results['verdict']}\n\n{results['interpretation']}",
                transform=ax4.transAxes, ha='center', va='center', fontsize=12,
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    plt.tight_layout()
    fig_path = os.path.join(FIGURE_DIR, 'test_l_radial_age_gradient.png')
    plt.savefig(fig_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"\nFigure saved to {fig_path}")
    
    return fig_path


def main():
    print("=" * 60)
    print("TEST L: Radial Age Gradient within Galaxies (MaNGA)")
    print("=" * 60)
    
    # Try to get gradient data
    df = download_gradient_data()
    
    # If that fails, use existing MaNGA data
    if df is None or len(df) < 50:
        print("Gradient query failed, using existing MaNGA data...")
        df = load_existing_manga_data()
    
    if df is None or len(df) < 50:
        print("ERROR: Insufficient data for analysis")
        results = {
            'verdict': 'Skipped',
            'interpretation': 'Could not retrieve MaNGA gradient data (HTTP errors).'
        }
        output_file = os.path.join(OUTPUT_DIR, 'sdss_test_l_radial_gradient_results.json')
        with open(output_file, 'w') as f:
            json.dump(results, f, indent=2)
        return results
    
    print(f"\nLoaded {len(df)} galaxies")
    
    # Analyze
    results, clean_df = analyze_radial_gradients(df)
    
    # Create figure
    fig_path = create_figure(clean_df, results)
    results['figure_path'] = fig_path
    
    # Save results
    output_file = os.path.join(OUTPUT_DIR, 'sdss_test_l_radial_gradient_results.json')
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\nResults saved to {output_file}")
    
    return results


if __name__ == '__main__':
    results = main()
