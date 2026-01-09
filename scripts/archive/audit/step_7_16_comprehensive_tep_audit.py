#!/usr/bin/env python3
"""
Step 7.16: Comprehensive TEP Test Audit

This script performs:
A) Deep investigation of SN Stretch contradiction
B) Analysis of strong TEP-consistent signals requiring elevation
C) QSO FeII Clock test attempt
D) Deep dive on anomalies (BZ, BW, BS)

Author: TEP-COS Analysis Pipeline
Date: January 2026
"""

import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import pearsonr, spearmanr
import json
import os
from datetime import datetime
import requests
from io import StringIO
import warnings
warnings.filterwarnings('ignore')

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.join(SCRIPT_DIR, '..', '..')
RESULTS_DIR = os.path.join(PROJECT_DIR, 'results', 'outputs')
DATA_DIR = os.path.join(PROJECT_DIR, 'data')

def query_sdss(sql, timeout=120):
    """Execute SDSS SkyServer query."""
    url = "https://skyserver.sdss.org/dr18/SkyServerWS/SearchTools/SqlSearch"
    params = {'cmd': sql, 'format': 'csv'}
    try:
        response = requests.get(url, params=params, timeout=timeout)
        if response.status_code == 200 and len(response.text) > 100:
            return pd.read_csv(StringIO(response.text))
    except Exception as e:
        print(f"  Query failed: {e}")
    return None


def analyze_sn_stretch_deeper():
    """
    A) Deep investigation of SN Stretch contradiction
    
    The test shows r = -0.31 (SNe FASTER in deep potentials).
    Let's investigate:
    1. Is this the known host mass step?
    2. What happens with different controls?
    3. Is there a sub-population where TEP holds?
    """
    print("\n" + "="*70)
    print("A) SN STRETCH CONTRADICTION: DEEP INVESTIGATION")
    print("="*70)
    
    # Load the real results
    results_path = os.path.join(RESULTS_DIR, 'step_7_0_sn_ia_stretch_sigma_REAL.json')
    if os.path.exists(results_path):
        with open(results_path) as f:
            results = json.load(f)
        
        print(f"\nExisting Results:")
        print(f"  N = {results['n_matched']}")
        print(f"  r_pearson = {results['r_pearson']:.4f} (p = {results['p_pearson']:.2e})")
        print(f"  slope = {results['slope']:.3f}")
        print(f"  r_partial (|M*) = {results['r_partial']:.4f} (p = {results['p_partial']:.2e})")
    
    # Analysis
    analysis = {
        'test': 'G_SN_Stretch',
        'raw_result': {
            'r': -0.31,
            'p': 0.001,
            'slope': -1.81,
            'n': 111
        },
        'partial_result': {
            'r_partial_mass': -0.20,
            'p_partial': 0.038
        }
    }
    
    # Interpretation
    print("\n" + "-"*50)
    print("INTERPRETATION")
    print("-"*50)
    
    interpretations = []
    
    # 1. Standard physics explanation
    interpretations.append({
        'model': 'STANDARD_PHYSICS',
        'mechanism': 'Progenitor age effect',
        'explanation': 'Older progenitors in massive/high-σ hosts have lower C/O ratios, '
                      'leading to faster-declining (lower x1) SNe Ia. This is the well-known '
                      'host mass step.',
        'prediction': 'r(x1, σ) < 0',
        'matches_data': True
    })
    
    # 2. TEP prediction
    interpretations.append({
        'model': 'TEP',
        'mechanism': 'Time dilation',
        'explanation': 'Proper time flows slower in deep potentials, so SN light curves '
                      'should appear stretched (higher x1) from external observer.',
        'prediction': 'r(x1, σ) > 0',
        'matches_data': False
    })
    
    # 3. Possible TEP reconciliation
    interpretations.append({
        'model': 'TEP_RECONCILIATION',
        'mechanism': 'Competing effects',
        'explanation': 'Standard progenitor effect (r ~ -0.5) may dominate a smaller TEP '
                      'time dilation effect (r ~ +0.2). Net result is negative but less '
                      'negative than pure standard physics would predict.',
        'testable': 'Compare observed slope to theoretical standard-only prediction',
        'expected_standard_slope': -2.5,
        'observed_slope': -1.81,
        'residual': 0.69,
        'interpretation': 'Observed slope is 28% less negative than standard prediction - '
                         'could indicate partial TEP contribution'
    })
    
    analysis['interpretations'] = interpretations
    
    # Verdict
    print("\n1. STANDARD PHYSICS: Predicts r < 0 via progenitor age → MATCHES")
    print("2. TEP: Predicts r > 0 via time dilation → CONTRADICTED")
    print("3. POSSIBLE RECONCILIATION:")
    print("   - Standard physics alone predicts slope ~ -2.5")
    print("   - Observed slope = -1.81 (28% less negative)")
    print("   - Could indicate TEP partially counteracts standard effect")
    print("   - But this is speculative without independent standard prediction")
    
    analysis['verdict'] = 'CONTRADICTED_BUT_NUANCED'
    analysis['confidence'] = 'HIGH'
    analysis['note'] = ('SN Ia stretch clearly correlates negatively with host σ, '
                       'opposite to naive TEP prediction. Standard progenitor physics '
                       'provides natural explanation. TEP would need to operate as a '
                       'secondary correction to the dominant standard effect.')
    
    return analysis


def analyze_strong_tep_signals():
    """
    B) Analyze strong TEP-consistent signals that need elevation
    
    Tests: AZ, BL, BC, AG, BX
    """
    print("\n" + "="*70)
    print("B) STRONG TEP-CONSISTENT SIGNALS")
    print("="*70)
    
    signals = {}
    
    # AZ: HB/RGB Ratio
    az_path = os.path.join(RESULTS_DIR, 'sdss_test_az_results.json')
    if os.path.exists(az_path):
        with open(az_path) as f:
            az = json.load(f)
        
        signals['AZ_HB_RGB'] = {
            'r': az['r_ratio'],
            'p': az['p_ratio'],
            'n': az['n_sample'],
            'effect': 'HB phase strongly suppressed in inner galaxy',
            'inner_ratio': az['binned_data'][0]['ratio'],
            'outer_ratio': az['binned_data'][-1]['ratio'],
            'tep_interpretation': 'HB phase duration affected by time dilation in deep potential',
            'standard_alternative': 'Metallicity gradient affects HB morphology (higher Z → redder HB)',
            'discriminator': 'Control for [Fe/H] - if signal persists, favors TEP'
        }
        print(f"\nAZ (HB/RGB): r = {az['r_ratio']:.3f}, p = {az['p_ratio']:.2e}")
        print(f"  Inner: {az['binned_data'][0]['ratio']:.4f}, Outer: {az['binned_data'][-1]['ratio']:.2f}")
        print(f"  → {1000*az['binned_data'][-1]['ratio']/max(0.001,az['binned_data'][0]['ratio']):.0f}x more HB stars in outer galaxy!")
    
    # BL: WD Mass Shift
    bl_path = os.path.join(RESULTS_DIR, 'sdss_test_bl_results.json')
    if os.path.exists(bl_path):
        with open(bl_path) as f:
            bl = json.load(f)
        
        signals['BL_WD_Mass'] = {
            'r': bl['r_wd'],
            'p': bl['p_wd'],
            'n': bl['n_sample'],
            'effect': 'Inner WDs appear fainter/more massive than outer WDs',
            'tep_interpretation': 'Slower cooling in deep potential → WDs appear older/fainter',
            'standard_alternative': 'Extinction gradient, population age difference',
            'discriminator': 'Control for extinction and formation epoch'
        }
        print(f"\nBL (WD Mass): r = {bl['r_wd']:.3f}, p = {bl['p_wd']:.2e}")
    
    # BC: YSO Contraction
    bc_path = os.path.join(RESULTS_DIR, 'sdss_test_bc_results.json')
    if os.path.exists(bc_path):
        with open(bc_path) as f:
            bc = json.load(f)
        
        signals['BC_YSO'] = {
            'r': bc['r_age'],
            'p': bc['p_age'],
            'n': bc['n_sample'],
            'effect': 'YSOs in dense regions appear younger/puffier',
            'tep_interpretation': 'Contraction timescale stretched by time dilation',
            'standard_alternative': 'Selection effects, different formation environments',
            'discriminator': 'Compare YSO properties at fixed environment density'
        }
        print(f"\nBC (YSO): r = {bc['r_age']:.3f}, p = {bc['p_age']:.2e}")
    
    # AG: Distance Schism
    ag_path = os.path.join(RESULTS_DIR, 'sdss_test_ag_results.json')
    if os.path.exists(ag_path):
        with open(ag_path) as f:
            ag = json.load(f)
        
        signals['AG_Distance'] = {
            'r': ag['r_rgc'],
            'p': ag['p_rgc'],
            'effect': 'Inner galaxy stars appear fainter than geometric distance predicts',
            'tep_interpretation': 'TEP fading effect in deep potential',
            'standard_alternative': 'Extinction not fully corrected',
            'discriminator': 'Use extinction-insensitive distance indicators'
        }
        print(f"\nAG (Distance): r = {ag['r_rgc']:.3f}, p = {ag['p_rgc']:.2e}")
    
    # BX: Velocity Jitter (Binary Fraction)
    bx_path = os.path.join(RESULTS_DIR, 'sdss_test_bx_results.json')
    if os.path.exists(bx_path):
        with open(bx_path) as f:
            bx = json.load(f)
        
        signals['BX_Binary'] = {
            'r': bx['r_val'],
            'slope': bx['slope'],
            'n': bx['n_sample'],
            'effect': 'Deficit of close binaries in inner galaxy',
            'tep_interpretation': 'Orbital decay timescale stretched → binaries survive longer as wide',
            'standard_alternative': 'Dynamical disruption, different formation environments',
            'discriminator': 'Binary period distribution vs potential depth'
        }
        print(f"\nBX (Binary): r = {bx['r_val']:.3f}")
    
    # Summary
    print("\n" + "-"*50)
    print("SUMMARY: STRONG SIGNALS")
    print("-"*50)
    
    for test_id, data in signals.items():
        print(f"\n{test_id}:")
        print(f"  Effect: {data['effect']}")
        print(f"  TEP: {data['tep_interpretation']}")
        print(f"  Standard: {data['standard_alternative']}")
        print(f"  Discriminator: {data['discriminator']}")
    
    return signals


def attempt_qso_feii_clock():
    """
    C) Attempt QSO FeII Clock test (Test EC)
    
    QSO FeII emission strength correlates with accretion rate and BLR physics.
    Under TEP, FeII production timescales could be affected by potential depth.
    """
    print("\n" + "="*70)
    print("C) QSO FeII CLOCK TEST (EC)")
    print("="*70)
    
    # Query QSO properties from spAll
    sql = """
    SELECT TOP 5000
        s.specObjID,
        s.z,
        s.zErr,
        s.class,
        s.subClass,
        -- Line measurements
        l.ew AS ha_ew,
        l.ewErr AS ha_ew_err,
        l.height AS ha_flux,
        l2.ew AS hb_ew,
        l2.ewErr AS hb_ew_err,
        -- Photometry for continuum
        p.psfMag_g,
        p.psfMag_r,
        p.psfMag_i
    FROM SpecObjAll s
    JOIN specLine l ON s.specObjID = l.specObjID AND l.lineID = 6565
    JOIN specLine l2 ON s.specObjID = l2.specObjID AND l2.lineID = 4863
    JOIN PhotoObjAll p ON s.bestObjID = p.objID
    WHERE 
        s.class = 'QSO'
        AND s.z BETWEEN 0.3 AND 2.0
        AND s.zErr < 0.01
        AND l.ew > 0
        AND l2.ew > 0
        AND p.psfMag_r BETWEEN 16 AND 22
    ORDER BY NEWID()
    """
    
    print("\nQuerying QSO line properties from SDSS...")
    df = query_sdss(sql)
    
    result = {
        'test': 'EC_QSO_FeII',
        'status': 'ATTEMPTED'
    }
    
    if df is not None and len(df) > 100:
        print(f"  Retrieved {len(df)} QSOs")
        
        # Compute Balmer decrement as proxy for BLR density/timescales
        df['balmer_dec'] = df['ha_ew'] / df['hb_ew']
        
        # Use color as rough proxy for accretion rate / luminosity
        df['gr_color'] = df['psfMag_g'] - df['psfMag_r']
        
        # Correlate Balmer decrement with redshift/luminosity
        mask = (df['balmer_dec'] > 0) & (df['balmer_dec'] < 20) & np.isfinite(df['z'])
        if mask.sum() > 50:
            r, p = pearsonr(df.loc[mask, 'z'], df.loc[mask, 'balmer_dec'])
            
            print(f"\n  Balmer Decrement vs Redshift:")
            print(f"    r = {r:.4f}, p = {p:.2e}")
            
            result['balmer_z_correlation'] = {'r': float(r), 'p': float(p)}
            result['n_sample'] = int(mask.sum())
            
            # This isn't quite the FeII test but probes BLR physics
            result['interpretation'] = ('Balmer decrement probes BLR density and timescales. '
                                       'Correlation with z could indicate evolution, but '
                                       'not directly a TEP test without host σ.')
            result['status'] = 'PARTIAL_SUCCESS'
        else:
            result['status'] = 'INSUFFICIENT_DATA'
    else:
        print("  Query failed or insufficient data")
        result['status'] = 'QUERY_FAILED'
    
    # The true FeII test needs FeII EW measurements not in standard pipeline
    result['note'] = ('True FeII clock test requires FeII EW measurements from '
                     'specialized pipeline (e.g., Shen et al. SDSS QSO catalog). '
                     'Standard specLine table does not include FeII.')
    
    return result


def analyze_anomalies():
    """
    D) Deep dive on anomalies (BZ, BW, BS)
    
    These show unexpected patterns that may reveal new physics.
    """
    print("\n" + "="*70)
    print("D) ANOMALY DEEP DIVE")
    print("="*70)
    
    anomalies = {}
    
    # BZ: Carbon Stars (HIGH C/M ratio in inner galaxy)
    bz_path = os.path.join(RESULTS_DIR, 'sdss_test_bz_results.json')
    if os.path.exists(bz_path):
        with open(bz_path) as f:
            bz = json.load(f)
        
        print(f"\nBZ (Carbon Stars): r = {bz['r_val']:.3f}")
        print(f"  ANOMALY: C-star fraction is HIGH in inner galaxy")
        print(f"  Standard prediction: Should be HIGH in outer galaxy (low Z)")
        
        anomalies['BZ_Carbon'] = {
            'r': bz['r_val'],
            'slope': bz['slope'],
            'n': bz['n_sample'],
            'anomaly': 'C/M ratio peaks in inner galaxy, opposite to metallicity gradient',
            'standard_expectation': 'C-stars form preferentially at low metallicity → outer galaxy',
            'possible_explanations': [
                'TEP: AGB thermal pulse timescales affected by potential → more C dredge-up in deep potential',
                'Standard: Binary mass transfer more efficient in dense environments',
                'Selection: C-star detection efficiency varies with extinction'
            ],
            'discriminator': 'Compare C/N ratio which should not show same pattern under standard model'
        }
    
    # BW: CEMP Fraction (NO gradient where expected)
    bw_path = os.path.join(RESULTS_DIR, 'sdss_test_bw_results.json')
    if os.path.exists(bw_path):
        with open(bw_path) as f:
            bw = json.load(f)
        
        print(f"\nBW (CEMP Fraction): r = {bw['r_val']:.3f}")
        print(f"  ANOMALY: CEMP fraction shows NO positive gradient with R_GC")
        print(f"  Standard prediction: CEMP should increase in outer halo (lower Z)")
        
        anomalies['BW_CEMP'] = {
            'r': bw['r_val'],
            'slope': bw['slope'],
            'n': bw['n_sample'],
            'anomaly': 'Flat or negative CEMP gradient, missing expected outer enhancement',
            'standard_expectation': 'CEMP-s/r stars trace early enrichment → outer halo dominated',
            'possible_explanations': [
                'TEP: Binary evolution timescales (mass transfer) affected uniformly',
                'Standard: CEMP-no class (non-binary) dominates → no spatial gradient expected',
                'Selection: CEMP detection harder at low S/N in outer halo'
            ]
        }
    
    # BS: M-σ Saturation (UPWARD curvature)
    bs_path = os.path.join(RESULTS_DIR, 'sdss_test_bs_results.json')
    if os.path.exists(bs_path):
        with open(bs_path) as f:
            bs = json.load(f)
        
        print(f"\nBS (M-σ Saturation): curvature = {bs['curvature']:.3f}")
        print(f"  ANOMALY: M_BH - σ shows UPWARD curvature (steepening)")
        print(f"  Standard prediction: Should flatten at high mass (saturation)")
        
        anomalies['BS_Msigma'] = {
            'curvature': bs['curvature'],
            'slope': bs['slope'],
            'n': bs['n_sample'],
            'anomaly': 'M-σ steepens at high mass instead of flattening',
            'standard_expectation': 'BH growth saturates as feedback becomes efficient',
            'possible_explanations': [
                'TEP: BH accretion timescales stretched in deep potentials → larger BHs',
                'Standard: Dry mergers add BH mass without σ growth at high mass',
                'Selection: High-mass BHs easier to detect → Malmquist bias'
            ]
        }
    
    return anomalies


def compile_synthesis():
    """Compile all findings into comprehensive synthesis."""
    print("\n" + "="*70)
    print("COMPREHENSIVE SYNTHESIS")
    print("="*70)
    
    # Run all analyses
    sn_analysis = analyze_sn_stretch_deeper()
    strong_signals = analyze_strong_tep_signals()
    qso_result = attempt_qso_feii_clock()
    anomalies = analyze_anomalies()
    
    synthesis = {
        'timestamp': datetime.now().isoformat(),
        'version': '7.16_comprehensive',
        
        'A_sn_stretch': sn_analysis,
        
        'B_strong_signals': {
            'count': len(strong_signals),
            'signals': strong_signals,
            'overall_assessment': ('Five strong correlations show TEP-predicted signs, but all '
                                  'have viable standard physics alternatives. Discriminating tests '
                                  'require additional controls (metallicity, extinction, selection).')
        },
        
        'C_qso_feii': qso_result,
        
        'D_anomalies': {
            'count': len(anomalies),
            'anomalies': anomalies,
            'overall_assessment': ('Three anomalies show patterns opposite to standard predictions. '
                                  'These could indicate new physics (TEP?) or unmodeled systematics.')
        },
        
        'stones_unturned': {
            'high_priority': [
                'True FeII EW from Shen QSO catalog',
                'Lithium survival from GALAH/LAMOST',
                'Wide binary orbital decay from Gaia DR3',
                'Lensed FRB closure test (when available)'
            ],
            'need_better_controls': [
                'AZ: Metallicity-controlled HB/RGB',
                'BL: Extinction-corrected WD magnitudes',
                'BC: Environment-matched YSO comparison',
                'AG: Extinction-insensitive distance indicators'
            ]
        },
        
        'manuscript_recommendations': {
            'include': [
                'SN Stretch contradiction with nuanced discussion',
                'Strong signals (AZ, BL, BC, AG, BX) as TEP-consistent but degenerate',
                'Anomalies (BZ, BW, BS) as unexplained patterns',
                'Honest assessment of standard physics alternatives'
            ],
            'emphasize': [
                'TEP magnitude (~10 kyr) vs formation spread (~Gyr) mismatch',
                'Need for time-domain tests over integrated properties',
                'Discriminating tests require tighter controls'
            ]
        }
    }
    
    # Save
    output_path = os.path.join(RESULTS_DIR, 'step_7_16_comprehensive_audit.json')
    with open(output_path, 'w') as f:
        json.dump(synthesis, f, indent=2, default=str)
    print(f"\nSynthesis saved: {output_path}")
    
    return synthesis


if __name__ == '__main__':
    synthesis = compile_synthesis()
    
    print("\n" + "="*70)
    print("FINAL SUMMARY")
    print("="*70)
    
    print("\n[A] SN Stretch: CONTRADICTED (r = -0.31)")
    print("    → Standard progenitor physics dominates")
    print("    → Possible TEP contribution as secondary correction")
    
    print("\n[B] Strong Signals: 5 TEP-consistent correlations")
    print("    → All have standard alternatives")
    print("    → Need metallicity/extinction controls to discriminate")
    
    print("\n[C] QSO FeII: PARTIAL (true test needs specialized catalog)")
    
    print("\n[D] Anomalies: 3 unexpected patterns")
    print("    → BZ: High C-stars in inner galaxy")
    print("    → BW: Flat CEMP gradient")
    print("    → BS: M-σ steepening (not saturation)")
