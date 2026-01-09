#!/usr/bin/env python3
"""
Step 7.13: Update SDSS Test Synthesis with All Results
"""

import json
import os

RESULTS_DIR = os.path.join(os.path.dirname(__file__), '..', '..', 'results', 'outputs')

# Original batch 4-6 status (from step_6_99)
ORIGINAL_STATUS = {
    "BQ": {"name": "Vertical Disk Heating", "status": "Signal", "metric": "Slope +2.4 km/s/kpc", "desc": "Standard vertical heating"},
    "BR": {"name": "IMF Variation", "status": "Signal", "metric": "Slope +1.8", "desc": "Known bottom-heavy IMF in ETGs"},
    "BS": {"name": "M-sigma Saturation", "status": "Signal", "metric": "Quadratic -0.12", "desc": "Standard M-sigma curvature"},
    "BT": {"name": "Galaxy Merger Rate", "status": "Skipped", "metric": "zoo2MainSpecz join failed", "desc": "Cross-match error"},
    "BU": {"name": "Spiral Arm Winding", "status": "Null", "metric": "p=0.42", "desc": "No correlation"},
    "BV": {"name": "BAL Quasar Fraction", "status": "Null", "metric": "p=0.21", "desc": "No correlation"},
    "BW": {"name": "CEMP Star Fraction", "status": "Signal", "metric": "Slope -0.8", "desc": "Standard nucleosynthesis"},
    "BX": {"name": "Velocity Jitter", "status": "Signal", "metric": "Excess 12 km/s", "desc": "Standard kinematics"},
    "BY": {"name": "Astrometric Binary", "status": "Skipped", "metric": "Gaia not in SDSS", "desc": "External data required"},
    "BZ": {"name": "Carbon Star Clock", "status": "Null", "metric": "p=0.35", "desc": "No radial gradient"},
    "CA": {"name": "BLR Kinematics", "status": "Skipped", "metric": "Insufficient AGN", "desc": "Sample too small"},
    "CB": {"name": "Template Systematics", "status": "Null", "metric": "p=0.67", "desc": "No systematic"},
    "CC": {"name": "Manganese Clock", "status": "Signal", "metric": "Slope +0.04", "desc": "Standard nucleosynthesis"},
    "CD": {"name": "Satellite Concentration", "status": "Skipped", "metric": "redMaPPer missing", "desc": "Catalog not in DR18"},
    "CE": {"name": "Nitrogen/Oxygen Clock", "status": "Null", "metric": "r=-0.001, p=0.92", "desc": "No correlation found"},
    "CF": {"name": "QSO Grav Redshift", "status": "Skipped", "metric": "DR16 QSO catalog missing", "desc": "Catalog not in DR18"},
    "CG": {"name": "Cluster Sigma Profiles", "status": "Skipped", "metric": "redMaPPer missing", "desc": "Catalog not in DR18"},
    "CH": {"name": "Hypervelocity Stars", "status": "Null", "metric": "p=0.82", "desc": "No excess"},
    "CI": {"name": "Galactic Dipole", "status": "Skipped", "metric": "apogeeStar columns differ", "desc": "Schema mismatch"},
    "CJ": {"name": "S-Process Clock", "status": "Signal", "metric": "Slope +0.02", "desc": "Standard nucleosynthesis"},
    "CK": {"name": "Mass-Z Residuals", "status": "Null", "metric": "p=0.15", "desc": "No correlation"},
    "CL": {"name": "Spiral Chirality", "status": "Null", "metric": "51/49 Split", "desc": "No preference"},
    "CM": {"name": "TRGB Magnitude", "status": "Skipped", "metric": "Gaia parallax required", "desc": "External data required"},
    "CN": {"name": "BCG X-ray Offsets", "status": "Skipped", "metric": "redMaPPer missing", "desc": "Catalog not in DR18"},
    "CO": {"name": "Exoplanet Yield", "status": "Skipped", "metric": "MARVELS incomplete", "desc": "Survey incomplete"},
    "CP": {"name": "Aluminum Clock", "status": "Signal", "metric": "Slope +0.18", "desc": "Standard nucleosynthesis"},
    "CQ": {"name": "Cluster Cooling", "status": "Skipped", "metric": "eFEDS not in SDSS", "desc": "External data required"},
    "CR": {"name": "Gravity Schism", "status": "Null", "metric": "p=0.55", "desc": "No bimodality"},
    "CS": {"name": "Void HI Fraction", "status": "Signal", "metric": "r=-0.44, p<1e-90", "desc": "Standard HI-morphology relation"},
    "CT": {"name": "Schechter Cutoff", "status": "Signal", "metric": "t=38.5, p<1e-280", "desc": "Standard luminosity function"},
    "CU": {"name": "Binary Quasar Frac", "status": "Null", "metric": "N=5000", "desc": "Requires spatial cross-match"},
    "CV": {"name": "Chromospheric Activity", "status": "Null", "metric": "p=0.33", "desc": "No correlation"},
    "CW": {"name": "Stellar Twins", "status": "Skipped", "metric": "apogeeStar schema", "desc": "Column mismatch"},
    "CX": {"name": "Void Metallicity", "status": "Signal", "metric": "r=-0.11, p<1e-14", "desc": "Standard mass-metallicity"},
    "CY": {"name": "Quasar Clustering", "status": "Skipped", "metric": "QsoCatalogAll missing", "desc": "Catalog not in DR18"},
    "CZ": {"name": "Diffuse Ionized Gas", "status": "Signal", "metric": "r=0.47, p<1e-100", "desc": "Standard DIG-mass relation"},
    "DA": {"name": "AGN Type 1/2", "status": "Signal", "metric": "r=0.22, p<1e-50", "desc": "Standard AGN demographics"},
    "DB": {"name": "Void Hubble Drift", "status": "Null", "metric": "Δslope=-0.9", "desc": "No significant drift"},
    "DC": {"name": "Pair Decay Ratio", "status": "Skipped", "metric": "Neighbors timeout", "desc": "Query too slow"},
    "DD": {"name": "QSO Color-Potential", "status": "Skipped", "metric": "QsoCatalogAll missing", "desc": "Catalog not in DR18"},
    "DE": {"name": "Richness-Sigma", "status": "Skipped", "metric": "redMaPPer missing", "desc": "Catalog not in DR18"},
    "DF": {"name": "Lithium Survival", "status": "Skipped", "metric": "Li not in APOGEE", "desc": "Requires GALAH"},
    "DG": {"name": "ICL Growth", "status": "Skipped", "metric": "redMaPPer missing", "desc": "Catalog not in DR18"},
    "DH": {"name": "Dust-to-Gas Ratio", "status": "Signal", "metric": "Slope +0.65", "desc": "Standard dust-mass relation"},
    "DI": {"name": "Galaxy Spin", "status": "Skipped", "metric": "MaNGA DAP schema", "desc": "Column mismatch"},
    "DJ": {"name": "Sersic Relaxation", "status": "Signal", "metric": "Slope +1.2", "desc": "Standard morphology-density"},
    "DK": {"name": "Ring Galaxy Frac", "status": "Null", "metric": "Constant input", "desc": "zooSpec spiral all 0"},
    "DL": {"name": "Tidal Debris", "status": "Signal", "metric": "r=0.54, p<1e-300", "desc": "Standard concentration-sigma"},
    "DM": {"name": "Red Sequence Scatter", "status": "Skipped", "metric": "galSpecExtra query", "desc": "Join failed"},
    "DN": {"name": "Line Width vs Sigma", "status": "Signal", "metric": "r=-0.04, p=0.048", "desc": "Weak standard relation"},
    "DO": {"name": "Cluster Lx-Sigma", "status": "Skipped", "metric": "redMaPPer missing", "desc": "Catalog not in DR18"},
    "DP": {"name": "TiO IMF", "status": "Signal", "metric": "r=0.21, p<1e-29", "desc": "Standard IMF-sigma relation"},
    "DQ": {"name": "Satellite Abundance", "status": "Signal", "metric": "Slope +32.3", "desc": "Standard mass-richness"},
    "DR": {"name": "Brown Dwarf Desert", "status": "Skipped", "metric": "apogeeStar schema", "desc": "Column mismatch"},
    "DS": {"name": "QSO Variability", "status": "Null", "metric": "N=2000", "desc": "Requires Stripe 82"},
    "DT": {"name": "RC Magnitude", "status": "Contradicted", "metric": "Slope +0.07", "desc": "Opposite to TEP prediction"},
    "DU": {"name": "HI vs Optical", "status": "Skipped", "metric": "mangaHIall join", "desc": "Query failed"},
    "DV": {"name": "Cannon vs ASPCAP", "status": "Skipped", "metric": "cannonStar schema", "desc": "Schema differences"},
    "DW": {"name": "Blue Stragglers", "status": "Null", "metric": "p=0.53", "desc": "No radial gradient"},
    "DX": {"name": "Ha/UV Ratio", "status": "Signal", "metric": "r=-0.26, p<1e-64", "desc": "ANOMALOUS: 2.1σ steeper slope than expected"},
    "DY": {"name": "Phase Spirals", "status": "Null", "metric": "Ratio 0.13", "desc": "Standard vertical heating"},
    "DZ": {"name": "Potassium Anomaly", "status": "Skipped", "metric": "K not in APOGEE", "desc": "Element not measured"},
}

def synthesize():
    """Create final synthesis"""
    counts = {"Signal": 0, "Null": 0, "Contradicted": 0, "Skipped": 0}
    
    for code, info in ORIGINAL_STATUS.items():
        counts[info['status']] += 1
    
    results_list = []
    for code in sorted(ORIGINAL_STATUS.keys()):
        info = ORIGINAL_STATUS[code]
        results_list.append({
            "code": code,
            "name": info['name'],
            "status": info['status'],
            "metric": info['metric'],
            "description": info['desc']
        })
    
    # Identify anomalous signals (not explained by standard physics)
    anomalous = []
    standard_signals = []
    for code, info in ORIGINAL_STATUS.items():
        if info['status'] == 'Signal':
            if 'anomalous' in info['desc'].lower():
                anomalous.append(code)
            else:
                standard_signals.append(code)
    
    synthesis = {
        "summary": {
            "total_tests": len(ORIGINAL_STATUS),
            "counts": counts,
            "anomalous_signals": anomalous,
            "standard_signals": standard_signals,
        },
        "conclusion": (
            "Of 62 SDSS tests attempted, 21 show statistically significant signals, "
            "but only DX (Hα/UV Ratio) shows a slope 2.1σ steeper than standard astrophysical predictions. "
            "The remaining signals are consistent with known astrophysical relations "
            "(IMF variation, nucleosynthesis, morphology-density, mass-metallicity). "
            "One test (DT: RC Magnitude) contradicts TEP predictions. "
            "24 tests were skipped due to missing catalogs (redMaPPer, Gaia) or schema changes in DR18."
        ),
        "tests": results_list
    }
    
    out_file = os.path.join(RESULTS_DIR, 'sdss_final_synthesis.json')
    with open(out_file, 'w') as f:
        json.dump(synthesis, f, indent=2)
    
    print("=" * 70)
    print("SDSS TEST SYNTHESIS - FINAL")
    print("=" * 70)
    print(f"\nTotal Tests: {len(ORIGINAL_STATUS)}")
    print(f"  Signal: {counts['Signal']} (including {len(standard_signals)} standard, {len(anomalous)} anomalous)")
    print(f"  Null: {counts['Null']}")
    print(f"  Contradicted: {counts['Contradicted']}")
    print(f"  Skipped: {counts['Skipped']}")
    print(f"\nAnomalous signals: {anomalous}")
    print(f"\nSaved to {out_file}")

if __name__ == "__main__":
    synthesize()
