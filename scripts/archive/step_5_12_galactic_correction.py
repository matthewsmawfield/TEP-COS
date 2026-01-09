#!/usr/bin/env python3
"""
Step 5.12: Galactic Potential Correction for Field MSPs

This script applies corrections to field pulsar Ṗ for:
1. Shklovskii effect (transverse velocity contribution)
2. Galactic vertical acceleration (z-component)
3. Galactic differential rotation (radial component)

These corrections make the field sample a cleaner "control" for comparison with GC pulsars.

Author: M. Smawfield
Date: 2026-01-03
"""

import json
import hashlib
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional
import numpy as np
from scipy import stats

# Paths
SCRIPT_DIR = Path(__file__).parent
PROJECT_ROOT = SCRIPT_DIR.parent.parent
RESULTS_DIR = PROJECT_ROOT / "results" / "outputs"
ATNF_DB_PATH = RESULTS_DIR / "atnf_psrcat.db"

# Constants
MSP_PERIOD_CUT_MS = 30.0  # P < 30 ms defines MSP
C_M_S = 299792458.0  # Speed of light in m/s
KPC_TO_M = 3.0857e19  # 1 kpc in meters
YR_TO_S = 365.25 * 24 * 3600  # 1 year in seconds
MAS_YR_TO_RAD_S = (1e-3 / 3600) * (np.pi / 180) / YR_TO_S  # mas/yr to rad/s

# Galactic parameters (Milky Way model)
R0_KPC = 8.34  # Solar galactocentric distance (Reid et al. 2014)
V0_KM_S = 240.0  # Circular velocity at Sun (Reid et al. 2014)
OORT_A = 15.3  # km/s/kpc (Oort constant A)
OORT_B = -11.9  # km/s/kpc (Oort constant B)
RHO_DISK = 0.1  # M_sun/pc^3 (local disk density)
G_SI = 6.674e-11  # Gravitational constant

# Regex for parsing
_NUM_RE = re.compile(r'^[+-]?\d+\.?\d*(?:[eE][+-]?\d+)?$')


def parse_atnf_extended(db_text: str) -> list[dict]:
    """
    Parse ATNF psrcat.db extracting extended parameters for Galactic corrections.
    
    Extracts: PSRJ, P0, P1, F0, F1, RAJ, DECJ, PMRA, PMDEC, DIST_DM, GL, GB, ASSOC
    """
    rows = []
    current = {}
    
    def parse_ra(ra_str: str) -> Optional[float]:
        """Parse RA in HH:MM:SS.sss format to degrees."""
        try:
            parts = ra_str.split(':')
            if len(parts) >= 2:
                h = float(parts[0])
                m = float(parts[1])
                s = float(parts[2]) if len(parts) > 2 else 0.0
                return 15.0 * (h + m/60 + s/3600)
        except:
            pass
        return None
    
    def parse_dec(dec_str: str) -> Optional[float]:
        """Parse Dec in DD:MM:SS.sss format to degrees."""
        try:
            sign = -1 if dec_str.startswith('-') else 1
            dec_str = dec_str.lstrip('+-')
            parts = dec_str.split(':')
            if len(parts) >= 2:
                d = float(parts[0])
                m = float(parts[1])
                s = float(parts[2]) if len(parts) > 2 else 0.0
                return sign * (d + m/60 + s/3600)
        except:
            pass
        return None
    
    def flush():
        nonlocal current
        if not current:
            return
        
        name = current.get("PSRJ") or current.get("PSRB")
        
        # Convert F0/F1 to P0/P1 if needed
        p0 = current.get("P0")
        p1 = current.get("P1")
        f0 = current.get("F0")
        f1 = current.get("F1")
        if p0 is None and f0 is not None and f0 != 0:
            p0 = 1.0 / f0
        if p1 is None and f0 is not None and f1 is not None and f0 != 0:
            p1 = -f1 / (f0 * f0)
        
        if name is None or p0 is None or p1 is None:
            current = {}
            return
        
        # Parse coordinates
        ra_deg = parse_ra(current.get("RAJ", ""))
        dec_deg = parse_dec(current.get("DECJ", ""))
        
        rows.append({
            "name": name,
            "P0_s": p0,
            "P_ms": p0 * 1000.0,
            "P1_sps": p1,
            "P1_e20": p1 / 1e-20,
            "ra_deg": ra_deg,
            "dec_deg": dec_deg,
            "pmra_mas_yr": current.get("PMRA"),
            "pmdec_mas_yr": current.get("PMDEC"),
            "dist_kpc": current.get("DIST_DM"),
            "gl_deg": current.get("GL"),
            "gb_deg": current.get("GB"),
            "assoc": current.get("ASSOC", ""),
        })
        current = {}
    
    for raw_line in db_text.splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        if line.startswith("@"):
            flush()
            continue
        
        parts = line.split()
        if len(parts) < 2:
            continue
        
        key = parts[0].strip()
        val0 = parts[1].strip()
        val_rest = " ".join(parts[1:]).strip()
        
        # Numeric fields
        if key in {"P0", "P1", "F0", "F1", "PMRA", "PMDEC", "DIST_DM", "GL", "GB"}:
            if _NUM_RE.match(val0):
                current[key] = float(val0)
            continue
        
        # String fields
        if key in {"PSRJ", "PSRB"}:
            current[key] = val0
        elif key == "ASSOC":
            current[key] = val_rest
        elif key == "RAJ":
            current[key] = val0
        elif key == "DECJ":
            current[key] = val0
    
    flush()
    return rows


def equatorial_to_galactic(ra_deg: float, dec_deg: float) -> tuple[float, float]:
    """Convert equatorial (J2000) to Galactic coordinates."""
    # Galactic pole (J2000): RA = 192.85948°, Dec = 27.12825°
    # Galactic center longitude: l = 0° at RA = 266.40510°
    ra_ngp = np.radians(192.85948)
    dec_ngp = np.radians(27.12825)
    l_ncp = np.radians(122.93192)  # l of North Celestial Pole
    
    ra = np.radians(ra_deg)
    dec = np.radians(dec_deg)
    
    sin_b = np.sin(dec) * np.sin(dec_ngp) + np.cos(dec) * np.cos(dec_ngp) * np.cos(ra - ra_ngp)
    b = np.arcsin(sin_b)
    
    cos_l_minus_l_ncp = (np.sin(dec) - sin_b * np.sin(dec_ngp)) / (np.cos(b) * np.cos(dec_ngp))
    sin_l_minus_l_ncp = np.cos(dec) * np.sin(ra - ra_ngp) / np.cos(b)
    l = l_ncp - np.arctan2(sin_l_minus_l_ncp, cos_l_minus_l_ncp)
    
    l_deg = np.degrees(l) % 360
    b_deg = np.degrees(b)
    
    return l_deg, b_deg


def compute_shklovskii_pdot(P_s: float, pm_total_mas_yr: float, dist_kpc: float) -> float:
    """
    Compute Shklovskii contribution to Ṗ.
    
    Ṗ_shk = P * μ² * d / c
    
    where μ is proper motion in rad/s, d is distance in meters.
    """
    mu_rad_s = pm_total_mas_yr * MAS_YR_TO_RAD_S
    d_m = dist_kpc * KPC_TO_M
    return P_s * mu_rad_s**2 * d_m / C_M_S


def compute_galactic_acceleration_pdot(P_s: float, dist_kpc: float, l_deg: float, b_deg: float) -> dict:
    """
    Compute Galactic acceleration contributions to Ṗ.
    
    Returns dict with:
    - pdot_z: vertical acceleration contribution
    - pdot_R: differential rotation contribution
    - pdot_gal_total: total Galactic contribution
    """
    l = np.radians(l_deg)
    b = np.radians(b_deg)
    d = dist_kpc  # kpc
    
    # Pulsar position in Galactic coordinates
    # z = d * sin(b)
    z_kpc = d * np.sin(b)
    
    # R = sqrt(R0^2 + (d*cos(b))^2 - 2*R0*d*cos(b)*cos(l))
    d_proj = d * np.cos(b)
    R_kpc = np.sqrt(R0_KPC**2 + d_proj**2 - 2 * R0_KPC * d_proj * np.cos(l))
    
    # Vertical acceleration (toward Galactic plane)
    # a_z ≈ -4πGρ_disk * z (for |z| < scale height)
    # Using Kuijken & Gilmore (1989) approximation
    z_pc = z_kpc * 1000
    rho_msun_pc3 = RHO_DISK
    a_z_m_s2 = -4 * np.pi * G_SI * (rho_msun_pc3 * 1.989e30 / (3.086e16)**3) * (z_pc * 3.086e16)
    
    # Line-of-sight component of vertical acceleration
    a_z_los = a_z_m_s2 * np.sin(b)
    
    # Differential rotation (Oort constants)
    # a_R ≈ (A^2 - B^2) * R0 * sin(l) * cos(b) / R  (simplified)
    # This is a rough approximation
    if R_kpc > 0.1:
        # Convert Oort constants to SI
        A_si = OORT_A * 1e3 / (KPC_TO_M / 1e3)  # 1/s
        B_si = OORT_B * 1e3 / (KPC_TO_M / 1e3)  # 1/s
        a_R_los = (A_si**2 - B_si**2) * (R0_KPC * KPC_TO_M) * np.sin(l) * np.cos(b) / (R_kpc * KPC_TO_M)
    else:
        a_R_los = 0.0
    
    # Total Galactic acceleration (line-of-sight)
    a_gal_los = a_z_los + a_R_los
    
    # Convert to Ṗ contribution: Ṗ_gal = P * a_los / c
    pdot_z = P_s * a_z_los / C_M_S
    pdot_R = P_s * a_R_los / C_M_S
    pdot_gal_total = P_s * a_gal_los / C_M_S
    
    return {
        "z_kpc": z_kpc,
        "R_kpc": R_kpc,
        "a_z_los_m_s2": a_z_los,
        "a_R_los_m_s2": a_R_los,
        "pdot_z": pdot_z,
        "pdot_R": pdot_R,
        "pdot_gal_total": pdot_gal_total,
    }


def is_gc_associated(assoc: str) -> bool:
    """Check if pulsar is associated with a globular cluster."""
    if not assoc:
        return False
    a = assoc.lower()
    # Common GC indicators
    gc_keywords = ["gc", "globular", "ngc", "47 tuc", "terzan", "m15", "m28", "m13", "m5", "m62"]
    return any(kw in a for kw in gc_keywords)


def main():
    # Load ATNF catalog
    if not ATNF_DB_PATH.exists():
        raise FileNotFoundError(f"ATNF catalog not found at {ATNF_DB_PATH}")
    
    atnf_text = ATNF_DB_PATH.read_text()
    atnf_sha256 = hashlib.sha256(atnf_text.encode()).hexdigest()
    
    # Parse catalog
    rows = parse_atnf_extended(atnf_text)
    print(f"Parsed {len(rows)} pulsars from ATNF catalog")
    
    # Filter to field MSPs with required parameters
    field_msps = []
    for r in rows:
        # MSP cut
        if r["P_ms"] is None or r["P_ms"] >= MSP_PERIOD_CUT_MS:
            continue
        
        # Exclude GC-associated
        if is_gc_associated(r["assoc"]):
            continue
        
        # Require distance and proper motion for corrections
        if r["dist_kpc"] is None or r["dist_kpc"] <= 0:
            continue
        
        # Compute Galactic coordinates if not provided
        if r["gl_deg"] is None or r["gb_deg"] is None:
            if r["ra_deg"] is not None and r["dec_deg"] is not None:
                r["gl_deg"], r["gb_deg"] = equatorial_to_galactic(r["ra_deg"], r["dec_deg"])
        
        if r["gl_deg"] is None or r["gb_deg"] is None:
            continue
        
        field_msps.append(r)
    
    print(f"Field MSPs with distance and coordinates: {len(field_msps)}")
    
    # Compute corrections
    corrected = []
    for r in field_msps:
        P_s = r["P0_s"]
        P1_obs = r["P1_sps"]
        dist_kpc = r["dist_kpc"]
        l_deg = r["gl_deg"]
        b_deg = r["gb_deg"]
        
        # Shklovskii correction
        pmra = r.get("pmra_mas_yr") or 0.0
        pmdec = r.get("pmdec_mas_yr") or 0.0
        pm_total = np.sqrt(pmra**2 + pmdec**2)
        
        if pm_total > 0:
            pdot_shk = compute_shklovskii_pdot(P_s, pm_total, dist_kpc)
        else:
            pdot_shk = 0.0
        
        # Galactic acceleration correction
        gal = compute_galactic_acceleration_pdot(P_s, dist_kpc, l_deg, b_deg)
        pdot_gal = gal["pdot_gal_total"]
        
        # Corrected intrinsic Ṗ
        # P1_int = P1_obs - P1_shk - P1_gal
        P1_int = P1_obs - pdot_shk - pdot_gal
        
        # Only keep positive intrinsic Ṗ (physically meaningful)
        if P1_int <= 0:
            P1_int = None
            log_P1_int = None
        else:
            log_P1_int = np.log10(P1_int)
        
        corrected.append({
            **r,
            "pm_total_mas_yr": pm_total,
            "pdot_shk": pdot_shk,
            "pdot_gal": pdot_gal,
            "P1_int": P1_int,
            "log_P1_obs": np.log10(abs(P1_obs)) if P1_obs != 0 else None,
            "log_P1_int": log_P1_int,
            "z_kpc": gal["z_kpc"],
        })
    
    # Statistics
    with_pm = [r for r in corrected if r["pm_total_mas_yr"] > 0]
    with_int = [r for r in corrected if r["log_P1_int"] is not None]
    
    print(f"\nField MSPs with proper motion: {len(with_pm)}")
    print(f"Field MSPs with positive intrinsic Ṗ after correction: {len(with_int)}")
    
    # Compare observed vs corrected
    log_obs = np.array([r["log_P1_obs"] for r in with_int if r["log_P1_obs"] is not None])
    log_int = np.array([r["log_P1_int"] for r in with_int])
    
    mean_obs = np.mean(log_obs)
    mean_int = np.mean(log_int)
    
    print(f"\nMean log|Ṗ_obs|: {mean_obs:.3f}")
    print(f"Mean log|Ṗ_int|: {mean_int:.3f}")
    print(f"Mean correction: {mean_obs - mean_int:.3f} dex")
    
    # Shklovskii vs Galactic contributions
    shk_contributions = [r["pdot_shk"] for r in with_int]
    gal_contributions = [r["pdot_gal"] for r in with_int]
    
    mean_shk = np.mean(shk_contributions)
    mean_gal = np.mean(gal_contributions)
    
    print(f"\nMean Shklovskii contribution: {mean_shk:.2e} s/s")
    print(f"Mean Galactic contribution: {mean_gal:.2e} s/s")
    
    # Build output
    output = {
        "meta": {
            "timestamp_utc": datetime.now(timezone.utc).isoformat(),
            "source": {
                "atnf_psrcat": {
                    "path": str(ATNF_DB_PATH),
                    "sha256": atnf_sha256,
                }
            },
            "selection": {
                "msp_cut": f"P < {MSP_PERIOD_CUT_MS} ms",
                "gc_excluded": True,
                "requires_distance": True,
                "requires_coordinates": True,
            },
            "galactic_model": {
                "R0_kpc": R0_KPC,
                "V0_km_s": V0_KM_S,
                "rho_disk_msun_pc3": RHO_DISK,
            },
        },
        "counts": {
            "total_parsed": len(rows),
            "field_msps_with_params": len(field_msps),
            "with_proper_motion": len(with_pm),
            "with_positive_intrinsic_pdot": len(with_int),
        },
        "statistics": {
            "mean_log_P1_obs": float(mean_obs),
            "mean_log_P1_int": float(mean_int),
            "mean_correction_dex": float(mean_obs - mean_int),
            "std_log_P1_obs": float(np.std(log_obs)),
            "std_log_P1_int": float(np.std(log_int)),
        },
        "correction_breakdown": {
            "mean_shklovskii_sps": float(mean_shk),
            "mean_galactic_sps": float(mean_gal),
            "median_shklovskii_sps": float(np.median(shk_contributions)),
            "median_galactic_sps": float(np.median(gal_contributions)),
        },
    }
    
    # Write outputs
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    
    json_path = RESULTS_DIR / "step_5_12_galactic_correction.json"
    with open(json_path, 'w') as f:
        json.dump(output, f, indent=2)
    print(f"\nWrote: {json_path}")
    
    # Write markdown summary
    md_path = RESULTS_DIR / "step_5_12_galactic_correction.md"
    md_lines = [
        "# Galactic Potential Correction for Field MSPs",
        f"**Source:** {ATNF_DB_PATH}",
        f"**SHA256:** `{atnf_sha256[:16]}...`",
        "",
        "## Overview",
        "",
        "This analysis applies corrections to field pulsar Ṗ for:",
        "1. **Shklovskii effect:** Transverse velocity contribution (Ṗ_shk = P × μ² × d / c)",
        "2. **Galactic vertical acceleration:** z-component toward the Galactic plane",
        "3. **Galactic differential rotation:** Radial component from Oort constants",
        "",
        "## Sample Sizes",
        f"- **Total ATNF pulsars parsed:** {len(rows)}",
        f"- **Field MSPs with distance + coordinates:** {len(field_msps)}",
        f"- **With proper motion:** {len(with_pm)}",
        f"- **With positive intrinsic Ṗ after correction:** {len(with_int)}",
        "",
        "## Results",
        "",
        "| Metric | Value |",
        "| --- | --- |",
        f"| Mean log|Ṗ_obs| | {mean_obs:.3f} |",
        f"| Mean log|Ṗ_int| | {mean_int:.3f} |",
        f"| **Mean correction** | **{mean_obs - mean_int:.3f} dex** |",
        "",
        "## Correction Breakdown",
        "",
        "| Component | Mean | Median |",
        "| --- | --- | --- |",
        f"| Shklovskii | {mean_shk:.2e} s/s | {np.median(shk_contributions):.2e} s/s |",
        f"| Galactic | {mean_gal:.2e} s/s | {np.median(gal_contributions):.2e} s/s |",
        "",
        "## Interpretation",
        "",
        "The Shklovskii and Galactic corrections are typically small compared to the observed Ṗ for MSPs,",
        "but they provide a cleaner 'intrinsic' estimate for population comparisons.",
        "",
        "The corrected field sample can be compared to GC MSPs (where cluster acceleration dominates)",
        "to better isolate environment-dependent effects.",
    ]
    
    with open(md_path, 'w') as f:
        f.write('\n'.join(md_lines))
    print(f"Wrote: {md_path}")


if __name__ == "__main__":
    main()
