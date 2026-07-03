#!/usr/bin/env python3
"""
TEP Universal Model Utilities (v0.8 Jakarta/Kos)

Shared functions for computing TEP quantities across all research papers:
COS (Paper 10), H0 (Paper 11), JWST (Paper 12), WB (Paper 13).

This module provides a unified point of truth for:
1. Universal Couplings (KAPPA_GAL, ALPHA_INT)
2. Chronological Enhancement (Gamma_t)
3. Screening Mechanisms (Temporal Topology)
4. Kinematic Profile Models (Wide Binaries)

Author: Matthew L. Smawfield
Date: April 2026
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import integrate

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))
from core import constants as tep_const

# =============================================================================
# 1. UNIVERSAL COUPLINGS & CONSTANTS
# =============================================================================

# CANONICAL OBSERVABLE RESPONSE COEFFICIENT (Paper 11)
# Cross-paper canonical theoretical prior, derived from geometric-factor
# and virial estimates (Appendix C of Paper 11).  Consistent with
# TEP-H0 (Paper 11) and TEP-JWST (Paper 12).
# Units: Magnitudes [mag]
KAPPA_GAL = tep_const.KAPPA_GAL
KAPPA_GAL_UNCERTAINTY = tep_const.KAPPA_GAL_UNCERTAINTY

# STELLAR EVOLUTION INDEX
# M/L ~ t^n from stellar isochrones.
ALPHA_NUCLEAR = tep_const.ALPHA_NUCLEAR

# POTENTIAL PARAMETERS
LOG_MH_REF = tep_const.LOG_MH_REF
PHI_REF_0 = tep_const.PHI_REF_0    # Dimensionless Phi/c^2 for 10^12 Msun halo at z=0
Z_REF = tep_const.Z_REF

# SCREENING SCALES (from core.constants)
RHO_CRIT_G_CM3 = tep_const.RHO_C

# PHYSICAL CONSTANTS
C_LIGHT_KM_S = 2.99792458e5
G_NEWTON_PC_MSUN = 4.30091e-3  # (pc/Msun) * (km/s)^2
G_AU = 887.1                   # (km/s)^2 * AU / M_sun

# =============================================================================
# 2. CHRONOLOGICAL ENHANCEMENT (Gamma_t)
# =============================================================================

def get_phi_from_log_mh(log_Mh):
    """Compute dimensionless virial potential Phi/c^2 at z=0."""
    return 1.6e-7 * (10**log_Mh / 1e12)**(2/3)

def compute_gamma_t(log_Mh, z, kappa=KAPPA_GAL, n=ALPHA_NUCLEAR):
    """
    Compute TEP chronological enhancement factor (Potential-Linear Form).
    
    Gamma_t = exp[ K * (Phi - Phi_ref) * sqrt(1+z) ]
    where K = kappa * ln(10) / (2.5 * n)
    """
    phi = get_phi_from_log_mh(log_Mh)
    k_exp = (kappa * np.log(10)) / (2.5 * n)
    z_scaling = np.sqrt(1 + z)
    argument = k_exp * (phi - PHI_REF_0) * z_scaling
    return np.exp(argument)

# =============================================================================
# 3. KINEMATIC & SCREENING MODELS
# =============================================================================

def tep_screening_model(s, r_s, alpha):
    """Canonical TEP exponential screening model for velocities."""
    return 1.0 + alpha * (1.0 - np.exp(-s / r_s))

def temporal_topology_suppression(rho, rho_c=RHO_CRIT_G_CM3, kappa_bare=KAPPA_GAL):
    """Continuous Temporal Topology suppression factor."""
    x = np.log10(rho / rho_c) / 0.5
    suppression = 1.0 / (1.0 + np.exp(x))
    return kappa_bare * suppression

# =============================================================================
# 4. MASS CORRECTIONS & BIAS
# =============================================================================

def isochrony_mass_bias(gamma_t, n=ALPHA_NUCLEAR):
    """M/L ratio bias: M_obs / M_true = Gamma_t^n."""
    return np.power(np.maximum(gamma_t, 0.01), n)

def correct_stellar_mass(log_Mstar, gamma_t, n=ALPHA_NUCLEAR):
    """Apply TEP correction to observed stellar mass."""
    return log_Mstar - np.log10(isochrony_mass_bias(gamma_t, n))

def stellar_to_halo_mass_behroozi_like(log_Mstar, z):
    """Empirical SMHM relation proxy for high-z."""
    log_ratio = -1.8 - 0.1 * (log_Mstar - 10) + 0.05 * (z - 5)
    return log_Mstar - log_ratio

# =============================================================================
# 5. COSMOLOGY UTILS
# =============================================================================

def cosmic_time_gyr(z, H0=67.4, Om=0.315):
    """Cosmic time at redshift z in Gyr (flat LCDM)."""
    z = np.atleast_1d(z)
    H0_s = H0 * 1e3 / 3.0857e22
    def integrand(zp):
        return 1.0 / ((1 + zp) * np.sqrt(Om * (1 + zp)**3 + (1-Om)))
    results = [integrate.quad(integrand, zi, np.inf)[0] / H0_s / 3.156e16 for zi in z]
    return np.array(results) if len(results) > 1 else results[0]
