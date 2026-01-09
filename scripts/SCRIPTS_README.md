# TEP-COS Scripts Documentation

**Generated:** 2026-01-07  
**Purpose:** Map manuscript sections to analysis scripts

---

## Manuscript Structure → Script Mapping

### Section 3: Gravitational Lensing
| Script | Purpose | Output | Manuscript Reference |
|--------|---------|--------|---------------------|
| `step_1_0_data_acquisition.py` | Download COSMOGRAIL light curves | `data/cosmograil/*.rdb` | §3.1 Data |
| `step_3_0_cosmograil_temporal_shear.py` | Core temporal shear analysis | `step_3_0_*.json` | §3.2 Results |
| `step_3_2_cosmograil_validation.py` | Injection-recovery, achromaticity | `validation_*.json` | §3.4 Validation |
| `step_3_5_advanced_lensing_analysis.py` | Redshift correlation | `step_3_5_*.json` | §3.6.1 Geometric Fingerprint |
| `step_3_6_high_z_predictions.py` | High-z predictions | `step_3_6_*.json` | §3.8 High-z Predictions |
| `step_3_7_chromaticity_simulation.py` | Required precision | `step_3_7_*.json` | §3.6.3 Chromaticity |
| `step_3_12_q2237_multiband_chromaticity.py` | Q2237 multi-band | `step_3_12_*.json` | §3.6.3 |
| `step_3_13_he0435_multiband_chromaticity.py` | HE0435 multi-band | `step_3_13_*.json` | §3.6.3 |
| `step_3_14_he1104_multiband_chromaticity.py` | HE1104 multi-band | `step_3_14_*.json` | §3.6.3 |
| `step_3_15_q2237_vakulik_multiband_chromaticity.py` | Q2237 Vakulik | `step_3_15_*.json` | §3.6.3 |

### Section 4: Pulsar Timing
| Script | Purpose | Output | Manuscript Reference |
|--------|---------|--------|---------------------|
| `step_5_9_freire_gcpsr_radial_analysis.py` | Radial analysis | `freire_gcpsr_radial_*.json` | §4.9 Radial Correlation |
| `step_5_10_pulsar_population_controls.py` | Population controls | `step_5_10_*.json` | §4.3 Results |
| `step_5_11_binary_pulsar_analysis.py` | Binary vs isolated (GC) | `step_5_11_*.json` | §4.6 Binary vs Isolated |
| `step_5_12_field_binary_analysis.py` | Field binary control | `step_5_12_*.json` | §4.7 Field Control |
| `step_5_13_cluster_acceleration_simulation.py` | Monte Carlo sim | `step_5_13_*.json` | §4.5 Simulation |

### Section 5: Galaxy Kinematics
| Script | Purpose | Output | Manuscript Reference |
|--------|---------|--------|---------------------|
| `step_2_0_cosmic_coriolis_analysis.py` | CMB dipole search | Log files | §5.2-5.4 |

### Section 5B: Stellar Archaeology
| Script | Purpose | Output | Manuscript Reference |
|--------|---------|--------|---------------------|
| `step_6_3_apogee_stellar_archaeology.py` | APOGEE ages | `apogee_*.json` | §5B.3 Results |
| `step_6_6_sdss_twin_galaxy_matched_pairs.py` | Twin galaxy pairs | `sdss_twin_*.json` | §5B.4 Matched Pairs |
| `step_6_6_sfr_holonomy.py` | SFR holonomy | `sdss_sfr_holonomy_*.json` | §6.1 SFR Holonomy |
| `step_6_6b_sfr_validation.py` | SFR validation | `sdss_sfr_holonomy_validation.json` | §6.1 |

### Section 6: Discussion (Bulletproof Tests)
| Script | Purpose | Output | Manuscript Reference |
|--------|---------|--------|---------------------|
| `step_6_12_sdss_test_h_chemical_clock.py` | Chemical Clock (H) | `sdss_test_h_*.json` | §5B.4.1, §6.1 |
| `step_6_94_sdss_test_dx_halpha_uv.py` | Timescale Ratios (DX) | `sdss_test_dx_*.json` | §5.9.1, §6.1 |
| `step_6_88_sdss_test_dq_satellite_abundance.py` | Satellite Abundance (DQ) | `sdss_test_dq_*.json` | §5.9.2, §6.1 |
| `step_6_99_sdss_test_l_radial_gradient.py` | LW-MW Age (L) | `sdss_test_l_*.json` | §6.1 |
| `step_6_101_sdss_test_m_mass_discrepancy.py` | Mass Discrepancy (M) | `sdss_test_m_*.json` | §6.1 |
| `step_6_96_sdss_test_g_sn_stretch.py` | SN Ia Stretch (G) ❌ | `sdss_test_g_*.json` | §6.1 (Contradiction) |
| `step_6_98_sdss_test_k_size_age.py` | Size-Age (K) ❌ | `sdss_test_k_*.json` | §6.1 (Contradiction) |
| `step_6_91_sdss_test_dt_red_clump.py` | Red Clump (DT) ❌ | `sdss_test_dt_*.json` | §5B.5 (Contradiction) |
| `step_6_95_sdss_test_dy_phase_spirals.py` | Phase Spirals (DY) | `sdss_test_dy_*.json` | §5.9.3 (Null) |

### Bulletproof Pipeline
| Script | Purpose | Output |
|--------|---------|--------|
| `step_7_00_bulletproof_tep_signals.py` | Test H, I, L | `bulletproof_tep_signals.json` |
| `step_7_01_bulletproof_timescale_ratios.py` | Test DX | `bulletproof_test_dx_timescale.json` |
| `step_7_02_bulletproof_additional_signals.py` | Additional | `bulletproof_additional_signals.json` |
| `step_7_03_bulletproof_sfr_holonomy.py` | SFR | `bulletproof_sfr_holonomy.json` |

### Utility Scripts
| Script | Purpose |
|--------|---------|
| `convert_cds_to_rdb.py` | Convert CDS format to RDB |
| `__init__.py` | Package init |

---

## Figure Scripts (scripts/figures/)
| Script | Purpose | Output |
|--------|---------|--------|
| `cosmograil_temporal_shear_figure.py` | Main temporal shear fig | `figures/manuscript/` |
| `cosmograil_comprehensive_figure.py` | Comprehensive fig | `figures/manuscript/` |
| `cosmograil_microlensing_figure.py` | Microlensing comparison | `figures/manuscript/` |

---

## Archived Scripts (scripts/archive/)
Scripts moved to archive because they were:
- Debug/diagnostic scripts (not used in final analysis)
- Experimental parameter sweeps (superseded by final analysis)
- SDSS tests that showed nulls/contradictions not featured in manuscript

---

## Key Outputs to Keep (results/outputs/)
### Lensing
- `step_3_0_cosmograil_temporal_shear_v3_expanded.json` - Original results
- `step_3_0_cosmograil_temporal_shear_opB_modejump.json` - Best operating point
- `validation_DESJ0408_full.json` - Definitive validation
- `step_3_2_validation_results.json` - Validation suite
- `step_3_12_*.json` through `step_3_15_*.json` - Chromaticity

### Pulsars
- `step_5_10_pulsar_population_controls.json`
- `step_5_11_binary_pulsar_analysis.json`
- `step_5_12_field_binary_analysis.json`
- `step_5_13_acceleration_sim.json`
- `freire_gcpsr_radial_summary.json`

### SDSS/Galaxy
- `bulletproof_tep_signals.json`
- `bulletproof_test_dx_timescale.json`
- `sdss_sfr_holonomy_results.json`
- `sdss_test_h_results.json`, `sdss_test_dx_results.json`, etc.

### Synthesis
- `TEP_EVIDENCE_SYNTHESIS_FINAL.md`
- `BULLETPROOF_TEP_SUMMARY.md`
