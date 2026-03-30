# TEP-COS Scripts Documentation

**Generated:** 2026-03-30  
**Purpose:** Map analysis pipeline to manuscript sections

---

## Pipeline Execution (run_pipeline.py)

The master pipeline executes scripts in the following phases:

| Phase | Scripts | Purpose |
|-------|---------|---------|
| **Data Acquisition** | utils/data_acquisition.py | Download Freire GCpsr, ATNF psrcat |
| **Population Controls** | step_5_10 | Match GC and Field MSPs by period, magnetic field |
| **Core Pulsar Analysis** | step_5_27, step_5_31, step_5_32, step_5_33, step_5_35 | Density scaling, hierarchical models, covariance validation |
| **Binary Analysis** | step_5_11, step_5_12, step_5_36 | Binary vs isolated MSP comparison |
| **Validation Suite** | step_5_33b, step_5_34, step_5_37-39, step_5_43 | Sensitivity, power analysis, Monte Carlo |
| **N-Body Pushback** | step_5_41-45, step_5_47-49 | Dynamical tests, core collapse, systematic ceilings |
| **Figure Generation** | step_5_13, step_5_32, step_5_40 | Publication figures |
| **Appendix** | step_7_0, step_7_1, step_7_2 | SN Ia tests, robustness, audit |

---

## Section 3: Pulsar Timing (Core Analysis)

| Script | Purpose | Output | Manuscript Ref |
|--------|---------|--------|----------------|
| `step_5_10_pulsar_population_controls.py` | Population controls, data ingestion | `step_5_10_*.json` | §4.3 |
| `step_5_27_hybrid_maximum_analysis.py` | Hybrid maximum Pdot analysis | `step_5_27_*.json` | §4.4 |
| `step_5_31_per_cluster_controlled_residuals.py` | Per-cluster controlled residuals | `step_5_31_*.json` | §4.4 |
| `step_5_32_full_density_scaling.py` | Full density scaling simulation | `step_5_32_*.json` | §4.5 |
| `step_5_33_hierarchical_density_scaling.py` | Hierarchical mixed-effects models | `step_5_33_*.json` | §4.5 |
| `step_5_35_covariance_validation.py` | Covariance-aware statistical validation | `step_5_35_*.json` | §4.4 |

### Binary Pulsar Analysis

| Script | Purpose | Output | Manuscript Ref |
|--------|---------|--------|----------------|
| `step_5_11_binary_pulsar_analysis.py` | GC binary vs isolated MSPs | `step_5_11_*.json` | §4.6 |
| `step_5_12_field_binary_analysis.py` | Field binary control sample | `step_5_12_*.json` | §4.7 |
| `step_5_36_integrated_binary_control.py` | Integrated binary control test | `step_5_36_*.json` | §4.6 |

### Validation & Sensitivity

| Script | Purpose | Output | Manuscript Ref |
|--------|---------|--------|----------------|
| `step_5_33b_outlier_exclusion_sensitivity.py` | Outlier exclusion sensitivity | `step_5_33b_*.json` | §4.4 |
| `step_5_34_shklovskii_sensitivity.py` | Shklovskii effect sensitivity | `step_5_34_*.json` | §4.4.2 |
| `step_5_37_rho_sensitivity.py` | Rho_intra sensitivity analysis | `step_5_37_*.json` | §4.4 |
| `step_5_38_power_analysis.py` | Statistical power validation | `step_5_38_*.json` | §4.4 |
| `step_5_39_monte_carlo_validation.py` | Monte Carlo Type I/II error validation | `step_5_39_*.json` | §4.4 |
| `step_5_43_sensitivity_cmc_report.py` | Sensitivity & CMC comparison report | `step_5_43_*.json` | §4.5 |

### N-Body Pushback & CMC Comparison

| Script | Purpose | Output | Manuscript Ref |
|--------|---------|--------|----------------|
| `step_5_41_pulsar_dynamical_calibration.py` | Pulsar dynamical calibration | `step_5_41_*.json` | §4.8 |
| `step_5_41b_sensitivity_analysis.py` | Dynamical calibration sensitivity | `step_5_41b_*.json` | §4.8 |
| `step_5_42_cmc_real_comparison.py` | CMC vs real cluster comparison | `step_5_42_*.json` | §4.5 |
| `step_5_44_theoretical_uncertainty.py` | Theoretical uncertainty quantification | `step_5_44_*.json` | §4.8 |
| `step_5_45_bayesian_posterior.py` | Bayesian posterior analysis | `step_5_45_*.json` | §4.4 |
| `step_5_47_core_collapse_test.py` | Core collapse test | `step_5_47_*.json` | §4.9 |
| `step_5_48_cmc_literature_comparison.py` | CMC literature comparison | `step_5_48_*.json` | §4.5 |
| `step_5_49_systematic_ceiling.py` | Systematic ceiling analysis | `step_5_49_*.json` | §4.8 |

### Legacy/Comparison Scripts

| Script | Purpose | Output |
|--------|---------|--------|
| `step_5_9_freire_gcpsr_radial_analysis.py` | Freire GCpsr radial analysis (legacy) | `freire_gcpsr_radial_*.json` |
| `step_5_46_spatial_gradient.py` | **REMOVED** — Spatial gradient analysis deprecated due to insufficient MSPs with radial positions (N=14) and Pdot measurements | — |

---

## Section 5: Galaxy Kinematics & Stellar Archaeology

| Script | Purpose | Output | Manuscript Ref |
|--------|---------|--------|----------------|
| `step_6_5_manga_spatially_resolved.py` | MaNGA spatially resolved analysis | `step_6_5_*.json` | §5.2 |
| `step_6_10_manga_test_e_age_discrepancy.py` | MaNGA age discrepancy test | `step_6_10_*.json` | §5.2 |
| `step_6_57_sdss_test_cc_manganese_clock.py` | Manganese clock test | `step_6_57_*.json` | §5B.3 |

---

## Section 6: Appendix - Supernova Ia Tests

| Script | Purpose | Output | Manuscript Ref |
|--------|---------|--------|----------------|
| `step_7_0_sn_ia_stretch_test.py` | SN Ia σ-mB correlation test | `step_7_0_*.json` | §6.1 |
| `step_7_1_sn_ia_robustness.py` | SN Ia robustness validation | `step_7_1_*.json` | §6.1 |
| `step_7_1_tep_vs_mass_step.py` | TEP vs mass-step analysis | `step_7_1b_*.json` | §6.1 |
| `step_7_2_sn_ia_audit.py` | SN Ia deep audit | `step_7_2_*.json` | §6.1 |

---

## Figure Generation Scripts

| Script | Purpose | Output |
|--------|---------|--------|
| `step_5_11_binary_spatial_figure.py` | Binary spatial figure | `figures/manuscript/` |
| `step_5_13_cluster_acceleration_figure.py` | Cluster acceleration figure | `figures/manuscript/` |
| `step_5_32_density_scaling_figure.py` | Density scaling figure | `figures/manuscript/` |
| `step_5_40_tep_summary_figure.py` | TEP summary figure | `figures/manuscript/` |

---

## Utility Scripts

| Script | Purpose |
|--------|---------|
| `convert_cds_to_rdb.py` | Convert CDS format to RDB |
| `__init__.py` | Package initialization |

### Utils (scripts/utils/)

| Script | Purpose |
|--------|---------|
| `check_manga_kinematics.py` | MaNGA kinematics utilities |
| `data_acquisition.py` | Data download and health checking |

---

## Usage

### Run Full Pipeline
```bash
python scripts/run_pipeline.py
```

### Run with Options
```bash
python scripts/run_pipeline.py --skip-validation    # Skip long validation steps
python scripts/run_pipeline.py --skip-figures       # Skip figure generation
python scripts/run_pipeline.py --only-core          # Fast mode: core analysis only
python scripts/run_pipeline.py --parallel           # Enable parallel processing
```

### Run Individual Steps
```bash
python scripts/steps/step_5_32_full_density_scaling.py
```

---

## Key Outputs (results/outputs/)

### Primary Results
- `step_5_10_pulsar_population_controls.json` - Population controls
- `step_5_32_full_density_scaling.json` - Density scaling
- `step_5_33_hierarchical_density_scaling.json` - Hierarchical models
- `step_5_35_covariance_validation.json` - Statistical validation

### Validation
- `step_5_38_power_analysis.json` - Power analysis
- `step_5_39_monte_carlo_validation.json` - Monte Carlo validation

### Figures
- `figures/manuscript/` - All publication figures

---

## Archive (scripts/archive/)

Scripts moved to archive include:
- Debug/diagnostic scripts (not used in final analysis)
- Experimental parameter sweeps (superseded)
- Deprecated SDSS tests (nulls/contradictions not featured)
- Old step versions (replaced by newer implementations)

