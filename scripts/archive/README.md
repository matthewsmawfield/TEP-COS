# Archived Scripts

These scripts were part of the exploratory analysis phase but are not used in the final published manuscript. They are preserved here for reference and reproducibility of the full research process.

## Categories

### Galaxy Kinematics Exploration (step_2_*)
Extensive exploration of MaNGA velocity fields, including:
- Grid scans, tomography, stratification
- Various normalization and selection criteria
- Ultimately yielded expected null result (signal below noise floor)

### Topology/Pattern Search (step_3_0_topology, step_3_1_pattern, etc.)
Early attempts to find spatial patterns in residuals. Superseded by focused temporal shear analysis.

### SDSS Analysis (step_4_*)
Alternative approaches using SDSS data. Not included in final paper.

### Early Pulsar Analysis (step_5_0 through step_5_8)
Exploratory pulsar analyses before settling on the final population control methodology (step_5_10, step_5_11).

## Final Analysis Scripts
The scripts used in the published paper are in `scripts/steps/`:
- `step_1_0_data_acquisition.py` - Data download
- `step_2_0_cosmic_coriolis_analysis.py` - Galaxy kinematics (null result)
- `step_3_0_cosmograil_temporal_shear.py` - Core lensing analysis
- `step_3_2_cosmograil_validation.py` - Validation suite
- `step_3_5_advanced_lensing_analysis.py` - Error budget, correlations
- `step_5_9_freire_gcpsr_radial_analysis.py` - Pulsar radial analysis
- `step_5_10_pulsar_population_controls.py` - Population controls
- `step_5_11_binary_pulsar_analysis.py` - Binary vs isolated MSPs
