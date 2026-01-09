# Advanced Lensing Analysis
**Generated:** 2026-01-03T22:51:37.689772+00:00

## A. Multi-Band Chromaticity Test
**Status:** NOT_FEASIBLE_WITH_CURRENT_DATA
**Reason:** COSMOGRAIL monitoring data is single-band (R-band)

### Falsifier Protocol
- **step_1:** Obtain multi-band light curves (e.g., g/r or V/I) for detection systems
- **step_2:** Compute Γ_band for each band using identical estimator settings
- **step_3:** Test ΔΓ = Γ_blue - Γ_red against zero
- **step_4:** TEP: ΔΓ ≈ 0; Microlensing: ΔΓ ≠ 0

- **TEP prediction:** Temporal shear should be achromatic (same Γ in all bands)
- **Microlensing prediction:** Microlensing produces chromatic effects (different Γ per band)

## C. Redshift Scaling Analysis
**Correlation (|Γ| vs geometric factor):** r = 0.545, p = 0.2055

### High-z Predictions
| z_source | z_lens (typical) | Geometric Factor | Predicted |Γ| | Exceeds 300? |
|----------|------------------|------------------|------------|--------------|
| 2.5 | 0.55 | 2.26 | 131 | ✗ |
| 3.0 | 0.60 | 2.50 | 165 | ✗ |
| 3.5 | 0.65 | 2.73 | 197 | ✗ |
| 4.0 | 0.70 | 2.94 | 227 | ✗ |

**Falsifiable prediction:** Systems with z_S > 2.5 should show |Γ| > 300 days/log(τ)

## E. Residual Cross-Correlation
**Pairs tested:** 39
**Mean correlation:** 0.284 ± 0.724
**t-test p-value:** 0.0204
**Significant at 0.05:** 18

- **TEP expectation:** Coherent residuals across pairs (positive correlation) if shared gravitational path
- **Noise expectation:** Independent residuals (correlation ~ 0)

## G. Error Budget Analysis
**Total pairs measured:** 23
**Detections (>3σ):** 5
**Combined p-value:** 0.00e+00

### Systematic Budget
| Source | Status | Notes |
|--------|--------|-------|
| statistical_bootstrap | QUANTIFIED | Bootstrap resampling uncertainty on Γ |
| filtering_sensitivity | PARTIALLY_QUANTIFIED | Sensitivity to Gaussian smoothing scale choice |
| microlensing_contamination | NOT_YET_TESTED | Stellar microlensing in lens galaxy |
| intrinsic_variability | MITIGATED | Non-stationary quasar variability structure |
| geometric_model | SECONDARY | Lens model uncertainties affecting path length |

**Dominant uncertainty:** Statistical (bootstrap)
**Key systematic:** Microlensing (requires multi-band test)
**Validation:** Null controls (HE0435, WFI2033) confirm no false positives