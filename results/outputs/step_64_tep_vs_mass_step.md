# SN Ia TEP vs. Mass Step Discrimination Analysis

**Analysis Date:** 2026-07-03T11:13:07.044891

## Summary

**Overall Verdict:** MASS_STEP_DOMINATED

**Assessment:** SN Ia mB-σ correlation is primarily the standard mass step effect

**TEP Contribution:** None detected above mass step baseline

## Individual Tests

### Partial Correlation Test

- **r_raw:** 0.21780141521973087
- **r_partial:** -0.04659393710077179
- **p_partial:** 0.4940104953811617
- **r_mass_mB:** 0.564285378825223
- **discrimination:** MASS_STEP_DOMINATED
- **interpretation:** Correlation fully explained by host mass; no residual TEP signal

### Collinearity Assessment

- **r_sigma_mass:** 0.44695869899755425
- **vif:** 1.2496439742598306
- **collinearity_level:** LOW
- **note:** σ and mass are separable; good discrimination potential

### Functional Form Test

- **unscreened_r:** 0.22217793042236308
- **screened_r:** -0.12755670065403865
- **r_difference:** 0.34973463107640174
- **verdict:** TEP_CONSISTENT
- **note:** Correlation present in unscreened, absent in screened - TEP signature

### Observable Comparison Test

- **mB_correlation_significant:** True
- **x1_correlation_significant:** True
- **x1_r:** -0.2645486552709711
- **x1_p:** 7.677889629587774e-05
- **verdict:** MASS_STEP_CONSISTENT
- **note:** Both observables correlate - suggests common driver (mass/metallicity)

## Recommendations

1. SN Ia channel cannot distinguish TEP from mass step. Consider this exploratory only; do not present as primary evidence.
2. Future work: Use independent mass estimates (e.g., SED fitting) to break σ-M degeneracy.
