# Step 53: CMC Parser Unit Tests

**Generated:** 2026-10-16T13:00:03.326075+00:00

## Tests

### deduplication

- n_before: 1439522
- n_after: 1213
- mean_logP_before: 1.1824967029207425
- mean_logP_after: 1.3703204380878775
- delta_mean_logP: 0.18782373516713502

### fixed_period

- fixed_period_s: 15.2061
- mean_pdot_fixed: 0.7943344608179587
- mean_pdot_variable: 0.062085938382168
- ratio_variable_to_fixed: 0.07816095290419048

### acceleration_conversion

- M_enc_analytic_msun: 89442.71909999158
- a_analytic_kms2: 1539.1303102726552
- P_s: 0.005
- pdot: 2.565217183787759e-05
- pdot_over_p: 0.005130434367575518
- a_from_pdot_over_p_ms2: 1539130.3102726554
- a_from_pdot_ms2: 1539130.3102726554
- consistency_pdot_over_p: 1.5127415924427084e-16
- consistency_pdot: 1.5127415924427084e-16

### projected_radius

- n: 1439522
- mean_r_3d: 6.438071406490977
- mean_r_proj: 5.056494236009217
- ratio_proj_to_3d: 0.785405118512847

### period_distribution_sensitivity

- n_gc_periods: 198
- mean_pdot_observed_dist: 0.00028695982015634766
- mean_pdot_narrow_dist: 0.0002286589825989367
- ratio_narrow_to_observed: 0.7968327498754137

## Interpretation

- **Deduplication** should not remove dynamically important states.
- **Fixed-period assumption** may under-predict amplitude if the real period
  distribution is broader.
- **Acceleration conversion** must be internally consistent (pdot/P == a/c).
- **Projected radius** sampling changes the mean radius and therefore the
  predicted acceleration.
- **Period-distribution sensitivity** quantifies how much the predicted
  amplitude depends on the assumed period draw.
