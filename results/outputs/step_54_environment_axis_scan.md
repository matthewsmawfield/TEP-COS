# Step 54: Environment Axis Scan

**Generated:** 2026-09-23T16:03:53.344855+00:00

## Results

| Variable | Gamma | Std Err | R^2 | p-value |
|----------|-------|---------|-----|---------|
| Central luminosity density (current) | +0.320 | 0.107 | 0.222 | 0.0056 |
| Mass density M/Rc^3 | +0.266 | 0.085 | 0.240 | 0.0038 |
| Potential depth M/Rc | +0.564 | 0.213 | 0.185 | 0.0125 |
| Acceleration scale M/Rc^2 | +0.398 | 0.125 | 0.245 | 0.0034 |
| Escape velocity v_esc | +1.052 | 0.537 | 0.110 | 0.0592 |
| Velocity dispersion sigma_v | +1.129 | 0.426 | 0.185 | 0.0125 |
| Relaxation time proxy Rc/sigma_v | -0.531 | 0.170 | 0.240 | 0.0038 |

## Interpretation

- If TEP is **density-driven**, the strongest correlation should be with
  **central luminosity density** or **mass density**.
- If TEP is **potential-driven** (coherence / potential depth), the signal
  should strengthen against **M/Rc** or **escape velocity**.
- If TEP is **dynamical-state-driven**, look for correlation with
  **velocity dispersion** or **relaxation time**.
- The axis with the highest R^2 and steepest positive slope is the preferred
  environmental variable for TEP.

## Channel Decomposition

Bivariate decomposition `excess ~ a·log(ρc) + b·log(Rc)` measured on the
cluster-mean log|Pdot| and on the per-cluster controlled residuals
(step_07). Channel predictions: acceleration a ∝ ρc·Rc → (a,b) = (1,1);
potential depth |Φ| ∝ ρc·Rc² → (a,b) = (1,2); ρ² variance bias → (a,b) = (2,0).

| Estimand | Fit | a (ρc) | b (Rc) | R² |
|----------|-----|--------|--------|----|
| cluster_mean_logPdot | ols | +0.282 ± 0.238 | -0.095 ± 0.525 | 0.223 |
| cluster_mean_logPdot | wls_n | +0.230 ± 0.135 | -0.469 ± 0.319 | 0.495 |
| controlled_residuals | wls_1_over_sem2 | +0.070 ± 0.025 | +0.030 ± 0.076 | 0.469 |
| controlled_residuals | ols | +0.093 ± 0.070 | -0.112 ± 0.154 | 0.445 |

| Estimand | Fit | accel (1,1) F/p | potential (1,2) F/p | ρ² (2,0) F/p |
|----------|-----|-----------------|---------------------|--------------|
| cluster_mean_logPdot | ols | 5.4 / 9.7e-03 | 8.6 / 1.1e-03 | 118.6 / 5.7e-15 |
| cluster_mean_logPdot | wls_n | 16.4 / 1.5e-05 | 30.6 / 5.6e-08 | 228.2 / 7.1e-19 |
| controlled_residuals | wls_1_over_sem2 | 1197.5 / 3.6e-23 | 710.3 / 1.0e-20 | 9614.3 / 4.3e-33 |
| controlled_residuals | ols | 121.2 / 1.3e-12 | 96.0 / 1.4e-11 | 1558.9 / 2.0e-24 |


The fitted Rc exponent is consistent with zero across estimands and weightings, while the unsuppressed channel predictions (b = +1 acceleration, b = +2 potential depth) are rejected by the joint tests. The measured excess is a saturated density response rather than a linear response to any unsuppressed field variable.
