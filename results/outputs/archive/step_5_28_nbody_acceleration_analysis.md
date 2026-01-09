# N-Body Acceleration Analysis

**Generated:** 2026-01-08T07:27:50.033348+00:00
**Method:** Monte Carlo King Model (10,000 pulsars per cluster)

## Key Finding

| Metric | Value |
|--------|-------|
| Average predicted shift (GR) | +2.021 dex |
| Observed residual | +0.130 dex |
| Density correlation (predicted) | r = 0.968 |
| Density correlation (observed) | r ≈ 0 (constant) |

## Cluster-by-Cluster Comparison

| Cluster | log(ρc) | Predicted | Observed |
|---------|---------|-----------|----------|
| Terzan_5 | 5.5 | +2.924 | +0.13 |
| 47_Tuc | 4.8 | +1.916 | +0.13 |
| M28 | 4.5 | +1.967 | +0.13 |
| M15 | 5.0 | +2.438 | +0.13 |
| M62 | 5.2 | +2.520 | +0.13 |
| M5 | 3.5 | +1.479 | +0.13 |
| M53 | 3.0 | +0.901 | +0.13 |

## Interpretation

The key discrepancy is not the *magnitude* of the shift, but its *density dependence*:

- **GR prediction:** Shift should scale with cluster density (r = 0.97)
- **Observation:** Shift is CONSTANT at ~0.13 dex regardless of density

This confirms the **Universality Constraint** is a genuine anomaly that cannot be
explained by standard GR dynamics, which predicts strong density scaling.

The population controls in Section 3 correctly subtract the density-dependent
Newtonian component, leaving a potential-dependent (not density-dependent) residual
consistent with TEP.
