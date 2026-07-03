# Step 55: Cluster Bootstrap

**Generated:** 2026-07-03T10:05:56.938497+00:00

Cluster-level bootstrap (resampling clusters with replacement) to account for
within-cluster dependence.

## Results

| Test | Mean | Std | 95% CI Lower | 95% CI Upper | p-value |
|------|------|-----|--------------|--------------|---------|
| Raw diff (dex) | +0.577 | 0.113 | +0.336 | +0.776 | 0.0000 |
| Period-only resid (dex) | +0.803 | 0.104 | +0.580 | +0.980 | 0.0000 |
| Period+B resid (dex) | +0.224 | 0.038 | +0.144 | +0.293 | 0.0000 |
| Gamma raw | 0.324 | 0.118 | 0.115 | 0.583 | 0.0008 |
| Gamma period-only | 0.398 | 0.092 | 0.232 | 0.591 | 0.0000 |

## Interpretation

- If cluster-bootstrap p-values are **larger** than pulsar-bootstrap p-values,
  the original analysis was inflating significance by treating within-cluster
  pulsars as independent.
- The 95% CI should be used for robust inference.
