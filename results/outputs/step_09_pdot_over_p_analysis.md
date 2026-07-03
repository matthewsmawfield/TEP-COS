# Step 09: Pdot/P Parallel Observable Analysis

**Timestamp:** 2026-07-03T10:04:43.113448+00:00

## Base Sample (step_02)

Counts: 198 GC MSPs, 202 field MSPs

### Base Comparison (log|Pdot/P|)

| Statistic | GC | Field | Difference |
|-----------|-----|-------|------------|
| Mean (dex) | -16.827 | -17.566 | 0.739 |
| Std (dex) | 0.731 | 0.557 | — |
| t-test | — | — | t = 11.36, p = 8.01e-26 |
| Mann-Whitney | — | — | p = 2.28e-26 |

### Controls on log|Pdot/P|

| Control | Diff (dex) | 16th–84th CI | p (two-sided) |
|---------|------------|--------------|---------------|
| Period-matched | 0.747 | [0.695, 0.800] | 0.0000 |
| Period+B-proxy | 0.746 | [0.697, 0.794] | 0.0000 |

## Hybrid Sample (step_06)

Counts: 199 GC MSPs, 351 field MSPs

### Base Comparison (log|Pdot/P|)

| Statistic | GC | Field | Difference |
|-----------|-----|-------|------------|
| Mean (dex) | -16.822 | -17.475 | 0.654 |
| Std (dex) | 0.733 | 0.664 | — |
| t-test | — | — | t = 10.39, p = 2.07e-22 |
| Mann-Whitney | — | — | p = 3.58e-30 |

### Hybrid Period-Matched Control

- Mean difference: 0.693 dex
- 16th–84th percentile CI: [0.636, 0.750]
- Two-sided p: 0.0000

## Manuscript Summary

| Observable | Raw diff (dex) | Controlled residual (dex) | p-value |
|------------|----------------|---------------------------|---------|
| log|Pdot| (existing) | 0.59 | 0.606 (period-matched) | < 10⁻¹³ |
| log|Pdot/P| (this work) | 0.74 | 0.75 (period-matched) | 0.0000 |
| log|Pdot/P| hybrid | 0.65 | 0.69 (period-matched) | 0.0000 |
