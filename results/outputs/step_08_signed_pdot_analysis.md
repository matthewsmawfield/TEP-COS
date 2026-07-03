# Step 08: Signed P-dot Latent-Mixture Analysis

**Timestamp:** 2026-07-03T10:04:31.819831+00:00

## Sample Counts

- GC MSPs: 198
- Field MSPs: 202

## Sign Analysis (Signed Pdot/P)

| Environment | N | Positive | Negative | Frac. Positive | Median signed Pdot/P (s⁻¹) |
|-------------|---|----------|----------|----------------|----------------------------|
| GC          | 198 | 109 | 89 | 0.551 | 2.929e-18 |
| Field       | 202 | 198 | 4 | 0.980 | 2.404e-18 |

**Difference in positive fraction:** -0.430  
**Fisher exact test p-value:** 9.49e-28

## Latent Mixture (2-component Gaussian on log|Pdot/P|)

### GC Mixture
- Means (dex): ['-17.237', '-16.581']
- Stds (dex): ['0.798', '0.554']
- Weights: ['0.375', '0.625']

### Field Mixture
- Means (dex): ['-17.664', '-17.419']
- Stds (dex): ['0.224', '0.813']
- Weights: ['0.600', '0.400']

### Inferred Acceleration Component
- Mean offset (dex): 0.182
- Std (dex): 0.000

## log|Pdot/P| Base Comparison

- GC mean: -16.827 dex
- Field mean: -17.566 dex
- Difference: 0.739 dex
- t-statistic: 11.36, p = 8.01e-26
- Mann-Whitney U p = 2.28e-26

## Controls on log|Pdot/P|

### Period-Matched Bootstrap (N_boot=2000)
- Mean difference: 0.747 dex
- 16th–84th percentile CI: [0.695, 0.800]
- Two-sided p: 0.0000

### Period+B-Proxy Matched Bootstrap (N_boot=2000)
- Mean difference: 0.746 dex
- 16th–84th percentile CI: [0.697, 0.794]
- Two-sided p: 0.0000
