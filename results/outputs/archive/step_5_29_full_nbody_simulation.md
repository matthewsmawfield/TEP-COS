# Full N-Body Simulation Results

**Generated:** 2026-01-08T07:41:28.553431+00:00
**Method:** N-body integration using gala (100 steps, dt=0.01 Myr)

## Key Finding

| Metric | N-body Prediction | Observation |
|--------|-------------------|-------------|
| Density correlation | r = 0.902 | r ≈ 0 |
| Shift pattern | Scales with ρ | CONSTANT |

## Cluster-by-Cluster Results

| Cluster | log(ρc) | N-body Shift | Observed |
|---------|---------|--------------|----------|
| Terzan_5 | 5.5 | +1.947 | +0.13 |
| NGC_6440 | 5.4 | +1.646 | +0.13 |
| M62 | 5.2 | +1.392 | +0.13 |
| M15 | 5.0 | +0.959 | +0.13 |
| 47_Tuc | 4.8 | +0.709 | +0.13 |
| M28 | 4.5 | +0.925 | +0.13 |
| NGC_6752 | 4.3 | +0.766 | +0.13 |
| M13 | 3.8 | +0.518 | +0.13 |
| M5 | 3.5 | +0.559 | +0.13 |
| M71 | 3.2 | +0.107 | +0.13 |
| M53 | 3.0 | +0.369 | +0.13 |

## Interpretation

The N-body simulation confirms that standard gravitational dynamics produces
acceleration noise that scales with cluster density. The observed constant
residual (~0.13 dex) across all clusters contradicts this prediction.
