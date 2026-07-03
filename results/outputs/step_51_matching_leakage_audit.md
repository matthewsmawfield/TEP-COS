# Step 51: Matching Leakage Audit

**Generated:** 2026-07-03T10:05:29.083472+00:00

## Results

| Variant | Residual (dex) | Gamma | Gamma err |
|---------|----------------|-------|-----------|
| Raw (no matching) | +0.591 | 0.320 | +/- 0.107 |
| Period-only matching (5 NN) | +0.817 | 0.395 | +/- 0.087 |
| Period + log(tau_c) matching (5 NN) | +0.198 | 0.090 | +/- 0.019 |
| Period + B-proxy matching (5 NN) | +0.233 | 0.118 | +/- 0.029 |
| Period + B-proxy matching (15 NN) | +0.357 | 0.184 | +/- 0.039 |
| Period + B-proxy cluster-residualized (5 NN) | +0.203 | 0.134 | +/- 0.029 |

## Interpretation

- **Raw** gives the unattenuated signal.
- **Period-only** should recover nearly the full amplitude if matching is clean.
- **Period + B-proxy** will attenuate if B-proxy leaks outcome information.
- **Expanded NN (15)** may over-smooth and further suppress amplitude.
- **Cluster-residualized** subtracts the field mean per cluster; if the field mean
carries density-correlated structure, this will attenuate the slope.
