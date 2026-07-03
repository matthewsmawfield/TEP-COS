# Step 52: Signed Observable Battery

**Generated:** 2026-07-03T10:05:29.475339+00:00

## Raw Comparisons

| Observable | GC | Field | Delta | N_GC | N_Field |
|------------|----|-------|-------|------|---------|
| log|Ṗ| | -19.170 | -19.761 | +0.591 | 198 | 202 |
| log|Ṗ/P| | -16.827 | -17.566 | +0.739 | 198 | 202 |
| signed Ṗ/P | 6.680e-18 | 8.288e-18 | -1.607e-18 | 198 | 202 |
| positive-only Ṗ | -19.074 | -19.757 | +0.683 | 109 | 198 |
| negative-only Ṗ | -19.288 | -19.940 | +0.651 | 89 | 4 |
| width / MAD | 0.499 | 0.336 | +0.163 | 198 | 202 |

## Matched Battery

| Variant | Residual | Gamma | Gamma err |
|---------|----------|-------|-----------|
| log|Ṗ|  raw | +0.591 | 0.320 | +/- 0.107 |
| log|Ṗ|  period-only | +0.817 | 0.395 | +/- 0.087 |
| log|Ṗ|  period+B-proxy | +0.233 | 0.118 | +/- 0.029 |
| log|Ṗ/P| raw | +0.739 | 0.343 | +/- 0.089 |
| log|Ṗ/P| period-only | +0.822 | 0.397 | +/- 0.088 |
| log|Ṗ/P| period+B-proxy | +0.281 | 0.134 | +/- 0.035 |
| signed Ṗ/P raw | -0.000 | 0.000 | +/- 0.000 |
| signed Ṗ/P period-only | +0.000 | 0.000 | +/- 0.000 |
| signed Ṗ/P period+B-proxy | -0.000 | 0.000 | +/- 0.000 |

## Interpretation

- **log|Ṗ/P|** is the direct acceleration observable; if TEP acts on Ṗ/P, this
  should show the cleanest signal.
- **signed Ṗ/P** preserves direction; an excess of negative values in GCs
  indicates acceleration-dominated line-of-sight components.
- **positive-only / negative-only** isolate the intrinsic and acceleration
  branches respectively.
- **width / MAD** tests whether the effect manifests as broadening rather than
  a mean shift (expected if TEP is a stochastic or multi-component process).
