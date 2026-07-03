# Binary Pulsar Analysis (Freire GCpsr)
**Source:** /Users/matthewsmawfield/www/Temporal Equivalence Principle/TEP-COS/data/freire_GCpsr.txt
**SHA256:** `a8520b1133eb7d16...`

## Sample Sizes
- **Total GC pulsars parsed:** 368
- **GC MSPs (P < 30 ms):** 344
- **Binary MSPs:** 117
- **Isolated MSPs:** 81

## Binary vs Isolated Comparison

| Metric | Binary MSPs | Isolated MSPs |
| --- | --- | --- |
| N | 117 | 81 |
| Mean log|Ṗ| | -19.296 | -18.967 |
| Std log|Ṗ| | 0.696 | 0.872 |
| Fraction negative Ṗ | 43.6% | 46.9% |

### Statistical Tests
- **Difference (binary - isolated):** -0.329 dex
- **Welch t-test p:** 0.005636
- **Mann-Whitney p:** 0.002964

## Interpretation

If the low |Ṗ| effect in GC pulsars were purely due to cluster acceleration, we would expect:
- Binary and isolated MSPs to be affected equally (same line-of-sight acceleration)
- No significant difference in log|Ṗ| between the two populations

Any significant difference would suggest:
- Population/selection effects (binary MSPs may have different intrinsic properties)
- Or a TEP-like effect that couples differently to binary vs isolated systems

**Result:** Significant difference detected (p = 0.005636). This warrants further investigation.