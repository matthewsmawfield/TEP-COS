# Binary Pulsar Analysis (Freire GCpsr)
**Source:** /Users/matthewsmawfield/www/TEP-COS/results/outputs/freire_GCpsr.txt
**SHA256:** `ad6a26cc3270d518...`

## Sample Sizes
- **Total GC pulsars parsed:** 349
- **GC MSPs (P < 30 ms):** 325
- **Binary MSPs:** 111
- **Isolated MSPs:** 81

## Binary vs Isolated Comparison

| Metric | Binary MSPs | Isolated MSPs |
| --- | --- | --- |
| N | 111 | 81 |
| Mean log|Ṗ| | -19.273 | -18.967 |
| Std log|Ṗ| | 0.705 | 0.872 |
| Fraction negative Ṗ | 43.2% | 46.9% |

### Statistical Tests
- **Difference (binary - isolated):** -0.305 dex
- **Welch t-test p:** 0.01092
- **Mann-Whitney p:** 0.006703

## Interpretation

If the low |Ṗ| effect in GC pulsars were purely due to cluster acceleration, we would expect:
- Binary and isolated MSPs to be affected equally (same line-of-sight acceleration)
- No significant difference in log|Ṗ| between the two populations

Any significant difference would suggest:
- Population/selection effects (binary MSPs may have different intrinsic properties)
- Or a TEP-like effect that couples differently to binary vs isolated systems

**Result:** Significant difference detected (p = 0.01092). This warrants further investigation.