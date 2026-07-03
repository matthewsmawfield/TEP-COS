# Pulsar Population Controls (Freire + ATNF)
**Freire GCpsr URL:** https://www3.mpifr-bonn.mpg.de/staff/pfreire/GCpsr.txt\
**Freire SHA256:** `a8520b1133eb7d16cec8f6ede01e9248141e4aa806915ade6a21e1eae6cecf3e`\
**ATNF psrcat_pkg URL:** https://www.atnf.csiro.au/research/pulsar/psrcat/downloads/psrcat_pkg.tar.gz\
**ATNF SHA256:** `2ea143fc2c9a01184f5a94715422351c49ecf97cae26004ea4897e56421d4cbd`\

## Sample sizes
- **GC MSPs (Freire, P<30 ms, measured Pdot):** 198
- **Field MSPs (ATNF, P<30 ms, Pdot present, non-GC ASSOC):** 202

## Base test (log10|Pdot|)
- **GC mean:** -19.170\
- **Field mean:** -19.761\
- **Difference (GC-Field):** 0.591 dex\
- **Welch t-test p:** 5.42e-14\
- **Mann-Whitney p:** 3.16e-16\

## Controls
### Period-matched bootstrap
- **Mean diff:** 0.612 dex (16–84%: 0.553 to 0.667)\
- **Two-sided p:** 0.0005\

### Period + B-proxy matched bootstrap
- **Mean diff:** 0.609 dex (16–84%: 0.560 to 0.657)\
- **Two-sided p:** 0.0005\

## Period Cut Sensitivity Analysis
Testing robustness of signal to MSP period boundary choice.

| Period Cut | GC N | Field N | Raw Diff (dex) | Period-Matched (dex) | 2D-Matched (dex) | p-value |
|------------|------|---------|----------------|----------------------|------------------|---------|
| P < 10 ms | 177 | 151 | 0.777 | 0.778 [0.717, 0.840] | 0.778 [0.717, 0.840] | 0.0005 |
| P < 30 ms | 198 | 202 | 0.591 | 0.612 [0.553, 0.667] | 0.609 [0.560, 0.657] | 0.0005 |
| P < 50 ms | 200 | 229 | 0.468 | 0.614 [0.556, 0.671] | 0.592 [0.549, 0.637] | 0.0005 |

**Interpretation:** The signal persists across period cut choices, demonstrating robustness to the P < 30 ms boundary definition.
