# Pulsar Population Controls (Freire + ATNF)
**Freire GCpsr URL:** https://www3.mpifr-bonn.mpg.de/staff/pfreire/GCpsr.txt\
**Freire SHA256:** `647ceabcd85221688ce4975cb270939186c3dd9667fdb80d99bd6623997f54c3`\
**ATNF psrcat_pkg URL:** https://www.atnf.csiro.au/research/pulsar/psrcat/downloads/psrcat_pkg.tar.gz\
**ATNF SHA256:** `c4330ee179cbec1e65c2d6e82c4890da107cd0f403ea488746a8a188220382f3`\

## Sample sizes
- **GC MSPs (Freire, P<30 ms, measured Pdot):** 196
- **Field MSPs (ATNF, P<30 ms, Pdot present, non-GC ASSOC):** 198

## Base test (log10|Pdot|)
- **GC mean:** -19.166\
- **Field mean:** -19.758\
- **Difference (GC-Field):** 0.592 dex\
- **Welch t-test p:** 8.99e-14\
- **Mann-Whitney p:** 6.73e-16\

## Controls
### Period-matched bootstrap
- **Mean diff:** 0.606 dex (16–84%: 0.551 to 0.663)\
- **Two-sided p:** 0\

### Period + B-proxy matched bootstrap
- **Mean diff:** 0.575 dex (16–84%: 0.524 to 0.630)\
- **Two-sided p:** 0\
