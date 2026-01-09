# Pulsar Population Controls (Freire + ATNF)
**Freire GCpsr URL:** https://www3.mpifr-bonn.mpg.de/staff/pfreire/GCpsr.txt\
**Freire SHA256:** `ad6a26cc3270d51840528bb59d1ef1c1310e5b20bf4cbb147e0659d7e75c2552`\
**ATNF psrcat_pkg URL:** https://www.atnf.csiro.au/research/pulsar/psrcat/downloads/psrcat_pkg.tar.gz\
**ATNF SHA256:** `c4330ee179cbec1e65c2d6e82c4890da107cd0f403ea488746a8a188220382f3`\

## Sample sizes
- **GC MSPs (Freire, P<30 ms, measured Pdot):** 181
- **Field MSPs (ATNF, P<30 ms, Pdot present, non-GC ASSOC):** 198

## Base test (log10|Pdot|)
- **GC mean:** -19.106\
- **Field mean:** -19.758\
- **Difference (GC-Field):** 0.652 dex\
- **Welch t-test p:** 2.92e-15\
- **Mann-Whitney p:** 3.47e-17\

## Controls
### Period-matched bootstrap
- **Mean diff:** 0.858 dex (16–84%: 0.796 to 0.919)\
- **Two-sided p:** 0\

### Period + B-proxy matched bootstrap
- **Mean diff:** 0.128 dex (16–84%: 0.111 to 0.145)\
- **Two-sided p:** 0\
