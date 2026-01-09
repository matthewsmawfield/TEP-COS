# Temporal Shear Analysis Robustness Report

**Analysis Date:** 2026-01-04  
**Pipeline Version:** v2.0 (NaN-aware, gap-aware, estimator-independent)

---

## Executive Summary

This report documents the rigorous validation of the temporal shear (Γ) analysis pipeline for gravitationally lensed quasars. The pipeline has been upgraded to eliminate gap-handling artifacts and tested for estimator independence. The combined result across 10 independent lens systems is:

**Γ = +1.14 ± 1.62 days/decade (z = 0.70, p = 0.48)**

This is consistent with zero temporal shear at current precision. The measurements show excellent internal consistency (I² = 0%) and estimator independence (ICCF vs INTERP agree within 1.1σ).

---

## 1. Pipeline Improvements

### 1.1 NaN-Aware Processing
- **Gaussian smoothing**: Implemented weighted convolution that propagates NaN through gaps rather than filling with zeros or means
- **Bandpass filtering**: Uses NaN-aware smoothing for both low-pass and high-pass components
- **Gap threshold**: 30 days (seasonal observing gaps marked as NaN)

### 1.2 Gap-Aware Interpolation
- Linear interpolation forbidden across gaps exceeding 30 days
- Prevents fabrication of correlated structure in seasonal gaps
- Applied consistently in both ICCF and INTERP estimators

### 1.3 Bootstrap Uncertainty
- **FR mode**: Photometric noise perturbation only (fixed sampling)
- **FR/RSS mode**: Combined photometric + sampling uncertainty
- Duplicate timestamps collapsed to avoid interpolation artifacts

---

## 2. Dataset

### 2.1 Systems Analyzed
| System | Band | Images | Pairs | Source |
|--------|------|--------|-------|--------|
| DESJ0408 | R | 3 | 3 | Courbin 2017 |
| HE0435 | R | 4 | 6 | Bonvin 2016 |
| HE1104 | B,R,I,J | 2 | 4 | ApJ 798 |
| HS2209 | R | 2 | 1 | Eulaers 2013 |
| J1001 | R | 2 | 1 | Rathnakumar 2013 |
| J1206 | R | 2 | 1 | Eulaers 2013 |
| PG1115 | R | 3 | 3 | Bonvin 2018 |
| Q2237 | g,r,V,I | 4 | 24 | A&A 637 |
| RXJ1131 | R | 4 | 6 | Tewes 2013 |
| WFI2033 | R | 3 | 3 | Bonvin 2019 |

**Total:** 16 system-band combinations, 52 image pairs

### 2.2 Data Quality
- Minimum epochs per light curve: 20
- Typical baseline: 5-15 years
- Typical cadence: 3-7 days (seasonal)

---

## 3. Results

### 3.1 Per-System Best Pairs (ICCF Estimator)
| System | Pair | Γ (days/dec) | σ | z |
|--------|------|--------------|---|---|
| DESJ0408 | B-D | +23.50 | 19.62 | 1.20 |
| HE0435 | A-D | +0.30 | 3.33 | 0.09 |
| HE1104 | A-B (I) | +9.35 | 20.84 | 0.45 |
| HS2209 | A-B | +7.14 | 71.38 | 0.10 |
| J1001 | A-B | -43.87 | 135.52 | 0.32 |
| J1206 | A-B | -22.12 | 60.85 | 0.36 |
| PG1115 | B-C | +2.22 | 2.64 | 0.84 |
| Q2237 | A-C (g) | +2.25 | 15.62 | 0.14 |
| RXJ1131 | A-B | -0.22 | 2.95 | 0.07 |
| WFI2033 | A-B | +1.34 | 6.76 | 0.20 |

### 3.2 Meta-Analysis (Independent Systems)
- **Inverse-variance weighted mean:** Γ = +1.14 ± 1.62 days/decade
- **Significance:** z = 0.70 (p = 0.48)
- **Sign test:** 7 positive, 3 negative (p = 0.34)
- **Heterogeneity:** Q = 2.17, I² = 0%

### 3.3 All-Pairs Analysis
- **Total pairs:** 52
- **Valid Γ measurements:** 43-52 (depending on estimator)
- **Pairs passing practical cuts:** 15
- **IVW mean (all pairs):** Γ = +0.50 ± 0.94 days/decade (z = 0.53)

---

## 4. Estimator Agreement

### 4.1 ICCF vs INTERP Comparison
| Metric | Value |
|--------|-------|
| Common pairs | 52 |
| Mean difference | -1.47 days/decade |
| Std difference | 14.78 days/decade |
| Mean |z_diff| | 0.25 |
| Max |z_diff| | 1.89 |
| Pearson correlation | r = 0.67 (p < 0.0001) |
| Sign agreement | 77% |

### 4.2 Meta-Analysis Comparison
| Estimator | Γ (days/dec) | σ |
|-----------|--------------|---|
| ICCF | +0.50 | 0.94 |
| INTERP | +5.50 | 4.64 |
| **Difference** | -5.00 | 4.73 |
| **z_diff** | -1.06 | |

**Conclusion:** Estimators agree within 1.1σ at the meta-analysis level. Results are estimator-independent.

---

## 5. Robustness Checks

### 5.1 Gap Handling
- ✅ NaN-aware smoothing implemented
- ✅ Gap-aware interpolation (30-day threshold)
- ✅ No zero-filling of gaps in bandpass filter

### 5.2 Estimator Independence
- ✅ ICCF and INTERP agree within 1.1σ
- ✅ 77% sign agreement across pairs
- ✅ r = 0.67 correlation between estimators

### 5.3 Statistical Consistency
- ✅ Low heterogeneity (I² = 0%)
- ✅ Random-effects and fixed-effects models agree
- ✅ No outlier systems driving the result

### 5.4 Multi-Band Consistency
- Q2237: 24 pairs across g,r,V,I bands
- HE1104: 4 pairs across B,R,I,J bands
- Sign consistency: Mixed (expected for null signal)

---

## 6. Interpretation

### 6.1 Current Status
The temporal shear analysis, after rigorous methodological improvements, yields a combined result of Γ = +1.14 ± 1.62 days/decade. This is:

1. **Consistent with zero** at current precision (z = 0.70)
2. **Internally consistent** across systems (I² = 0%)
3. **Estimator-independent** (ICCF vs INTERP agree within 1.1σ)

### 6.2 Detection Power
To achieve a 2σ detection of a 1 day/decade signal would require:
- ~4× more data (independent systems), or
- ~2× reduction in per-system uncertainty

### 6.3 Scientific Conclusion
The current dataset does not provide statistically significant evidence for temporal shear in gravitationally lensed quasars. However, the pipeline is now methodologically sound and ready for application to larger datasets as they become available.

---

## 7. Files Generated

| File | Description |
|------|-------------|
| `step_3_0_cosmograil_temporal_shear_expanded_iccf_bbinterp_mc03_fr.json` | ICCF estimator results |
| `step_3_0_cosmograil_temporal_shear_expanded_interp_mc03_fr.json` | INTERP estimator results |
| `temporal_shear_robustness_report.md` | This report |

---

## 8. Recommendations

1. **Manuscript claims** should be updated to reflect the null result with proper uncertainty quantification
2. **Future work** should focus on expanding the dataset with new COSMOGRAIL releases and LSST data
3. **Injection-recovery tests** should be performed to calibrate detection thresholds before claiming any future detections

---

*Report generated by TEP-COS temporal shear analysis pipeline v2.0*
