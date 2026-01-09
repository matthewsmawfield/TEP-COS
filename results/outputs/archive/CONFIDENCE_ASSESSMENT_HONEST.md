# Honest Confidence Assessment: TEP-COS Results

**Date**: 2026-01-05  
**Purpose**: Rigorous audit of all claimed results and methods  
**Updated**: 2026-01-05 (post-investigation)

---

## Executive Summary

**The lensing "discovery" is not robust.** The claimed Temporal Shear detections are dominated by mode-jumping artifacts in the cross-correlation estimator. The pulsar channel shows a real environmental signal but is confounded by known dynamical effects. The SDSS stellar archaeology channel is sensitive to control-set choices and does not provide independent confirmation.

**Confidence Level: LOW for discovery claims.**

---

## 1. Gravitational Lensing Channel

### What was claimed
- "Temporal Shear" Γ detected at 8.2σ in multiple lens systems
- Γ = −333 days/decade in DESJ0408 A-D (strongest detection)
- Effect correlates with source redshift as TEP predicts

### What the audit found

**CRITICAL ISSUE: Mode Jumping**

The Temporal Shear estimator measures time delay at multiple smoothing scales (τ = 5, 10, 20, 40, 80, 160 days) and fits a linear slope Γ = dΔt/d(log τ).

However, the delays at different scales are **not stable**:

| System | Pair | τ=10 delay | τ=20 delay | JUMP |
|--------|------|-----------|-----------|------|
| DESJ0408 | A-D | +196 days | −156 days | **−352 days** |
| DESJ0408 | B-D | −191 days | −50 days | +141 days |
| PG1115 | B-C | −153 days | +15 days | +168 days |
| J1206 | A-B | +2 days | −115 days | −118 days |

The large jumps (100-350 days) occur because the cross-correlation estimator finds **different alias peaks** at different smoothing scales. These aliases correspond to:
- Seasonal observing gaps (~400 days)
- Physical time-delay degeneracies
- Noise-driven spurious peaks

**Stability Test Result**

When excluding small-scale points (τ < 20) to avoid mode-jumping:
- **4 out of 5 claimed detections become UNSTABLE**
- Signs flip, magnitudes collapse, or both
- Only PG1115 A-B shows any stability

### Extended Investigation (2026-01-05)

We tested whether the mode-jumping could be resolved by:
1. Using lowpass smoothing instead of bandpass
2. Restricting to stable scale ranges
3. Analyzing long-baseline systems (5+ years)

**Results**:
- Long-baseline systems (HE0435, RXJ1131, WFI2033) show **STABLE** delays
- But their Γ values are **small**: typically |Γ| < 35 days/decade
- The large claimed signals (|Γ| > 100) come exclusively from short-baseline systems
- **NO STABLE DETECTION** exists with |Γ| > 50 days/decade

| System | Baseline | Stability | Max |Γ| |
|--------|----------|-----------|--------|
| HE0435 | 4583 days | STABLE | 11 days/dec |
| RXJ1131 | 3087 days | STABLE | 34 days/dec |
| WFI2033 | 5141 days | STABLE | 10 days/dec |
| DESJ0408 | 262 days | UNSTABLE | (artifact) |
| PG1115 | 178 days | UNSTABLE | (artifact) |

### Verdict: NOT A DISCOVERY

The "Temporal Shear" is an artifact of the estimator failing at small smoothing scales in short-baseline data. When stable long-baseline systems are analyzed, the signal is consistent with **zero** (|Γ| < 35 days/decade). This is a **mode-jumping artifact**, not a physical detection.

**The lensing discovery cannot be salvaged with existing data.**

---

## 2. Pulsar Timing Channel

### What was claimed
- Cluster MSPs spin down more slowly than field MSPs
- Residual offset of 0.13 dex after population controls
- Consistent with TEP enhancement factor ~10^6

### What the audit found

**Real environmental signal**: The 0.13 dex residual after matching on period and magnetic-field proxy is statistically robust (p < 0.01).

**Critical confound**: Globular cluster Ṗ measurements are contaminated by **line-of-sight acceleration** from the cluster potential. This is a well-known effect (Freire et al., Prager et al.) that produces both positive and negative apparent spin-down contributions.

**Field Binary Control**: The null result (p = 0.70) between binary and isolated field MSPs is valuable. It shows the cluster signal is environmental, not intrinsic to binarity.

**However**: This does not prove TEP. The environmental signal could be:
1. Standard GR dynamical acceleration (established mechanism)
2. TEP-enhanced time dilation (speculative)
3. Some combination

Without independent calibration of the cluster acceleration field, these cannot be distinguished.

### Verdict: DIAGNOSTIC, NOT DISCOVERY

The pulsar channel shows a real environmental effect but is confounded by known astrophysics. It is consistent with TEP but equally consistent with standard dynamical explanations. Treated as a diagnostic cross-check, not a standalone detection.

---

## 3. SDSS Stellar Archaeology Channel

### What was claimed
- At fixed [Mg/Fe], higher-σ galaxies appear younger (r = −0.27)
- Consistent with TEP time dilation

### What the audit found

**Global correlations**: Yes, the correlation exists (r = −0.04 to −0.27 depending on controls).

**Critical test**: The twin-galaxy matched-pair experiment controls for:
- Redshift (same coordinate epoch)
- [Mg/Fe] proxy (nucleosynthesis clock)
- Local environment (kNN density)
- Geometric size (half-light radius in kpc)

**Result**: Under geometric controls (no stellar mass), the sign **flips** — higher-σ galaxies appear slightly **older**, not younger.

**Sensitivity**: When stellar mass is included in matching (despite potential TEP contamination), the sign flips back to TEP-consistent. This means the "signal" is sensitive to whether you use potentially contaminated observables in the control set.

### Verdict: CONSTRAINT, NOT CONFIRMATION

The SDSS analysis provides a constraint on the conceptual framework but does not independently confirm TEP. The sign-sensitivity to control-set choice indicates vulnerability to confounding.

---

## 4. Methods Assessment

### Data Provenance
- **COSMOGRAIL**: Public .rdb files from VizieR/COSMOGRAIL collaboration ✓
- **ATNF PSRCAT**: Standard pulsar catalog ✓
- **Freire GCpsr**: Maintained globular cluster pulsar database ✓
- **SDSS DR18**: Standard SQL queries to SkyServer ✓

Data provenance is solid.

### Analysis Pipeline
- **Lensing**: The mode-locking was intended to prevent alias jumping, but it fails at small scales. The estimator is fundamentally vulnerable to seasonal gaps.
- **Pulsars**: Population matching is appropriate, but cannot disentangle TEP from dynamical acceleration.
- **SDSS**: Control-set sensitivity makes conclusions fragile.

### Reproducibility
- All scripts are in `scripts/steps/`
- Results are in `results/outputs/`
- The analysis is reproducible, but reproducibility of an artifact is not validation

---

## 5. What Would Be Needed for a Defensible Discovery

### Lensing
1. **Stable delays**: The delay at each τ should be within ~10-20 days of the broadband value, not jumping by 100-400 days
2. **Multi-band achromaticity**: Same Γ in g, r, i bands (currently untested for key systems)
3. **Independent estimator**: Different delay algorithm (e.g., PyCS, JAVELIN) producing consistent Γ

### Pulsars
1. **Dynamical modeling**: Independent calibration of cluster acceleration field from proper motions + mass models
2. **Radial gradient within clusters**: If TEP is real, Ṗ should correlate with cluster-centric radius (currently heterogeneous results)

### SDSS
1. **Coeval populations**: Stars that truly formed together but are now in different potentials (e.g., stellar streams)
2. **Geometric-only controls**: No reliance on stellar-population-derived quantities

---

## 6. Bottom Line

| Channel | Claimed | Actual Status |
|---------|---------|---------------|
| Lensing | 8.2σ detection | **Artifact** (mode-jumping) |
| Pulsars | Diagnostic signal | **Confounded** (dynamical acceleration) |
| SDSS | TEP-consistent r = −0.27 | **Fragile** (sign-sensitive to controls) |

**Honest assessment**: There is no defensible discovery-level result in the current TEP-COS analysis. The lensing "detection" is the most serious problem — it is likely an artifact of the analysis method, not a physical signal.

The pulsar and SDSS channels show patterns that are *consistent with* TEP but are equally consistent with standard astrophysics or confounded by methodological choices.

---

## 7. Recommendations

1. **Do not publish lensing claims** until the mode-jumping issue is resolved or acknowledged
2. **Reframe pulsars** as "environmental phenomenology" requiring future disambiguation
3. **Reframe SDSS** as a constraint/stress-test, not confirmation
4. **Consider retraction** of strong claims in manuscript if they cannot be substantiated

This is the honest assessment.
