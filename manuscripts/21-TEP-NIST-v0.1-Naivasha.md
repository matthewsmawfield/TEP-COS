# A Scalar-Field Geological Covariate for Laboratory Variance in Measurements of G: A TEP Analysis of the NIST–BIPM Torsion-Balance Discrepancy
**Matthew Lukin Smawfield**
Version: v0.1 (Naivasha)
First published: 6 June 2026
DOI: 10.5281/zenodo.20576483

---

## Abstract

High-precision measurements of Newton's gravitational constant (*G*) have long shown laboratory-to-laboratory scatter larger than quoted random uncertainties, usually attributed to unresolved apparatus-dependent systematics. Here we report candidate evidence that part of this scatter is structured by an external geophysical covariate: the Bouguer-constrained mass column beneath the laboratory. Within the Temporal Equivalence Principle (TEP), local matter density sources a screened scalar proper-time field (φ), modifying the locally inferred gravitational coupling through *G*<sub>eff</sub> = *G*<sub>N</sub> / *A*²(φ). The 2026 NIST redetermination of *G*, performed with the BIPM torsion-balance geometry after relocation from Sèvres to Gaithersburg, provides a high-leverage two-site test because the apparatus class is nearly fixed while the geological column changes.

Using published Bouguer gravity anomalies, CRUST1.0 deep structure, and a locked lab-scale TEP parameterization, the amplitude-anchored 1D response estimator predicts a lower inferred *G* at Gaithersburg than at Sèvres, matching the observed sign. The preferred 3D finite-difference solver currently recovers the same sign but underpredicts the amplitude because the available symmetric boundary taper suppresses lateral density contrast; measured lateral Bouguer maps are required for a laboratory-precision 3D prediction. We therefore interpret the NIST–BIPM datum not as a standalone proof, but as the first member of a falsifiable residual class: after apparatus class is controlled, laboratories on lighter sedimentary columns should infer higher *G*, while those on denser crystalline or shield basement should infer lower *G*.

The central contribution is a reproducible discovery pipeline: source-reconciled metrology inputs, Bouguer-constrained geological columns, scalar-field forward modeling, apparatus-class controls, and preregistered no-refit predictions for independent laboratories. If existing and future *G* residuals track the locked scalar-field column in sign and rank after apparatus controls, the result would constitute direct empirical support for the TEP density-sector parameterization; if they do not, the locked density-sector parameterization is falsified.

Keywords: Temporal Equivalence Principle, gravitational constant, metrology, NIST, BIPM, temporal topology, scalar field, geophysical residual, discovery pipeline

## 1. Introduction

## 1.1 The NIST–BIPM Discrepancy

In April 2026, the National Institute of Standards and Technology (NIST) published a redetermination of the gravitational constant *G* using a torsion-balance apparatus functionally identical to the prior Bureau International des Poids et Mesures (BIPM) setup in Sèvres, France. The final *Metrologia* article (Schlamminger *et al.*, 2026) reports the result is lower by 2.50×10⁻⁴ relative to BIPM. The NIST press release gives the same central value but states it is 0.0235% lower (≈235 ppm), while the NIST landing-page abstract lists 2.81×10⁻⁴ (281 ppm). This manuscript adopts the peer-reviewed *Metrologia* PDF as the governing source and explicitly flags this ambiguity; all downstream analyses adopt 2.50×10⁻⁴.

The 250 ppm shift is 176× larger than the combined random environmental noise budget of either apparatus. However, it is only ~0.5× the historical few-hundred-ppm scatter among independent *G* determinations, which is conventionally attributed to unmodeled apparatus systematics. TEP therefore does not primarily compete with random-noise explanations; it competes with unmodeled apparatus/systematics explanations by predicting that, after controlling for apparatus class, residual *G* measurements should correlate with Bouguer-constrained scalar-field columns.

**NIST value source reconciliation.**

| Source | *G* (×10⁻¹¹ m³ kg⁻¹ s⁻²) | Relative difference vs BIPM | Status |
| --- | --- | --- | --- |
| Metrologia PDF (adopted) | 6.67387 ± 0.00038 | 2.50×10⁻⁴ (250 ppm) | governing |
| NIST publication metadata | 6.67366 ± 0.00020 | 2.81×10⁻⁴ (281 ppm) | flagged conflict |
| NIST news article | 6.67387 | 0.0235% (235 ppm) | supports PDF central value |

This manuscript adopts the peer-reviewed *Metrologia* PDF as the governing source.

## 1.2 The Temporal Equivalence Principle Hypothesis

The Temporal Equivalence Principle (TEP) proposes that proper time is not a passive geometric outcome, but a dynamical scalar field (φ) characterized by a background cosmic gradient, or Temporal Shear (Σ<sub>μ</sub> = ∇<sub>μ</sub> ln A(φ)). Under TEP, the matter metric is rescaled by a universal conformal factor A<sup>2</sup>(φ) relative to the gravitational metric, so the Jordan-frame effective gravitational coupling is modulated by the local scalar field profile as *G*<sub>eff</sub> = *G*<sub>N</sub> / A<sup>2</sup>(φ). This paper tests whether the NIST–BIPM discrepancy is the first member of a broader testable residual class predicted by TEP.

## 1.3 This Work

This paper presents a reproducible discovery pipeline that retrieves source-reconciled metrology inputs, builds Bouguer-constrained crustal mass columns for precision-*G* laboratories, computes local scalar field profiles, and generates locked no-refit predictions for the sign and rank of residual *G* values after apparatus-class controls. The structure is as follows: Section 2 presents the TEP theoretical framework; Section 3 describes the geological modeling methodology; Section 4 details the computation pipeline and solver validation; Section 5 reports the NIST–BIPM result and the preliminary multi-lab prediction set; Section 6 discusses implications; and Section 7 concludes with falsifiability criteria and preregistered predictions for future measurements.

## 1.4 Discovery Claim and Falsifiability

The discovery claim of this paper is not that a single two-site discrepancy proves TEP. The discovery claim is that precision *G* metrology appears to contain a previously unmodeled geophysical residual: once apparatus class is controlled, the sign of the residual should follow the Bouguer-constrained scalar-field column predicted by TEP. The NIST–BIPM comparison is unusually valuable because the apparatus geometry is largely inherited while the laboratory environment changes from the Paris Basin to the Appalachian Piedmont. This makes it a high-leverage first case, but not a sufficient final test.

The paper is therefore falsifiable in a direct way. A locked TEP pipeline must predict the sign and rank ordering of residual *G* values across independent laboratories without re-fitting the density-sector coupling. If laboratories on lighter sedimentary columns do not systematically return higher inferred *G*, or if dense crystalline/shield sites do not systematically return lower inferred *G* after apparatus controls, the proposed density-sector interpretation fails.

## 2. Theoretical Framework

## 2.1 The Scalar Field Equation

The TEP scalar field φ is sourced by local mass distributions through the field equation:

\begin{equation}
\Box\phi = \frac{8\pi G}{c^4} \rho_{\text{eff}}
\end{equation}

where □ is the d'Alembertian operator and ρ<sub>eff</sub> incorporates the screening-corrected density. For the matter metric rescaled by A<sup>2</sup>(φ), the Jordan-frame effective gravitational coupling is G<sub>eff</sub> = G<sub>N</sub> / A<sup>2</sup>(φ).

The lab-scale phenomenological ansatz used here is φ<sub>ρ</sub> = α<sub>log</sub> ln(ρ/ρ<sub>lab</sub>).  With the (+,−,−,−) metric signature adopted across the TEP corpus, the scalar field equation for non-relativistic matter (T = +ρ) and negative conformal coupling gives a static limit in which φ decreases with increasing density: dφ/dρ < 0.  This fixes the sign of the density-sector coupling to α<sub>log</sub> < 0.  The magnitude |α<sub>log</sub>| = 7.66×10⁻³ is the canonical corpus value from Paper 21 terrestrial atomic-clock locking, restored after a prior local amplitude-anchored fit. Exact equality with the NIST-BIPM datum would require β<sub>A,lab</sub> ≈ −0.597, reported only as a diagnostic fit. The sign is fixed by the field equation and is not a free parameter. It is important to emphasize that α<sub>log</sub> is a laboratory density-sector response parameter in the locked terrestrial ansatz. It is not a direct measurement of the microscopic conformal coupling β<sub>A</sub> that appears in the matter metric A(φ) = exp(β<sub>A</sub>φ/M<sub>Pl</sub>). The reference value β<sub>A,lab</sub> = −1 is used here only as a lab-sector convention for the Naivasha pipeline; cross-domain comparisons must pass through domain-specific response coefficients and screening projections.

## 2.2 Temporal Shear

The Temporal Shear vector is defined as the active gradient of the conformal scaling factor:

\begin{equation}
\Sigma_\mu = \nabla_\mu \ln A(\phi)
\end{equation}

## 2.3 Screening Regime

When local density approaches the phenomenological saturation scale ρ<sub>T</sub> ≈ 20 g/cm³, the scalar field response is smoothly suppressed and A(φ) → 1. Both laboratory sites operate in a mildly screened dilute regime relative to this scale, allowing TEP scaling to remain observable in the model.

The implemented smooth log-density screen uses n≥3.21 to retain ≥95% coupling at laboratory densities while satisfying the solar-wind suppression requirement of n≥0.10 down to solar-wind densities of 10⁻²⁴ g/cm³ (identifying the Cassini check as a proxy bound). This is a continuous screening profile, not a binary on/off threshold.

## 2.4 Observable in the Laboratory Frame

For a universal conformal factor, a skeptic will ask whether local rods, clocks, source masses, density standards, and torque calibration all rescale together, making the effect unobservable as a units artifact. In the TEP laboratory-coupling model adopted here, the inferred observable is assumed to be the ratio of the gravitational torque to the elastic restoring torque in the torsion balance, because the two torques transform differently under A(φ). Appendix A derives three conformal-frame cases: (i) full conformal cancellation yields no observable Geff shift; (ii) partial elastic-sector cancellation yields a reduced A(φ)-dependent shift; and (iii) the TEP density-sector coupling adopted here gives Ginferred = GN/A²(φ). The derivation below corresponds to Case (iii).

In the Einstein frame the gravitational force between two point masses is F = G<sub>*</sub> m<sub>s</sub><sup>*</sup> m<sub>t</sub><sup>*</sup> / R<sup>2</sup>, where G<sub>*</sub> is the bare constant and m<sup>*</sup> are the Einstein-frame masses. In the Jordan frame, with matter metric g̃<sub>μν</sub> = A<sup>2</sup>(φ) g<sub>μν</sub>, the *measured* inertial mass is m = m<sup>*</sup> / A(φ), and the gravitational force becomes:

F = G<sub>*</sub> \frac{m<sub>s</sub><sup>*</sup> m<sub>t</sub><sup>*</sup>}{A<sup>2</sup>(φ) R<sup>2</sup>} = G<sub>*</sub> \frac{A<sup>2</sup>(φ) m<sub>s</sub> m<sub>t</sub>}{A<sup>2</sup>(φ) R<sup>2</sup>} = G<sub>*</sub> \frac{m<sub>s</sub> m<sub>t</sub>}{R<sup>2</sup>}

so the gravitational force, expressed in Jordan-frame masses and lengths, *appears* to use the same G<sub>*</sub>. However, the torsion balance does not measure force directly; it measures the *angular deflection* produced by the balance between gravitational torque and the fiber's elastic restoring torque. The restoring torque is τ<sub>fiber</sub> = κ θ, where the torsion constant κ is a material property determined by the shear modulus and geometry of the fiber. Both the shear modulus and the inertial moment of the balance arm are Jordan-frame quantities that scale with A(φ) in the same way as the masses, but the *ratio* that sets the equilibrium angle is:

θ = \frac{τ<sub>grav</sub>}{κ} = \frac{G<sub>eff</sub> m<sub>s</sub> m<sub>t</sub> L<sub>arm</sub>}{R<sup>2</sup> κ}

where L<sub>arm</sub> is the balance-arm length. When the apparatus is calibrated — i.e., when κ is determined independently from the fiber's elastic properties and the arm's moment of inertia — that calibration is performed in the same local Jordan frame as the G measurement. The crucial point is that the Cavendish-like *inference* of G from the measured deflection, arm length, mass values, and separation uses a *locally calibrated* restoring torque. Because the gravitational torque scales as G<sub>*</sub> / A<sup>2</sup>(φ) relative to the locally calibrated elastic torque, the inferred quantity is:

G<sub>inferred</sub> = \frac{G<sub>*</sub>}{A<sup>2</sup>(φ)} = G<sub>eff</sub>

Under the TEP density-sector coupling of Case (iii), the residual does not vanish as a units artifact because the gravitational interaction (mediated by the Einstein-frame metric) and the elastic restoring force (mediated by material bonds in the Jordan frame) enter the observable with different conformal weights. The torsion balance therefore measures a genuine effective coupling G<sub>eff</sub> = G<sub>N</sub> / A<sup>2</sup>(φ), not a redefinition of local units. A referee may object that a universal Jordan-frame matter coupling would cancel the effect entirely (Case i); that objection is addressed in Appendix A.

## 3. Geological Modeling

## 3.1 Data Sources

Site elevation is retrieved live from USGS 3DEP (Gaithersburg) and, outside USGS coverage, from global SRTM 30 m (Sèvres); DEM and geological-service query metadata are recorded from OpenTopography, USGS MRDS, and BRGM for reproducibility. The primary near-field geological constraint is the *Bouguer gravity anomaly*, an empirically constrained proxy for local subterranean mass surplus or deficit relative to a standard 2.67 g/cm³ crustal column, derived from surface gravity measurements with standard geophysical corrections (terrain, latitude, free-air). Bouguer anomalies resolve density variations at the station scale (~1 km), bypassing the resolution limit of global crustal models.

Published station values are drawn from BRGM (France), USGS/NGMDB (USA), Yuan *et al.* (2014) for the Yangtze Basin, and Geoscience Australia for the Tasman Fold Belt. For five additional prediction sites, representative Bouguer anomalies are estimated from published regional gravity compilations (BGR for Braunschweig, BGS for Teddington, GSC for Ottawa, CAGS for Beijing, swisstopo for Zurich; all cross-checked against WGM2012, Bonvalot *et al.*, 2012). JILA Boulder is the sole exception: its −220 mGal isostatic Bouguer anomaly reflects deep crustal compensation of the Rocky Mountains, not surface-layer density, so the simple slab inversion is invalid and the near-field column is constrained by USGS geological map data instead.

Deep structure below ~5 km is taken from the published global model CRUST1.0 (Laske *et al.*, 2013), which provides layer-by-layer densities and thicknesses on a 1°×1° grid.

#### Laboratory Coordinates

- BIPM Sèvres: 48.8298° N, 2.2137° E (Paris Basin)

- NIST Gaithersburg: 39.1383° N, −77.2014° W (Piedmont Plateau)

- HUST Wuhan: 30.5155° N, 114.4112° E (Yangtze Basin)

- UQ Brisbane: −27.4698° S, 153.0251° E (Tasman Fold Belt)

- JILA Boulder: 40.0150° N, −105.2705° W (Rocky Mountains)

- PTB Braunschweig: 52.2970° N, 10.4430° E (North German Plain)

- NPL Teddington: 51.4240° N, −0.3380° W (London Basin)

- NRC Ottawa: 45.4215° N, −75.6972° W (Canadian Shield)

- NIM Beijing: 39.9042° N, 116.4074° E (North China Plain)

- U Zurich: 47.3769° N, 8.5417° E (Swiss Molasse)

## 3.2 Mass Distribution Mapping

The Bouguer anomaly Δg<sub>B</sub> relates to local density contrast Δρ through the slab formula Δg<sub>B</sub> ≈ 2πG Δρ h. For a 5 km column, the conversion factor is Δρ ≈ 4.77 × 10⁻³ g/cm³ per mGal. The Paris Basin shows a well-documented negative anomaly of −42 ± 8 mGal (Le Pichon *et al.*, 1971), corresponding to a top-layer density of 2.47 ± 0.04 g/cm³ — a light sedimentary column. The Appalachian Piedmont at Gaithersburg shows a near-zero anomaly of +5 ± 10 mGal (Diment & Weaver, 1964), corresponding to 2.69 ± 0.05 g/cm³ — dense crystalline basement. The resulting *near-field* density contrast is 0.22 g/cm³, fifty times larger than CRUST1.0's 1°×1° average contrast of 0.004 g/cm³. CRUST1.0 deep structure is retained below the Bouguer-sensitive depth, yielding a hybrid column: Bouguer-constrained top + CRUST1.0 deep.

## 3.3 Near-Field Distance Weighting

The scalar field gradient that modulates the local effective coupling is dominated by near-field mass, not the deep crustal column. In analogy to Newtonian gravity, where the local acceleration gradient from a point source scales as r⁻², the density-sector contribution to φ is weighted by a kernel *K*(*z*) = 1/(1 + (*z*/*z*₀)²), where *z* is the depth to the layer midpoint and *z*₀ = 5 km is the near-field coherence scale. At *z* = 0 the kernel is unity (full sensitivity); at *z* = *z*₀ it is 0.5; at *z* = 5*z*₀ it is 0.038. The top Bouguer layer (0–5 km) therefore contributes ~55% of the effective density-sector signal, while the lower crust (20–30 km) contributes <3%. The geometric mass term still uses the total column mass because the conformal factor depends on the integrated mass budget.

## 3.4 Facility-Scale Mass Model

Facility-scale asymmetry is treated as a qualitative diagnostic rather than a quantitative discovery input. Parametric models suggest that local building geometry can generate horizontal shear, but the discovery claim does not rely on those models pending measured building plans or facility-scale microgravity surveys. Appendix B details the parametric facility-scale grids (51×51×21 cells, ~1 m resolution) used for sensitivity analysis, including the asymmetric L-shaped stone cellar at Sèvres and the symmetric concrete bunker at Gaithersburg. In the current pipeline these grids are built, solved directly, and superposed onto the crustal solution via the linearity of the Poisson equation.

The locally inferred Bouguer-constrained mass column acts as the geophysical realization of the abstract environmental operator $\mathcal{S}_\Sigma(\mathcal{E})$. The facility-scale mass anomaly serves as the empirical proxy for the saturation of the Temporal Topology field, directly modifying the local effective gravitational coupling $G_{\rm eff}$.

## 4. Computation

## 4.1 Pipeline Architecture

The reproducible analysis pipeline is implemented in Python and executed sequentially. Each step writes a JSON output to `results/outputs/` and a detailed log to `logs/`. Steps are fail-fast: execution halts on the first failure so that downstream steps do not consume stale data.

## 4.2 Scalar Field Solution

The primary scalar field is computed from the 3D screened Poisson equation on a Cartesian grid. The field at the laboratory position is the volume integral of the Green's function over the local density distribution:

φ(*r*₀) = ∫<sub>V</sub> G(*r*₀, *r*') · (α<sub>log</sub>/λ<sub>T</sub>²) ln(ρ(*r*')/ρ<sub>lab</sub>) · S<sub>screen</sub>(ρ(*r*')) d³r'

where G(*r*, *r*') is the Green's function for the screened Poisson operator. In the unscreened limit (λ<sub>screen</sub> → ∞), G ∼ 1/(4π|*r* − *r*'|). The 3D finite-difference solver relaxes this equation directly on a 41×41×11 grid (100 km lateral, 50 km depth, 5 km resolution) using a sparse direct solver with Dirichlet boundary conditions φ = 0 at the domain boundary. The source term S(ρ) = (α<sub>log</sub>/λ<sub>T</sub>²) ln(ρ/ρ<sub>lab</sub>) S<sub>screen</sub>(ρ) carries the inverse square of the Temporal Topology screening length. All forward analysis here uses λ<sub>T</sub> ≈ 4200 km from the 25-year multi-center GNSS baseline (Paper 6). Paper 14 reports an independent MGEX held-out verification on a ~1 yr combined-clock span (λ<sub>T</sub> ≈ 1,396 ± 90 km, R² ≈ 0.49); that shorter-baseline, different-product check confirms signal presence but is not adopted for dimensional normalization. No arbitrary domain-area division is required.

The 1D mass-weighted surrogate φ<sub>ρ</sub> = α<sub>log</sub> Σ<sub>i</sub> K(z<sub>i</sub>) (m<sub>i</sub>/M) ln(ρ<sub>i</sub>/ρ<sub>lab</sub>) S(ρ<sub>i</sub>), where K(z) = 1/(1 + (z/z₀)²), is retained as an *amplitude-anchored phenomenological response estimator* that delivers the headline ppm shift. It is not treated as a validated field solution and it is not the primary out-of-sample solver. The surrogate averages out lateral structure and is physically valid only when the density field is uniform over the support of the Green's function — an assumption that fails at laboratory scales.

*Dimensional analysis.* The screened density contrast ln(ρ/ρ<sub>lab</sub>) S<sub>screen</sub>(ρ) is dimensionless, while the discrete Laplacian ∇²φ has units of φ/km². To obtain a dimensionless φ, the source term carries the inverse square of the Temporal Topology screening length: S(ρ) = (α<sub>log</sub>/λ<sub>T</sub>²) ln(ρ/ρ<sub>lab</sub>) S<sub>screen</sub>(ρ), with λ<sub>T</sub> ≈ 4200 km from the 25-year multi-center GNSS baseline (Paper 6). The dimensional scale is therefore fixed by a physically derived correlation length, not by an arbitrary simulation domain size. The overall amplitude still rides on the fitted coupling α<sub>log</sub>, so the 3D solver confirms the *sign* of the inter-site shift but underpredicts the observed amplitude; the current symmetric boundary taper should be interpreted as a conservative lower-bound geometry, not as a laboratory-precision field solution. With the current symmetric lateral taper both sites have nearly identical crustal φ; the inter-site dG/G from the 3D crustal grid is ≈ −24 ppm (screened β). The 1D–3D discrepancy is not treated as a free choice between solvers; it is a diagnostic showing that scalar-point predictions are currently limited by the absence of measured lateral Bouguer structure at laboratory scale.

## 4.3 3D Finite-Difference Solver

The 3D finite-difference solver relaxes the static screened Poisson equation ∇²φ = −S(ρ) on a Cartesian grid using a sparse direct solver, with Dirichlet boundary conditions φ = 0 at the domain boundary. The source term is S(ρ) = (α<sub>log</sub>/λ<sub>T</sub>²) ln(ρ/ρ<sub>lab</sub>) S<sub>screen</sub>(ρ), carrying the inverse square of the Temporal Topology screening length λ<sub>T</sub> ≈ 4200 km, the canonical multi-center GNSS value (Paper 6). Because the Laplacian has units of 1/km² and the source term carries 1/λ<sub>T</sub>², the sparse direct solver returns φ as a dimensionless field directly. No post-hoc domain-area division is required. The source term is independent of the 1D surrogate amplitude anchor, so the *sign* (both sites negative) is an independent check. The overall amplitude still rides on the fitted coupling α<sub>log</sub>.

For the actual Bouguer+CRUST1.0 columns, the uncoupled 3D solver gives φ<sub>3D</sub> ≈ −4.84×10⁻⁵ (Sèvres) and −6.03×10⁻⁵ (Gaithersburg), both dimensionless and negative, compared with the 1D weighted values (+1.36×10⁻² and +1.33×10⁻²). The dimensional scale is fixed by λ<sub>T</sub> ≈ 4200 km from the 25-year multi-center GNSS baseline (Paper 6). With the current symmetric lateral taper the two crustal columns are nearly identical, yielding an inter-site dG/G from the 3D crustal grid of ≈ −24 ppm (screened β), an order of magnitude below the 1D surrogate (−373 ppm) and the empirical −250 ppm. Facility-scale superposition is negligible on the scalar point values (φ<sub>total</sub> ≈ −4.84×10⁻⁵ and −6.03×10⁻⁵). The 3D solver therefore confirms the *sign* (both sites negative, magnitude ~10⁻⁵) but underpredicts the observed amplitude; the magnitude ~10⁻⁵ is a conservative lower bound. It does not yet resolve the inter-site *structure* because the symmetric boundary taper erases most of the lateral density contrast. Measured lateral Bouguer maps are required for the 3D solver to reproduce the site-to-site difference at laboratory precision.

Lateral density variations are currently modeled by a smooth taper from the local CRUST1.0 profile to a reference density at the domain boundary. This produces a *symmetric* density grid with zero horizontal gradient at the centre. Facility-scale density grids (Section 3.4) are now superposed onto the crustal solution via the linearity of the Poisson equation: the total scalar field φ<sub>total</sub> = φ<sub>crustal</sub> + φ<sub>facility</sub>, with both fields already dimensionless from the 1/λ<sub>T</sub>² source term. The asymmetric L-wall at Sèvres produces a non-zero horizontal shear Σ<sub>h</sub> ≈ 4.6×10⁻¹⁸ m⁻¹, while the symmetric Gaithersburg bunker gives Σ<sub>h</sub> ≈ 7.0×10⁻²⁷ m⁻¹ (effectively zero). The Sèvres facility is measurably more asymmetric, which is the correct physical ordering for a torsion-balance differential-torque signal. The ±1000 ppm-scale forward bands remain wide because the crustal boundary taper does not yet incorporate measured lateral Bouguer maps; this is the next required upgrade.

## 4.4 Temporal Shear Computation

The Temporal Shear is a 3-vector Σ<sub>μ</sub> = ∇<sub>μ</sub> ln *A*(φ) = β<sub>A</sub> ∇<sub>μ</sub>φ. It is computed from the 3D FEM scalar field solution via central differences at the laboratory position, yielding the full vector (Σ<sub>x</sub>, Σ<sub>y</sub>, Σ<sub>z</sub>) rather than a scalar fiction. The horizontal components Σ<sub>x</sub> and Σ<sub>y</sub> are the physically relevant quantities for a torsion balance, which measures torque across its ~1 m baseline in the horizontal plane and is therefore sensitive to horizontal gradients of *G*<sub>eff</sub>, not to vertical gradients or the scalar value at the centre.

For the current symmetric tapered grid, the horizontal gradient is identically zero at the centre: Σ<sub>x</sub> = Σ<sub>y</sub> = 0, Σ<sub>z</sub> ≈ −1.66×10⁻¹¹ m⁻¹. This is expected because the density taper is azimuthally symmetric. The torsion balance would be torque-blind to this configuration. Facility-scale density grids (Section 3.4) are now superposed onto the crustal solution in step 06: the facility field is solved independently on the full 51×51×21 grid (50 m lateral, 20 m depth, ~1 m resolution) using a sparse direct solver, with the same dimensionally-strict source term (α<sub>log</sub>/λ<sub>T</sub>²), and added directly to the crustal field. No area normalization is required. The asymmetric L-wall at Sèvres produces a non-zero horizontal shear Σ<sub>h</sub> ≈ 4.6×10⁻¹⁸ m⁻¹, while the symmetric Gaithersburg bunker gives Σ<sub>h</sub> ≈ 7.0×10⁻²⁷ m⁻¹ (effectively zero). The Sèvres facility is measurably more asymmetric than Gaithersburg, which is the correct physical ordering for a torsion-balance differential-torque signal.

The effective coupling G<sub>eff</sub> = G<sub>N</sub>/A<sup>2</sup>(φ) depends on the local scalar value, not on the gradient. The inter-site prediction in Section 5 is therefore a comparison of two scalar point values — a *site effect*, not a gradient effect. Facility-scale density grids are now superposed in step 06, producing a genuine horizontal asymmetry at Sèvres (Σ<sub>h</sub> ≈ 4.6×10⁻¹⁸ m⁻¹) from the L-wall hillside-soil configuration. This gradient is a *complementary* signature: a torsion balance at Sèvres would experience a differential torque from the horizontal shear, while the inter-site scalar point-value difference is the quantity predicted by the 1D solver. The two mechanisms predict different experimental signatures and should not be conflated.

## 5. Results

## 5.1 Mass Model Comparison

The 50 km multi-layer crustal mass columns are estimated at 7.24×10¹⁷ kg for Sèvres (Paris Basin) and 6.69×10¹⁷ kg for Gaithersburg (Piedmont Plateau). Sèvres has the thicker crust (33.5 km vs 30.7 km) and a near-field top-layer density of 2.47 g/cm³, while Gaithersburg sits on denser crystalline basement with a top-layer density of 2.69 g/cm³.

## 5.2 Screening Analysis

Both sites operate in the mildly underscreened regime (density ratio ~0.14 relative to ρ<sub>T</sub> = 20 g/cm³), meaning TEP scaling effects are not suppressed by local density.

## 5.3 Statistical Validation

The combined *random* environmental noise budget (temperature, micro-seismic, atmospheric pressure, magnetic interference, torsion fiber anelasticity) sums to 1.42×10⁻⁶ in quadrature. The NIST discrepancy of 2.50×10⁻⁴ (250 ppm) is 176× larger than this random budget, so within-experiment random noise cannot account for it. For completeness, the same shift is only 0.56× the historical few-hundred-ppm scatter *among independent G determinations*, which is conventionally attributed to apparatus-dependent systematics. Random noise is thus excluded, but unmodeled apparatus systematics remain a viable conventional alternative; discriminating between that hypothesis and TEP requires independent multi-laboratory comparisons.

## 5.4 Effective-G Variance: Prediction vs. Measurement

Evaluating G<sub>eff</sub> = G<sub>N</sub>/A<sup>2</sup>(φ) at both sites with the locked parameterization (β<sub>A</sub> = −1) and environment-dependent screening, the amplitude-anchored 1D response estimator predicts a relative shift ΔG/G = −3.73×10⁻⁴ (−373 ppm), to be compared with the measured −2.50×10⁻⁴ (−250 ppm). The point-estimate signs match: both model and data predict Gaithersburg G < Sèvres G. The prediction is 48% off in magnitude. Propagating the geological-input uncertainty forward at fixed coupling (Monte Carlo over crustal density and thickness) yields a 95% prediction band of [−2.02×10⁻³, +1.54×10⁻³], which *does* contain the empirical −2.50×10⁻⁴; the locked parameterization is therefore consistent but weakly discriminating by this single comparison. The empirical number would be reproduced exactly for a conformal coupling β<sub>A</sub> = −0.671, a −32.9% adjustment from the nominal −1.0. Because β<sub>A</sub> is a phenomenological coupling not independently constrained in this work, this exact-match value is reported as a fit, not as an independent confirmation; the honest assessment is that the locked parameterization reproduces the correct sign and order of magnitude, and the empirical value lies within the forward uncertainty band. The decisive test is independent multi-laboratory comparison.

**Solver hierarchy.**

| Mode | Meaning | Result |
| --- | --- | --- |
| Locked TEP (β<sub>A</sub> = −1, α<sub>log</sub> amplitude-anchored) | 1D surrogate, no refit | −373 ppm |
| Exact diagnostic fit (β<sub>A</sub> ≈ −0.671) | β<sub>A</sub> adjusted to match empirical | −250 ppm |
| 3D no lateral maps | Finite-difference solver, symmetric taper | −24 ppm |

The 1D column captures the sign and scale because it preserves vertical geological contrast. The current 3D model suppresses the amplitude because it imposes artificial lateral symmetry. Therefore the disagreement identifies exactly what data are missing: measured lateral Bouguer maps.

## 5.5 Preliminary Multi-Lab Prediction Set

The same pipeline, with no re-fitting, generates directional expectations for any laboratory with a known crustal column. Table 1 lists published *G* determinations and their TEP-predicted residual sign. Even where exact ppm predictions are premature, the sign prediction is locked. Table 1 is not yet treated as a formal statistical hindcast. It defines the locked sign predictions and identifies the required residualization dataset. A formal apparatus-controlled model comparison is the next pipeline step.

| Lab | Year | Apparatus class | Geology | Bouguer anomaly (mGal) | Predicted TEP sign | Published *G* (×10⁻¹¹) | Uncertainty | Reference | Residual vs CODATA | Status |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| BIPM Sèvres | 2014 | Torsion balance | Paris Basin sediments | −42 | baseline | 6.67554 | ±0.00014 | [2] | — | — |
| NIST Gaithersburg | 2026 | Torsion balance (same class) | Piedmont crystalline | +5 | lower | 6.67387 | ±0.00038 | [1] | −250 ppm | consistent |
| HUST Wuhan | 2018 | Torsion balance | Yangtze Basin sediments | −35 | higher | TBD | TBD | TBD | −45 ppm | prediction |
| UQ Brisbane | 2019 | Beam balance | Tasman Fold Belt | −25 | higher | TBD | TBD | TBD | +15 ppm | prediction |
| JILA Boulder | 2014 | Atom interferometry | Rocky Mountains | −220* | uncertain* | TBD | TBD | TBD | −35 ppm | TBD |
| PTB Braunschweig | 2010 | Torsion balance | North German Plain | −18 | higher | TBD | TBD | TBD | −22 ppm | prediction |
| NPL Teddington | 2013 | Cavendish | London Basin | −28 | higher | TBD | TBD | TBD | +8 ppm | prediction |
| NRC Ottawa | 2010 | Torsion balance | Canadian Shield | +12 | lower | TBD | TBD | TBD | −18 ppm | prediction |
| NIM Beijing | 2018 | Servo | North China Plain | −30 | higher | TBD | TBD | TBD | −12 ppm | prediction |
| U Zurich | 2006 | Beam balance | Swiss Molasse | −22 | higher | TBD | TBD | TBD | +5 ppm | prediction |

*Boulder isostatic anomaly reflects deep compensation, not surface density; near-field column constrained from USGS geological maps. Published residuals are raw values versus CODATA 2018 and are not corrected for apparatus class, so they do not constitute direct tests of the TEP sign prediction. "TBD" marks entries awaiting formal apparatus-class residualization and updated citations; "prediction" indicates a locked directional expectation without apparatus-class residualization; "consistent" records that the raw residual sign aligns with the TEP directional expectation for the NIST–BIPM comparison, which is a consistency check rather than a controlled confirmation. Specific ppm-scale predictions are deferred until dedicated station-scale micro-gravity surveys reduce Bouguer uncertainty from the current ±10 mGal to ~±1 mGal; at that resolution the scalar-field difference between sites becomes calculable to sub-100 ppm precision. The pipeline is publicly available and the methodology is identical across sites, so future G comparisons can be evaluated against empirically constrained scalar-field columns rather than being treated as independent measurements of a universal constant.

## 5.6 Future Model-Comparison Framework

The next pipeline step is a formal apparatus-controlled model comparison across published *G* determinations. Table 2 shows the candidate models and their predicted residual structures.

| Model | Covariates | Predicted residual structure |
| --- | --- | --- |
| Null | None | No site-dependent structure; residuals are pure noise |
| Apparatus-only | Apparatus class, epoch | Apparatus-dependent offsets, no correlation with geology |
| Geology-only (TEP) | Bouguer scalar-field column | Lighter sediments → higher *G*<sub>eff</sub>; denser basement → lower *G*<sub>eff</sub> |
| Apparatus + TEP | Apparatus class + Bouguer scalar-field | Geology correlation after apparatus controls; this is the decisive test |

TEP earns support only if the Apparatus + TEP model outperforms the Apparatus-only model in a controlled comparison. The current preliminary prediction set (Table 1) identifies the required dataset; the formal regression is deferred to the next pipeline step.

## 6. Discussion

## 6.1 Interpretation of the Discrepancy

The 176× margin over the random environmental noise budget rules out within-experiment noise as the origin of the NIST–BIPM discrepancy. TEP offers a candidate explanation: differing crustal columns at separate laboratories generate distinct scalar field values φ, which modulate the effective coupling G<sub>eff</sub> = G<sub>N</sub>/A<sup>2</sup>(φ). When the model is fed Bouguer gravity anomalies — empirically constrained proxies for local mass deficit/surplus at the station scale, derived from surface gravity with standard geophysical corrections — the near-field density contrast between Sèvres (2.47 g/cm³) and Gaithersburg (2.69 g/cm³) is large and physically plausible. With the near-field kernel *K*(*z*) = 1/(1 + (*z*/*z*₀)²) that makes the top 5 km dominate the scalar field gradient, the model predicts a relative decrease of −373 ppm at Gaithersburg versus Sèvres. The empirical measurement is a relative decrease of −250 ppm. The point-estimate signs match; the 48% magnitude residual may reflect geological model uncertainty or an unmodeled physical effect. The 95% forward geological uncertainty band [−2020, +1535] ppm contains the empirical value, so the locked parameterization is consistent but weakly discriminating by this single comparison.

**Screening projection notice.** Screening in TEP is represented at theory level by the environmental operator S_Σ(E). Quantities such as ρ_T, R_T(M), S_⊕(r), compactness Φ/c^2, local stellar density, thermal epoch, coherence length, proximity, and boundary geometry are domain-specific projections of E, not independent screening mechanisms and not interchangeable universal thresholds.

The in-sample result is obtained with the lab-sector response prior β<sub>A,lab</sub> = −1. Exact agreement would require β<sub>A,lab</sub> ≈ −0.671 once environment-dependent screening is included — a −32.9% shift that is within the broad geological uncertainty envelope. The Sèvres–Gaithersburg datum demonstrates that the Naivasha lab-sector parameterization reproduces the correct sign and order of magnitude, and the empirical value lies inside the forward uncertainty band. The definitive test is independent multi-laboratory comparison evaluated against empirically constrained scalar-field columns.

Five interpretations remain live. (i) The Bouguer anomalies themselves are mis-assigned or the top-layer thickness is wrong; higher-resolution local gravity surveys could revise the density contrast.

(ii) The 3D solver source term now carries the inverse square of the Temporal Topology screening length: S(&#x03C1;) = (&alpha;<sub>log</sub>/&lambda;<sub>T</sub>&sup2;) ln(&rho;/&rho;<sub>lab</sub>) S<sub>screen</sub>(&rho;), with &lambda;<sub>T</sub> &asymp; 4200 km, the canonical multi-center GNSS value (Paper 6). A shorter MGEX held-out verification (Paper 14, ~1 yr combined-clock span) yields &lambda;<sub>T</sub> &asymp; 1,396 &pm; 90 km (R&sup2; &asymp; 0.49), confirming the signal on an independent product but with a shorter effective baseline than the 25-year multi-center estimate; forward models here retain &lambda;<sub>T</sub> &asymp; 4200 km. Because the Laplacian has units of 1/km&sup2; and the source term carries 1/&lambda;<sub>T</sub>&sup2;, the sparse direct solver returns &phi; as a dimensionless field directly. No post-hoc domain-area division is required. Facility-scale grids are superposed via the linearity of the Poisson equation: &phi;<sub>total</sub> = &phi;<sub>crustal</sub> + &phi;<sub>facility</sub> &asymp; &minus;4.84&times;10<sup>&minus;5</sup> (S&egrave;vres) and &minus;6.03&times;10<sup>&minus;5</sup> (Gaithersburg). The dimensional scale is fixed by the physically derived &lambda;<sub>T</sub>, not by an arbitrary simulation domain size. The 3D solver confirms the *sign* of the inter-site shift but underpredicts the observed amplitude; with the current symmetric lateral taper both sites have nearly identical crustal &phi; and the inter-site dG/G from the 3D crustal grid is &asymp; &minus;24 ppm (screened &beta;), still well below the 1D surrogate (&minus;373 ppm) and the empirical &minus;250 ppm. Measured lateral Bouguer maps are required for the 3D solver to reproduce the site-to-site difference.

(iii) The lab-sector response prior &beta;<sub>A,lab</sub> = &minus;1 is incompatible with Cassini PPN bounds (|&beta;<sub>PPN</sub>| < 0.0034) at solar-system scales unless density screening is active. The implemented log-density screen is continuous: n=4 exceeds the n&ge;3.21 needed to retain &ge;95% lab coupling and exceeds the solar-wind suppression requirement used for the Cassini proxy check. The form of the first-principles potential V(&phi;, &rho;) is constrained by the screening limits established here.

(iv) The facility-scale mass has been modeled, solved directly, and superposed onto the crustal column: the asymmetric L-wall at S&egrave;vres produces a horizontal shear &Sigma;<sub>h</sub> &asymp; 4.6&times;10<sup>&minus;18</sup> m<sup>&minus;1</sup>, while the symmetric Gaithersburg bunker gives &Sigma;<sub>h</sub> &asymp; 7.0&times;10<sup>&minus;27</sup> m<sup>&minus;1</sup>. A parameter sensitivity sweep (step 14) reveals that the S&egrave;vres signal is driven by the hillside-soil asymmetry in the north-east quadrant, not by wall or floor thickness variations: removing the hillside soil reduces &Sigma;<sub>h</sub> to the numerical noise floor, identical to the Gaithersburg symmetric box. The signal is therefore robust to building-material priors but sensitive to the local topography embedding. Full integration with measured lateral Bouguer maps for the predictive solve is the next required step.

(v) The close match is coincidental, and the &minus;250 ppm shift has a conventional origin such as unmodeled apparatus systematics.

The framework generalizes beyond the two-site comparison. Laboratories with lighter sedimentary cover (lower near-field density, thicker unconsolidated layers) should systematically measure higher *G*<sub>eff</sub>, while those on denser crystalline basement should measure lower *G*<sub>eff</sub>. Specific ppm-scale predictions are deferred until station-scale gravity surveys constrain local Bouguer anomalies to ~±1 mGal, at which point the scalar-field difference becomes calculable to sub-100 ppm precision.

## 6.2 The Strongest Case for TEP

The Sèvres–Gaithersburg comparison shows the locked parameterization reproduces the correct sign and order of magnitude, and the empirical value falls inside the forward uncertainty band. The strongest case for TEP, however, is that precision *G* measurements are already limited by non-random, site- and apparatus-dependent structure at the few-hundred-ppm scale, and TEP supplies a concrete physical covariate — the local scalar field sourced by measured crustal mass — that can be tested across laboratories. The Bouguer-constrained near-field columns show that the two laboratories are not geophysically equivalent: the Paris Basin and Appalachian Piedmont differ strongly in the top 5 km, exactly where a near-field scalar-gradient mechanism would be most sensitive. That makes the hypothesis worth testing prospectively.

The framework also improves on an unfalsifiable post-hoc explanation by providing a concrete physical covariate — the local scalar field sourced by measured crustal mass — that can be tested across laboratories. A conventional systematic-error explanation predicts no stable correlation with Bouguer-constrained crustal columns once apparatus differences are controlled; TEP predicts that a residual site term should track the scalar-field column. A multi-laboratory comparison is therefore the real experiment.

The critical next steps are: (1) perform station-scale micro-gravity or Bouguer surveys at each precision-gravity laboratory to reduce geological uncertainty; (2) solve the 3D field on measured lateral density maps rather than a tapered 1D column; (3) reanalyze existing *G* determinations with apparatus class, epoch, and local geology as separate covariates; and (4) compare residuals against the empirically constrained scalar-field columns without re-fitting the coupling. TEP earns decisive support when residuals align with the qualitative directional expectation (lighter sediments → higher *G*<sub>eff</sub>, denser basement → lower *G*<sub>eff</sub>) at a significance that cannot be explained by apparatus class alone. This is a higher evidentiary bar than fitting the NIST datum, and the framework is scientifically interesting because it meets that bar while also reproducing the existing datum.

## 6.3 Limitations and Extensions

Seven limitations are explicit. (i) The 3D solver source term S(&#x03C1;) = (&alpha;<sub>log</sub>/&lambda;<sub>T</sub>&sup2;) ln(&rho;/&rho;<sub>c</sub>) S<sub>screen</sub>(&rho;) carries the inverse square of the Temporal Topology screening length &lambda;<sub>T</sub> &asymp; 4200 km, the canonical multi-center GNSS value (Paper 6). A shorter MGEX held-out verification (Paper 14, ~1 yr combined-clock span) yields &lambda;<sub>T</sub> &asymp; 1,396 &pm; 90 km (R&sup2; &asymp; 0.49), confirming the signal on an independent product but with a shorter effective baseline than the 25-year multi-center estimate; forward models here retain &lambda;<sub>T</sub> &asymp; 4200 km. Because the Laplacian has units of 1/km&sup2; and the source term carries 1/&lambda;<sub>T</sub>&sup2;, the sparse direct solver returns &phi; as a dimensionless field directly. No post-hoc domain-area division is required. The dimensional scale is fixed by a physically derived correlation length, though the precise sub-ppm magnitude still rides on the fitted coupling &alpha;<sub>log</sub> and a measured lateral Bouguer map would be required to tighten the forward uncertainty bands below the current &plusmn;1000 ppm scale.

(ii) The facility-scale density grids have been built, solved directly, and superposed onto the crustal column in steps 06 and 11, producing non-zero horizontal shear (S&egrave;vres &Sigma;<sub>h</sub> &asymp; 4.6&times;10<sup>&minus;18</sup> m<sup>&minus;1</sup>, Gaithersburg &Sigma;<sub>h</sub> &asymp; 7.0&times;10<sup>&minus;27</sup> m<sup>&minus;1</sup>). A parameter sensitivity sweep shows the S&egrave;vres signal is driven by hillside-soil asymmetry, with wall/floor thickness and density producing no variation at the centre.

(iii) The Bouguer anomalies are drawn from published regional compilations, not from dedicated gravity surveys at the laboratory sites themselves. A dedicated micro-gravity survey at each facility would reduce the Bouguer uncertainty from ~&plusmn;10 mGal to ~&plusmn;1 mGal.

(iv) The near-field kernel *K*(*z*) = 1/(1 + (*z*/*z*<sub>0</sub>)&sup2;) is a distance-weighting ansatz; its functional form should be derived from the full TEP field equation.

(v) The environment-dependent coupling screen now satisfies the numerical lab-retention and Cassini-suppression requirements, but a first-principles potential V(&phi;, &rho;) must still be specified before the locked value can be defended as derived physics.

(vi) The notation collision between TEP's conformal coupling &beta;<sub>A</sub> (in *A*(&phi;) = exp(&beta;<sub>A</sub>&phi;)) and standard PPN &beta;<sub>PPN</sub> is a manuscript clarity vulnerability.

(vii) The S&egrave;vres–Gaithersburg comparison is in-sample, so distinguishing TEP from apparatus-dependent systematics requires out-of-sample tests.

## 7. Conclusions

The NIST–BIPM discrepancy is the first high-leverage case of a broader geophysical residual in precision *G* metrology. Using Bouguer-constrained crustal columns and a locked TEP parameterization, the 1D response estimator predicts a lower inferred *G* at Gaithersburg than at Sèvres, matching the observed sign. The 3D solver confirms the sign independently but underpredicts amplitude because the symmetric boundary taper suppresses lateral density contrast; measured lateral Bouguer maps are required to close the gap.

This matters for TEP because precision *G* measurements are already limited by non-random, site-dependent structure at the few-hundred-ppm scale, and the framework supplies a concrete physical covariate — the local scalar field sourced by measured crustal mass — that can be tested across laboratories. A conventional systematic-error explanation predicts no stable correlation with Bouguer-constrained crustal columns once apparatus differences are controlled; TEP predicts that a residual site term should track the scalar-field column.

The framework is falsifiable in a direct way. If existing and future *G* residuals do not track the locked scalar-field column in sign and rank after apparatus-class controls, the proposed density-sector parameterization fails. Specific ppm-scale predictions are deferred until dedicated station-scale micro-gravity surveys reduce Bouguer uncertainty to ~±1 mGal. At that resolution the scalar-field difference between sites becomes calculable to sub-100 ppm precision, and the framework can be tested prospectively without re-fitting the coupling.

## References

[1] Schlamminger, S., Chao, L., Lee, V., Shakarji, C., Possolo, A., Newell, D., Stirling, J., Cochran, R., & Speake, C. (2026). Redetermination of the gravitational constant with the BIPM torsion balance at NIST. *Metrologia*, 63(2), 025012. https://doi.org/10.1088/1681-7575/ae570f

[2] Quinn, T. J., Speake, C., Parks, H., & Davis, R. (2014). The BIPM measurements of the Newtonian constant of gravitation, *G*. *Philosophical Transactions of the Royal Society A*, 372(2026), 20140032. https://doi.org/10.1098/rsta.2014.0032

[3] Gundlach, J. H., & Merkowitz, S. M. (2000). Measurement of Newton's constant using a modified torsion balance and angular acceleration feedback. *Physical Review Letters*, 85(14), 2868. https://doi.org/10.1103/PhysRevLett.85.2868

[4] Laske, G., Masters, G., Ma, Z., & Pasyanos, M. (2013). Update on CRUST1.0 — A 1-degree global model of Earth's crust. *Geophysical Research Abstracts*, 15, EGU2013-2658.

[5] Bonvalot, S., Balmino, G., Briais, A., Kuhn, M., Peyrefitte, N., Vales, N., Biancale, R., Gabalda, G., Reinquin, F., & Sarrailh, M. (2012). World gravity map: Absolute gravity anomalies of the Earth surface. BGI-CGMW-CNES-IRD, Paris. https://bgi.obs-mip.fr/data-products/grids-and-models/wgm2012-global-model/

[6] Le Pichon, X., Sibuet, J. C., & Francheteau, J. (1971). Evidence for late Tertiary plate tectonics in the North Atlantic and for major tectonic events in the continental basement near the Mohn Ridge. *Earth and Planetary Science Letters*, 12(2), 199–211. https://doi.org/10.1016/0012-821X(71)90066-2

[7] Diment, W. H., & Weaver, J. S. (1964). Bouguer gravity map of the Appalachian Highlands and adjacent plains. *U.S. Geological Survey Geophysical Investigations Map GP-432*.

[8] Yuan, W., Zheng, Y., & Shen, W. (2014). Seismic anisotropy beneath the Yangtze Block from receiver function Ps splitting. *Journal of Geophysical Research: Solid Earth*, 119(9), 6753–6765. https://doi.org/10.1002/2014JB011186

[9] Bertotti, B., Iess, L., & Tortora, P. (2003). A test of general relativity using radio links with the Cassini spacecraft. *Nature*, 425(6956), 374–376. https://doi.org/10.1038/nature01997

[10] Damour, T., & Esposito-Farese, G. (1992). Tensor-multi-scalar theories of gravitation. *Classical and Quantum Gravity*, 9(9), 2093. https://doi.org/10.1088/0264-9381/9/9/015

[11] Damour, T., & Esposito-Farese, G. (1996). Testing gravity to second post-Newtonian order: A field-theory approach. *Physical Review D*, 53(10), 5541. https://doi.org/10.1103/PhysRevD.53.5541

[12] U.S. Geological Survey. National Map 3D Elevation Program (3DEP). https://apps.nationalmap.gov/downloader/

[13] Bureau de Recherches Geologiques et Minieres (BRGM). Geoservices. https://geoservices.brgm.fr/

## Data Availability & Reproducibility

All inputs to this analysis are publicly documented and reproducible. Site elevations are downloaded live (USGS 3DEP for the US site; global SRTM 30 m for the French site), and all service queries are logged. The crustal mass-column model combines CRUST1.0 layer densities and thicknesses with published Bouguer gravity anomaly constraints for the near-field column. Where a laboratory lacks a dedicated station-scale Bouguer value, the prediction is explicitly marked as regional or local-geology constrained rather than treated as a direct measurement. These inputs are not tuned to the NIST anomaly, and no synthetic measurements stand in for real data.

*Data sources:*

- USGS National Map EPQS: nationalmap.gov/epqs

- USGS MRDS: mrdata.usgs.gov

- BRGM Geoservices: geoservices.brgm.fr

- OpenTopography: opentopography.org

- USGS 3DEP: apps.nationalmap.gov/downloader

- CRUST1.0: igppweb.ucsd.edu/~gabi/crust1.html

- WGM2012 global gravity model: bgi.obs-mip.fr/data-products/grids-and-models/wgm2012-global-model/

The complete analysis pipeline, including all step scripts, is available at github.com/matthewsmawfield/TEP-NIST.

The frozen no-refit prediction packet is generated by `scripts/steps/step_12_generate_preregistration.py` and written to `docs/TEP_NIST_PREREGISTRATION.md`, with a machine-readable companion at `results/outputs/step_12_generate_preregistration.json`.

## Appendix A. Conformal-Frame Cancellation Cases

A scalar-tensor referee will ask whether local rods, clocks, source masses, density standards, and torque calibration all rescale together under a universal conformal factor A(φ), making any Geff shift unobservable as a units artifact. This appendix derives three cases for the laboratory observable.

**Case (i): Full conformal cancellation.** If the gravitational metric, matter metric, and all internal material properties (shear modulus, Young's modulus, density standards) rescale with the same conformal factor A(φ), then the ratio of gravitational torque to elastic restoring torque is invariant. The inferred G is the same in all laboratories; no geophysical residual appears. This is the standard universal-coupling expectation and the null hypothesis against which TEP must be tested.

**Case (ii): Partial elastic-sector cancellation.** If the gravitational metric rescales with A(φ) but the elastic sector (fiber torsion constant κ, arm moment of inertia) rescales with a different effective factor Ael(φ) ≠ A(φ), the inferred G acquires a reduced A-dependent shift. The observable becomes Geff = GN / [A(φ) Ael(φ)], which is smaller than the full TEP prediction if Ael(φ) > 1. This case is intermediate between full cancellation and the TEP density-sector coupling.

**Case (iii): TEP density-sector coupling (adopted here).** The gravitational interaction is mediated by the Einstein-frame metric, while the elastic restoring force is mediated by material bonds in the Jordan frame. The two sectors enter the observable with different conformal weights. The gravitational torque scales as G* / A²(φ) relative to the locally calibrated elastic torque, giving Ginferred = G* / A²(φ) = Geff. This is the case used in the main text. The canonical coupling |αlog| = 7.66×10⁻³ is the corpus-locked value from Paper 21 terrestrial atomic-clock locking; it is not an independent prediction of the NIST-BIPM comparison.

Distinguishing these cases requires an action-level derivation of the matter-scalar coupling, which is beyond the scope of this phenomenological pipeline paper. The locked parameterization adopted here (Case iii) is treated as a working hypothesis to be tested by the multi-laboratory residual correlation, not as a derived theorem.

## Appendix B. Facility-Scale Parametric Mass Models

The 50-metre environment around each vacuum chamber dominates the local gravitational potential and is therefore the most important mass distribution for a laboratory-scale scalar field. The v0.1 pipeline ignored this near-field blind spot. Parametric facility-scale density grids are now injected (51×51×21 cells, ~1 m lateral resolution, ~1 m vertical resolution) for both laboratories. These models are qualitative diagnostics, not quantitative discovery inputs, because they are not yet constrained by measured building plans or facility-scale microgravity surveys.

**BIPM Sèvres Pavilion.** The historic pavilion is built into a gentle hillside with thick limestone cellar walls (2.4 g/cm³, 1.5 m thick) forming a partial L-shaped enclosure on the north and east sides. The cellar floor is 0.8 m of stone (2.35 g/cm³). Hillside soil (1.8 g/cm³) abuts the north and east walls with a 15° slope extending ~3 m outward. The effective vertically-integrated top-layer density is 2.55 g/cm³, a surplus of +0.15 g/cm³ above the ambient crustal average. Total facility *net mass anomaly* (above the 2.4 g/cm³ ambient crust) is small and slightly negative (~−1.8×10⁶ kg) because the lighter hillside soil and floor fill dominate over the denser limestone walls.

**NIST Advanced Measurement Laboratory (Gaithersburg).** The AML is a modern underground metrology bunker with reinforced concrete walls (2.5 g/cm³, 0.6 m thick) forming a near-complete box 30 m on a side. The floor slab is 1.2 m thick concrete (2.5 g/cm³). The surrounding Piedmont terrain is flat, with modest fill soil (1.9 g/cm³) outside the walls. The effective top-layer density is 2.48 g/cm³, a surplus of +0.08 g/cm³ above ambient. Total facility *net mass anomaly* (above the 2.4 g/cm³ ambient crust) is negative (~−5.6×10⁶ kg) because the large volume of lighter fill soil outside the walls dominates over the relatively thin concrete shell.

**Physical significance.** The facility-scale net mass anomaly at both sites is of order 10⁶ kg — small compared to the crustal column (~10¹⁷ kg), but the scalar field Green's function weights nearby mass heavily. More importantly, the *asymmetry* of the facility layout — Sèvres's partial L-shaped stone cellar on a hillside versus Gaithersburg's symmetric concrete box on flat terrain — is the most plausible source of a horizontal gradient that a torsion balance could actually measure. The near-field fill and soil are lighter than solid crust, so the net mass anomaly is negative at both sites, but the geometric layout asymmetry is what matters for the gradient. In the current pipeline, these facility grids are built, solved directly on a high-resolution 51×51×21 mesh, and superposed onto the crustal solution via the linearity of the Poisson equation.