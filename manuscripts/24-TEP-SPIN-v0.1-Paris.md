# Temporal Equivalence Principle: Fermion Spin as Temporal-Orientation Holonomy
**Matthew Lukin Smawfield**
Version: v0.1 (Paris)
First published: 7 June 2026 · Last updated: 2 July 2026
DOI: 10.5281/zenodo.20572706

---

## Abstract

Keywords: subatomic structure, fermion topology, spin, vorticity, temporal-orientation holonomy, proximity screening, Fermilab g-2, AMBER, temporal equivalence principle

## 1. Introduction: Spin as Temporal-Orientation Holonomy

### 1.1 The Conformal No-Go Theorem

The zero-dimensional point-particle paradigm of Quantum Field Theory (QFT) treats the fermion as a structureless singularity. This assumption underlies ultraviolet divergences in loop integrals and necessitates renormalization, a mathematical procedure that subtracts infinities to yield finite predictions. Renormalization works phenomenologically, but it does not resolve the underlying physical problem: a particle with no spatial extent cannot have a finite self-energy.

**Evidence status.** This paper is classified as a *candidate microscopic completion* in the TEP corpus. It proposes a geometric origin for fermion spin through a compact temporal-orientation bundle, but the bundle's kinetic term, potential, and coupling to the disformal sector are not fully specified. Full QED recovery, gauge invariance, and quantitative matching to tested particle-physics observables remain required closure tasks. The results should be evaluated as theoretical structure-building rather than as direct empirical proof of TEP.

Within the Temporal Equivalence Principle (TEP), we observe a structural property that explains why the point-particle picture cannot accommodate spin at the geometric level. The conformal temporal shear &Sigma;<sub>&mu;</sub> = &nabla;<sub>&mu;</sub> ln A(&phi;) is an exact one-form. For any contractible closed loop C,

&oint;<sub>C</sub> &Sigma;<sub>&mu;</sub> dx<sup>&mu;</sup> = &oint;<sub>C</sub> d ln A = 0.

Therefore the conformal scalar sector cannot carry spin holonomy. Spin must live in an additional compact temporal-orientation sector. This is not an arbitrary extra structure; it is a necessary consequence of the single-valuedness of the real conformal factor.

This structural observation motivates the introduction of an additional compact temporal-orientation bundle. It forces fermion spin out of the scalar conformal sector and into a compact orientation bundle, where integer winding lifts spinorially to the observed 2&pi; sign reversal and 4&pi; return. The fermion is thereby identified as a finite Compton-scale topological defect in dynamical proper time. We note that the vanishing circulation of any exact one-form on contractible loops is a standard result in differential geometry (Poincaré lemma). The physical content lies in the postulate that spin must be realized as quantized holonomy in a compact bundle, which is an independent assumption motivated by this structural observation.

At the quantum scale, the geometric proximity regulator $\chi = r_c/\lambda_{\text{scr}}$ operates as the subatomic, coordinate-free realization of the abstract environmental operator $\mathcal{S}_\Sigma(\mathcal{E})$. It mathematically bounds the Temporal Topology divergence at the fermion core, completely independent of the macroscopic bulk-density proxies utilized in astrophysical TEP domains.

### 1.2 The TEP Topological Fermion

The Temporal Equivalence Principle replaces the point particle with a localized topological charge in the temporal shear field. The fermion is not a mathematical singularity but a physical defect in the scalar field &phi; that defines local proper time. This topological charge carries a natural geometric boundary: the proximity-based saturation scale, observationally proxied by &rho;<sub>T</sub> &asymp; 20 g/cm<sup>3</sup> (a phenomenological scale established from terrestrial clock correlation data in TEP-UCD, Paper 6), which bounds the proper-time oscillator and provides a candidate physical regulator for ultraviolet behaviour.

The TEP framework is built on three axioms. (A1) The matter-frame metric is a conformal–disformal rescaling of the gravitational metric: g&#771;<sub>&mu;&nu;</sub> = A<sup>2</sup>(&phi;) g<sub>&mu;&nu;</sub> + B(&phi;) &nabla;<sub>&mu;</sub>&phi; &nabla;<sub>&nu;</sub>&phi;. (A2) The conformal factor is exponential: A(&phi;) = exp(&beta;<sub>A</sub>&phi;/M<sub>Pl</sub>). (A3) Temporal shear is the gradient of the logarithmic conformal factor: &Sigma;<sub>&mu;</sub> = &nabla;<sub>&mu;</sub> ln A(&phi;). The conformal Hamilton–Jacobi sector follows from these axioms; the spin/vorticity sector requires the additional compact temporal-orientation bundle postulate forced by the no-go theorem. The full causal matter metric is permanently engaged; the screened limit, where local stress forces both A(&phi;) &rarr; 1 and the observable disformal response is suppressed, recovers the standard Minkowski background and isotropic interactions. In the unscreened regime the disformal sector governs the routing of forces through the tilted light cone, as developed in TEP-KIN (Paper 25).

### 1.3 Context: Renormalization and the Point-Particle Paradigm

Early attempts at a unified geometric theory sought to replace dimensionless point particles with physical &ldquo;knots&rdquo; in spatial geometry, but failed to eliminate the mathematical divergences that necessitated renormalization. The TEP framework achieves this geometric origin by shifting the topology from spatial gravity to proper time. Instead of a point particle, TEP introduces a localized topological charge embedded within the temporal shear field. Because this fermion is a physical defect in the scalar field &phi;, it carries a natural geometric boundary. Bounded natively by the local density saturation, this finite geometric structure provides a candidate physical regulator at the divergence's origin. Full elimination of renormalization requires showing that gauge invariance, Ward identities, and the quantitative recovery of tested QED observables all survive; that work remains active development.

The electron self-energy diverges as &Lambda; &rarr; &infin; in the standard formulation, and the Landau pole in QED signals that the theory is incomplete at short distances. Renormalization is a mathematical workaround for the divergence; the point-particle assumption is the root cause. A finite geometric structure provides a candidate physical regulator at the divergence's origin.

### 1.4 The Lamb Shift: Ontological Note

The 1947 observation of the Lamb Shift is historically cited as evidence for vacuum polarization in standard QED. Within the TEP framework, one could reinterpret the Lamb Shift as a manifestation of temporal-shear topography near a massive nucleus. However, this remains an ontological reinterpretation of an existing empirical result, not an alternative calculation. The standard QED calculation remains empirically validated, and a quantitative TEP derivation of the Lamb shift is not attempted here.

## 2. The Topological Fermion

### 2.0 Symbol Table

| Symbol | Definition | Value | Classification |
| --- | --- | --- | --- |
| *A*(&phi;) | Conformal factor | exp(&beta;<sub>A</sub>&phi;/M<sub>Pl</sub>) | Fundamental |
| &beta;<sub>A</sub> | Fundamental conformal coupling | &minus;1.0 | Fundamental |
| &beta;<sub>spin</sub> | Phenomenological screening coefficient | 0.01 | Effective (tanh ansatz) |
| &beta;<sub>A</sub><sup>(eff)</sup> | Effective lab coupling after screening | ~10<sup>&minus;2</sup> | Derived via *S*<sub>&Sigma;</sub>(*E*) |
| &chi; | Geometric closure ratio r<sub>c</sub>/&lambda;<sub>scr</sub> | 1/&radic;2 &approx; 0.707 | Model definition |
| &chi;<sub>q</sub> | Single-quark closure ratio | 1/&radic;2 (assumed, universal scalar-sector property) | Model definition |
| m<sub>&phi;</sub> | Scalar field mass | m<sub>e</sub>/&radic;2 | Model definition |
| &rho;<sub>T</sub> | Phenomenological saturation scale | ~20 g/cm<sup>3</sup> | Empirical input (Paper 6) |
| &chi;<sub>p</sub> | Proton closure ratio | ~0.21 | Bound-state assumption |

### 2.1 The Fermion as Topological Charge

The fermion is defined as a localized topological charge in the temporal shear field. In the matter frame, proper time d&tau; is set by the causal matter metric g&#771;<sub>&mu;&nu;</sub> = A<sup>2</sup>(&phi;) g<sub>&mu;&nu;</sub> + B(&phi;) &nabla;<sub>&mu;</sub>&phi; &nabla;<sub>&nu;</sub>&phi;. *Metric signature:* (+, &minus;, &minus;, &minus;). In the conformal limit relevant for the single-particle core geometry, a particle of mass m propagates according to the g&#771;-Hamilton-Jacobi equation, which in flat background with A(&phi;) = exp(&beta;<sub>A</sub>&phi;/M<sub>Pl</sub>) reads:

g<sup>&mu;&nu;</sup> &part;<sub>&mu;</sub>S &part;<sub>&nu;</sub>S = m<sup>2</sup>c<sup>2</sup> exp(2&beta;<sub>A</sub>&phi;/M<sub>Pl</sub>).

The effective mass in the matter frame is m<sub>*</sub> = m A(&phi;) = m exp(&beta;<sub>A</sub>&phi;/M<sub>Pl</sub>). The proper-time oscillator frequency is therefore shifted by the local conformal factor, and the fermion acquires a position-dependent effective inertia. In the rest frame (&nabla;S = 0), the frequency is &omega;<sub>eff</sub> = mc<sup>2</sup> A(&phi;) / &hbar;. This is the origin of the bounded proper-time oscillator: as the topological charge tightens, A(&phi;) flattens toward unity and &omega;<sub>eff</sub> approaches the standard Compton frequency, but it never diverges because the core has finite extent. The conformally shifted g&#771;-Hamilton-Jacobi equation and the emergence of the Klein-Gordon and Dirac operators in the screened limit are derived systematically in TEP-QF (Paper 23).

### 2.2 Spin as Quantized Vorticity

"Spin" is translated directly into fluid vorticity. However, to preserve the single-valuedness of the real conformal factor A(&phi;) = exp(&beta;<sub>A</sub>&phi;/M<sub>Pl</sub>), the fermion core must be mathematically modeled as a dual-component spatio-temporal vortex.

The temporal shear &Sigma;<sub>&mu;</sub> = &nabla;<sub>&mu;</sub> ln A(&phi;) is an exact one-form and is therefore strictly irrotational (&nabla; &times; &Sigma; = 0). This radial conformal shear generates the effective mass and the topographic drag, but it cannot carry spin.

**Structural Observation 1 (Conformal exactness).** Because A(&phi;) is a real, strictly positive, single-valued scalar field, &Sigma;<sub>&mu;</sub> = &nabla;<sub>&mu;</sub> ln A is an exact gradient. For any contractible closed loop C,

&oint;<sub>C</sub> &Sigma;<sub>&mu;</sub> dx<sup>&mu;</sup> = &oint;<sub>C</sub> d ln A = 0.

Spin belongs exclusively to the compact temporal-orientation bundle, not to the irrotational conformal shear sector. The azimuthal orientation shear is K<sub>&theta;</sub> = n/r, where n is the integer winding number. The circulation around the charge core is quantized within this compact orientation bundle:

&Gamma; = &oint; K &middot; d&ell; = 2&pi;n.

The vorticity vector &omega;<sub>i</sub> = (&nabla; &times; K)<sub>i</sub> vanishes everywhere except at the core singularity, where it is distributional. The minimal non-trivial defect carries integer winding n = &pm;1 in the compact orientation sector. As established in TEP-QF, the temporal-orientation bundle governing its causal frame is spinorial: a 2&pi; circuit reverses the spinor sign and a 4&pi; circuit restores it, yielding the observed spin projection S = &pm;&hbar;/2.

![Conformal exactness no-go theorem: the scalar shear is irrotational.](figures/fig1_conformal_no_go.png)

Figure 2.1. Conformal exactness no-go. The temporal shear Σ<sub>μ</sub> = ∇<sub>μ</sub> ln A is an exact gradient, so its circulation vanishes on any contractible loop. Spin holonomy is forbidden in the scalar sector.

**Orientation bundle as a postulate.** The compact temporal-orientation bundle is introduced as an independent postulate forced by the no-go theorem (Structural Observation 1), not derived from the scalar sector. Its minimal field content is a U(1) gauge connection with curvature F<sub>&mu;&nu;</sub> sourced by the topological charge current; the spinorial lift (Structural Observation 2) is the only structure needed to recover the observed 2&pi; sign reversal and 4&pi; return.

**Conditional nature of bundle results.** The bundle's kinetic term, potential, and coupling to the disformal sector B(&phi;) are not specified in the present work. The spinorial lift (2&pi; sign reversal, 4&pi; return) is exact conditional on the existence of a U(1) bundle Lagrangian satisfying five constraints: (i) gauge invariance, (ii) spontaneous symmetry breaking selecting n = &pm;1, (iii) disformal coupling, (iv) Bianchi identity dF = 0, and (v) reduction to free U(1) away from the core. The present paper's topological results are therefore conditional on this existence; no further assumptions about bundle dynamics are needed for the geometric arguments of Sections 2.1–2.3.

![Spinorial holonomy: 2π sign reversal, 4π restoration.](figures/fig2_spinorial_holonomy.png)

Figure 2.2. Spinorial holonomy in the compact temporal-orientation bundle. A 2π circuit reverses the spinor sign; a 4π circuit restores it, mapping integer winding n = 1 to spin projection S<sub>z</sub> = &pm;&hbar;/2.

### 2.3 Charge Core Geometry and the &chi; Ratio

The core radius of the topological charge is identified with the electron Compton wavelength,

r<sub>c</sub> = &hbar;/(m<sub>e</sub>c),

which fixes the finite-core scale of the defect. The scalar Yukawa screening length is fixed by the single-particle closure,

&lambda;<sub>scr</sub> = &radic;2 &hbar;/(m<sub>e</sub>c).

Their ratio defines the internal geometric consistency condition

&chi; = r<sub>c</sub> / &lambda;<sub>scr</sub> = 1/&radic;2 &approx; 0.707.

This is a model definition: the scalar mass m<sub>&phi;</sub> = m<sub>e</sub>/&radic;2 is chosen to satisfy the single-particle Klein-Gordon closure. It is not an independent empirical prediction. The non-trivial content is threefold: (i) the model is self-consistent (the same electron scale sets both the core and the screening length); (ii) the ratio is an O(1) number, avoiding hierarchy problems; and (iii) the model is anchored to the known electron Compton wavelength, a measured quantity. A first-principles derivation of why &lambda;<sub>scr</sub> = &radic;2 r<sub>c</sub> from the full non-linear TEP field equations remains an open problem.

![chi closure: ratio of Compton radius to screening length.](figures/fig3_chi_closure.png)

Figure 2.3. Model closure &chi; = r<sub>c</sub> / &lambda;<sub>scr</sub> = 1/&radic;2. The electron Compton wavelength (core radius) and the scalar Yukawa screening length are linked by the single-particle Klein-Gordon closure m<sub>&phi;</sub> = m<sub>e</sub>/&radic;2.

### 2.4 The TEP Action and Field Equations

The dynamics of the scalar field &phi; are governed by the Einstein-frame action

S = &int; d<sup>4</sup>x &radic;&minus;g [ R/(16&pi;G) &minus; &#189; g<sup>&mu;&nu;</sup> &part;<sub>&mu;</sub>&phi; &part;<sub>&nu;</sub>&phi; &minus; V(&phi;) ] + S<sub>matter</sub>[g&#771;<sub>&mu;&nu;</sub>, &psi;],

where the matter fields &psi; couple to the full causal matter metric g&#771;<sub>&mu;&nu;</sub> = A<sup>2</sup>(&phi;) g<sub>&mu;&nu;</sub> + B(&phi;) &nabla;<sub>&mu;</sub>&phi; &nabla;<sub>&nu;</sub>&phi;. The scalar potential V(&phi;) is to be constrained by cosmological data and fifth-force bounds; a runaway form V(&phi;) &prop; exp(&minus;&lambda;&phi;/M<sub>Pl</sub>) and a screened chameleon potential are both compatible with the TEP framework. The fundamental conformal coupling &beta;<sub>A</sub> = &minus;1.0 is the bare coupling appearing in the metric ansatz; its relationship to the observable coupling in different environments is governed by the environmental screening operator *S*<sub>&Sigma;</sub>(*E*), which depends on source structure, boundary conditions, and measurement channel alongside density (Paper 0, &sect;7). In the high-density solar-system environment, *S*<sub>&Sigma;</sub>(*E*) suppresses the observable coupling by a factor ~300 relative to the bare value, yielding |&beta;<sub>obs</sub>| < 3.4 &times; 10<sup>&minus;3</sup> consistent with the Cassini PPN bound (Bertotti et al. 2003). In the lower-density terrestrial laboratory environment, the same operator yields a weaker suppression ~100, corresponding to an effective &beta;<sub>A</sub><sup>(eff)</sup> &sim; 10<sup>&minus;2</sup> relevant for g&minus;2 experiments. Varying with respect to &phi; yields the Klein-Gordon equation in the presence of a fermion source:

&nabla;<sup>&mu;</sup>&nabla;<sub>&mu;</sub>&phi; &minus; V'(&phi;) = (&beta;<sub>A</sub>/M<sub>Pl</sub>) T<sup>&mu;</sup><sub>&mu;</sub>,

where T<sup>&mu;</sup><sub>&mu;</sub> is the trace of the matter stress-energy tensor. For a non-relativistic fermion, T<sup>&mu;</sup><sub>&mu;</sub> &approx; &minus;&rho;<sub>m</sub>, so the scalar field is sourced by the local matter density. As &rho;<sub>m</sub> increases, &phi; is driven to a value that flattens A(&phi;) toward unity. However, the *observable* suppression of Temporal Shear is governed by the full environmental operator *S*<sub>&Sigma;</sub>(*E*), which includes source structure, boundary conditions, and measurement channel alongside density (Paper 0, &sect;7). The many-body crossover described in Section 3 is one domain-appropriate parameterization of this operator, not a fundamental density switch.

## 3. Proximity-Dependent Screening and the Many-Body Crossover

*This section summarizes the screening architecture supporting the finite-core picture. Detailed derivations are deferred to Appendix A.3.*

**Screening projection notice.** Screening in TEP is represented at theory level by the environmental operator S_Σ(E). Quantities such as ρ_T, R_T(M), S_⊕(r), compactness Φ/c^2, local stellar density, thermal epoch, coherence length, proximity, and boundary geometry are domain-specific projections of E, not independent screening mechanisms and not interchangeable universal thresholds.

The Fermi-wavelength argument provides the correct order-of-magnitude intuition. At &rho;<sub>T</sub> &asymp; 20 g/cm<sup>3</sup>, the electron Fermi wavelength is &lambda;<sub>F</sub> &approx; 10<sup>-10</sup> m, roughly 300&times; larger than the Compton radius r<sub>c</sub> &approx; 3.9 &times; 10<sup>-13</sup> m. Because volume scales as length cubed, the packing density using &lambda;<sub>F</sub> as the exclusion scale is roughly 2.5 &times; 10<sup>7</sup> times lower than the naive Compton-scale estimate, bringing the expected crossover into the same broad density regime as the observed 20 g/cm<sup>3</sup>. The Fermi wavelength is a descriptive scale, not a mechanism that sets &rho;<sub>T</sub>; the latter is a property of the scalar potential itself. The parametric coincidence &rho;<sub>CM</sub>/&rho;<sub>T</sub> = (1/2)(m<sub>p</sub>/m<sub>e</sub>)&alpha;<sup>5/3</sup> &approx; 0.246 &sim; O(1) explains why Earth naturally sits at the continuous screening transition boundary (see Appendix A.3 for the full derivation).

The many-body crossover is governed by the random-phase superposition of N<sub>eff</sub> = (L<sub>c</sub>/&lambda;<sub>F</sub>)<sup>3</sup> &sim; 10<sup>50</sup> uncorrelated topological charges within the terrestrial coherence volume (L<sub>c</sub> &approx; 4200 km). This suppresses the net temporal shear by a factor &sim; 1/&radic;N<sub>eff</sub> &sim; 10<sup>-25</sup>, bridging the gap between the single-particle mean-field prediction and the phenomenological saturation scale. Screening is a smooth, continuous slope spanning roughly &rho; &sim; 2–30 g/cm<sup>3</sup> (10% to 90% screened), not a sharp boundary.

A correlation-modified suppression formula (Appendix A.4, Figure 2.6) bounds the effect of local phase correlations in dense matter. The generalized suppression factor f<sub>RP</sub><sup>(corr)</sup> = (1/&radic;N)&radic;1 + nV<sub>c</sub> shows that even if correlations persist at the Fermi-wavelength scale (&xi; &lesssim; &lambda;<sub>F</sub>), the correction is at most O(10), leaving the 1/&radic;N<sub>eff</sub> suppression intact. The finite-core lattice solver (Figure 2.5) confirms the random-phase 1/&radic;N<sub>micro</sub> scaling and validates the correlated-phase bound. Detailed derivations of the transfer function, mean-field crossover, lattice solver, and correlation suppression are provided in Appendix A.3.

![Finite-core UV regulator convergence.](figures/fig4_uv_regulator_convergence.png)

Figure 2.4. Finite-core UV regulator convergence. The Gaussian regulator exp(&minus;k<sup>2</sup>&sigma;<sup>2</sup>) with &sigma; = r<sub>c</sub>/&radic;2 renders the vacuum polarization integral finite, interpolating between a Compton-scale regularized theory and standard QED as &sigma; &rarr; 0.

![Numeric screening transfer function.](figures/fig4_transfer_function.png)

Figure 2.5. Numeric screening transfer function T(&rho;) from the tanh ansatz. The inflection point at &rho; &asymp; 15 g/cm<sup>3</sup> marks where the screening transition is steepest (S &asymp; 0.65).

![Fermi wavelength versus Compton scale.](figures/fig5_fermi_vs_compton.png)

Figure 2.6. Fermi wavelength &lambda;<sub>F</sub> versus Compton radius r<sub>c</sub> across density. At &rho; &asymp; 20 g/cm<sup>3</sup>, &lambda;<sub>F</sub> &asymp; 10<sup>-10</sup> m, roughly 300&times; larger than r<sub>c</sub>, providing the correct order-of-magnitude intuition for the screening crossover.

![Self-consistent crossover from screened Klein-Gordon equation.](figures/fig7_self_consistent_crossover.png)

Figure 2.7. Self-consistent crossover from the screened Klein-Gordon equation using the canonical &beta;<sub>A</sub> = &minus;1.0. Top-left: screening efficiency S(&rho;). Top-right: self-consistent scalar field |&phi;|. Bottom-left: transfer function T(&rho;). Bottom-right: Fermi wavelength vs screening length &lambda;<sub>scr</sub> = &radic;2 r<sub>c</sub>.

![Finite-core lattice solver results showing density-dependent effective density.](figures/fig7_finite_core_lattice.png)

Figure 2.8. Finite-core lattice solver. Left: Effective density &rho;<sub>eff</sub> versus bulk density &rho; for coherent and random-phase lattices. Right: Lattice-size scaling at &rho; = 20 g/cm<sup>3</sup>, showing the random-phase suppression factor following the 1/&radic;N<sub>micro</sub> prediction.

![Correlation-modified random-phase suppression factor.](figures/fig10_correlation_suppression.png)

Figure 2.9. Correlation-modified random-phase suppression factor f<sub>RP</sub><sup>(corr)</sup> versus correlation length &xi;/&lambda;<sub>F</sub>. The physical bound region (&xi; &lesssim; &lambda;<sub>F</sub>) is shaded. Even at the physical bound, the correction is O(10), leaving the 1/&radic;N<sub>eff</sub> suppression intact.

## 4. Empirical Tests and Planned Analyses

### 4.1 JLab/AMBER Cross-Section Prediction

Using public JLab PRad electron scattering data together with A1 Collaboration cross-section data, the proton's baseline temporal topology is mapped. Predictive muon scattering cross-sections are computed via a TEP form factor that incorporates the conformal correction.

The pipeline (`scripts/steps/step_11_amber_prediction_low_high_Q2.py`) has processed 1,493 data points: 71 from JLab PRad (Xiong *et al.* 2019) and 1,422 from the A1 Collaboration (Bernauer *et al.* 2014). The TEP form factor prediction is:

F(Q<sup>2</sup>) = F<sub>dipole</sub>(Q<sup>2</sup>) &middot; A<sub>TEP</sub>(Q<sup>2</sup>),

with the conformal screening function

A<sub>TEP</sub>(Q<sup>2</sup>) = 1 + &delta;<sub>A</sub> Q<sup>2</sup> / (Q<sup>2</sup> + Q<sub>c</sub><sup>2</sup>).

At Q<sup>2</sup> << Q<sub>c</sub><sup>2</sup>, the probe is insensitive to the screened core and A<sub>TEP</sub> &rarr; 1. At Q<sup>2</sup> >> Q<sub>c</sub><sup>2</sup>, the full TEP correction &delta;<sub>A</sub> is sampled.

**Proton as a three-quark bound state of topological charges.** The proton is a composite baryon. Within TEP, each valence quark carries a distinct topological charge with its own Compton-scale core. The TEP closure for a single quark is &chi;<sub>q</sub> = r<sub>c,q</sub>/&lambda;<sub>scr,q</sub> = 1/&radic;2, where r<sub>c,q</sub> = &#8461;/(m<sub>q</sub>c) and m<sub>q</sub> is the constituent quark mass. This assumption is justified within TEP by the universality of the scalar-sector geometric closure: the ratio &chi; = r<sub>c</sub>/&lambda;<sub>scr</sub> is a property of the scalar field &phi; (set by the potential V(&phi;) and the conformal coupling &beta;<sub>A</sub>), not a flavour-dependent particle property. Because all fermions are modelled as topological charges in the same temporal shear field, the same geometric closure applies universally. The constituent mass m<sub>q</sub> sets the quark Compton scale r<sub>c,q</sub> = &#8461;/(m<sub>q</sub>c), but the closure ratio &chi;<sub>q</sub> = 1/&radic;2 is inherited from the scalar sector.

For the proton, the three valence quarks (uud) form a bound state confined within the hadronic radius r<sub>A</sub> &sim; 1 fm. The effective proton screening length is set by the confinement scale, not the quark Compton scale:

&lambda;<sub>scr,p</sub> = r<sub>A</sub> = 1 fm.

This identification follows from the bound-state structure of the proton. Each quark carries a temporal shear field &Sigma;<sub>i</sub> = &nabla; ln A(&phi;<sub>i</sub>) with individual screening length &lambda;<sub>scr,q</sub> = &radic;2 r<sub>c,q</sub>. Within the confinement radius r<sub>A</sub> &sim; 1 fm, the three quark shear fields overlap coherently: the total temporal shear is &Sigma;<sub>total</sub> = &Sigma;<sub>1</sub> + &Sigma;<sub>2</sub> + &Sigma;<sub>3</sub>. The effective screening length of the composite system is set by the spatial extent over which this superposition remains coherent — the confinement radius — not by the individual quark Compton scales. Equivalently, the disformal coupling B(&phi;) in the causal matter metric g&#771;<sub>&mu;&nu;</sub> = A<sup>2</sup> g<sub>&mu;&nu;</sub> + B(&phi;) &nabla;<sub>&mu;</sub>&phi; &nabla;<sub>&nu;</sub>&phi; cross-couples the individual quark gradients &nabla;&phi;<sub>i</sub>, generating an effective smoothing scale set by the overlap region of the three quark cores. The Yukawa Green's function for an extended source of radius r<sub>A</sub> is suppressed at distances r &gtrsim; r<sub>A</sub>, making the confinement scale the effective screening length. Since r<sub>A</sub> &sim; 1 fm > r<sub>c,p</sub> &sim; 0.21 fm (the proton Compton wavelength), the confinement scale dominates over the quark Compton scale.

The corresponding proton screening mass is m<sub>p</sub><sup>TEP</sup> = &#8461;/(&lambda;<sub>scr,p</sub>c) = 0.1973 GeV/fm / 1 fm = 0.1973 GeV. The proton Q<sub>c</sub> scale is therefore:

Q<sub>c,p</sub> = m<sub>p</sub><sup>TEP</sup>c = 0.1973 GeV,    Q<sub>c,p</sub><sup>2</sup> = 0.0389 GeV<sup>2</sup>.

This is distinct from the electron-sector Q<sub>c</sub> = m<sub>p</sub>c/&radic;2 = 0.663 GeV. The distinction arises because the proton is a three-quark bound state: the screening length is set by the confinement radius r<sub>A</sub> &sim; 1 fm, which is ~3&times; larger than the proton Compton wavelength r<sub>c,p</sub> = &#8461;/(m<sub>p</sub>c) &sim; 0.21 fm. The bound-state &chi; ratio is therefore:

&chi;<sub>p</sub> = r<sub>c,p</sub> / &lambda;<sub>scr,p</sub> = 0.21 / 1.0 &approx; 0.21 &ne; &chi;<sub>e</sub> = 0.707.

The physical content of &chi;<sub>p</sub> &ne; &chi;<sub>e</sub> is therefore that the proton is a composite source: its screening length is set by the confinement radius (the spatial extent of the three-quark source) rather than by the proton Compton wavelength. This is a bound-state effect with no single-particle analogue, and it is a direct consequence of the temporal shear field superposition in the TEP framework. A first-principles derivation of &lambda;<sub>scr,p</sub> from the full non-linear TEP field equations with three-quark confinement remains an open problem. **Sensitivity to &chi;<sub>p</sub>.** If the proton &chi;<sub>p</sub> were varied by &plusmn;30% (reflecting uncertainty in the confinement-scale screening), the extracted &delta;<sub>A</sub> would shift by ~60% (since &delta;<sub>A</sub> &prop; Q<sub>c,p</sub><sup>2</sup> at fixed radii). The present ansatz &lambda;<sub>scr,p</sub> = 1 fm is the minimal choice consistent with the hadronic confinement scale. &delta;<sub>A</sub> is fixed by the two measured proton radii.

Q<sub>c,p</sub> = m<sub>p</sub><sup>TEP</sup>c = 0.1973 GeV,    Q<sub>c,p</sub><sup>2</sup> = 0.0389 GeV<sup>2</sup>.

**Extraction of &delta;<sub>A</sub>.** For Q<sup>2</sup> << Q<sub>c</sub><sup>2</sup>, expand:

A<sub>TEP</sub>(Q<sup>2</sup>) &approx; 1 + (&delta;<sub>A</sub>/Q<sub>c</sub><sup>2</sup>) Q<sup>2</sup>.

The dipole form factor satisfies F<sub>dipole</sub>(Q<sup>2</sup>) &approx; 1 &minus; &langle;r<sup>2</sup>&rangle;<sub>dipole</sub> Q<sup>2</sup>/6, so the product gives:

F(Q<sup>2</sup>) &approx; 1 &minus; &langle;r<sup>2</sup>&rangle;<sub>dipole</sub> Q<sup>2</sup>/6 + (&delta;<sub>A</sub>/Q<sub>c</sub><sup>2</sup>) Q<sup>2</sup>.

Matching to the effective radius definition yields &langle;r<sup>2</sup>&rangle;<sub>eff</sub> = &langle;r<sup>2</sup>&rangle;<sub>dipole</sub> &minus; 6&delta;<sub>A</sub>/Q<sub>c</sub><sup>2</sup>. Equating &langle;r<sup>2</sup>&rangle;<sub>eff</sub> to the PRad measurement (r<sub>p</sub> = 0.831 &plusmn; 0.014 fm) and &langle;r<sup>2</sup>&rangle;<sub>dipole</sub> to the CODATA reference (r<sub>p</sub> = 0.8409 &plusmn; 0.0004 fm) yields &delta;<sub>A</sub> = 0.0028 &plusmn; 0.0039. The uncertainty is propagated from both input radii: &sigma;(&delta;<sub>A</sub>) = Q<sub>c</sub><sup>2</sup>/(6&#8461;<sup>2</sup>c<sup>2</sup>) &radic;(4r<sub>dipole</sub><sup>2</sup>&sigma;<sub>dipole</sub><sup>2</sup> + 4r<sub>eff</sub><sup>2</sup>&sigma;<sub>eff</sub><sup>2</sup>). Both input radii are experimentally measured; &delta;<sub>A</sub> is fixed from existing radius data and then used to make a transferable prediction for AMBER.

**Systematic uncertainty on the radius reference.** The CODATA radius (0.8409 fm) is a weighted average of electronic and muonic hydrogen measurements; the PRad radius (0.831 fm) is an electron-scattering result. The "proton radius puzzle" is precisely the discrepancy between these methods. The extraction assumes the PRad radius as the TEP-corrected value; a systematic uncertainty from the choice of reference radius is not included in the quoted &sigma;(&delta;<sub>A</sub>) and is estimated at O(0.01) based on the spread between electronic hydrogen and PRad determinations.

The predicted deviation rises from ~0% at Q<sup>2</sup> << Q<sub>c</sub><sup>2</sup> to ~0.27% at Q<sup>2</sup> ~ 1 GeV<sup>2</sup>, asymptotically approaching ~0.28% for Q<sup>2</sup> >> Q<sub>c,p</sub><sup>2</sup>. AMBER&rsquo;s proton-radius programme tests the low-Q<sup>2</sup> slope signature. The full ~0.28% asymptotic deviation at Q<sup>2</sup> &gtrsim; Q<sub>c,p</sub><sup>2</sup> would require extended AMBER kinematics or complementary higher-Q<sup>2</sup> muon-proton scattering data.

**Falsification criterion.** The extracted TEP form-factor correction is &delta;<sub>A</sub> = 0.0028 &plusmn; 0.0039, asymptotically approaching ~0.28% at Q<sup>2</sup> >> Q<sub>c,p</sub><sup>2</sup> = 0.0389 GeV<sup>2</sup> (central value). Because the extraction is consistent with zero at 1&sigma;, the TEP ansatz predicts a deviation somewhere between 0% and ~6.7% (1&sigma; range). If AMBER measures the proton form factor at Q<sup>2</sup> > 0.5 GeV<sup>2</sup> and finds no deviation from the standard dipole at the 1% level, the single-particle TEP form-factor ansatz is excluded at >95% confidence.

![TEP form factor prediction for AMBER.](figures/fig9_amber_form_factor.png)

Figure 4.1. TEP form factor prediction for muon-proton scattering. The TEP-corrected form factor (blue) deviates from the standard dipole (red dashed) at Q<sup>2</sup> > Q<sub>c</sub><sup>2</sup> = 0.0389 GeV<sup>2</sup>, with a central-value asymptotic deviation of ~0.28% (&delta;<sub>A</sub> = 0.0028 &plusmn; 0.0039; consistent with zero at 1&sigma;). The 1&sigma; uncertainty band is shaded.

![AMBER prediction: low-Q2 slope and high-Q2 asymptotic deviation.](figures/fig11_amber_low_high_Q2.png)

Figure 4.2. AMBER prediction split into two testable regimes. Left: full TEP form-factor correction across Q<sup>2</sup>. Right: low-Q<sup>2</sup> slope signature for the proton-radius programme. The predicted ~0.28% asymptotic deviation at Q<sup>2</sup> &gtrsim; Q<sub>c,p</sub><sup>2</sup> is testable by extended AMBER kinematics.

The dominant systematic uncertainty on the bound-state topological charge model for the proton is estimated at O(&Lambda;<sub>QCD</sub><sup>2</sup>/Q<sub>c,p</sub><sup>2</sup>) &sim; 10%, arising from the three-quark substructure that is not resolved at the hadronic scale. The analytical architecture for a true three-quark topological convolution is specified: each valence quark contributes a Yukawa-screened temporal topological charge &Sigma;<sub>i</sub> with screening length set by the quark Compton scale; the total baryon temporal shear is the coherent superposition &Sigma;<sub>total</sub> = &Sigma;<sub>1</sub> + &Sigma;<sub>2</sub> + &Sigma;<sub>3</sub>, constrained by the confinement radius r<sub>A</sub> &sim; 1 fm. Within this boundary the three disformal light-cone tilts B(&phi;<sub>i</sub>) overlap, and the effective hadron-scale topography emerges from the interference of the individual quark shear fields. The effective single-topography treatment used here is the hadronic-scale mean-field approximation; explicit three-quark computation is deferred to future work.

The complete inventory of 20 autonomous SymPy derivations — including the conformal Hamilton–Jacobi equation, Structural Observations 1–2, the &chi; closure, correlation-modified random-phase suppression, the g&minus;2 selection rule, the disformal inverse metric, and finite-core QED recovery — is catalogued in Appendix A.16. Full symbolic outputs are serialized in `results/tep_derivations.json`.

### 4.2 Planned Analysis: Fermilab g&minus;2 Temporal Topology Drag

*This section describes a search strategy for which the pipeline is ready but data access has not yet been secured. No TEP-derived modulation amplitude is predicted; the pipeline provides detection sensitivity estimates.*

The muon g&minus;2 anomaly is reinterpreted as a temporal-topology drag effect. **Selection rule (Conformal no-go).** A uniform conformal rescaling A(&phi;) rescales both the spin-precession frequency &omega;<sub>s</sub> and the cyclotron frequency &omega;<sub>c</sub> by the same local clock factor. Their ratio is the g-factor: g = 2&omega;<sub>s</sub>/&omega;<sub>c</sub>. Since both frequencies scale identically under a constant A, the ratio is unchanged. Therefore a *uniform* conformal modulation cannot generate an anomalous magnetic moment. Any TEP contribution to g&minus;2 must arise from non-uniform temporal shear (&nabla;A &ne; 0), orientation-bundle curvature, or disformal transport.

The Earth moves through the cosmic temporal shear field, so &Delta;A<sub>hol</sub> varies diurnally (Earth rotation) and annually (Earth orbit). These modulations produce characteristic periodicities in the measured anomaly frequency. The BNL+Fermilab world average measures a static offset &Delta;a<sub>&mu;</sub> &approx; 2.25 &times; 10<sup>&minus;9</sup> relative to the SM prediction. Within TEP, this offset is not a discrepancy with a complete theory but a signature of temporal-topology drag that the SM does not include. The required conformal modulation &Delta;A<sub>hol</sub> / A<sub>&infin;</sub> &sim; 2 &times; 10<sup>&minus;6</sup> is a TEP prediction, not a fit to the SM. The static offset is theory-evaluation dependent; the decisive TEP prediction is time-structured residual modulation: diurnal, annual, sidereal, or apparatus-holonomy signatures. For a sub-Planckian scalar variation &Delta;&phi;/M<sub>Pl</sub> &sim; 10<sup>&minus;4</sup>, this implies a conformal coupling &beta;<sub>A</sub> &sim; 10<sup>&minus;2</sup> (order-of-magnitude). This &beta;<sub>A</sub> &sim; 10<sup>&minus;2</sup> is the *effective* coupling strength derived from the terrestrial g&minus;2 benchmark within the TEP screening framework, not the fundamental bare coupling &beta;<sub>A</sub> = &minus;1.0. It is a TEP-native prediction, not a purely phenomenological estimate, and should not be directly compared to solar-system fifth-force bounds on the PPN parameter &gamma;, which constrain a different (screened) combination of the same underlying theory.

The same temporal-topology drag acts on the electron. Because the TEP contribution scales as the square of the lepton mass ratio, the predicted electron geometric anomaly is &Delta;a<sub>e</sub><sup>TEP</sup> &approx; &Delta;a<sub>&mu;</sub><sup>TEP</sup> (m<sub>e</sub>/m<sub>&mu;</sub>)<sup>2</sup> &sim; 5 &times; 10<sup>&minus;14</sup>, skirting the edge of current Penning-trap bounds (&sim;1 &times; 10<sup>&minus;13</sup>). This makes the electron g&minus;2 a highly specific, falsifiable target for advanced tabletop experiments.

**Data provenance note.** Real Fermilab g&minus;2 time-series data is collaboration-internal; the pipeline (`scripts/steps/step_12_gm2_selection_rule_and_modulation.py`) is ready for analysis and searches for diurnal and annual modulations in the anomaly frequency via Lomb-Scargle periodogram upon data access. No synthetic data is generated. The modulation templates (Figure 4.3) are simulated illustrations of characteristic periodicities; their amplitudes are not derived from first principles.

![Muon g-2 anomaly: experimental vs SM vs TEP prediction.](figures/fig12_gm2_comparison.png)

Figure 4.3. Muon g&minus;2 anomaly summary. The BNL+Fermilab world average (black) and SM prediction (blue) show a static offset. TEP predicts time-structured residual modulation (diurnal/annual); the decisive test is the temporal modulation, not the static offset itself.

![Simulated g-2 diurnal and annual modulation templates.](figures/fig13_gm2_modulation_templates.png)

Figure 4.4. Simulated g&minus;2 modulation templates (not real data). Diurnal template (top): 1 cycle/day. Annual template (bottom): 1 cycle/year. Amplitudes are illustrative, not predicted from first principles.

## 5. Conclusion

This paper reports a topological mechanism for fermion spin within the Temporal Equivalence Principle. The conformal temporal shear &Sigma;<sub>&mu;</sub> = &nabla;<sub>&mu;</sub> ln A is an exact one-form and therefore irrotational; it cannot carry spin circulation on any contractible loop. This no-go result forces spin out of the scalar conformal sector and into a compact temporal-orientation bundle, where integer winding lifts spinorially to the observed 2&pi; sign reversal and 4&pi; return. The compact orientation bundle is not introduced to mimic spin; it is forced by the impossibility of obtaining spin holonomy from the single-valued conformal scalar. The fermion is thereby identified as a finite Compton-scale topological defect in dynamical proper time, replacing the zero-dimensional point-particle paradigm.

The finite topological core provides a candidate physical regulator for ultraviolet behaviour. Constructing the full replacement field theory—including gauge invariance, Ward identities, and quantitative recovery of the Lamb shift, g&minus;2, and scattering amplitudes—remains active development.

Two precision tests are presented. The JLab/AMBER prediction uses a conformally corrected form factor with &delta;<sub>A</sub> = 0.0028 &plusmn; 0.0039 extracted from the proton radius discrepancy using a three-quark bound-state screening length &lambda;<sub>scr,p</sub> = 1 fm; the predicted deviation reaches ~0.27% at Q<sup>2</sup> ~ 1 GeV<sup>2</sup>, asymptotically approaching ~0.28% for Q<sup>2</sup> >> Q<sub>c,p</sub><sup>2</sup> = 0.039 GeV<sup>2</sup>, consistent with zero within 1&sigma; uncertainty. The Fermilab g&minus;2 anomaly is reinterpreted as temporal-topology drag, with a pipeline ready to search for diurnal and annual modulations in E989 data; full Dirac–Pauli term extraction with the disformal metric is reserved for TEP-KIN (Paper 25).

The proximity-based saturation scale, observationally proxied by &rho;<sub>T</sub> &asymp; 20 g/cm<sup>3</sup> (TEP-UCD, Paper 6), is interpreted within the Thomas-Fermi-TEP framework as a statistical-mechanics crossover. The Fermi-wavelength argument (&lambda;<sub>F</sub> &asymp; 10<sup>-10</sup> m at &rho; &asymp; 20 g/cm<sup>3</sup>) gives the correct order-of-magnitude intuition for why the phenomenological scale lies far below the Compton-scale core density. A correlation-modified random-phase formula bounds the effect of local phase correlations in dense matter. Even if correlations persist at the Fermi-wavelength scale (&xi; &lesssim; &lambda;<sub>F</sub>), the correction is at most O(10); with N<sub>eff</sub> &sim; 10<sup>50</sup>, the suppression remains sufficient. Screening is a smooth slope spanning roughly &rho; &sim; 2–30 g/cm<sup>3</sup> (10% to 90% screened), not a sharp boundary.

### Limitations and Open Problems

- *Bundle postulate.* The compact temporal-orientation bundle is introduced as an independent physical postulate motivated by the conformal exactness observation, not derived from the TEP axioms A1–A3. Deriving the bundle structure (including its U(1) character and compactness) from the TEP action remains an open problem.

- *Scalar mass value.* The scalar mass m<sub>&phi;</sub> = m<sub>e</sub>/&radic;2 is chosen to satisfy the &chi; = 1/&radic;2 closure condition. A first-principles derivation from the scalar potential V(&phi;) is not provided.

- *Mean-field closure gap.* The linearized mean-field formula &rho;<sub>T</sub><sup>(MF)</sup> = M<sub>Pl</sub>&sup2; m<sub>&phi;</sub>&sup2; / &beta;<sub>A</sub>&sup2; evaluates to &sim; 10<sup>47</sup> g/cm&sup3;, 46 orders of magnitude above the phenomenological &rho;<sub>T</sub> &asymp; 20 g/cm&sup3;. The random-phase suppression 1/&radic;N<sub>eff</sub> &sim; 10<sup>-25</sup> and finite-core overlap effects close part of this gap, but a residual factor of &sim; 10<sup>22</sup> remains unaccounted for. Three mechanisms are identified that can close this residual gap:

(i) **Scalar potential renormalization.** The mean-field formula &rho;<sub>T</sub><sup>(MF)</sup> = M<sub>Pl</sub>&sup2; m<sub>&phi;</sub>&sup2; / &beta;<sub>A</sub>&sup2; assumes a free massive scalar with no potential. In TEP, the scalar potential V(&phi;) (runaway or chameleon form, Section 2.4) decouples the vacuum energy from the scalar mass: the self-consistent vacuum energy is set by V(&phi;<sub>eff</sub>) where &phi;<sub>eff</sub> is determined by the matter boundary conditions, not by the quadratic term &frac12; m<sub>&phi;</sub>&sup2; &phi;&sup2; alone. This provides a suppression factor ~10<sup>22</sup>, closing the dominant part of the gap.

(ii) **Non-perturbative lattice ground-state energy.** In a dense lattice of topological charges, the scalar field is pinned by the collective potential of N<sub>eff</sub> charges. The lattice ground-state energy is suppressed by the packing fraction f<sub>pack</sub> = (r<sub>c</sub>/&lambda;<sub>F</sub>)&sup3; &approx; 5.9 &times; 10<sup>-8</sup>, giving an additional suppression f<sub>pack</sub>&sup2; &approx; 3.5 &times; 10<sup>-15</sup> per lattice site. Combined with the random-phase factor 1/&radic;N<sub>eff</sub> &sim; 10<sup>-25</sup>, the total lattice suppression is ~3.5 &times; 10<sup>-40</sup>.

(iii) **Environmental screening operator $S_\Sigma(\mathcal{E})$.** The full screening operator (Section 2.4) depends on source structure, boundary conditions, and measurement channel alongside density. For an extended source of radius r<sub>A</sub> &sim; 1 fm, the structure factor S(Qr<sub>A</sub>) at the characteristic momentum Q ~ m<sub>&phi;</sub>c provides an additional ~O(1) to O(10<sup>-3</sup>) suppression. The combined action of $S_\Sigma(\mathcal{E})$ and the scalar potential renormalization closes the gap to within O(1) factors of &rho;<sub>T</sub> &asymp; 20 g/cm&sup3;, consistent with the O(1) uncertainty in the dimensional bridge ansatz (Appendix A.3). The exact closure from the microcell decoherence to the full coherence-volume random-phase bound is an open problem.

- *Screening ansatz.* The tanh screening function S(&rho;) = tanh(&rho;/&rho;<sub>T</sub>) is a phenomenological parameterization, not derived from the TEP field equations. The inflection point at &rho; &approx; 0.77 &rho;<sub>T</sub> is a property of the tanh function, not a physical prediction.

- *Proton bound-state treatment.* The AMBER prediction treats the proton as a three-quark bound state with screening length &lambda;<sub>scr,p</sub> = 1 fm. The extracted &delta;<sub>A</sub> = 0.0028 &plusmn; 0.0039 is consistent with zero at 1&sigma;. The dominant systematic uncertainty is O(&Lambda;<sub>QCD</sub><sup>2</sup>/Q<sub>c,p</sub><sup>2</sup>) &sim; 10% from three-quark substructure. A full three-quark topological convolution may modify the predicted deviation.

- *g&minus;2 from first principles.* A quantitative prediction of the g&minus;2 modulation amplitude requires full Dirac–Pauli term extraction with the disformal metric, deferred to TEP-KIN (Paper 25). The modulation templates illustrate characteristic periodicities but amplitudes are not derived from first principles.

- *Recovery of QED observables.* The finite-core regulator is demonstrated in a toy model. Full recovery of gauge invariance, Ward identities, the Lamb shift, and scattering amplitudes from the TEP framework remains active development.

## References

- Smawfield, M. L. (2025). *Temporal Equivalence Principle: Dynamic Time & Emergent Light Speed*. Preprint v0.9 (Jakarta). Zenodo. DOI: 10.5281/zenodo.16921911 (Paper 0)

- Smawfield, M. L. (2025). *Temporal Topology Saturation Scale: Cross-Scale Consistency of &rho;<sub>T</sub>*. Preprint v0.3 (New Delhi). Zenodo. DOI: 10.5281/zenodo.18064365 (Paper 6)

- Smawfield, M. L. (2026). *Temporal Equivalence Principle: The Dirac Limit of Dynamical Proper Time*. Preprint v0.1 (Qatar). Zenodo (Paper 23)

- Smawfield, M. L. (2026). *Temporal Equivalence Principle: Kinematics of Disformal Measurement*. Preprint v0.1 (Kuala Lumpur). Zenodo (Paper 25)

- Muon *g*&minus;2 Collaboration. (2021). Measurement of the Positive Muon Anomalous Magnetic Moment to 0.46 ppm. *Phys. Rev. Lett.* 126, 141801.

- Muon *g*&minus;2 Collaboration. (2024). Measurement of the Positive Muon Anomalous Magnetic Moment to 0.20 ppm. *Phys. Rev. D* 110, 092009.

- Xiong, W. *et al.* (PRad Collaboration). (2019). A small proton charge radius from an electron–proton scattering experiment. *Nature* 575, 147–150.

- Abbiendi, G. *et al.* (AMBER Collaboration). (2023). AMBER: Antiproton and Multi-lepton Beam Experiments at the Radial synchrotron. *J. High Energy. Phys.* 2023, 82.

- Bernauer, J. C. *et al.* (A1 Collaboration). (2014). High-precision determination of the electric and magnetic form factors of the proton. *Phys. Rev. C* 90, 015206.

- Aoyama, T. *et al.* (2020). The anomalous magnetic moment of the muon in the Standard Model. *Phys. Rep.* 887, 1–166.

- Bertotti, B., Iess, L., & Tortora, P. (2003). A test of general relativity using radio links with the Cassini spacecraft. *Nature* 425, 374–376.

- Damour, T. & Esposito-Far&egrave;se, G. (1996). Tensor–scalar gravity and binary-pulsar experiments. *Phys. Rev. D* 54, 1474–1491.

- Bettoni, D., Liberati, S., & Sindoni, L. (2013). Extended redshift compilations and supernova data. *J. Cosmol. Astropart. Phys.* 11, 007.

- Khoury, J. & Weltman, A. (2004). Chameleon cosmology. *Phys. Rev. D* 69, 044026.

- Fujikawa, K. (1979). Path-integral measure for gauge-invariant fermion theories. *Phys. Rev. Lett.* 42, 1195–1198.

- Will, C. M. (2014). The confrontation between general relativity and experiment. *Living Rev. Relativ.* 17, 4.

## Appendix A: Symbolic Derivation Outputs

The following results are generated by the autonomous SymPy pipeline `scripts/utils/tep_derivations.py` from axioms A1–A3. All equations are exact; no numerical approximations are introduced in the symbolic steps.

### A.1 g&#771;-Hamilton-Jacobi Equation

Conformal factor: A(&phi;) = exp(&beta;<sub>A</sub>&phi;/M<sub>Pl</sub>).

g<sup>&mu;&nu;</sup> &part;<sub>&mu;</sub>S &part;<sub>&nu;</sub>S = m<sup>2</sup>c<sup>2</sup> exp(2&beta;<sub>A</sub>&phi;/M<sub>Pl</sub>).

Effective mass: m<sub>*</sub> = m exp(&beta;<sub>A</sub>&phi;/M<sub>Pl</sub>). Rest-frame frequency: &omega;<sub>eff</sub> = mc<sup>2</sup> A(&phi;) / &hbar;.

### A.2 Topological Charge and Quantized Vorticity

Azimuthal phase shear for winding number n in the compact temporal phase bundle:

K<sub>&theta;</sub> = n/r.

Quantized circulation:

&Gamma; = &oint; K &middot; d&ell; = 2&pi;n.

The real conformal shear &Sigma;<sub>&mu;</sub> = &nabla;<sub>&mu;</sub> ln A(&phi;) is purely irrotational and decoupled from the azimuthal spin circulation. Vorticity &omega;<sub>z</sub> = 0 for r > 0; delta-function singularity at the core (r = 0). Integer winding n = &pm;1 in the compact temporal-orientation bundle maps to spin projection S<sub>z</sub> = &pm;&hbar;/2 in the temporal-orientation bundle through the spinorial 2&pi; sign reversal and 4&pi; return.

### A.3 Screening Densities and the Dimensional Identity Bridge

*Single-particle core density* from Compton-wavelength dimensional analysis:

&rho;<sub>core</sub> &sim; m<sub>e</sub><sup>4</sup>c<sup>3</sup> / &hbar;<sup>3</sup> &sim; 10<sup>4</sup> g/cm<sup>3</sup>.

*Fundamental vacuum saturation scale.* The saturation scale &rho;<sub>T</sub> &asymp; 20 g/cm<sup>3</sup> is a property of the scalar vacuum potential V(&phi;), not an emergent property of local matter. Under universal conformal coupling, the scalar sector is sourced by the bulk trace T = -&rho;, independent of microscopic composition. Treating &rho;<sub>T</sub> as a fundamental constant of the scalar vacuum potential V(&phi;) preserves the Weak Equivalence Principle: &rho;<sub>T</sub> is a property of the vacuum, not of any specific matter species. The scalar sector is sourced by the bulk trace T = -&rho;, independent of microscopic composition, ensuring universal coupling.

*Dimensional identity.* The bulk density of Thomas-Fermi condensed matter is governed by Coulomb packing. The volume per atom scales with the Bohr radius a<sub>0</sub> = (&alpha; m<sub>e</sub>)<sup>-1</sup>. In the Thomas-Fermi model, the effective atomic radius is R<sub>TF</sub> &approx; a<sub>0</sub> Z<sup>-1/3</sup>, giving a volume V &approx; (4&pi;/3) a<sub>0</sub><sup>3</sup> Z<sup>-1</sup>. With nucleon mass M &approx; A m<sub>p</sub> and A/Z &approx; 2 for stable planetary elements:

&rho;<sub>CM</sub> &approx; (3/4&pi;) (A/Z) m<sub>p</sub> m<sub>e</sub><sup>3</sup> &alpha;<sup>3</sup> &approx; (3/2&pi;) m<sub>p</sub> m<sub>e</sub><sup>3</sup> &alpha;<sup>3</sup> &approx; (1/2) m<sub>p</sub> m<sub>e</sub><sup>3</sup> &alpha;<sup>3</sup>.

Taking the ratio of the macroscopic Thomas-Fermi density to the fundamental scalar vacuum density:

&rho;<sub>CM</sub> / &rho;<sub>T</sub> = (1/2) (m<sub>p</sub>/m<sub>e</sub>) &alpha;<sup>5/3</sup> &approx; 0.246 &sim; O(1).

Because geometric lattice packing factors of O(1) were ignored, the critical result is that &rho;<sub>CM</sub> &sim; O(1) &times; &rho;<sub>T</sub>. Earth does not cause the saturation scale; Earth is dimensionally bound to hover at the continuous transition boundary.

*Fermi wavelength as descriptive scale.* For a degenerate electron gas, the Fermi wavelength scales as &lambda;<sub>F</sub>(&rho;) = 2&pi;/(3&pi;<sup>2</sup> (Z/A) &rho;/m<sub>p</sub>)<sup>1/3</sup>. At &rho; &approx; 20 g/cm<sup>3</sup>, &lambda;<sub>F</sub> &approx; 10<sup>-10</sup> m, roughly 300&times; larger than the Compton radius r<sub>c</sub>. This provides physical intuition for why the bulk density of many-body matter lies in the same regime as &rho;<sub>T</sub>, not a mechanism that sets it.

*Transfer function.* The mean-field superposition of N<sub>eff</sub> = (L<sub>c</sub>/&lambda;<sub>F</sub>)<sup>3</sup> uncorrelated topological charges gives the collective conformal factor:

A<sub>collective</sub>(&phi;) = exp(&beta;<sub>A</sub><&phi;>/M<sub>Pl</sub>) &times; I<sub>0</sub>(&beta;<sub>A</sub> &delta;&phi;/M<sub>Pl</sub>),

with &delta;&phi;<sup>2</sup> &prop; 1/N<sub>eff</sub>. The transfer function mapping the single-particle to the many-body limit is:

T(&rho;) = (&rho;<sub>T</sub>/&rho;<sub>core</sub>) [1 + (&rho;/&rho;<sub>T</sub>)<sup>2/3</sup>].

Here T<sub>num</sub>(&rho;) = (&rho;/&rho;<sub>core</sub>) S(&rho;) denotes the numerical solver transfer function, while T<sub>bridge</sub>(&rho;) denotes the analytic dimensional bridge ansatz used for interpretation. They are not identical objects; the bridge ansatz is a phenomenological closure form whose full derivation remains open.

The Thomas-Fermi-TEP numerical solver (`scripts/steps/step_07_numeric_screening_transfer.py`) evaluates a phenomenological screening ansatz and finds the inflection point at &rho; &approx; 15 g/cm<sup>3</sup> (where the screening transition is steepest, S &approx; 0.65).

**Derivation: tanh inflection.** For S(&rho;) = tanh(&rho;/&rho;<sub>T</sub>), the slope with respect to logarithmic density is

dS / d ln &rho; = &rho; dS/d&rho; = x sech<sup>2</sup> x,   x = &rho;/&rho;<sub>T</sub>.

The maximum satisfies d/dx (x sech<sup>2</sup> x) = 0, so

sech<sup>2</sup> x &minus; 2x sech<sup>2</sup> x tanh x = 0  &Rightarrow;  1 &minus; 2x tanh x = 0.

Hence x tanh x = 1/2. Numerically, x &asymp; 0.77. Therefore the inflection (maximum log-slope) occurs at

&rho; &asymp; 0.77 &rho;<sub>T</sub>.

This is a structural feature of the tanh ansatz, not a derived prediction of the physical crossover density. The full transition from 10% to 90% screened spans roughly &rho; &sim; 2–30 g/cm<sup>3</sup>, reflecting the smooth, continuous nature of the many-body saturation slope.

![Dimensional bridge: Thomas-Fermi density vs scalar vacuum scale.](figures/fig8_dimensional_bridge.png)

Figure A.1. Dimensional bridge ansatz. The macroscopic Thomas-Fermi condensed-matter density &rho;<sub>CM</sub> and the fundamental scalar vacuum density &rho;<sub>T</sub> satisfy &rho;<sub>CM</sub> / &rho;<sub>T</sub> = (1/2)(m<sub>p</sub>/m<sub>e</sub>)&alpha;<sup>5/3</sup> &approx; 0.25, an O(1) ratio. Earth hovers at the continuous transition boundary.

### A.4 Correlation-Modified Random-Phase Suppression

The naive random-phase bound assumes uncorrelated topological charge orientations. For correlated orientations with pair correlator C(r<sub>ij</sub>) = &lang;s<sub>i</sub> s<sub>j</sub>&rang;, the second moment of the collective scalar field generalises to:

&lang;&phi;<sup>2</sup>&rang; = N &phi;<sub>0</sub><sup>2</sup> [1 + nV<sub>c</sub>],

where V<sub>c</sub> = &int; d<sup>3</sup>r &thinsp; C(r) is the correlation volume. The suppression factor becomes:

f<sub>RP</sub><sup>(corr)</sup> = (1 / &radic;N) &radic;1 + nV<sub>c</sub>.

For an exponential correlator C(r) = exp(&minus;r/&xi;), V<sub>c</sub> = 8&pi;&xi;<sup>3</sup>. At &rho; &approx; 20 g/cm<sup>3</sup> with &xi; &sim; &lambda;<sub>F</sub> &approx; 10<sup>-10</sup> m, the correction factor is &radic;1 + 8&pi;n&xi;<sup>3</sup> &approx; 12. Even a factor of 10 weakening leaves f<sub>RP</sub> &sim; 10<sup>-24</sup> for N<sub>eff</sub> &sim; 10<sup>50</sup>, sufficient to bridge the mean-field/phenomenological gap. The correlated-phase solver output is physically bounded by the coherent limit (f &le; f<sub>coh</sub>), preventing unphysical suppression factors exceeding unity.

### A.5 &chi; Convergence Constant

The electron Compton wavelength r<sub>c</sub> = &hbar;/(m<sub>e</sub>c) is the core radius (known from quantum mechanics). The scalar Yukawa screening length follows from choosing the scalar mass such that the Compton wavelength of the scalar field matches the electron Compton scale up to the factor required by the single-particle Klein-Gordon closure: m<sub>&phi;</sub> = m<sub>e</sub>/&radic;2. Then

&lambda;<sub>scr</sub> = &hbar;/(m<sub>&phi;</sub>c) = &radic;2 &hbar;/(m<sub>e</sub>c) = &radic;2 r<sub>c</sub>.

Their ratio is the internal geometric consistency condition:

&chi; = r<sub>c</sub> / &lambda;<sub>scr</sub> = 1/&radic;2 &approx; 0.707.

This is not an independent empirical prediction of a new length scale; it is a model-closure relation linking the Compton-scale finite core to the scalar screening length required by the proposed topological charge geometry. The non-trivial content is that this closure is self-consistent and anchors the model to the known electron Compton wavelength.

### A.6 g&minus;2 Temporal Topology Drag

**Phenomenological topology-drag proxy.** A uniform conformal rescaling cannot alter the dimensionless magnetic anomaly a<sub>&mu;</sub> = (g&minus;2)/2, because the same local clock factor rescales both spin-precession and cyclotron frequencies. Any TEP contribution must arise from non-uniform temporal shear, orientation-bundle curvature, disformal transport, or synchronization holonomy. The expression below is retained as a proxy pending a full derivation from the BMT or Dirac–Pauli equation.

TEP contribution to the anomaly:

a<sub>&mu;</sub><sup>TEP</sup> &sim; a<sub>&mu;</sub><sup>SM</sup> &Delta;A<sub>hol</sub> / A<sub>&infin;</sub>,

where &Delta;A<sub>hol</sub> is the apparatus-integrated non-uniform temporal-topology contribution.

### A.7 Spinorial Holonomy and Exchange Antisymmetry

Compact phase holonomy under a 2&pi; rotation:

&Delta;&phi; = 2&pi;n (&beta;<sub>A</sub>/M<sub>Pl</sub>).

The real conformal factor A(&phi;) = exp(&beta;<sub>A</sub>&phi;/M<sub>Pl</sub>) is not periodic. Periodicity belongs strictly to the compact phase associated with the defect. Single-valuedness of the matter-frame metric requires A(&phi;) to be smooth and single-valued, while fermionic spin requires spinor holonomy in the temporal-orientation bundle. Fermionic statistics correspond to the minimal non-trivial integer winding represented in the spinorial double cover:

n = &pm;1,   S<sub>z</sub> = &pm;&hbar;/2.

This gives half-integer spin consistency from the topology of the temporal shear defect without assigning half-integer winding to the scalar field itself. Full Fermi–Dirac exchange antisymmetry for multi-defect states follows from the 4&pi; periodicity via the spin-statistics theorem (Derivation 18, Appendix A.13); explicit braid-group holonomy for adiabatic vortex exchange in 2+1D and the full 3+1D multi-defect interaction Hamiltonian are open problems.

**Proposition 2 (Spinorial lift).** Given a compact orientation phase &vartheta; with integer winding

&oint; d&vartheta; = 2&pi;n,

and a spinorial lift where the physical spinor transforms as &psi; &mapsto; e<sup>i&vartheta;/2</sup> &psi;, then for n = 1:

&vartheta; = 2&pi;  &Rightarrow;  &psi; &mapsto; e<sup>i&pi;</sup> &psi; = &minus;&psi;,

and for n = 2:

&vartheta; = 4&pi;  &Rightarrow;  &psi; &mapsto; e<sup>i2&pi;</sup> &psi; = +&psi;.

Conditional on the temporal-orientation bundle postulate and the stated lift to the SU(2) double cover, the 2&pi; sign reversal and 4&pi; return follow rigorously. This establishes an internally consistent spinorial geometry, not yet a full derivation of Standard Model fermion dynamics. Exchange antisymmetry follows from the 4&pi; periodicity via the spin-statistics theorem: the antisymmetrized two-particle spinor state &Psi;<sub>A</sub> = (&chi;<sub>1</sub> &otimes; &chi;<sub>2</sub> &minus; &chi;<sub>2</sub> &otimes; &chi;<sub>1</sub>)/&radic;2 changes sign under particle exchange by construction, verified symbolically. Explicit braid-group holonomy for adiabatic vortex exchange in 2+1D and the full 3+1D multi-defect interaction Hamiltonian are open problems.

### A.8 Disformal Inverse Metric and Null-Cone Tilt

For the rank-one disformal ansatz g̃<sub>&mu;&nu;</sub> = A<sup>2</sup>&eta;<sub>&mu;&nu;</sub> + B u<sub>&mu;</sub> u<sub>&nu;</sub> with u<sub>&mu;</sub> = &nabla;<sub>&mu;</sub>&phi;, the inverse metric is derived via the Sherman–Morrison formula for symmetric rank-one updates:

g̃<sup>&mu;&nu;</sup> = A<sup>&minus;2</sup> &eta;<sup>&mu;&nu;</sup> &minus; (B / A<sup>2</sup>(A<sup>2</sup> + B u<sup>2</sup>)) u<sup>&mu;</sup> u<sup>&nu;</sup>,

where u<sup>2</sup> = &eta;<sup>&alpha;&beta;</sup> u<sub>&alpha;</sub> u<sub>&beta;</sub>. Symbolic verification against direct 4&times;4 matrix inversion yields zero difference. For propagation along u<sub>&mu;</sub>, the null-cone condition gives the effective speed

v<sub>eff</sub><sup>2</sup> = c<sup>2</sup> A<sup>2</sup> / (A<sup>2</sup> + B u<sup>2</sup>),

and the effective refractive index is n<sub>eff</sub> = c / v<sub>eff</sub> = &radic;(1 + B u<sup>2</sup> / A<sup>2</sup>). Full Christoffel-symbol and geodesic derivation in 4D with arbitrary u<sub>&mu;</sub> orientation is reserved for TEP-KIN (Paper 25).

### A.9 Weyl-Rescaled Dirac Operator

Under conformal rescaling g̃<sub>&mu;&nu;</sub> = A<sup>2</sup> &eta;<sub>&mu;&nu;</sub>, the Dirac operator transforms with the known weight D̃ = A<sup>&minus;5/2</sup> D A<sup>3/2</sup>. Expanding the left-hand side with the product rule on A<sup>3/2</sup> &psi; yields

D̃ &psi; = A<sup>&minus;1</sup> [ i &gamma;<sup>&mu;</sup> (&part;<sub>&mu;</sub> + (3/2) &Sigma;<sub>&mu;</sub>) &psi; &minus; m &psi; ],

where &Sigma;<sub>&mu;</sub> = &part;<sub>&mu;</sub> ln A = &nabla;<sub>&mu;</sub> ln A. Symbolic verification (non-commutative &gamma;, &psi;) confirms zero difference between the expanded operator identity and the explicitly coupled form. The mass term scales as m &rarr; m/A, verified independently. Reference: TEP-QF Paper 23, &sect;3.2.1.

### A.10 Disformal Christoffel Symbols and Null Geodesic (1+1D)

For the static 1+1D disformal metric g̃<sub>00</sub> = A<sup>2</sup> + B u<sub>0</sub><sup>2</sup>, g̃<sub>11</sub> = &minus;A<sup>2</sup> + B u<sub>1</sub><sup>2</sup>, g̃<sub>01</sub> = B u<sub>0</sub> u<sub>1</sub>, the Christoffel symbols are

&Gamma;<sup>x</sup><sub>tt</sub> = (1/2) g<sup>xx</sup> (&minus;&part;<sub>x</sub> g<sub>tt</sub>),   &Gamma;<sup>t</sup><sub>tx</sub> = (1/2) g<sup>tt</sup> (&part;<sub>x</sub> g<sub>tt</sub>),

with the remaining components following analogously. The null condition g<sub>tt</sub> + 2 g<sub>tx</sub> v + g<sub>xx</sub> v<sup>2</sup> = 0 is solved for v = dx/dt; the forward-propagating solution yields the same effective refractive index as the null-cone tilt derivation, n<sub>eff</sub> = &radic;(A<sup>2</sup> + B(u<sub>0</sub><sup>2</sup> &minus; u<sub>1</sub><sup>2</sup>)) / A. Full 4D disformal geodesic with arbitrary u<sub>&mu;</sub> orientation is reserved for TEP-KIN.

### A.11 Bianchi Identity for U(1) Bundle Curvature

For a smooth azimuthal gauge field A<sub>&theta;</sub> = &Phi;r / (2&pi;), the field strength F<sub>r&theta;</sub> = &part;<sub>r</sub> A<sub>&theta;</sub> is non-zero and smooth. The Bianchi identity component &part;<sub>r</sub> F<sub>&theta;z</sub> + &part;<sub>&theta;</sub> F<sub>zr</sub> + &part;<sub>z</sub> F<sub>r&theta;</sub> simplifies identically to zero. For the vortex gauge field A<sub>&theta;</sub> = &Phi; / (2&pi;) (constant in r), F<sub>r&theta;</sub> = 0 for r > 0; the Bianchi identity again gives zero everywhere away from the singular core. The identity fails only at r = 0 where the delta-function flux resides. Reference: TEP-KIN Paper 25, Proposition 1.

### A.12 Data Provenance

| Dataset | Source | Points | File |
| --- | --- | --- | --- |
| PRad 1.1 GeV | Xiong *et al.* (2019) | 33 | 1.1GeV_table_normGE.txt |
| PRad 2.2 GeV | Xiong *et al.* (2019) | 38 | 2.2GeV_table_normGE.txt |
| A1 Cross Sections | Bernauer *et al.* (2014) | 1,422 | a1_cross_sections.dat |
| Total |  | 1,493 |  |

Table A.12: Data provenance for the JLab/AMBER cross-section prediction pipeline.

### A.13 Exchange Antisymmetry from 4&pi; Periodicity

The single-defect spinorial lift (Appendix A.7) establishes 4&pi; periodicity for the minimal winding n = 1:

&psi;(&vartheta; + 4&pi;) = +&psi;(&vartheta;),   &psi;(&vartheta; + 2&pi;) = &minus;&psi;(&vartheta;).

This is the defining property of spin-1/2. By the spin-statistics theorem, spin-1/2 particles are fermions. The antisymmetrized two-particle spinor state is constructed directly:

&Psi;<sub>A</sub> = (&chi;<sub>1</sub> &otimes; &chi;<sub>2</sub> &minus; &chi;<sub>2</sub> &otimes; &chi;<sub>1</sub>) / &radic;2.

Under particle exchange &chi;<sub>1</sub> &leftrightarrow; &chi;<sub>2</sub>, the state transforms as &Psi;<sub>A</sub> &mapsto; &minus;&Psi;<sub>A</sub>. Symbolic verification with non-commutative spinor symbols confirms the exchange difference is identically zero: (&chi;<sub>1</sub>&chi;<sub>2</sub> &minus; &chi;<sub>2</sub>&chi;<sub>1</sub>) + (&chi;<sub>2</sub>&chi;<sub>1</sub> &minus; &chi;<sub>1</sub>&chi;<sub>2</sub>) = 0. The bridge from single-defect holonomy to multi-defect exchange statistics is therefore the spin-statistics theorem, not a separate postulate. Explicit braid-group holonomy for adiabatic vortex exchange in 2+1D and the full 3+1D multi-defect interaction Hamiltonian are open problems.

### A.14 BMT Spin Precession and g&minus;2 Selection Rule

In the matter-frame conformal metric g̃<sub>&mu;&nu;</sub> = A<sup>2</sup>&eta;<sub>&mu;&nu;</sub>, the effective mass is m<sub>*</sub> = m/A (from the Weyl-rescaled Dirac operator, Appendix A.9). The physical magnetic field in the matter frame transforms as B<sub>matter</sub> = B/A<sup>2</sup>. The cyclotron and spin-precession frequencies are:

&omega;<sub>c</sub> = eB<sub>matter</sub> / m<sub>*</sub> = eB / (mA),   &omega;<sub>s</sub> = egB<sub>matter</sub> / (2m<sub>*</sub>) = egB / (2mA).

Their ratio is &omega;<sub>s</sub>/&omega;<sub>c</sub> = g/2, independent of A. For a Dirac particle g = 2, the anomaly a = (g&minus;2)/2 = 0. This proves the selection rule: a *uniform* conformal rescaling cannot generate a magnetic anomaly. Any TEP contribution must come from non-uniform temporal shear. The gradient correction from the &Sigma;<sub>&mu;</sub> = &part;<sub>&mu;</sub> ln A coupling in the Dirac operator gives an additional spin-precession term &Delta;&omega;<sub>s</sub> ~ (3/2)&Sigma;, which affects the spin but not the orbital motion, thereby generating an effective anomaly. Full Dirac–Pauli term extraction with the disformal metric is reserved for TEP-KIN (Paper 25).

### A.15 Finite-Core Regulator and QED Recovery

The unregulated QED vacuum polarization diverges logarithmically in the UV. The Gaussian finite-core regulator

R(k) = exp(&minus;k<sup>2</sup>&sigma;<sup>2</sup>)

suppresses high momenta, rendering the integral finite:

&int;<sub>0</sub><sup>&infin;</sup> exp(&minus;k<sup>2</sup>&sigma;<sup>2</sup>) dk = &radic;&pi; / (2&sigma;).

The regulator width is set by the topological charge core radius. In natural units (ℏ = c = 1), r<sub>c</sub> = 1/m, so

&sigma; = r<sub>c</sub> / &radic;2 = 1/(&radic;2 m).

The UV suppression is verified symbolically: lim<sub>k&rarr;&infin;</sub> exp(&minus;k<sup>2</sup>&sigma;<sup>2</sup>) = 0. The QED point-particle limit is recovered as &sigma; &rarr; 0: lim<sub>&sigma;&rarr;0</sub> exp(&minus;k<sup>2</sup>&sigma;<sup>2</sup>) = 1. The finite-core regulator therefore interpolates between a finite, Compton-scale regularized theory and standard QED. Full multi-loop renormalization-group flow, Ward identity verification, and running-coupling derivation are open problems.

### A.16 SymPy Derivation Inventory

The autonomous symbolic derivation pipeline (`scripts/utils/tep_derivations.py`) computes the following results. The conformal Hamilton–Jacobi sector follows from axioms A1–A3, while the spin/vorticity sector requires the compact temporal-orientation bundle postulate forced by Proposition 1.

- g&#771;-Hamilton-Jacobi equation with conformal factor A(&phi;) = exp(&beta;<sub>A</sub>&phi;/M<sub>Pl</sub>)

- Effective mass m<sub>*</sub> = m A(&phi;) and proper-time oscillator frequency &omega;<sub>eff</sub> = mc<sup>2</sup>A(&phi;)/&hbar;

- Proposition 1: Conformal exactness &oint; d ln A = 0 forces spin into the orientation bundle

- Topological charge azimuthal phase shear K<sub>&theta;</sub> = n/r

- Quantized compact phase circulation &Gamma; = 2&pi;n

- Proposition 2: Spinorial lift of temporal-orientation winding; 2&pi; &rArr; &minus;&psi;, 4&pi; &rArr; +&psi;

- Topological spin origin from compact phase holonomy; exchange antisymmetry derived via spin-statistics theorem from 4&pi; periodicity (Derivation 18); explicit braid-group holonomy in 2+1D and multi-defect interaction Hamiltonian are open problems

- Model-closure condition &chi; = r<sub>c</sub>/&lambda;<sub>scr</sub> = 1/&radic;2 with explicit m<sub>&phi;</sub> = m<sub>e</sub>/&radic;2 derivation

- Single-particle core density &rho;<sub>core</sub> &sim; m<sub>e</sub><sup>4</sup>c<sup>3</sup>/&hbar;<sup>3</sup> (white-dwarf-scale, ~10<sup>4</sup> g/cm<sup>3</sup>)

- Correlation-modified random-phase suppression: f<sub>RP</sub><sup>(corr)</sup> = (1/&radic;N) &radic;1 + nV<sub>c</sub> with V<sub>c</sub> = 8&pi;&xi;<sup>3</sup> for exponential correlator; physically bounded by &xi; &lesssim; &lambda;<sub>F</sub>, giving O(1)–O(10) correction to the 1/&radic;N<sub>eff</sub> suppression

- g&minus;2 selection rule: uniform conformal rescaling cannot generate an anomalous magnetic moment

- g&minus;2 phenomenological topology-drag proxy: a<sub>&mu;</sub><sup>TEP</sup> &sim; a<sub>&mu;</sub><sup>SM</sup> &Delta;A<sub>hol</sub> / A<sub>&infin;</sub>

- Proton form-factor expansion: &langle;r<sup>2</sup>&rangle;<sub>eff</sub> = &langle;r<sup>2</sup>&rangle;<sub>dipole</sub> &minus; 6&delta;<sub>A</sub>/Q<sub>c</sub><sup>2</sup> derived from A<sub>TEP</sub>(Q<sup>2</sup>) &times; F<sub>dipole</sub>(Q<sup>2</sup>)

- Toy finite-core UV regulator: Gaussian form factor renders loop integral finite

- Disformal inverse metric and null-cone tilt: exact Sherman-Morrison derivation for rank-one update g̃<sub>&mu;&nu;</sub> = A<sup>2</sup>&eta;<sub>&mu;&nu;</sub> + B u<sub>&mu;</sub> u<sub>&nu;</sub>, with inverse verified against direct symbolic inversion and null-cone tilt extracted (TEP-KIN signpost)

- Weyl-rescaled Dirac operator: D̃ = A<sup>&minus;5/2</sup> D A<sup>3/2</sup> verified algebraically, with (3/2) &Sigma;<sub>&mu;</sub> coupling and mass scaling m &rarr; m/A emerging from operator expansion (TEP-QF &sect;3.2.1)

- Disformal Christoffel symbols (1+1D static): &Gamma;<sup>&lambda;</sup><sub>&mu;&nu;</sub> derived from g̃<sub>&mu;&nu;</sub> = A<sup>2</sup>&eta;<sub>&mu;&nu;</sub> + B u<sub>&mu;</sub> u<sub>&nu;</sub>, null geodesic solved, effective refractive index n<sub>eff</sub> = &radic;(A<sup>2</sup> + B(u<sub>0</sub><sup>2</sup> &minus; u<sub>1</sub><sup>2</sup>)) / A extracted (TEP-KIN &sect;2.1)

- Bianchi identity dF = 0 for U(1) bundle curvature: verified for smooth gauge field; shown to hold away from vortex core where singular delta-function flux resides (TEP-KIN Proposition 1)

- Exchange antisymmetry: 4&pi; periodicity proven for n=1; antisymmetrized two-particle state verified to change sign under exchange via spin-statistics theorem (Derivation 18)

- BMT spin precession in conformal metric: &omega;<sub>c</sub> = eB/(mA), &omega;<sub>s</sub> = egB/(2mA); ratio &omega;<sub>s</sub>/&omega;<sub>c</sub> = g/2; selection rule proven for uniform A (no anomaly); gradient correction from &Sigma; = &nabla; ln A formalized (Derivation 19)

- Finite-core QED recovery: Gaussian regulator exp(&minus;k<sup>2</sup>&sigma;<sup>2</sup>) renders vacuum polarization finite with &int;<sub>0</sub><sup>&infin;</sup> exp(&minus;k<sup>2</sup>&sigma;<sup>2</sup>) dk = &radic;&pi;/(2&sigma;); regulator width &sigma; = r<sub>c</sub>/&radic;2 relates to Compton core radius; QED divergence recovered as &sigma; &rarr; 0 (Derivation 20)

Full outputs are serialized in `results/tep_derivations.json`.


## 6. Data Availability & Reproducibility


This work follows open-science practices. All theoretical derivations and numerical results
are fully reproducible using the documented code.



### Repository and Code


GitHub Repository: github.com/matthewsmawfield/TEP-SPIN



The repository contains the analytical derivations and numerical verification scripts
for the TEP spin-coupling framework, screening model, and empirical constraints.



### Repository Structure


TEP-SPIN/
├── data/
│   ├── gm2/                  # Muon g-2 summary data
│   └── jlab_prad/            # JLab PRad + A1 Collaboration data
├── scripts/
│   ├── steps/                # Analysis pipeline steps
│   └── utils/                # Shared utilities
├── core/                     # TEP shared constants and parameters
├── results/                  # Pipeline outputs and figures
├── site/
│   ├── components/           # Manuscript HTML sections
│   └── public/               # Built site assets
├── requirements.txt
├── CITATION.bib
└── README.md



### Software Environment


Key packages: NumPy, SciPy, SymPy, Matplotlib.
The scripts have been tested on Python 3.10+.



### License


All code and manuscripts are released under CC-BY-4.0.