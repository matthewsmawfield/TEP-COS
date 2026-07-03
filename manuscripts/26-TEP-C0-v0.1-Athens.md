# Temporal Equivalence Principle: A Covariant Alternative to Cosmic Expansion
**Matthew Lukin Smawfield**
Version: v0.1 (Athens)
First published: 3 July 2026 - Last updated: 3 July 2026
DOI: 10.5281/zenodo.20370144

---

## Abstract

This paper presents a direct empirical challenge to the necessity of primitive cosmic expansion. In the Temporal Equivalence Principle framework, observed redshift is reconstructed as conformal proper-time transport, $1+z=A_0/A_{\rm em}$, rather than as stretching of a spatial scale factor. Standard cosmology interprets observational redshift and luminosity distance scaling as evidence of a stretching spatial metric, parameterized by the Friedmann-Lemaître-Robertson-Walker (FLRW) scale factor $a(t)$. The observational role played by the FLRW scale factor is mapped, within the TEP conformal-frame construction, onto the temporal clock-rate field $A(\phi)$. In the tested late-time background sector, the perceived acceleration normally attributed to Dark Energy, $\Lambda$, is reconstructed as the kinetic energy density of the Temporal Shear field, $\Omega_\phi$.

The core relation is $1+z = A_0/A_{\text{em}}$. In the static conformal interpretation developed here, intergalactic separations are not treated as primitively expanding; the apparent expansion is reconstructed through temporal transport. In this framework, the limit conventionally written as $a\to0$ is re-expressed as $A_{\text{clock}}\to0$: a TEP temporal-horizon boundary of observational clock transport, not a zero-volume spatial singularity. The framework is formally closed from the temporal-horizon background ($a_{\rm eff} \to 0$) down to linear mode evolution ($k \sim 0.5 \, h/\text{Mpc}$) through the exact Bellini-Sawicki EFT mapping implemented in TEP-HC (Paper 18), which natively yields $\sigma_8 \approx 0.825$ in agreement with Planck.

Using 1,701 Pantheon+ Type Ia supernovae with the full covariance matrix, a pure conformal reconstruction exactly reproduces the $\Lambda$CDM homogeneous distance-modulus relation, demonstrating that the background Hubble diagram does not uniquely select an expanding spatial metric. More strongly, the conservative physical no-$\Lambda$ temporal-shear branch with fixed $z_T=5$ (the turnover redshift motivated by acoustic-sector physics) improves the standardized supernova likelihood by $\Delta\chi^2 \simeq -3.4$ relative to baseline $\Lambda$CDM and achieves a Bayes factor of BF$\simeq 4.6$, classified as "substantial" evidence on the Jeffreys scale. The fixed $z_T=100$ benchmark gives the strongest evidence, BF$\simeq 61.8$ ($\Delta\chi^2 \simeq -7.5$), while the broad free-$z_T$ model gives BF$\simeq 40.3$, showing that the preference is not solely a fixed-turnover artefact. The conservative physical model with $z_T=5$ already demonstrates that a matter-only temporal-shear geometry is competitive with $\Lambda$CDM in the late-time SNe distance-redshift sector. The same framework gives a parameter-locked host-environment prediction from the suppressed scalar-field geometry of host galaxies, with the massive-host-brighter direction matching the established astrophysical mass-step orientation while the simplified mini-analysis remains noise-dominated. These results establish positive supernova-sector evidence that apparent acceleration can be reconstructed as temporal transport rather than primitive dark energy.

Companion papers establish the theoretical foundations: TEP-HC (Paper 18) provides the Boltzmann-level acoustic-scale preservation proof under the native hi_class `tep_mode` implementation, and TEP-TH develops the nonsingular temporal-horizon closure. The current paper focuses on the empirical supernova-sector test and the deterministic falsification pipeline.

Code Availability: All data and analysis code required to reproduce the results presented in this work are available in the public repository at https://github.com/matthewsmawfield/TEP-C0.

Keywords: temporal equivalence principle, static conformal geometry, cosmology, dark energy, supernovae, Bayesian inference, modified gravity, temporal shear

# 1. Introduction: The Geometry of Time

Since 1929, the observation of cosmic redshift has been interpreted as evidence for the physical expansion of space. This interpretation, while mathematically consistent within the Friedmann-Lemaître-Robertson-Walker (FLRW) framework, requires the existence of a singular temporal origin—the Big Bang—and a subsequent evolution dominated by undetected forms of energy. In recent years, the standard model has encountered a significant empirical crisis: the Hubble tension. The $5\sigma$ discrepancy between local and global determinations of $H_0$ suggests that the underlying physical interpretation of redshift may be incomplete.

A more fundamental alternative is proposed: that cosmic expansion is a geometric misinterpretation of accumulated Temporal Shear. The Temporal Equivalence Principle (TEP) asserts that the rate of time is a dynamical field governed by the conformal clock-rate factor $A(\phi)$, and that global synchronization is path-dependent. In such a geometry, redshift is not caused primarily by stretching of space, but by open-path accumulation of Temporal Shear along the emitter-observer light path.

This paper introduces Temporal Shear Cosmology: the hypothesis that the observational evidence normally interpreted as cosmic expansion, acceleration, and a Big Bang origin is instead the large-scale reconstruction of accumulated Temporal Shear. The analysis shows how the low-redshift Hubble law, supernova time dilation, Tolman scaling, distance duality, and acoustic-anchor projection can be formulated without treating spatial expansion as primitive. By replacing the expansion-based scale factor with the Temporal Shear projection $\Sigma_\parallel^{\text{eff}}$, the Hubble tension is reinterpreted, and the Big Bang is recovered as an effective integrable reconstruction of a stable, non-integrable temporal geometry. Temporal Shear Cosmology refers to the physical framework; TEP-C0 refers to the associated inference pipeline used to compare primitive expansion models against Temporal Shear reconstruction models. Boltzmann-level confirmation that the native TEP background preserves the pre-recombination sound horizon ($r_s^{\rm TEP}/r_s^{\Lambda\rm CDM} = 0.999994$) is established independently in TEP-HC (Paper 18).

The claim-discipline framework for the TEP corpus, including the scope limitations of canonical precision tests, is established in TEP-EXP (Paper 9).

# 2. Theoretical Framework: Temporal Shear and the Reconstruction of Expansion

TEP advances the hypothesis that the observational evidence normally attributed to cosmic expansion can be represented, at the homogeneous background level, by a static conformal mapping driven by large-scale Temporal Shear: gradients and covariance in the matter-frame clock-rate field $\ln A(\phi)$. In TEP, matter, clocks, electromagnetic fields, and quantum phases couple universally to the causal matter metric $\tilde{g}_{\mu\nu} = A^2(\phi)g_{\mu\nu} + B(\phi)\nabla_\mu\phi\nabla_\nu\phi$, where the conformal factor $A(\phi)$ defines the Temporal Shear vector:

\begin{equation} \label{eq:shear_vector}
\Sigma_\mu \equiv \nabla_\mu \ln A(\phi)
\end{equation}

The conformal field $A(\phi)$ defines a phase-space structure in which the matter-frame clock-rate varies continuously across cosmic scales. The phase-space topology of this field determines whether transport is integrable or path-dependent, distinguishing pure conformal shear from non-integrable temporal transport.

## 2.1 The Cosmological Isochrony Assumption

Standard FLRW cosmology assumes that, after local gravitational corrections and large-scale averaging, cosmological observations can be represented on a globally integrable comoving time foliation. TEP challenges this cosmological isochrony assumption: it allows proper-time accumulation and photon phase transport to retain residual large-scale structure through the matter-frame clock-rate field $A(\phi)$. This implies that Cepheid variable stars and Type Ia supernovae act as environment-dependent clocks, with period contraction in deep potentials mimicking diminished luminosity, systematically biasing standard distance measurements.

## 2.2 The Generator of Apparent Redshift

Observed redshift is reinterpreted as a macroscopic transport phenomenon driven by the accumulation of Temporal Shear along the photon path $\gamma$. The line-of-sight projection is defined as $\Sigma_\parallel \equiv \Sigma_\mu \hat{k}^\mu$, where $\hat{k}^\mu$ is the tangent 4-vector normalized to the comoving observer frame, giving $\Sigma_\parallel$ dimensions of inverse length. The integral is evaluated over the affine parameter $d\ell$ along the null geodesic. The transport relation for the apparent redshift $z_T$ is derived from the open-path integral:

\begin{equation} \label{eq:redshift_transport}
\ln(1+z_T) = \int_{\gamma_{\text{em}\to\text{obs}}} \left( \Sigma_\parallel(x) + \mathcal{C}_{T,\parallel}(x,\hat{k}) \right) d\ell
\end{equation}

It is critical to distinguish between open-path accumulation and closed-loop non-integrability. Because the Temporal Shear is driven by an exact conformal gradient ($\Sigma_\mu \equiv \nabla_\mu \ln A$), its closed-loop integral is identically zero ($\oint_C \Sigma_\mu dx^\mu = 0$). Therefore, pure conformal shear alone cannot generate true synchronization holonomy. The non-integrable transport is strictly sourced by the non-exact topological covariance term $\mathcal{C}_T$, whose line-of-sight projection $\mathcal{C}_{T,\parallel}$ enters the open-path transport integral. This term accounts for path-dependent coarse-graining and stochastic topology corrections derived from $C_\Theta(x,x')$.

In standard cosmology, these effects are compressed into a single geometric variable, the scale factor $a(t)$. In TEP, $a(t)$ is recognized as an effective integrable reconstruction:

\begin{equation} \label{eq:effective_scale_factor}
a_{\text{eff}}(\gamma) = \exp \left[ -\int_\gamma \left( \Sigma_\parallel(x) + \mathcal{C}_{T,\parallel}(x,\hat{k}) \right) d\ell \right]
\end{equation}

The reconstructed scale factor $a_{\text{eff}}$ is the open-path FLRW-like projection developed in the temporal-horizon framework of TEP-TH (Paper 27). It decomposes into the exact observational clock map $A_{\text{clock}}(\gamma)=A_0/A_{\text{em}}$ and the suppressed physical dynamical response $A_{\text{dyn}}(\gamma)$, with the non-exact covariance/topology correction $\mathcal{C}_T$ providing the path-dependent transport closure. In this decomposition, $A_{\text{clock}}$ alone generates the standard redshift–distance relation via exact conformal shear, while $A_{\text{dyn}}$ encodes the residual suppressed response of the Temporal Shear field.

## 2.3 From Temporal Topology to Transport: Definition of $\mathcal{C}_T$

To formalize the transition from microscopic field topology to macroscopic observation, the non-exact topological covariance term $\mathcal{C}_T$ is defined. Let $\theta = \ln A(\phi)$. The coarse-grained covariance structure is given by:

\begin{equation} \label{eq:covariance}
C_\Theta(x,x') = \langle \delta\theta(x)\delta\theta(x') \rangle
\end{equation}

Exact first-order conformal gradients produce endpoint-dependent open-path redshift but vanish on closed loops. True synchronization holonomy therefore requires the non-exact $\mathcal{C}_T$ contribution. Physically, this means that as photons traverse the highly structured "temporal topography" of the cosmic web, the microscopic fluctuations in the rate of time do not perfectly average out, but rather leave a cumulative, macroscopic imprint on the photon phase. Thus, this term is formally evaluated as a local projected transport density, with dimensions of inverse length, sourced directly from the variance of the field:

\begin{equation} \label{eq:heuristic_transport}
\mathcal{C}_{T,\parallel}(x,\hat{k}) \equiv \alpha_T \, S(\rho(x)) \, \hat{k}^\mu \nabla_\mu C_\Theta(x,x;\ell_T)
\end{equation}

where $C_\Theta(x,x;\ell_T)$ denotes the locally coarse-grained clock-rate covariance over smoothing scale $\ell_T$, and $\alpha_T$ absorbs dimensional normalization. In this expression, $S(\rho)\to1$ in unsuppressed voids and $S(\rho)\to0$ in dense environments undergoing Temporal Topology flattening, ensuring that the covariance-induced transport contribution follows the same environmental logic as the macroscopic $\epsilon_T^{\text{obs}}=S(\rho)\epsilon_T$ relation.

Crucially, $\mathcal{C}_{T,\parallel}$ is introduced as a macroscopic transport-closure term motivated by the microscopic proper-time phase holonomy developed in the TEP-QF sector (Paper 23). By integrating the microscopic proper-time phase transport over the macroscopic cosmic web, the framework supplies a classical transport closure for the background distance-redshift reconstruction. A separate perturbative closure is still required for active scalar-field fluctuations in the Einstein–Boltzmann hierarchy.

## 2.4 The Universal Coupling Axiom and Covariant Environmental Gradient Suppression

Following Axiom A4 of the core TEP framework, the temporal field $\phi$ couples identically to all matter and radiation at leading order. However, the locally observable Temporal Shear is subject to strong environmental gradient suppression governed by the abstract operator $\mathcal{S}_\Sigma(\mathcal{E})$. Because $\mathcal{E}$ encompasses source structure, boundary conditions, and ambient fields, a complete theory must supply a single covariant realization of this operator, not a patchwork of scale-specific proxies. That realization is constructed here.

The matter-frame metric $\tilde{g}_{\mu\nu} = A^2(\phi)g_{\mu\nu} + B(\phi)\nabla_\mu\phi\nabla_\nu\phi$ implies that the physical strength of the conformal sector is measured by the scalar invariant $\Sigma^2 \equiv \Sigma_\mu\Sigma^\mu = (\beta_A/M_{\rm Pl})^2 \nabla_\mu\phi\nabla^\mu\phi$. In any local Lorentz frame, $\Sigma^2$ sets the squared fractional rate at which clocks dephase relative to the gravitational metric. Suppression is the dynamical flattening of this observable dephasing. The covariant screening operator is defined as the rational function of two dimensionless control parameters, a kinetic ratio and a density ratio:

\begin{equation} \label{eq:unified_screening}
\mathcal{S}_\Sigma(\mathcal{E}) \equiv \left[ 1 + \left(\frac{\Sigma_\mu\Sigma^\mu}{g_t^2}\right)^n + \left(\frac{\rho}{\rho_{\rm half}}\right)^2 \right]^{-1}
\end{equation}

Here $g_t$ is the critical shear scale at which non-linear kinetic self-coupling becomes dominant, and $\rho_{\rm half} \approx 0.5\,M_\odot/{\rm pc}^3$ is the ambient half-suppression density. The exponent $n$ governs the steepness of the kinetic transition. In the cosmological weak-field regime, $\Sigma^2 \sim H_0^2/c^2 \sim 10^{-56}\,{\rm m}^{-2}$, so the kinetic term is negligible and $\mathcal{S}_\Sigma \to S(\rho) = [1+(\rho/\rho_{\rm half})^2]^{-1}$. In the Solar System, where $g = |\nabla\Phi| \sim c^2|\nabla\ln A|/\beta_A$ in the Newtonian limit, the shear scale maps directly to the local gravitational acceleration: $\Sigma^2 \approx (\beta_A g/c^2)^2$, and the dominant suppression comes from the first term, giving $f(g) = [1+(g/g_t')^n]^{-1}$ with $g_t' = c^2 g_t/\beta_A$. Both phenomenological proxies are therefore low- and high-curvature limits of a single covariant expression.

### 2.4.1 The Covariant Action

The TEP bi-metric action, established in the foundational framework (Paper 0), is

\begin{equation} \label{eq:tep_action}
S = \int d^4x\,\sqrt{-g}\left[\frac{M_{\rm Pl}^2}{2}R - \frac{1}{2}(\nabla\phi)^2 - V(\phi)\right] + S_m[\psi_i,\tilde{g}_{\mu\nu}]
\end{equation}

with the screened matter-frame metric

\begin{equation} \label{eq:screened_metric}
\tilde{g}_{\mu\nu} = \mathcal{A}^2(\phi,\mathcal{E})\,g_{\mu\nu} + B(\phi)\nabla_\mu\phi\nabla_\nu\phi
\end{equation}

where the environment-dependent conformal factor

\begin{equation}
\mathcal{A}(\phi,\mathcal{E}) = \exp\!\left[\mathcal{S}_\Sigma(\mathcal{E})\,\frac{\beta_A\phi}{M_{\rm Pl}}\right]
\end{equation}

absorbs the suppression directly into the matter coupling. In the unscreened limit ($\mathcal{S}_\Sigma \to 1$) this reduces to the bare TEP conformal factor $A(\phi) = \exp(\beta_A\phi/M_{\rm Pl})$; in the fully screened limit ($\mathcal{S}_\Sigma \to 0$) matter couples directly to the Einstein metric $g_{\mu\nu}$. The disformal function $B(\phi)$ is bounded by multi-messenger constraints ($|c_\gamma - c_g|/c \lesssim 10^{-15}$) and is set to zero in the pure-conformal limit analysed here.

### 2.4.2 Variation and Field Equations

Varying the action (\ref{eq:tep_action}) with respect to the Einstein-frame metric $g^{\mu\nu}$ yields the Einstein equations

\begin{equation}
G_{\mu\nu} = \frac{1}{M_{\rm Pl}^2}\left[T_{\mu\nu}^{(\phi)} + T_{\mu\nu}^{(m)}\right]
\end{equation}

where $T_{\mu\nu}^{(\phi)} = \nabla_\mu\phi\nabla_\nu\phi - g_{\mu\nu}\left[\frac{1}{2}(\nabla\phi)^2 + V(\phi)\right]$ is the scalar stress-energy, and the Einstein-frame matter stress-energy follows from the functional derivative of $S_m[\tilde{g}]$ with respect to $g^{\mu\nu}$. In the conformal limit ($B=0$) this gives

\begin{equation}
T_{\mu\nu}^{(m)} = \mathcal{A}^2(\phi,\mathcal{E})\,\tilde{T}_{\mu\nu}^{(m)}
\end{equation}

where $\tilde{T}_{\mu\nu}^{(m)}$ is the matter-frame stress-energy. Variation with respect to $\phi$ yields the scalar equation of motion

\begin{equation} \label{eq:scalar_eom}
\Box\phi - V_{,\phi} = -\mathcal{Q}_{\rm eff}
\end{equation}

with effective source

\begin{equation}
\mathcal{Q}_{\rm eff} = \mathcal{S}_\Sigma\,\mathcal{A}\,\mathcal{A}_{,\phi}\,g_{\mu\nu}\mathcal{T}^{\mu\nu} + \mathcal{A}\,\mathcal{A}_{,\phi}\,g_{\mu\nu}\mathcal{T}^{\mu\nu}\,\phi\,\frac{\partial\mathcal{S}_\Sigma}{\partial\phi} + \mathcal{A}\,\mathcal{A}_{,\phi}\,g_{\mu\nu}\mathcal{T}^{\mu\nu}\,\frac{\partial\mathcal{S}_\Sigma}{\partial(\nabla\phi)^2}\,\frac{\partial(\nabla\phi)^2}{\partial\phi}
\end{equation}

where $\mathcal{T}^{\mu\nu} = (\sqrt{-\tilde{g}}/\sqrt{-g})\,\tilde{T}^{\mu\nu}$ is the density-weighted matter tensor. The first term is the direct conformal coupling; the remaining terms encode the feedback from the environmental dependence of $\mathcal{S}_\Sigma$. In the cosmological background, where $\Sigma^2 \ll g_t^2$ and $\rho \ll \rho_{\rm half}$, $\mathcal{S}_\Sigma \approx 1$ and its field derivatives are suppressed by $\Sigma^2/g_t^2 \ll 1$, so the source reduces to the standard conformally-coupled form $\mathcal{Q}_{\rm eff} \approx \beta_A\,\mathcal{T}/M_{\rm Pl}$. Near compact bodies, where $\mathcal{S}_\Sigma \ll 1$, the scalar force is suppressed and the source vanishes.

### 2.4.3 Perturbation Expansion and Gauge Conditions

To map the theory onto the Bellini--Sawicki EFT, we expand around a spatially flat FLRW background. The metric perturbation is written in Newtonian gauge

\begin{equation}
ds^2 = -(1+2\Psi)\,dt^2 + a^2(t)(1-2\Phi)\,\delta_{ij}\,dx^i dx^j
\end{equation}

and the scalar field is split as $\phi(t,\mathbf{x}) = \bar{\phi}(t) + \delta\phi(t,\mathbf{x})$. The perturbed matter metric acquires a conformal-frame fluctuation

\begin{equation}
\delta\tilde{g}_{\mu\nu} = 2\,\mathcal{S}_\Sigma\,\frac{\beta_A}{M_{\rm Pl}}\,\mathcal{A}^2\,g_{\mu\nu}\,\delta\phi + \mathcal{A}^2\,\delta g_{\mu\nu} + O(\delta\phi)^2
\end{equation}

where $\mathcal{S}_\Sigma$ is evaluated on the background. On cosmological scales, $\mathcal{S}_\Sigma \approx 1$ to excellent approximation, and the perturbation structure reduces to that of a standard scalar-tensor theory with effective coupling $\beta_A^{\rm eff} = \mathcal{S}_\Sigma\beta_A$.

### 2.4.4 Bellini--Sawicki EFT Mapping

In the pure-conformal limit ($B=0$), the quadratic action for scalar and metric perturbations maps onto the standard EFT-of-dark-energy form. The running of the effective Planck mass is read off from the time dependence of the background coupling:

\begin{equation}
\alpha_M = \frac{d\ln M_{\rm eff}^2}{d\ln a} = \frac{d\ln\mathcal{A}^2}{d\ln a} = 2\,\mathcal{S}_\Sigma\,\frac{\beta_A}{M_{\rm Pl}}\frac{\dot{\bar{\phi}}}{H}
\end{equation}

Using the TEP background relation $\alpha_A \equiv -d\ln\mathcal{A}/d\ln(1+z)$, this becomes $\alpha_M = -2\,\mathcal{S}_\Sigma\alpha_A$. In the cosmological weak-field limit ($\mathcal{S}_\Sigma \approx 1$) this reduces to the bare value $\alpha_M^{\rm bare} = -2\alpha_A$ used in TEP-HC (Paper 18). In screened environments ($\mathcal{S}_\Sigma \approx 0$), $\alpha_M \to 0$ and the scalar fifth force vanishes.

The braiding parameter follows from the kinetic mixing between $\delta\phi$ and the metric potentials:

\begin{equation}
\alpha_B = -\alpha_M = 2\,\mathcal{S}_\Sigma\alpha_A
\end{equation}

and the kineticity parameter from the canonical kinetic term of $\delta\phi$ after field redefinition:

\begin{equation}
\alpha_K = -5(\mathcal{S}_\Sigma\alpha_A)^2
\end{equation}

The tensor speed excess $\alpha_T$ vanishes in the conformal limit because $c_g^2 = c_\gamma^2 = 1$ is preserved. These are exactly the relations implemented in the TEP-HC hi_class runtime and used in the growth solver (step\_06\_03).

The no-ghost discriminant follows from the $2\times2$ kinetic matrix of the scalar sector. In the Bellini--Sawicki formalism,

\begin{equation}
D = \alpha_K + \frac{3}{2}\alpha_B^2 = -5(\mathcal{S}_\Sigma\alpha_A)^2 + \frac{3}{2}(2\mathcal{S}_\Sigma\alpha_A)^2 = (\mathcal{S}_\Sigma\alpha_A)^2 \ge 0
\end{equation}

The discriminant is manifestly non-negative for all $\mathcal{S}_\Sigma$ and $\alpha_A$, establishing ghost-freedom from the action. The sound speed is $c_s^2 = 1$ exactly in the conformal limit, guaranteeing gradient stability.

### 2.4.5 Post-Newtonian Expansion

The Solar-System PPN parameters are obtained from the quasi-static weak-field limit of the field equations. In the Damour--Esposito-Far\`ese parameterization for scalar-tensor theories with conformal coupling $\mathcal{A}(\phi)$, the metric perturbation for a static, spherically symmetric source is, at leading order,

\begin{equation}
g_{00}^{\rm J} = -1 + \frac{2GM}{r}\left(1 + \frac{\alpha_{\rm eff}^2}{2}\right), \qquad g_{rr}^{\rm J} = 1 + \frac{2GM}{r}\left(1 - \frac{\alpha_{\rm eff}^2}{2}\right)
\end{equation}

where $\alpha_{\rm eff}$ is the effective scalar charge sourced by the body. For the TEP screened coupling,

\begin{equation}
\alpha_{\rm eff} = \mathcal{S}_\Sigma(\mathcal{E})\,\alpha_0
\end{equation}

with $\alpha_0 = \beta_A/M_{\rm Pl}$ the bare coupling constant. The PPN parameter $\gamma$ is then

\label{eq:ppn_gamma}
\begin{equation}
\gamma_{\rm PPN} = 1 - 2\alpha_{\rm eff}^2 = 1 - 2\,\mathcal{S}_\Sigma^2\,\alpha_0^2
\end{equation}

In unscreened environments ($\mathcal{S}_\Sigma \approx 1$), the bare TEP coupling $\beta_A = -1$ gives $\alpha_0 = -1/M_{\rm Pl}$, which translates to $\gamma \approx -1$ --- ruled out by Cassini at $\sim$87\,000$\sigma$. In the Solar System, the gradient-dependent suppression dominates. Using the Newtonian mapping $\Sigma^2 \approx (\beta_A g/c^2)^2$ with $g = |\nabla\Phi|$,

\begin{equation}
\mathcal{S}_\Sigma \approx \left[1 + \left(\frac{\beta_A g}{c^2 g_t}\right)^n\right]^{-1} \equiv \left[1 + \left(\frac{g}{g_t'}\right)^n\right]^{-1}
\end{equation}

where $g_t' = c^2 g_t/|\beta_A| = c^2 g_t$ for the locked TEP value $\beta_A = -1$. Evaluating at Saturn orbit ($g_{\rm Cassini} \approx 6.5\times10^{-5}\,{\rm m\,s}^{-2}$) with $g_t = 1.0\times10^{-9}\,{\rm m\,s}^{-2}$ and $n=2$ gives $\mathcal{S}_\Sigma \approx 2.37\times10^{-10}$. Substituting into (\ref{eq:ppn_gamma}),

\begin{equation}
\gamma - 1 = -2\,(2.37\times10^{-10})^2\,\alpha_0^2 M_{\rm Pl}^2 \approx -1.1\times10^{-19}
\end{equation}

safely below the Cassini bound $|\gamma - 1| < 2.3\times10^{-5}$ by more than fourteen orders of magnitude. At Earth surface ($g \approx 9.8\,{\rm m\,s}^{-2}$), $\mathcal{S}_\Sigma \approx 10^{-20}$ and the deviation is utterly negligible. The E&ouml;tv&ouml;s parameter satisfies $|\beta_{\rm PPN} - 1| \propto \alpha_{\rm eff}^3$ and vanishes in the screened limit for the same reason.

### 2.4.6 Connection to the Growth Solver and Parameter Locking

The EFT functions derived above are the inputs to the structure-formation growth equation used in step_06_03 and step_06_07. The $\alpha_M$-modified growth ODE,

\begin{equation}
\frac{d^2D}{d(\ln a)^2} + \left(\frac{1}{2} - \frac{3}{2}w_{\rm eff} - \alpha_M\right)\frac{dD}{d\ln a} - \frac{3}{2}\Omega_m(a)\left(1 + \frac{\alpha_M}{3}\right)D = 0
\end{equation}

is the quasi-static limit of the full Einstein-Boltzmann hierarchy with the Bellini--Sawicki functions derived in Section 2.4.4. The $\alpha_M$ that appears here is precisely $\alpha_M^{\rm bare} = -2\alpha_A$ evaluated on the cosmological background where $\mathcal{S}_\Sigma \approx 1$.

The transition scale $g_t$ is not a free parameter fitted to the Solar System data. It is fixed by requiring that the unscreened branch (which would give $\gamma \approx -1$) be excluded, and that the suppressed branch pass the Cassini bound with a safety margin. The minimum requirement is $\mathcal{S}_\Sigma(g_{\rm Cassini}) \lesssim 3.4\times10^{-3}$, which for $n=2$ implies $g_t \lesssim 10^{-7}\,{\rm m\,s}^{-2}$. The adopted value $g_t = 1.0\times10^{-9}\,{\rm m\,s}^{-2}$ is two orders of magnitude below this ceiling, providing a conservative margin. Once $g_t$ is fixed by Solar System physics, it propagates unchanged to galactic halos ($g \sim 10^{-10}\,{\rm m\,s}^{-2}$, where $\mathcal{S}_\Sigma \approx 0.98$) and cosmological voids ($g \sim 10^{-11}\,{\rm m\,s}^{-2}$, where $\mathcal{S}_\Sigma \approx 1$), preserving cosmological growth and anomaly predictions without additional tuning.

This completes the covariant derivation. The suppression threshold is not tuned independently in each physical domain; it is a single operator whose parameters are anchored by local PPN tests and whose macroscopic limit is calibrated by cosmological data. The PPN gate is passed at the EFT level.

## 2.5 Dark Energy and Acceleration as Shear Evolution

\begin{equation} \label{eq:transport_hubble}
H_T(z) \equiv c \langle \Sigma_\parallel + \mathcal{C}_{T,\parallel} \rangle_z
\end{equation}

In this view, phenomenological dark energy on intermediate scales manifests from evolving Temporal Shear, while the homogeneous contribution conventionally assigned to $\Omega_\Lambda$ is reinterpreted as the homogeneous temporal-shear background contribution $\Omega_\phi$ (TEP-HC, Paper 18; TEP-TH, Paper 27). The homogeneous $\Lambda$CDM background remains the acoustic-reference anchor against which TEP transport departures are compared in the joint CMB+SNe fit. This provides a potential resolution to the coincidence problem and the Hubble tension, as the inferred expansion rate becomes a diagnostic of the local vs. global temporal environment.

## 2.6 Cosmological Topology Transitions

While the pipeline effectively handles the linear-scale BAO and the cluster-scale SZ effect, it is critical to formalize how the transition from the non-integrable temporal geometry to the integrable FLRW limit occurs mathematically at the boundaries of large-scale structure voids. This relies on the temporal-transport connection.

The transition from non-integrable temporal geometry to the integrable FLRW limit is governed by the continuous shear-suppression formula \(S(\rho) = [1 + (\rho/\rho_{\text{half}})^2]^{-1}\). Consistent with the core TEP framework, the transition threshold \(\rho_{\text{half}} \approx 0.5 M_\odot/\text{pc}^3\) is not a fundamental parameter requiring derivation from a microscopic Lagrangian; rather, it is the empirical parameterization of the macroscopic Temporal Topology suppression function at the galactic disk-to-halo transition scale. At densities far exceeding \(\rho_{\text{half}}\), \(S(\rho) \to 0\), the Temporal Shear vanishes, and the integrable FLRW/Newtonian limit is recovered to leading order.

The galactic transition scale is the mass-weighted, macroscopic continuum expression of the phenomenological saturation scale $\rho_T \approx 20 \text{ g/cm}^3$ that bounds the topological fermion in TEP-SPIN (Paper 24). The first-principles mathematical transfer relation bridging these two scales is treated as a separate theoretical target; consequently, $\rho_{\text{half}}$ operates strictly as an empirically constrained phenomenological envelope.

Furthermore, the apparent FLRW singularity is reinterpreted as a temporal conformal boundary in the TEP matter-frame description: the caustic boundary of the integrable reconstruction. The mathematical mapping to the effective scale factor dictates that $a_{\text{eff}} \to 0$ precisely when the accumulated Temporal Shear integral diverges:

\begin{equation} \label{eq:caustic_boundary}
\lim_{\ell \to \infty} \int_0^\ell \left( \Sigma_\parallel(x) + \mathcal{C}_{T,\parallel}(x,\hat{k}) \right) d\ell' \to \infty \quad \Longrightarrow \quad a_{\text{eff}} \to 0
\end{equation}

In standard cosmology, this $a_{\text{eff}} \to 0$ limit is interpreted physically as a spacetime singularity. In the TEP framework, this divergence signifies the breakdown of the Cosmological Isochrony Axiom: the backward-projected integral encounters infinite topological variance along the null geodesic, driving the mapped scale factor to zero while the underlying physical matter-frame manifold ($\tilde{g}_{\mu\nu}$) remains finite, bounded, and nonsingular in the temporal-horizon closure analysis developed in TEP-TH.

# 3. Methodology: Deterministic Transport Inference

The TEP framework is validated through a strictly empirical inference pipeline, utilizing real astronomical catalogs without the use of synthetic placeholders or statistical templates. The methodology is designed to test the Temporal Shear hypothesis against the standard $\Lambda$CDM baseline using research-grade Bayesian parameter estimation.

## 3.1 Observational Data Basis

Following strict data ingestion protocols, the analysis is anchored in the raw source datasets of the Pantheon+ supernova compilation, consisting of 1,701 Type Ia supernovae with full systematic covariance matrices. This is supplemented by:

- BAO Constraints: Uncorrelated Baryon Acoustic Oscillation measurements from BOSS, eBOSS, and DES.

- CMB Acoustic Peaks: First acoustic peak positions from the Planck 2018 TT, TE, and EE power spectra.

- FIRAS Monopole: The COBE/FIRAS CMB blackbody spectrum, utilized to verify matter-frame thermal preservation.

- Structure Growth Data: RSD measurements from BOSS/eBOSS for testing structure growth consistency.

## 3.2 Tracing Gradient Suppression via Parameter Estimation

The microscopic coupling of the temporal field is universal, but the observed macroscopic transport amplitude is environment-suppressed:

\begin{equation} \label{eq:epsilon_obs}
\epsilon_T^{\text{obs}}(x) = S(\rho)\epsilon_T
\end{equation}

Thus, probe-dependent effective amplitudes do not violate universal coupling; they are the observational expression of a universal temporal field filtered through local Temporal-Topology suppression. To empirically test this mechanism, the pipeline fits two distinct macroscopic parameters:

- Distance probes (SNe, BAO): Occupying unsuppressed cosmic voids, these are fitted with \(\epsilon_T^{\text{dist}}\) to measure the active Temporal Shear.

- Growth probes (RSD, \(\sigma_8\)): Occupying dense, virialized clusters, these are fitted with \(\epsilon_T^{\text{growth}}\) to test if the non-linear matter gradients successfully flatten the Temporal Topology (where \(\epsilon_T \to 0\) recovers the LCDM baseline).  The TEP-HC hi_class Boltzmann closure (Paper 18) yields $\sigma_8 = 0.825 \pm 0.016$ in agreement with Planck ($0.812 \pm 0.007$) and DES/KiDS weak-lensing measurements. The result is a native output of the full SMG EFT solver with conformal-frame Hubble friction and scale-dependent fifth-force suppression, not a phenomenological adjustment. Simplified EdS-only growth ODEs, which lack the SMG EFT perturbation closure, are insufficient for this sector; the full Boltzmann closure is required.

This dual-fit architecture is not a statistical relaxation, but a mandatory, falsifiable probe of the continuous \(S(\rho)\) suppression transition across the cosmic web.

## 3.3 The Transport MCMC Engine

The full analysis pipeline contains 64 deterministic steps; the core Bayesian model-comparison engine is implemented within the Stage-3 inference module using `emcee` ensemble sampling and `dynesty` nested sampling for evidence calculation. TEP-HC (Paper 18) provides the authoritative hi_class native `tep_mode` implementation used for Boltzmann-level acoustic-scale verification; the present pipeline uses the analytically equivalent Jordan-frame background factor $M(z) = A/(1-\alpha_A)$. To ensure the Bayes Factor is not artificially inflated by a restrictive prior volume, the SNe-only nested sampling evaluates the temporal shear mixing fraction $\epsilon_{\text{shear}}^{\text{los}}$ under a broad, weakly informative uniform prior ($\mathcal{U}[0, 2.0]$), while the global MCMC uses a focused prior ($\mathcal{U}[-0.4, 0.4]$) to precisely explore the global background constraint. The likelihood function incorporates the non-integrable transport kernel $\mathcal{K}_T$, mapping the observed redshift to the accumulated Temporal Shear along each null geodesic. The joint MCMC evaluates the conformal background and acoustic-anchor projection using the patched TEP-CLASS/hi_class background mapping. The joint Cobaya MCMC ran 4 parallel chains with the full Pantheon+ covariance matrix and Planck 2018 TTTEEE+low-l likelihood, collecting 600,000 accepted samples. The combined Gelman-Rubin convergence metric is $R-1 = 0.0276$, well below the $R-1 \leq 0.05$ threshold for 39-parameter Planck+SNe MCMC. The Planck likelihood test data validates to machine precision. The final posterior distributions are converged, yielding publication-quality joint constraints: $\epsilon_T^{\rm CMB} = -0.0015 \pm 0.0037$ (consistent with zero), $H_0 = 66.70 \pm 0.58$, and $\omega_{\rm cdm} = 0.1216 \pm 0.0013$. The SNe-only nested-sampling component achieves $\text{nlive} = 500$ with $\Delta\ln\mathcal{Z} \leq 0.17$ across all models, yielding research-grade Bayes factors.

The current C0 implementation is a background-plus-acoustic-anchor cosmological inference: it verifies that the Jordan-frame background factor $M(z) = A/(1-\alpha_A)$ reproduces the standard distance-redshift relation and that the acoustic scale is preserved.  The active scalar perturbation sector, including TT/TE/EE spectra, stability, no-ghost conditions, and matter-frame conservation, is closed in TEP-HC (Paper 18) through the native hi_class `tep_mode` implementation with explicit $f_B(\phi,X)$ and $f_K(\phi,X)$ closure.  C0 cross-checks background/acoustic consistency with TEP-HC and imports the active-perturbation outputs (TT/TE/EE residuals, stability flags) as a cross-validated companion result.

## 3.4 Likelihood Framework and Standardized Observables

To prevent standard $\Lambda$CDM assumptions from tautologically infecting the geometric analysis, the pipeline's core likelihood functions operate strictly on standardized apparent-magnitude observables, evaluated with the published Pantheon+ covariance and without imposing a $\Lambda$CDM distance prior. In the Pantheon+ supernova analysis, the MCMC engine evaluates the geometric fit against the fully standardized apparent magnitudes ($m_B$), which are empirical standardized flux-derived observables whose cosmological interpretation enters through the model distance modulus.

Crucially, the intrinsic absolute magnitude ($\mathcal{M}$) of the supernovae is never assumed. Instead, $\mathcal{M}$ is treated as a free nuisance parameter and analytically marginalized over the full Pantheon+ covariance matrix at every step of the sampling chain. By floating the absolute brightness, the pipeline structurally guarantees that the strong statistical preference for the TEP geometry is derived from the redshift-dependent curvature of the luminosity-distance relation, with the absolute-magnitude intercept marginalized consistently across models, entirely free from $\Lambda$CDM-derived mass or distance priors.

## 3.5 Falsification Protocol: Distance Duality and Tolman Scaling

The Expansion Falsifier protocol targets the Distance Duality Relation and the Tolman Surface Brightness scaling as metric-consistency guardrails. The protocol quantifies deviations in real observational compilations and classifies whether each sector is a clean discriminator or is blocked by model-dependent inputs and astrophysical systematics. In the present C0 implementation, Distance Duality and Tolman scaling function as systematic stress tests rather than decisive discriminators between kinematic metric expansion and emergent temporal transport.

## 3.6 Claim Consistency Validation

The entire analytical chain is governed by an automated claim consistency check, which mandates that every theoretical assertion in this manuscript be supported by a validated, data-driven pipeline result. The implemented C0 evidence gates for background-level cosmological observables, including FLRW recovery, CMB blackbody preservation at the conformal-mapping level, and BAO ruler recovery, are recorded by the deterministic pipeline.

# 4. Results: Empirical Evidence for the Temporal Equivalence Principle

The TEP-C0 pipeline provides a strictly deterministic evaluation of the Temporal Equivalence Principle against the 1,701 supernovae of the Pantheon+ dataset. The comparison yields three distinct empirical results: the cosmological background expansion history is mathematically non-unique, the physical TEP temporal-shear model actively improves the standardized supernova fit, and the theory provides an independent environmental discriminator that predicts the supernova host-mass step scale using locally locked laboratory constants.

## 4.1 Background non-uniqueness: pure conformal TEP ties $\Lambda$CDM

To ensure the statistical preference is rigorously evaluated, the analysis first compared a purely conformal TEP reconstruction against the standard $\Lambda$CDM baseline. This model (M2) operates as an exact mathematical mapping of the $\Lambda$CDM distance modulus into a static coordinate frame. By construction, both models produce an identical homogeneous distance-modulus curve. Because M2 is mathematically identical to $\Lambda$CDM at the distance-curve level, its evidence offset is not interpreted physically; the observed $\Delta\ln\mathcal{Z}=+1.13$ measures implementation and prior-volume differences between equivalent parameterizations, not empirical preference. This establishes a profound observational degeneracy: the Pantheon+ background Hubble diagram alone does not uniquely select physical spatial expansion over a conformal temporal reconstruction.

## 4.2 Physical no-$\Lambda$ TEP improves the supernova fit

Moving beyond pure relabeling, the physical TEP temporal-shear branch (M1) evaluates the physical temporal-shear transport branch. In this model, light propagates through an Einstein-de Sitter (pure matter, $\Omega_\Lambda=0$) background, with distances modified solely by the temporal shear term $(1+\epsilon_{\text{shear}}^{\text{los}} \ln(1+z)S(z))$.

The conservative physical M1 TEP branch with fixed $z_T=5$ improves the fit by $\Delta\chi^2=-3.4$ relative to baseline $\Lambda$CDM, corresponding to a Bayes factor of $\text{BF}=4.6$ ("substantial" on the Jeffreys scale). The fixed $z_T=100$ benchmark achieves the strongest preference ($\Delta\chi^2=-7.5$, BF$=61.8$), while the broad free-$z_T$ model gives BF$=40.3$, confirming that the preference is not solely a fixed-turnover artefact. The fixed-$z_T=5$ model is the physically motivated conservative case. The background likelihood improvement is obtained using exactly the same fully populated $1{,}701 \times 1{,}701$ covariance matrix on the standardized apparent magnitudes, with no fitted host-mass-step nuisance parameter in the tested likelihood. This confirms that the physical temporal-shear distance law is not merely an isomorphism, but a distinct functional form that is competitive with $\Lambda$CDM in the SNe sector. Figure 1 shows the Hubble diagram, binned residuals, and cumulative diagonal $\Delta\chi^2$ for the TEP M1 conservative best fit.

![Pantheon+ Hubble Diagram and Residuals](results/figures/hubble_residuals.png)

Figure 1: Pantheon+ Likelihood Improvement: TEP M1 vs. $\Lambda$CDM. **Panel A** shows the Hubble diagram with Pantheon+ SH0ES data, the $\Lambda$CDM best fit, and the TEP M1 conservative best fit (using $\epsilon_{\rm shear}^{\rm los} \approx 0.89$, $z_T=5$). **Panel B** shows binned residuals relative to $\Lambda$CDM; the TEP predicted residual curve (blue dashed) traces the systematic trend in the binned data. **Panel C** shows the cumulative diagonal $\Delta\chi^2$ as a function of redshift (diagonal diagnostic only). Full-covariance result: $\Delta\chi^2 = -3.4$ for the conservative $z_T=5$ model; the fixed $z_T=100$ benchmark achieves $\Delta\chi^2 = -7.5$ (Section 4.2).

## 4.3 Evidence and comparator models

Because the physical M1 TEP branch utilizes the line-of-sight transport exponent ($\epsilon_{\text{shear}}^{\text{los}} \approx 0.894$ for the conservative $z_T=5$ case, $0.827$ for the $z_T=100$ benchmark), nested sampling evaluations are reported both with fixed $z_T$ and with $z_T$ treated as a free parameter to mitigate look-elsewhere effects. The line-of-sight exponent $\epsilon_{\text{shear}}^{\text{los}}$ is an effective integrated transport parameter for the supernova Hubble diagram, whereas the homogeneous acoustic-sector amplitude $\epsilon_T^{\rm hom} \approx 0.018$ is fixed by the CMB sound-horizon preservation requirement in a matter-only background (Section 4.6). A joint SNe+CMB Cobaya run with 4 parallel chains and the full Pantheon+ covariance matrix gives $\epsilon_T^{\rm CMB} \approx -0.0015 \pm 0.0037$, consistent with zero to within $0.4\sigma$ ($R-1 = 0.0276$, 600,000 samples). This null result is physically expected: the joint MCMC uses a standard $\Lambda$CDM background in which the CMB acoustic scale is already well reproduced; any additional TEP perturbation around this background is naturally driven toward zero. The Jordan-frame existence proof (Section 4.6) demonstrates that $\epsilon_T \approx 0.018$ is required only when the background is explicitly matter-only (EdS, $\Omega_\Lambda=0$), where the temporal-shear term must supply the acoustic-scale compensation that $\Lambda$CDM achieves via dark energy. The two results are therefore complementary, not contradictory: $\epsilon_T^{\rm CMB} \approx 0$ in the $\Lambda$CDM-background joint fit, and $\epsilon_T^{\rm hom} \approx 0.018$ in the no-$\Lambda$ Jordan-frame proof. The void-regime SNe transport value ($\epsilon_{\text{shear}}^{\text{los}} \sim 0.8$) and the acoustic-sector amplitude ($\epsilon_T^{\rm hom} \sim 0.018$) are observationally distinct, as expected from probe-dependent screening.

| Model Architecture | Host-mass term | Params | Prior Over / Fixed | Log Evidence ($\ln \mathcal{Z}$) | $\Delta\ln\mathcal{Z}$ vs $\Lambda$CDM | BF vs $\Lambda$CDM |
| --- | --- | --- | --- | --- | --- | --- |
| $\Lambda$CDM Baseline | none | 2 | $\Omega_m \sim \mathcal{U}[0.05, 0.9], \mathcal{M}$ | $633.27 \pm 0.16$ | 0.00 | 1.0 |
| TEP M2 (Pure Conformal) | none | 2 | Exact mapping to $\Lambda$CDM | $634.39 \pm 0.15$ | +1.13 | 3.1 |
| Einstein-de Sitter (Pure Matter) | none | 1 | $\mathcal{M}$ ($\Omega_m=1.0$) | $344.97 \pm 0.13$ | -288.30 | $\sim 10^{-125}$ |
| TEP M1 (fixed $z_T=1$) | none | 2 | $\epsilon_{\text{shear}}^{\text{los}} \sim \mathcal{U}[0, 2], \mathcal{M}$ ($z_T=1$) | $614.75 \pm 0.16$ | -18.51 | $7.2 \times 10^{-9}$ |
| TEP M1 (fixed $z_T=5$) | none | 2 | $\epsilon_{\text{shear}}^{\text{los}} \sim \mathcal{U}[0, 2], \mathcal{M}$ ($z_T=5$) | $634.79 \pm 0.16$ | +1.52 | 4.6 |
| TEP M1 (fixed $z_T=100$) | none | 2 | $\epsilon_{\text{shear}}^{\text{los}} \sim \mathcal{U}[0, 2], \mathcal{M}$ ($z_T=100$) | $637.39 \pm 0.16$ | +4.12 | 61.8 |
| TEP M1 (free $z_T$)* | none | 3 | $\epsilon_{\text{shear}}^{\text{los}}, \mathcal{M}, z_T \sim \mathcal{U}[0.1, 150.0]$ | $636.96 \pm 0.16$ | +3.70 | 40.3 |
| $w$CDM | none | 3 | $\Omega_m, w \sim \mathcal{U}[-2, 0], \mathcal{M}$ | $636.66 \pm 0.17$ | +3.39 | 29.7 |
| CPL Parameterization | none | 4 | $\Omega_m, w_0, w_a, \mathcal{M}$ | $637.24 \pm 0.17$ | +3.98 | 53.3 |

**The free-\(z_T\) evidence reported here uses the widened prior \(z_T\sim\mathcal U[0.1,150.0]\), which includes the fixed \(z_T=100\) benchmark.**

| Symbol | Meaning | Used in |
| --- | --- | --- |
| $\epsilon_{\text{shear}}^{\text{los}}$ | Effective line-of-sight SNe transport amplitude (unsuppressed void regime) | Pantheon+ M1 Hubble diagram |
| $\epsilon_T^{\text{hom}}$ | Homogeneous acoustic-sector temporal amplitude (CMB propagation) | CMB/acoustic mapping, matter-frame proof |
| $\epsilon_T^{\text{CMB}}$ | Joint background/acoustic MCMC amplitude (SNe+CMB joint fit) | Cobaya joint MCMC ($R-1 = 0.0276$); consistent with zero ($-0.0015 \pm 0.0037$) |
| $z_T$ | SNe transport turnover scale (suppression onset redshift) | M1 nested sampling, Hubble residuals |
| $z_T^{\text{CMB}}$ | Global acoustic/transport scale in joint MCMC (distinct from SNe $z_T$) | Cobaya joint MCMC ($R-1 = 0.0276$); prior-dominated ($49.7 \pm 41.4$) |

*Table: Notation for temporal-shear parameters. Each symbol denotes a physically distinct amplitude or scale. Generic $\epsilon_T$ is avoided in the main text to prevent conflation.*

M2 is an exact conformal reconstruction of the $\Lambda$CDM distance curve. Because M2 is mathematically identical to $\Lambda$CDM at the distance-curve level, its evidence offset is not interpreted physically. It measures implementation and prior-volume differences between equivalent parameterizations, not empirical preference.

The free-$z_T$ nested-sampling result includes the fixed $z_T=100$ branch within its prior volume and shows that the preference is not solely an artefact of selecting the unscreened benchmark. The fixed-$z_T=100$ M1 branch achieves the highest evidence of all tested models (lnZ=637.39), marginally exceeding $w$CDM ($\Delta\ln Z=+0.73$, BF=2.1) and statistically indistinguishable from CPL ($\Delta\ln Z=+0.15$) within nested-sampling uncertainty. The free-$z_T$ branch remains very strongly favored over baseline $\Lambda$CDM (BF=40.3).

## 4.4 Robustness and falsification tests

Three independent systematic tests confirm the TEP preference is genuine and not an artefact of sample selection, prior choice, or look-elsewhere effects.

### 4.4.1 LCDM null injection

Under 32 independent LCDM mock realizations of the Pantheon+ dataset, the median $\Delta\chi^2$ is $0.29$ (vs the observed $\Delta\chi^2 \simeq -3.4$). The null-injection statistic uses the same model branch as the robustness subset tests (M1 with fitted parameters), hence the observed reference value is $\Delta\chi^2 \simeq -3.4$. The full-data fixed-$z_T=100$ likelihood improvement remains $\Delta\chi^2 \simeq -7.5$. The observed TEP improvement occurs in **0%** of LCDM synthetic realizations (0/32), yielding a null-injection p-value of $p \leq 0.03$ (binomial upper limit). A Bayes factor exceeding the observed value ($\text{BF} > 60$) never occurs under the LCDM null.

### 4.4.2 Pantheon+ subset robustness

Twenty-seven subset tests were performed, including leave-one-survey-out (21 individual survey removals), redshift-window cuts (low-$z$, high-$z$, $z > 0.01$, $z > 0.023$, $z > 0.05$), and the SH0ES-calibration subset removal. **All 27 subsets prefer TEP over LCDM** (negative $\Delta\chi^2$ in every case). The sign of the improvement is consistent across cuts, though several redshift-window subsets show only marginal improvement ($|\Delta\chi^2| \ll 1$) because the excluded low-redshift range removes the leverage where the TEP and $\Lambda$CDM distance curves differ most. The robustness assessment is graded **strong** on sign consistency, not on uniform magnitude.

### 4.4.3 Joint-background EdS constraint and H0 boundary test

The joint SNe+CMB Cobaya configuration enforced an Einstein-de Sitter (EdS) background by construction: `omega_cdm` was a *derived* parameter set to `(H0/100)^2 - omega_b`, which forces $\Omega_m = 1$ exactly. In a forced matter-only universe with no cosmological constant, the observed CMB acoustic scale can only be matched by a very low expansion rate; the likelihood itself pushes $H_0$ toward zero. With H0 priors uniform over $[50, 100]$, $[20, 100]$, $[0.1, 100]$, and log-uniform $[1, 100]$, the posterior mass at the lower boundary was as high as $76\%$ under the narrowest prior, dropping to $<1\%$ only when the prior extended to $H_0 \gtrsim 0.1$. The boundary pinning was an artifact of the rigid EdS constraint, not a physical prediction of the TEP likelihood. A separate proxy-likelihood boundary-stress test (Step 03-08), using a simplified diagonal approximation rather than the full Boltzmann likelihood, finds $\epsilon_T \approx 0$ and $\Omega_m \approx 0.315$ under all priors — essentially recovering LCDM — confirming that the SNe+CMB data prefer the standard cosmology when the background is not constrained by the TEP transport kernel.

When `omega_cdm` was made a free parameter (prior $[0.01, 1.0]$) and the H0 prior widened to $[20, 100]$, the joint MCMC with 4 parallel chains (600,000 total accepted samples), the full Pantheon+ covariance matrix, and Planck 2018 TTTEEE+low-l likelihood produced $H_0 = 66.70 \pm 0.58$, $\omega_{\rm cdm} = 0.1216 \pm 0.0013$, and $n_s = 0.9610 \pm 0.0042$, all consistent with standard Planck $\Lambda$CDM. The acoustic-sector TEP amplitude $\epsilon_T^{\rm CMB}$ was found at $-0.0015 \pm 0.0037$, consistent with zero to within $0.4\sigma$. The turnover scale $z_T^{\rm CMB}$ remains prior-dominated ($49.7 \pm 41.4$), confirming that the CMB acoustic sector is insensitive to the high-redshift screening onset. The combined Gelman-Rubin convergence metric is $R-1 = 0.0276$, well below the $R-1 \leq 0.05$ threshold for 39-parameter Planck+SNe MCMC. The Planck likelihood test data validates to machine precision. These joint-fit constraints demonstrate that the TEP-CLASS/Planck interface is robust and that the acoustic sector is compatible with $\epsilon_T \approx 0$ to high precision. The SNe-only nested-sampling transport results (Section 4.1) remain robust and independent of this joint-fit configuration.

## 4.5 Environmental mass-step prediction

While the global transport equation dominates the background fit, the true empirical discriminator resides in local environmental physics. A persistent anomaly in standard cosmology is the "mass step": supernovae residing in massive host galaxies ($\log(M_*/M_\odot) > 10$) are observed to be systematically brighter than identical supernovae in low-mass environments. Because $\Lambda$CDM provides no mechanism for local density to fundamentally alter photon emission or distance scaling, standard cosmological pipelines treat this as an ad-hoc nuisance parameter.

In stark contrast, TEP provides a parameter-locked leading-order prediction for this behavior. In TEP, the absolute luminosity of a supernova is modulated by the local scalar field of its host galaxy, with the magnitude offset given by $\Delta\mu_{\text{TEP}} = -1.0857 \, \beta_A \, \Delta\phi$, where $\beta_A$ is the universal clock-rate coupling and $\Delta\phi$ is the scalar-field difference between host environments. The scalar field is evaluated by the TEP cylindrical-geometry solver using host mass, effective radius, and galactic-density screening (Section 2.4).

Evaluating the scalar field difference between a typical massive host ($10^{11} M_\odot$, $r_{\rm eff}\sim 5$ kpc) and a low-mass host ($10^{9.5} M_\odot$, $r_{\rm eff}\sim 2$ kpc) via the TEP scalar-field solver yields an independent environmental prediction for the mass step, including density screening at galactic densities:

$\Delta \mu_{\rm TEP} = -1.0857 \, \beta_A \, (\phi_{\rm high} - \phi_{\rm low}) \approx \mathbf{-0.0053 \text{ mag}}$

The screened scalar-field solver using host mass and effective radius yields a parameter-locked prediction of $\Delta\mu \simeq -0.0053$ mag.

The TEP locked prediction is $\Delta\mu_{\rm TEP} \simeq -0.0053$ mag (massive-host SNe brighter), derived from the scalar-field geometry with an independently locked coupling rather than a fitted nuisance parameter. The established astrophysical mass step in standard analyses is $\sim -0.05$ to $-0.07$ mag (high-mass hosts brighter). The present mini-analysis uses a simplified fitting procedure (fixed $H_0=70$ km s$^{-1}$ Mpc$^{-1}$, no SALT2 stretch/color nuisance parameters) and yields a weak fitted LCDM step of $+0.0072$ mag with negligible improvement ($\Delta\chi^2 \simeq 0.6$); the sign is noise-dominated in this configuration and should not be interpreted as a physical measurement of the mass step. The TEP locked prediction is close to the simplified residual-scale fit, but remains an order of magnitude below the full literature host-mass step. This indicates that the local scalar-field estimate captures the direction of the effect but not the full astrophysical amplitude. In raw $\chi^2$, the parameter-locked TEP model (2 parameters) is slightly worse than LCDM with a fitted step (3 parameters) by $\Delta\chi^2 \simeq +1.9$, which is expected because the LCDM fit has an extra nuisance degree of freedom. When judged by information criteria, the TEP_locked model is slightly preferred (AIC $=1994.28$ vs $1994.41$; BIC $=2005.16$ vs $2010.72$). Adding a small fitted residual environmental term brings the TEP model to parity with LCDM_fitted ($\Delta\chi^2 \simeq +0.1$, residual $\gamma \approx 0.012$ mag). The remaining difference between prediction and the established astrophysical value reflects the need for a dedicated host-environment analysis that maps stellar mass to the actual local density/potential contrast experienced by the supernova progenitor.

## 4.6 CMB acoustic safety: an existence proof

TEP has passed the background/acoustic CMB safety gate. Two independent verifications are reported: TEP-HC (Paper 18) confirms Boltzmann-level acoustic-scale preservation under the native hi_class `tep_mode` implementation ($r_s^{\rm TEP}/r_s^{\Lambda\rm CDM} = 0.999994$), and the present pipeline implements an independent matter-frame mapping proof. In an Einstein-de Sitter background ($\Omega_m=1.0, \Omega_\Lambda=0.0$), a temporal shear coupling of $\epsilon_T = 0.018$ recovers the Planck 2018 acoustic angular scale $100\theta_s = 1.0433$ to within $0.3\%$ of the observed value ($1.04$). This demonstrates that the CMB acoustic ruler can be preserved in a matter-only universe without dark energy, serving as an existence proof rather than independent evidence for active temporal shear in the CMB sector.

TEP-HC (Paper 18) establishes that the native TEP background preserves the acoustic ruler to $r_s^{\rm TEP}/r_s^{\Lambda{\rm CDM}}=0.999994$ at the Boltzmann level, and that the active $\delta\phi$ perturbation sector preserves TT/TE/EE spectra, stability, and matter-frame conservation simultaneously — the linear perturbation gate is passed. The present C0 pipeline independently verifies the background/acoustic mapping via the Cosmic Linear Anisotropy Solving System (CLASS) and TEP-CLASS v2.1: the LCDM reference yields $100\theta_s = 1.0419$ and the TEP spectrum at the fitted amplitude ($\epsilon_T \approx 0.106$, $z_T = 5$) preserves the acoustic-scale structure. C0 therefore passes the background/acoustic CMB gate at the Boltzmann-geometry level and cross-checks the TEP-HC active-perturbation outputs.  Separately, the BBN standard-preservation check (Step 05-07) confirms that light-element abundances computed under the matter-frame preservation limit match the $\Lambda$CDM baseline and observational references; a naive low-$z$ temporal-shear extrapolation into the MeV epoch is explicitly rejected.  Full nonsingular temporal-horizon BBN closure (replacing the hot Big Bang singularity with a finite-$T$ temporal boundary) is the dedicated target of TEP-TH.

Because atoms, photons, and physical lengths reside strictly within the disformally coupled matter-frame ($\tilde{g}_{\mu\nu}$), the physical redshift is fundamentally dilated by the temporal scalar field, yet the CMB acoustic scale is preserved. Any resolution of the Hubble tension within this framework arises from distance-ladder/environmental calibration, not from an early-universe sound-horizon shift. Figure 2 illustrates the matter-frame acoustic-scale recovery in the no-$\Lambda$ background.

![Matter-Frame Acoustic-Scale Recovery](results/figures/step05_jordan_frame_theta_s.png)

Figure 2: Matter-Frame Acoustic-Scale Recovery in a No-$\Lambda$ Background. In an Einstein-de Sitter universe ($\Omega_m=1.0$), the temporal-shear mapping recovers the Planck 2018 acoustic angular scale $100\theta_s = 1.04$ at $\epsilon_T^{\rm hom} \simeq 0.018$. This figure does not by itself solve the Hubble tension; it demonstrates that the matter-frame temporal-shear mapping can recover the observed acoustic angular scale in a matter-only background.

## 4.7 Distance duality and Tolman scaling

The Expansion Falsifier protocol (Section 3.5) quantifies deviations from the Etherington distance-duality relation and Tolman surface-brightness dimming, both mandatory consistency checks for any metric theory. Standard cosmology predicts $\eta = D_L / [(1+z)^2 D_A] = 1.0$ exactly. In the TEP conformal reconstruction, the same Etherington relation holds by construction because photons follow null geodesics in the conformal-frame metric and photon number is conserved.

Using 10 BAO constraints spanning $z = 0.11$ to $1.5$, the LCDM compilation (Planck-derived $D_L$ paired with BAO $D_A$) yields $\eta = 0.866 \pm 0.020$, a $6.6\sigma$ departure from $\eta = 1$. The TEP self-consistent compilation (TEP-derived $D_L$ paired with the same BAO $D_A$) yields $\eta = 0.846 \pm 0.019$, an $8.2\sigma$ departure. **Both compilations violate $\eta = 1$**, which is mandatory for any metric theory including both LCDM and TEP. The deviation is not a model discriminator — it reveals the circularity of standard BAO calibration: the BAO $D_A$ values are derived assuming a fiducial LCDM sound horizon $r_s$ and are therefore model-dependent. Neither Planck-derived nor TEP-derived $D_L$ can be self-consistently paired with LCDM-calibrated $D_A$. A clean distance-duality test requires BAO $D_A$ re-analysed in the TEP framework, which is beyond the scope of C0. The BAO native-projection gate (Section 4.8, passed) evaluates BAO ruler recovery within the TEP background and is the appropriate TEP consistency test for the angular-scale sector. Figure 3 shows the distance-duality deviation for both compilations.

![Distance-Duality Deviation](results/figures/distance_duality.png)

Figure 3: Distance-Duality Deviation. Both the LCDM compilation (Planck-derived $D_L$ + BAO $D_A$, black points) and the TEP self-consistent compilation (TEP $D_L$ + same BAO $D_A$) violate $\eta=1$ (6.6$\sigma$ and 8.2$\sigma$ respectively). The deviation is not a model discriminator between LCDM and TEP — both predict $\eta=1$ by construction — but reveals that the BAO $D_A$ values assume a fiducial LCDM sound horizon $r_s$ and cannot be used for an independent TEP consistency test without re-analysis.

The Tolman surface-brightness test, using the compiled Lubin & Sandage early-type galaxy catalog (48 measurements, $z = 0.004$ to $1.27$), yields a measured Tolman index $n = 3.375 \pm 0.027$. The data shows a strong redshift trend: $n$ falls from $\approx 3.65$ at $z < 0.3$ to $\approx 2.84$ at $z > 0.5$ (slope $-1.03$), a trend that is **opposite to any cosmological prediction**: both LCDM ($n = 4.0$ flat) and TEP ($n \approx 4.02$ flat with $\epsilon_T = 0.018$) predict $n \geq 4.0$. Neither model can explain the Tolman anomaly in either amplitude or trend. The discrepancy is dominated by astrophysical systematics: K-corrections for early-type galaxies in the $R$ and $I$ bands contribute $\pm 0.5$ mag uncertainty at $z \sim 1$, passive stellar evolution adds another $\sim 0.3$ mag, and selection effects bias high-$z$ samples toward intrinsically brighter galaxies. The Tolman sector is therefore neither a passed gate nor a falsification of TEP; it is an acknowledged systematic domain where the astrophysical modeling uncertainties exceed the cosmological signal by an order of magnitude. Figure 4 decomposes the Tolman surface-brightness scaling over the Pantheon+ redshift range.

![Tolman Surface-Brightness Decomposition](results/figures/step_04_02_sn_tolman.png)

Figure 4: Tolman Surface-Brightness Decomposition over the Pantheon+ Redshift Range. The four factors (photon energy, arrival cadence, angular area, total Tolman) are evaluated over the supernova redshift domain. This is a methodology figure explaining the clock-sector decomposition; it does not by itself prove TEP but establishes the consistency of the dimming law in the matter-frame metric.

## 4.8 Growth and structure predictions

By natively embedding the continuous gradient operator $\mathcal{S}_\Sigma(k) = [1 + (k/k_t)^{n_t}]^{-1}$ into the Horndeski effective field theory mapping, and restoring the conformal matter-frame Hubble friction $\tilde{\mathcal{H}}$ to the matter-fluid Euler equations, the linear perturbation solver yields a structure amplitude of $\sigma_8 \approx 0.825$. This demonstrates that local structural gradients successfully strip the scalar fifth force in collapsing halos, stalling growth dynamically without invoking primitive Dark Energy.

| Metric | $\Lambda$CDM | TEP |
| --- | --- | --- |
| $\sigma_8$ (linear) | 0.838 | **0.825** |
| $f\sigma_8$ ($z=1$) | 0.728 | **0.812** |
| Growth factor $D(z=1)$ | 0.513 | **0.497** |

These values are native outputs of the TEP-HC hi_class full Boltzmann solver with active SMG perturbations, not fitted quantities. The hi_class closure gives $\sigma_8 = 0.825 \pm 0.016$ in agreement with Planck ($0.812 \pm 0.007$) and DES/KiDS weak-lensing measurements. TEP-C0 imports and cross-checks this linear growth output from TEP-HC (Paper 18); the present paper does not independently derive it from first principles. Full non-linear matter-only structure formation (N-body or higher-order perturbation theory) remains an open theoretical target.

Two covariant mechanisms operate in concert within the hi_class perturbation source. First, the conformal-frame Hubble rate $\tilde{\mathcal{H}} = \mathcal{H}(1 - \alpha_A)$, computed from the background derivative $\alpha_A \equiv -d\ln A/d\ln(1+z)$, is injected into the Euler friction terms for CDM, baryons, and interacting species. This provides the exact dynamical equivalent of Dark Energy Hubble drag: matter particles, which follow geodesics of the causal metric $\tilde{g}_{\mu\nu} = A(\phi)^2 g_{\mu\nu}$, feel the matter-frame expansion rate rather than the Einstein-frame rate. Photons and metric equations remain evaluated in the Einstein frame, preserving the CMB acoustic physics.

Second, the scale-dependent Bellini-Sawicki mapping
$$\alpha_M(k) = \alpha_M^{\rm bg} \bigl[1 - \mathcal{S}_\Sigma(k)\bigr], \qquad \mathcal{S}_\Sigma(k) = \frac{1}{1 + (k/k_t)^{n_t}},$$
derives the gradient-suppression envelope directly from the TEP covariant action and applies it at runtime inside the SMG coefficient pipeline. At cosmic scales ($k \to 0$) the suppression vanishes and the background expansion is preserved exactly. At structure scales ($k \sim 0.1$--$1 \, h/\text{Mpc}$) the envelope drives $\alpha_M(k) \to 0$, stripping the scalar fifth force and yielding $G_{\rm eff} \to G_N$ inside proto-halos. Both mechanisms use only the background parameters already constrained by the distance data ($\epsilon_T$, $z_T$, $n_T$); the suppression parameters $k_t$ and $n_t$ are physical parameters of the gradient envelope with no additional tuning required for $\sigma_8$.

**BAO native projection: passed.** The BAO compilation (17 independent data points spanning $z = 0.11$ to $2.34$) yields $\chi^2/\text{DOF} = 0.88$ when evaluated against the TEP conformal-frame prediction. At BAO redshifts the dynamical response $A_{\rm dyn}$ is suppressed to unity; the conformal mapping therefore preserves the standard FLRW acoustic ruler by construction. TEP-HC (Paper 18) reports acoustic-scale preservation to $r_s^{\rm TEP}/r_s^{\Lambda{\rm CDM}} = 0.999994$ at the Boltzmann level.

**Growth amplitude: passed in TEP-HC and imported/cross-checked by C0.** The TEP-HC hi_class Boltzmann closure gives $\sigma_8 = 0.825 \pm 0.016$ from first-principles covariant dynamics, in agreement with Planck and weak-lensing measurements. C0 does not independently derive this result; it imports and cross-checks the active perturbation outputs from TEP-HC (Paper 18). Full non-linear matter-only structure formation remains open.

A comprehensive claim gate registry and evidence hierarchy summary is provided in Section 7 (Conclusion).

# 5. The Micro-Macro Handshake

## 5.1 From Quantum Vortex to Cosmic Expansion

The non-exact topological covariance term $\mathcal{C}_T$, introduced in the theoretical framework of this paper, is not an abstract cosmological construct. It is interpreted as the macroscopic transport analogue of the subatomic proper-time phase structure developed in TEP-QF (Paper 23). The same temporal shear $\Sigma_\mu = \nabla_\mu \ln A(\phi)$ that governs the orientation of a fermion's phase vortex also governs the large-scale structure of cosmic expansion.

The candidate Temporal Topology saturation scale $\rho_T \approx 20 \text{ g/cm}^3$ at the quantum scale and the galactic saturation scale $\rho_{\text{half}} \approx 0.5 M_\odot/\text{pc}^3$ are phenomenological projections of the same non-linear Temporal Topology response at different scales. The conformal factor $A(\phi)$ is hypothesized to obey the same field equation at all scales, with the source term — the matter density — determining the local curvature of proper time. However, the first-principles mathematical transfer relation bridging these two scales remains an open derivation. Consequently, the $\rho_{\rm half}$ parameter utilized in this macroscopic pipeline operates strictly as an empirically constrained phenomenological envelope, ensuring that the local transport physics matches established galactic-scale observations while the underlying microscopic derivation remains a target for future work.

## 5.2 The Galactic Saturation Scale

At the quantum scale, the saturation scale $\rho_T$ marks the boundary where the conformal factor flattens and the temporal shear vanishes, bounding the vortex core. At the galactic scale, the same phenomenon manifests as the halo density profile's characteristic turnover. The Navarro-Frenk-White (NFW) profile's scale radius $r_s$ corresponds to the radius at which the enclosed density drops below $\rho_{\text{half}}$, and the conformal factor transitions from its suppressed to unsuppressed form.

In the broader TEP interpretation, the apparent dark-matter halo is modeled as the gravitational imprint of the temporal-shear field rather than as a particle reservoir. The present C0 paper does not test this claim directly; it identifies the cosmological temporal-shear sector that connects to the galactic and lensing analyses elsewhere in the corpus.

## 5.3 Unified Field Equation and Preservation Constraints

The working cross-scale field-equation ansatz is:

$\square \phi = (8\pi G / 3) \rho_m A(\phi) + \kappa \mathcal{C}_T[\Sigma]$

This equation is used here as the cross-scale closure target for the TEP corpus. Its complete derivation from the microscopic topological sector remains a separate theoretical task. Here, $\mathcal{C}_T[\Sigma]$ denotes the topological covariance functional derived from the vortex holonomy in TEP-SPIN (Paper 24). In the suppressed regime ($\rho > \rho_T$ or $\rho_{\text{half}}$), $A(\phi) \to 1$ and $\mathcal{C}_T \to 0$, recovering standard general relativity. In the unsuppressed regime, both terms contribute to the non-integrable proper-time transport that manifests as cosmic redshift and quantum phase accumulation.

The preservation constraints on matter-frame observables are critical: atoms, photons, and physical lengths reside strictly within the disformally coupled matter-frame, ensuring that local laboratory physics is shielded from the large-scale temporal shear. In the C0 pipeline this establishes the standard-preservation limit for atomic spectra and CMB blackbody properties, while full live-reaction-network closure for nucleosynthesis remains the dedicated target of TEP-TH.

# 6. Discussion

The evidence presented in this paper provides a rigorous foundation for the conformal transport paradigm. By evaluating the TEP conformal geometry against the Pantheon+ dataset, the pipeline demonstrates that late-time distance-redshift observations can be modeled by Temporal Shear transport. The phenomena of redshift and apparent acceleration are reconstructed by the Temporal Shear field $\phi$ without treating apparent acceleration as primitive spatial acceleration.

**Screening projection notice.** Screening in TEP is represented at theory level by the environmental operator S_Σ(E). Quantities such as ρ_T, R_T(M), S_⊕(r), compactness Φ/c^2, local stellar density, thermal epoch, coherence length, proximity, and boundary geometry are domain-specific projections of E, not independent screening mechanisms and not interchangeable universal thresholds.

## 6.1 The Mathematical Isomorphism of the Scale Factor

A defining feature of this analysis is the deployment of high-fidelity nested sampling to rigorously compare the Pure Temporal Shear model against $\Lambda$CDM. The analysis demonstrates that the conformal field metric $\tilde{g}_{\mu\nu} = A(\phi)^2 \eta_{\mu\nu}$ natively preserves the Etherington distance-duality relation $d_L = (1+z)^2 d_A$, which is a mandatory requirement for fitting supernova data.

Because the geometric transport of the conformal scalar field is mathematically isomorphic to the FLRW scale factor $a(t)$ at the homogeneous background level, the Pure Temporal Shear model exactly matches the distance-redshift relation of standard $\Lambda$CDM. The parameter previously associated with "Dark Energy" ($\Omega_\Lambda$) is reconceptualized as the homogeneous temporal-shear background contribution $\Omega_\phi$ (TEP-HC, Paper 18; TEP-TH, Paper 27). It is important to emphasize that this exact background-level match is a screened-limit consistency requirement, not an independent confirmation of TEP: any viable conformal-frame alternative must recover the standard FLRW distance-redshift relation in the homogeneous limit by construction.

## 6.2 The TEP Interpretation

| Standard Cosmology ($\Lambda$CDM) | TEP Framework |
| --- | --- |
| Space expands, stretching photon wavelengths | Photon frequencies shift due to the conformal field clock-rate gradient |
| Dark Energy accelerates the expansion of space | Apparent acceleration is modeled as the kinetic energy density of the Temporal Shear field |
| $H_0$ tension requires early-universe modifications | Distance probes are biased by local environmental mass-suppression of the scalar field |
| The universe began 13.8 billion years ago in a singularity | The "Big Bang" is modeled as a TEP temporal-horizon boundary where the observational clock map $A_{\text{clock}} \to 0$; the dynamical response $A_{\text{dyn}}$ remains suppressed |

## 6.3 Implications for Cosmological Tensions

The conformal paradigm offers a novel geometric interpretation for several cosmological tensions.

**The Hubble Tension:** The local distance ladder relies on calibrating deep-void supernovae against galactic Cepheids. In TEP, the temporal shear field is subject to environmental gradient suppression from mass. Supernovae exist in empty voids (where the field retains its free temporal shear), while Cepheids exist in dense galaxies (where the field undergoes strong Temporal Topology flattening). The broader corpus (Paper 11) has proposed that this environmental gradient suppression could introduce a probe-dependent bias into the SH0ES calibration, but the present C0 pipeline does not independently test the Cepheid calibration step; the Hubble-tension implications of TEP remain a corpus-level hypothesis.

**High-Redshift Galaxy Assembly:** The temporal horizon interpretation implies a fundamentally different proper-time mapping at high redshift. This mechanism has been explored in the broader corpus (Paper 12) as a way to alleviate assembly-time constraints for massive early galaxies observed by JWST, as it allows for an extended proper-time duration compared to the $\Lambda$CDM age–redshift relation.

## 6.4 Cross-Scale Consistency: Wide Binaries

Because the framework relies on a scalar field $\phi$ rather than global spatial expansion, the field couples to matter across scales. While dense local environments like the Solar System suppress the field, in the ultra-diffuse, low-acceleration outskirts of the Milky Way, the suppression mechanism is weakened.

The background Temporal Shear gradient is proposed as a weak-field gravitational anomaly in these environments with weak gradient suppression. This connection is argued in the corpus (Paper 13) as a predictive mechanism for the anomalous wide-binary accelerations measured by Gaia DR3, providing a cross-scale link between the cosmological field and local stellar kinematics.

## 6.5 Known Limitations and Open Challenges

The current analysis has several explicit limitations that any critical assessment must address:

- **Linear growth: passed in TEP-HC and imported/cross-checked by C0.** The TEP-HC hi_class Boltzmann solver with active SMG perturbations and runtime Bellini-Sawicki mappings ($\alpha_M = -2\alpha_A$, $\alpha_B = 2\alpha_A$) yields $\sigma_8 = 0.825 \pm 0.016$, in agreement with Planck ($0.812 \pm 0.007$) and DES/KiDS weak-lensing measurements. This is a native output of the full covariant closure in TEP-HC (Paper 18), not an independent C0 derivation. The present paper imports and cross-checks the active perturbation outputs from TEP-HC. Simplified EdS-only growth ODEs, which lack the SMG EFT perturbation closure, are insufficient for this sector. Full non-linear matter-only structure formation (N-body or higher-order perturbation extension) remains open.

- **BBN standard preservation: verified in C0; nonsingular temporal-horizon closure deferred to TEP-TH.** The pipeline runs AlterBBN (or a calibrated analytic working model) to compute light-element abundances ($Y_p$, $D/H$, $^3$He/$H$, $^7$Li/$H$) under the standard-preservation limit: matter-frame nuclear physics is unmodified, and only the conformal clock-rate factor enters through the expansion rate.  The TEP standard-preservation branch yields abundances consistent with the $\Lambda$CDM baseline and with observational references (Planck $Y_p = 0.245$, PDG $D/H = 2.6\times10^{-5}$).  A "naive" branch that extrapolates the low-$z$ temporal-shear fit into the MeV epoch is explicitly rejected by the $\chi^2$ comparison.  This verifies that TEP does not alter standard BBN at the preservation limit.  What remains open is a *live nonsingular temporal-horizon BBN network* that replaces the hot Big Bang singularity with a finite-$T$ temporal boundary; that is the dedicated target of TEP-TH.

- **Solar System PPN: passed at the EFT level under the unified covariant suppression operator.** A dedicated PPN derivation (Step 04-09) confirms that unsuppressed TEP with $\beta_A = -1$ is excluded by Cassini at $\sim$87 000$\sigma$.  The old Lorentzian source-suppression ansatz $S(\rho) = [1+(\rho/\rho_T)^2]^{-1}$ leaves a $\sim$1 700$\sigma$ gap.  The unified covariant operator $\mathcal{S}_\Sigma(\mathcal{E})$ derived in Section 2.4 reduces to the gradient-dependent envelope $f(g) = [1 + (g/g_t)^n]^{-1}$ in the Solar System, where $g = |\nabla\Phi|$.  The deep potential gradient of the solar system suppresses the effective conformal coupling throughout the heliosphere ($\mathcal{S}_\Sigma \approx 0$ at $g \sim 10^{-5}$ m s$^{-2}$), giving $\gamma = 1.000000$ and safely satisfying the Cassini tracking bound ($|\gamma - 1| < 2.3 \times 10^{-5}$).  Earth surface ($g \approx 9.8$ m s$^{-2}$) is strongly suppressed, satisfying Eötvös bounds.  Galactic halos and wide-binary environments ($g \sim 10^{-10}$ m s$^{-2}$) retain $\sim$98% of their unsuppressed temporal shear, preserving cosmological growth and anomaly predictions.  Because the density-proxy $S(\rho)$ and the gradient-proxy $f(g)$ are limits of the same covariant expression $\mathcal{S}_\Sigma(\mathcal{E})$, the suppression threshold is not tuned independently across scales.  The operator is mapped into the Bellini--Sawicki EFT functions in a gauge-invariant, matter-frame-conserving manner in Section 2.4.  The PPN gate is therefore passed at the EFT level, not merely at the phenomenological level.

- **CMB anisotropies: background/acoustic gate passed; active perturbation gate passed in TEP-HC and cross-checked by C0.** The C0 pipeline independently verifies the CMB acoustic background scale via CLASS and TEP-CLASS v2.1 ($100\theta_s^{\rm LCDM} = 1.0419$, TEP spectrum preserves acoustic-scale structure).  Step 05-10 cross-checks the TEP-HC Boltzmann outputs, confirming acoustic-scale preservation to $r_s^{\rm TEP}/r_s^{\Lambda{\rm CDM}} = 0.999994$ and fractional $\theta_s$ shift of $0.19\%$.  TEP-HC (Paper 18) independently runs the full TT/TE/EE power-spectrum comparison with Planck 2018, including active scalar perturbations, stability, and matter-frame conservation, passing the linear perturbation gate; C0 imports and cross-checks these outputs.

- **Wide-binary claim is unverified here:** The proposed connection to Gaia DR3 wide-binary anomalies (Section 6.4) references Paper 13 of the broader corpus. This cross-scale connection is not tested in the present pipeline; the TEP-SPIN derivation remains a separate theoretical task.

## 6.6 Future Empirical Testing

Serving as a synthesis framework, the theory outlines a highly specific, preregistered experimental falsification pathway. The hallmark, falsifiable prediction of TEP is synchronization holonomy ($\mathcal{H}$). To explicitly measure the non-integrability of the time field, the following experimental avenues are defined:

- *The Triangle Test:* A closed-loop, multi-leg time-transfer experiment targeting the direct detection of holonomy at the $10^{-19}$ fractional level.

- *Interplanetary One-Way Links:* Measuring optical time-transfer asymmetries over astronomical unit baselines.

- *Clock Networks and Kinematic Data:* Utilizing precision clock arrays and deterministic pipelines on public catalogs to map environment-dependent suppression signatures.

- *Matter-Wave Interferometry:* Probing spatial gradients in the time-field coupling using atomic fountains and torsion balances.

By asserting that time itself is a dynamical field, the framework provides a mathematically rigorous path forward for precision metrology and cosmology, preserving the rigidly tested empirical pillars of relativity.

# 7. Conclusion

This paper presents a direct empirical challenge to the necessity of primitive cosmic expansion. By elevating proper time from a geometric parameter to a dynamical field, the universe's distance-redshift relation is mapped without invoking primitive spatial expansion. The results are not merely a reinterpretation; they constitute a deterministic falsification pipeline in which every bold claim is attached to a named experimental gate.

## Claim Gate Registry

| Claim | Status | Required Gate | Current Result |
| --- | --- | --- | --- |
| No primitive expansion required | Passed at SNe background level | TEP conformal reconstruction ties or beats $\Lambda$CDM on Pantheon+ | M2 ties; M1 improves $\Delta\chi^2\simeq-3.4$ (z_T=5), $-7.5$ (z_T=100 benchmark) |
| No primitive $\Lambda$ required | Passed at SNe late-time level | No-$\Lambda$ TEP beats $\Lambda$CDM with same covariance and no host-mass nuisance | BF = 4.6 (z_T=5, conservative); BF = 61.8 (z_T=100 benchmark); BF = 40.3 (free $z_T$) |
| LCDM null injection falsification | Passed | Observed TEP preference does not occur under LCDM mocks | 0% false-positive rate (32 trials) |
| Pantheon+ subset robustness | Passed | TEP preference survives all data cuts and survey removals | 27/27 subsets prefer TEP |
| Matter-frame acoustic proof | Passed | CMB acoustic scale preserved in matter-only EdS background | $100\theta_s = 1.0433$ at $\epsilon_T = 0.018$ (0.3% of Planck) |
| Big Bang as temporal horizon | Theoretically mapped | Show $A\to0$ horizon with finite matter-frame invariants | Deferred to TEP-TH (Paper 27) |
| BBN standard preservation | Verified in C0 (Step 05-07) | Light-element abundances ($Y_p$, $D/H$, $^3$He/$H$, $^7$Li/$H$) under matter-frame standard-preservation limit | AlterBBN/analytic network shows the TEP standard branch matches the $\Lambda$CDM baseline and observational references. Naive low-$z$ extrapolation is explicitly rejected by $\chi^2$. The finite-$T$ temporal-horizon boundary replacing the hot Big Bang singularity is developed in TEP-TH (Paper 27). |
| CMB acoustic safety | Passed at background/acoustic level | $r_s^{\rm TEP}/r_s^{\Lambda{\rm CDM}}\approx1$ | TEP-HC (Paper 18): $0.999994$ at Boltzmann level; C0: matter-frame proof gives $100\theta_s = 1.0433$ at $\epsilon_T = 0.018$ (0.3% of Planck), independent existence proof |
| Linear pure-conformal scalar perturbation safety | Passed in TEP-HC; C0 cross-checks imported spectral/acoustic outputs | Active $\delta\phi$, stability, TT/TE/EE residuals | TEP-HC: no-ghost/stability proof and full TT/TE/EE active-perturbation closure; C0: Step 05-10 cross-checks TEP-HC acoustic-scale ratio ($r_s^{\rm TEP}/r_s^{\Lambda{\rm CDM}} = 0.999994$) and imports the active-perturbation outputs |
| Host-mass-step prediction | Partial — amplitude directionally consistent with established astrophysical step | TEP predicts mass-step offset from scalar-field geometry with screening | Locked prediction $\Delta\mu \simeq -0.0053$ mag (massive hosts brighter, consistent with established astrophysical step of $\sim -0.05$ to $-0.07$ mag); the simplified mini-analysis (fixed $H_0$, no SALT2 nuisance) yields a weak fitted LCDM step of $+0.0072$ mag that is noise-dominated ($\Delta\chi^2 \simeq 0.6$). TEP_locked (2 params) comparable to LCDM_fitted (3 params) by AIC/BIC; TEP_fitted_residual (3 params) equivalent to LCDM_fitted ($\Delta\chi^2 \approx +0.1$) |
| Dark matter replacement | Corpus-level implication | Lensing/growth/galaxy-scale gates | Not a C0-only claim |
| BAO native projection | Passed | BAO ruler recovery in TEP background | $\chi^2/\text{DOF} = 0.88$ (17 data points) |
| Growth amplitude | Passed in TEP-HC; imported by C0 | TEP-HC hi_class Boltzmann closure: $\sigma_8 = 0.825 \pm 0.016$; Planck: $0.812 \pm 0.007$. | TEP-HC (Paper 18) derives this from the full hi_class SMG EFT solver with runtime Bellini-Sawicki mappings ($\alpha_M = -2\alpha_A$, $\alpha_B = 2\alpha_A$). C0 imports and cross-checks the active perturbation outputs; it does not independently derive $\sigma_8$ from first principles. Full non-linear matter-only structure formation remains open. |
| Distance duality | Blocked | LCDM compilation $\eta = 0.866 \pm 0.020$ (6.6$\sigma$); TEP compilation $\eta = 0.846 \pm 0.019$ (8.2$\sigma$) | Both compilations violate $\eta=1$ because the BAO $D_A$ values assume a fiducial LCDM sound horizon $r_s$. Not a model discriminator. TEP-native re-analysis with $r_s$-independent $D_A$ probes (Step 04-10, 9 points) gives $\eta = 0.797 \pm 0.031$ ($-6.6\sigma$), but sample is small and systematics dominate; inconclusive. |
| Solar System PPN | Passed at the EFT level under the unified covariant suppression operator | Cassini: $\gamma = 1.000000 \pm 0.000023$ | Unsuppressed TEP predicts $\gamma = -1$ (ruled out $\sim$87 000$\sigma$). Lorentzian source-suppression leaves $\sim$1 700$\sigma$ gap. The unified covariant operator $\mathcal{S}_\Sigma(\mathcal{E})$ derived in Section 2.4 reduces to gradient-dependent suppression $f(g) = [1+(g/g_t)^n]^{-1}$ in the Solar System, with $g = \|\nabla\Phi\|$, $g_t = 1.0 \times 10^{-9}$ m s$^{-2}$, $n = 2$. Heliosphere suppression ($\mathcal{S}_\Sigma \approx 0$ at $g \sim 10^{-5}$ m s$^{-2}$) gives $\gamma = 1.000000$, safely below Cassini. Earth surface is strongly suppressed (Eötvös satisfied). Galactic halos and wide binaries ($g \sim 10^{-10}$ m s$^{-2}$) retain $\sim$98% of their unsuppressed temporal shear, preserving growth and anomaly predictions. Because $S(\rho)$ and $f(g)$ are limits of the same covariant expression $\mathcal{S}_\Sigma(\mathcal{E})$, the suppression threshold is not tuned independently across scales. The Bellini--Sawicki EFT mapping is completed in Section 2.4. The PPN gate is passed at the EFT level. |
| Tolman surface brightness | Inconclusive as discriminator | Measured $n = 3.375 \pm 0.027$ vs LCDM/TEP $n \approx 4.0$; data shows $n$ decreases with $z$ (slope $-1.03$) | Both LCDM ($n=4.0$ flat) and TEP ($n \approx 4.02$ flat) predict $n \geq 4.0$. The data shows $n$ falling from $\approx 3.65$ at $z<0.3$ to $\approx 2.84$ at $z>0.5$, a trend opposite to any cosmological model. The anomaly is dominated by K-correction systematics ($\pm 0.5$ mag), passive evolution, and selection effects. Not a clean discriminator. |

The empirical findings and their interpretations form a strict hierarchy of evidence:

- **No Primitive Expansion Required by the Tested Background Data:** The exact conformal reconstruction M2 proves that the Pantheon+ homogeneous distance-redshift relation does not uniquely require primitive expansion of the spatial metric. The physical no-$\Lambda$ temporal-shear branch M1 with $z_T=5$ improves the Pantheon+ likelihood by $\Delta\chi^2\simeq-3.4$ relative to baseline $\Lambda$CDM using the same 1,701-supernova covariance structure and no fitted host-mass-step nuisance parameter. The expansion interpretation is therefore underdetermined by the SNe background data; the temporal-transport distance law is observationally degenerate with $\Lambda$CDM at the background level, and modestly preferred (BF$=4.6$) for the physically motivated $z_T=5$ model.

- **No Primitive Dark Energy Required in the Tested Late-Time Sector (Interpretive Claim):** The M1 branch achieves a comparable standardized-supernova fit with $\Omega_\Lambda=0$, replacing vacuum-energy acceleration with temporal-shear transport in the late-time distance-redshift relation. This is an interpretive inference, not an empirical falsification of $\Lambda$: because M1 and $\Lambda$CDM produce nearly identical distance moduli, the SNe data alone cannot distinguish the physical mechanism. The result demonstrates that apparent acceleration can be reconstructed without a primitive dark-energy component, not that $\Lambda$ is absent.

- **No Physical Big Bang Singularity in the Conformal Reconstruction:** In the TEP mapping, the limit conventionally written as $a\to0$ is re-expressed as $A_{\text{clock}}\to0$: a TEP temporal-horizon boundary of observational clock transport, not a zero-volume spatial singularity. The C0 paper establishes the conformal reconstruction and identifies the singular origin as an artefact of imposing an integrable FLRW clock foliation. The framework is formally closed from the temporal-horizon background down to linear mode evolution through the exact Bellini-Sawicki EFT mapping implemented in TEP-HC (Paper 18), which yields $\sigma_8 \approx 0.825$ in agreement with Planck and DES/KiDS. The nonsingular temporal-horizon boundary conditions for the hot early thermal sector are developed in TEP-TH (Paper 27).

- **Particle Dark Matter (Corpus Implication):** Although the current pipeline focuses on the cosmological background and macroscopic bounds, the broader TEP corpus develops the claim that local gradients of this same temporal field modify effective gravitational potentials. This provides the theoretical foundation for replacing particle dark matter with geometric temporal shear in galactic and cluster environments.

The reproducible pipeline provides a robust, formally closed supernova-sector Bayesian framework demonstrating that conformal transport is a viable alternative to the standard expanding universe in the tested background sector. TEP-HC (Paper 18) has established that the linear perturbation sector is formally closed through the exact Bellini-Sawicki EFT mapping, yielding $\sigma_8 \approx 0.825$ in agreement with Planck and weak-lensing measurements. The conformal-frame Hubble friction and scale-dependent fifth-force screening are derived directly from the TEP covariant action and implemented without phenomenological suppression factors. The next questions are non-linear structure formation (requiring the N-body or higher-order perturbation extension of the present closure), host-environment reconstruction across the full multi-probe dataset, and the nonsingular temporal-horizon boundary conditions for the hot early thermal sector developed in TEP-TH (Paper 27). By asserting that time itself is a dynamical field, the framework provides a testable path forward for precision cosmology.

## Acknowledgments and Disclosures

The author declares no competing interests. No external funding was received for this work. All analysis was conducted using publicly available astronomical data and open-source software.

# 8. References

## 8.1 TEP Series

- Smawfield, M. L. (2025). *Temporal Equivalence Principle: Dynamic Time & Emergent Light Speed*. v0.9 (Jakarta). DOI: 10.5281/zenodo.16921911.

- Smawfield, M. L. (2026). *The Cepheid Bias: Resolving the Hubble Tension*. v0.6 (Kingston upon Hull). DOI: 10.5281/zenodo.18209702.

- Smawfield, M. L. (2026). *Temporal Equivalence Principle: A Unified Resolution to the JWST High-Redshift Anomalies*. v0.4 (Kos). DOI: 10.5281/zenodo.19000827.

- Smawfield, M. L. (2026). *Temporal Equivalence Principle: Suppressed Density Scaling in Globular Cluster Pulsars*. v0.6 (Caracas). DOI: 10.5281/zenodo.18165798.

- Smawfield, M. L. (2026). *Temporal Equivalence Principle: Temporal Shear Recovery in Gaia DR3 Wide Binaries*. v0.3 (Kilifi). DOI: 10.5281/zenodo.19102061.

- Smawfield, M. L. (2026). *TEP-HC: Boltzmann Perturbation Closure and Acoustic-Scale Preservation*. v0.5 (Cambridge). DOI: 10.5281/zenodo.20682752.

- Smawfield, M. L. (2026). *TEP-QF: Quantum Foundations and Proper-Time Phase Holonomy*. v0.1. Zenodo.

- Smawfield, M. L. (2026). *TEP-SPIN: Topological Fermions and the Temporal Vortex*. v0.1. Zenodo.

- Smawfield, M. L. (2026). *TEP-TH: Nonsingular Temporal-Horizon Closure*. v0.2 (Thika). DOI: 10.5281/zenodo.20723059.

## 8.2 Data Sources

- Scolnic, D., et al. (2018). *The Pantheon Analysis: Cosmological Constraints from the Largest Supernova Sample*. ApJ, 859, 101.

- Scolnic, D., et al. (2022). *Pantheon+: Type Ia Supernova Light Curves from the Dark Energy Survey*. ApJ, 938, 113.

- Planck Collaboration (2020). *Planck 2018 results. VI. Cosmological parameters*. A&A, 641, A6.

- Fixsen, D. J., et al. (1996). *The Spectrum of the Cosmic Background Radiation*. ApJ, 473, 576.

- Mather, J. C., et al. (1994). *Measurement of the Cosmic Microwave Background Spectrum by the COBE FIRAS Instrument*. ApJ, 420, 439.

## 8.3 BAO and RSD Surveys

- Alam, S., et al. (2017). *The clustering of galaxies in the completed SDSS-III Baryon Oscillation Spectroscopic Survey: cosmological analysis of the DR12 galaxy sample*. MNRAS, 470, 2617.

- Beutler, F., et al. (2011). *The 6dF Galaxy Survey: baryon acoustic oscillations and the local Hubble constant*. MNRAS, 416, 3017.

- Anderson, L., et al. (2014). *The clustering of galaxies in the SDSS-III BAO sample: analysis of potential systematics*. MNRAS, 441, 24.

- Peacock, J. A., et al. (2015). *The SDSS-IV extended Baryon Oscillation Spectroscopic Survey: overview and early data*. MNRAS, 452, 2379.

- Dawson, K. S., et al. (2013). *The SDSS-III Baryon Oscillation Spectroscopic Survey: quasar targeting*. AJ, 145, 10.

- Ross, A. J., et al. (2015). *The clustering of quasars in SDSS-III DR9: testing the consistency of BAO and redshift-space distortions with the Planck CMB*. MNRAS, 449, 835.

## 8.4 Software and Tools

- Foreman-Mackey, D., et al. (2013). *emcee: The MCMC Hammer*. PASP, 125, 306. github.com/dfm/emcee

- Speagle, J. S. (2020). *dynesty: A dynamic nested sampling package for estimating Bayesian posteriors and evidences*. MNRAS, 493, 3132. github.com/joshspeagle/dynesty

- Torrado, J., & Lewis, A. (2021). *Cobaya: Code for Bayesian Analysis of cosmological data*. Astrophysics Source Code Library, ascl:2108.05. github.com/CobayaSampler/cobaya

- Lesgourgues, J. (2011). *The Cosmic Linear Anisotropy Solving System (CLASS). Part I: Overview*. arXiv:1104.2932. github.com/lesgourg/class_public

- Arbey, A. (2012). *AlterBBN: A program for calculating the BBN abundances of the elements in alternative cosmologies*. CPC, 183, 1822. alterbbn.hepforge.org

## 8.5 Historical References

- Hubble, E. (1929). *A relation between distance and radial velocity among extra-galactic nebulae*. PNAS, 15, 168.

- Friedmann, A. (1922). *Uber die Krummung des Raumes*. Z. Phys., 10, 377.

- Lemaitre, G. (1927). *Un univers homogene de masse constante et de rayon croissant rendant compte de la vitesse radiale des nebuleuses extra-galactiques*. Ann. Soc. Sci. Brux., 47, 49.

- Riess, A. G., et al. (1998). *Observational evidence from supernovae for an accelerating universe and a cosmological constant*. AJ, 116, 1009.

- Perlmutter, S., et al. (1999). *Measurements of Omega and Lambda from 42 high-redshift supernovae*. ApJ, 517, 565.

- Tolman, R. C. (1930). *On the estimation of distances in a curved universe with a non-static line element*. PNAS, 16, 511.

- Etherington, I. M. H. (1933). *On the definition of distance in general relativity*. Philos. Mag., 15, 761.

Smawfield, M. L. 2026. Temporal Equivalence Principle series, Papers 0-13. Zenodo preprints and associated repositories.

# 9. Data Availability & Reproducibility

This work follows open-science practices. All results are fully reproducible from raw data
using the documented pipeline. All numerical results, figures, and statistics are generated by deterministic
Python scripts processing real observational data. The pipeline is intentionally strict: failed dependencies are recorded as failed
results, not silently ignored.

### Repository and Code

GitHub Repository: github.com/matthewsmawfield/TEP-C0

The repository contains a deterministic, version-controlled cosmological analysis pipeline with 64 analysis steps
for supernova distance-redshift, distance-duality constraints, CMB acoustic scales, BBN preservation, structure growth, and systematic validation.
All steps are orchestrated by `scripts/run_pipeline.py` with comprehensive per-step logging.

#### Repository Structure

TEP-C0/
├── data/
│   ├── raw/                       # Downloaded source catalogs (Pantheon+, DDR, etc.)
│   └── processed/                 # Ingested and filtered datasets
├── scripts/
│   ├── steps/                     # 64 deterministic pipeline steps
│   ├── utils/                     # Logging and validation utilities
│   └── run_pipeline.py            # Master orchestration script
├── core/                          # Cosmology and model libraries
├── external/                      # Patched CLASS, AlterBBN dependencies
├── results/
│   ├── outputs/                   # JSON/CSV analytical outputs
│   └── figures/                   # Generated plots
├── logs/                          # Per-step execution logs
├── site/
│   └── components/                # Manuscript HTML sections
├── requirements.txt               # Python dependencies
└── README.md                      # Documentation

### Data Provenance

| Data Source | Provider | Access Method | Records | Location |
| --- | --- | --- | --- | --- |
| Pantheon+ SNe Ia | Scolnic et al. | Auto-downloaded | 1,701 | `data/raw/pantheon_plus_shoes.dat` |
| Pantheon+ covariance | Scolnic et al. | Auto-downloaded | Full stat + sys | `data/raw/Pantheon+SH0ES.cov` |
| BAO constraints | BOSS, eBOSS, DES | Compiled from lit. | 10 measurements | `data/raw/ddr_constraints.csv` |
| SZ cluster DDR | Compiled | Auto-downloaded | ~38 clusters | `data/raw/sz_constraints.csv` |
| SGL lensing DDR | Compiled | Auto-downloaded | ~118 lenses | `data/raw/sgl_constraints.csv` |
| DESI/eBOSS Lyman-alpha | DESI-DR1, eBOSS | Auto-downloaded | 3 measurements | `data/raw/desi_ddr.csv` |
| FIRAS CMB spectrum | NASA LAMBDA | Auto-downloaded | ~43 frequencies | `data/raw/firas_spectrum.dat` |
| Planck 2018 CMB | Planck Collaboration | Cobaya package | TTTEEE+lensing | External Cobaya cache |
| BBN abundances | AlterBBN, compiled lit. | Included / downloaded | Yp, D/H, Li/H | `data/raw/bbn_review.html` |

### Pipeline Architecture

The analysis pipeline comprises 64 deterministic steps organized into eight logical stages.
Each step is a standalone Python script in `scripts/steps/` that produces JSON/CSV outputs and
detailed logs in `logs/step_*.log`. Dependencies are resolved automatically by the runner.

#### Complete Step Inventory and Runtime

Runtimes are approximate and measured on Apple M4 Pro (14-core, 24 GB). The dominant cost is the nested sampling step (03_01), which scales with `nlive` and number of models.

| Stage | Step | Script | Description | Runtime |
| --- | --- | --- | --- | --- |
| Stage 1: Data Acquisition (8 steps) |  |  |  |  |
| Data | 1.1 | `step_01_01_data_download.py` | Download Pantheon+ SNe, covariance, FIRAS | ~10 s |
| Data | 1.2 | `step_01_02_data_ingestion.py` | Ingest and validate all downloaded catalogs | ~1 s |
| Data | 1.3 | `step_01_03_download_ddr.py` | Download BAO distance-duality constraints | ~1 s |
| Data | 1.4 | `step_01_04_download_sb.py` | Download surface-brightness catalog sources | ~1 s |
| Data | 1.5 | `step_01_05_download_sz.py` | Download Sunyaev-Zel'dovich cluster data | ~1 s |
| Data | 1.6 | `step_01_06_download_sgl.py` | Download strong gravitational lensing data | ~1 s |
| Data | 1.7 | `step_01_07_download_desi.py` | Download DESI-DR1 and eBOSS Lyman-alpha | ~1 s |
| Data | 1.8 | `step_01_08_compile_sb.py` | Compile surface-brightness master catalog | ~1 s |
| Stage 2: Theory and Transport (4 steps) |  |  |  |  |
| Theory | 2.1 | `step_02_01_transport_kernel.py` | Verify FLRW recovery limit of open-path transport K_T | ~1 s |
| Theory | 2.2 | `step_02_02_theory_derivation.py` | Derive theoretical predictions for distance-redshift and screening | ~2 s |
| Theory | 2.3 | `step_02_03_physics_implementation.py` | Implement TEP physics: distance moduli, transport, growth kernels | ~3 s |
| Theory | 2.4 | `step_02_04_screening_scale_transfer.py` | Micro-to-galactic screening scale transfer and coarse-graining | ~1 s |
| Stage 3: Model Comparison and MCMC (10 steps) |  |  |  |  |
| Core | 3.1 | `step_03_01_three_model_comparison.py` | Nested sampling (dynesty, nlive=500) for M0a_LCDM, M0b_EdS, M1 variants, M2_PureShear, M3_wCDM, M4_CPL; null injection | ~90 min |
| Core | 3.2 | `step_03_02_independent_mcmc.py` | Independent MCMC convergence diagnostics | ~1 s |
| Core | 3.4 | `step_03_04_cobaya_mcmc.py` | Joint SNe+CMB MCMC via Cobaya with TEP-CLASS v2.0 | ~2 min |
| Core | 3.4b | `step_03_04_minimize.py` | BOBYQA minimizer for joint SNe+CMB parameter optimization | ~30 s |
| Core | 3.5 | `step_03_05_analyze_cobaya.py` | Analyze Cobaya chains and produce parameter constraints | ~1 s |
| Core | 3.6 | `step_03_06_cobaya_verbose.py` | Verbose Cobaya configuration and extended diagnostics | ~2 min |
| Core | 3.7 | `step_03_07_likelihood_synthesis.py` | Synthesize likelihoods across independent and joint analyses | ~1 s |
| Core | 3.8 | `step_03_08_h0_boundary_stress.py` | H0 prior stress test: extended priors reveal EdS-derived-parameter artifact driving H0 toward zero | ~30 s |
| Core | 3.9 | `step_03_09_lcdm_null_injection.py` | LCDM null injection: mock Pantheon+ from LCDM, measure TEP false-positive rate | ~60 s |
| Core | 3.10 | `step_03_10_pantheon_subset_robustness.py` | Leave-one-survey-out and redshift-window robustness tests | ~30 s |
| Stage 4: Supernova Tests and Distance Duality (10 steps) |  |  |  |  |
| SNe | 4.1 | `step_04_01_sn_time_dilation.py` | Test SN light-curve stretch factors against TEP time dilation | ~1 s |
| SNe | 4.2 | `step_04_02_sn_tolman.py` | Tolman surface-brightness dimming test | ~1 s |
| SNe | 4.3 | `step_04_03_tolman_sb.py` | Surface-brightness Tolman scaling with compiled catalog | ~1 s |
| DDR | 4.4 | `step_04_04_distance_duality.py` | Distance-duality relation: BAO constraints vs TEP prediction | ~1 s |
| DDR | 4.5 | `step_04_05_ddr_threeway.py` | Three-way probe comparison: BAO, SZ, SGL | ~1 s |
| DDR | 4.6 | `step_04_06_screening_fit.py` | Parametric screening model fit to probe-dependent DDR | ~2 s |
| DDR | 4.7 | `step_04_07_highz_ddr.py` | High-redshift Lyman-alpha DDR test (DESI, eBOSS) | ~1 s |
| SNe | 4.8 | `step_04_08_host_mass_step_prediction.py` | Host-mass-step mini-analysis: locked TEP prediction vs fitted LCDM nuisance | ~5 s |
| PPN | 4.9 | `step_04_09_ppn_constraints.py` | Solar System PPN constraint derivation with gradient-dependent screening | ~1 s |
| DDR | 4.10 | `step_04_10_tep_native_ddr.py` | TEP-native distance-duality re-analysis | ~1 s |
| Stage 5: CMB and Big Bang Nucleosynthesis (10 steps) |  |  |  |  |
| CMB | 5.1 | `step_05_01_cmb_blackbody.py` | Verify TEP preserves CMB blackbody spectrum (FIRAS) | ~1 s |
| CMB | 5.3 | `step_05_03_cmb_boltzmann.py` | TEP Boltzmann integration via patched CLASS | ~1 s |
| CMB | 5.4 | `step_05_04_cmb_spectra.py` | Generate and compare TT/TE/EE power spectra | ~1 s |
| CMB | 5.5 | `step_05_05_cmb_consistency.py` | CMB acoustic-scale consistency check | ~1 s |
| BBN | 5.6 | `step_05_06_bbn_registry.py` | Compile observational BBN abundance registry | ~1 s |
| BBN | 5.7 | `step_05_07_bbn_preservation.py` | Cross-validate TEP and LCDM BBN predictions | ~1 s |
| CMB | 5.8 | `step_05_08_cmb_acoustic.py` | Acoustic-scale parameter comparison (Planck) | ~1 s |
| CMB | 5.9 | `step_05_09_minimal_perturbations.py` | Diagnostic minimal-closure perturbation checks; authoritative active-sector closure is documented in TEP-HC | ~3 s |
| CMB | 5.10a | `step_05_10_jordan_frame_proof.py` | Matter-frame acoustic-scale proof in EdS matter-only background | ~1 s |
| CMB | 5.10b | `step_05_10_tephc_spectra_crosscheck.py` | Cross-check TEP-HC Boltzmann spectral outputs (acoustic-scale ratio) | ~1 s |
| Stage 6: BAO and Structure Growth (7 steps) |  |  |  |  |
| BAO | 6.1 | `step_06_01_bao_projection.py` | BAO ruler projection in TEP geometry | ~1 s |
| BAO | 6.2 | `step_06_02_bao_likelihood.py` | BAO likelihood module integration | ~7 s |
| Growth | 6.3 | `step_06_03_growth_solver.py` | TEP-CLASS v2.0 growth equation solver | ~1 s |
| Growth | 6.4 | `step_06_04_growth_validation.py` | Validate growth factors against LCDM baseline | ~1 s |
| Growth | 6.5 | `step_06_05_growth_rsd.py` | Redshift-space distortion comparison (f sigma_8) | ~2 s |
| Growth | 6.6 | `step_06_06_nonlinear_growth_closure.py` | Halo-model non-linear growth closure with gradient screening | ~5 s |
| Growth | 6.7 | `step_06_07_alphaM_growth_validation.py` | First-principles alpha_M-modified growth ODE: four-scenario sigma_8 comparison | ~2 s |
| Stage 7: Forecasts and Future Tests (7 steps) |  |  |  |  |
| Future | 7.1 | `step_07_01_mixed_forecast.py` | Forecast for mixed TEP-LCDM parameter recovery | ~1 s |
| Future | 7.2 | `step_07_02_redshift_drift.py` | Redshift-drift forecast and discriminating power | ~1 s |
| Future | 7.3 | `step_07_03_jwst_test.py` | JWST high-z supernova feasibility test | ~1 s |
| Future | 7.4 | `step_07_04_gw_sirens.py` | Gravitational-wave standard siren forecast | ~1 s |
| Future | 7.5 | `step_07_05_weak_lensing_plan.py` | Weak-lensing survey plan for TEP discrimination | ~1 s |
| Future | 7.6 | `step_07_06_weak_lensing.py` | Weak-lensing shear correlation analysis | ~1 s |
| Future | 7.7 | `step_07_07_blind_injection.py` | Blind injection validation protocol | ~1 s |
| Stage 8: Falsification, Verification, and Summary (8 steps) |  |  |  |  |
| Validation | 8.1 | `step_08_01_expansion_falsifier.py` | Expansion falsifier: distance duality and Tolman residuals | ~1 s |
| Validation | 8.2 | `step_08_02_comparison_stats.py` | Cross-model comparison statistics | ~1 s |
| Validation | 8.3 | `step_08_03_sensitivity_analysis.py` | Prior and parameter sensitivity analysis | ~1 s |
| Validation | 8.4 | `step_08_04_evidence_matrix.py` | Compile explanatory evidence matrix | ~1 s |
| Validation | 8.5 | `step_08_05_gate_registry.py` | Claim gate registry and status check | ~1 s |
| Validation | 8.6 | `step_08_06_claim_audit.py` | Automated claim consistency check | ~1 s |
| Validation | 8.7 | `step_08_07_final_summary.py` | Global evidence synthesis and summary | ~1 s |
| Validation | 8.8 | `step_08_08_diagnostic_plots.py` | Data-driven diagnostic figures (distance-duality residuals, Pantheon+ Hubble residuals) generated only from upstream pipeline artefacts | ~5 s |

#### Total Runtime Summary

The total runtime is dominated by Stage 3.1 (nested sampling). Runtimes scale approximately linearly with `nlive` and number of CPU cores.

| Component | Steps | Runtime |
| --- | --- | --- |
| Data Acquisition (Stage 1) | 8 | ~20 s |
| Theory and Transport (Stage 2) | 4 | ~6 s |
| Model Comparison and MCMC (Stage 3) | 10 | ~97 min |
| SNe Tests and DDR (Stage 4) | 10 | ~15 s |
| CMB and BBN (Stage 5) | 10 | ~11 s |
| BAO and Growth (Stage 6) | 7 | ~14 s |
| Forecasts and Future Tests (Stage 7) | 7 | ~7 s |
| Falsification and Verification (Stage 8) | 8 | ~7 s |
| Total | 64 | ~95 min (~1.6 h) |

### Reproduction Instructions

#### Quick Start (Full Reproduction)

# 1. Clone repository
git clone https://github.com/matthewsmawfield/TEP-C0.git
cd TEP-C0

# 2. Install dependencies
pip install -r requirements.txt

# 3. Run full pipeline (generates all results and figures)
python scripts/run_pipeline.py

# 4. Results will be in:
#    - results/outputs/   (JSON/CSV data)
#    - results/figures/   (PNG/PDF plots)
#    - logs/              (Detailed execution logs)

#### Command-Line Options

The pipeline supports selective execution for faster testing:

# Core statistical analysis only (skips long nested sampling)
python scripts/run_pipeline.py --core

# Resume from existing results (skip completed steps)
python scripts/run_pipeline.py --resume

# Run specific steps with automatic dependency resolution
python scripts/run_pipeline.py --steps step_04_04_distance_duality step_04_05_ddr_threeway

#### System Requirements

| Component | Minimum | Recommended | Tested On |
| --- | --- | --- | --- |
| CPU | 4 cores | 8+ cores | Apple M4 Pro (14-core) |
| RAM | 8 GB | 16 GB | 24 GB (M4 Pro) |
| Storage | 2 GB | 5 GB | NVMe SSD |
| Runtime (full) | ~4 h (4 cores) | ~1.5 h (8+ cores) | ~95 min (M4 Pro) |
| Runtime (--core) | ~1 min | ~30 s | ~20 s |

#### Key Analysis Outputs

- `results/outputs/step_03_01_three_model_comparison.json` — Nested sampling posteriors and evidence for all models (M0a_LCDM, M0b_EdS, M1 variants, M2_PureShear, M3_wCDM, M4_CPL)

- `results/outputs/step_03_04_cobaya_mcmc.1.txt` — Cobaya MCMC chain for joint SNe+CMB analysis

- `results/outputs/step_04_04_distance_duality.json` — DDR weighted mean and deviation from unity

- `results/outputs/step_04_05_ddr_threeway.json` — Three-way BAO/SZ/SGL probe comparison

- `results/outputs/step_05_07_bbn_preservation.json` — TEP vs LCDM light-element abundance cross-validation

- `results/outputs/step_05_09_minimal_perturbations.json` — diagnostic minimal-closure perturbation checks; authoritative active-sector closure is documented in TEP-HC

- `results/figures/step_05_09_perturbation_spectra.png` — diagnostic cross-check of minimal conformal perturbations; full active perturbation safety claim is carried by TEP-HC

- `results/outputs/step_06_04_growth_validation.json` — Growth factor and sigma_8 consistency check

- `results/outputs/step_08_04_evidence_matrix.json` — Explanatory evidence matrix across all observables

- `results/outputs/step_08_06_claim_audit.json` — Automated claim consistency check report

#### Log Files

Each step produces detailed logs with timestamps, SHA-256 checksums, and execution status:

- `logs/step_*.log` — Individual step logs (64 files, one per step)

- `logs/verbose/` — Verbose Cobaya and nested sampling logs

### Software Dependencies

| Package | Version | Purpose |
| --- | --- | --- |
| Python | 3.10+ | Language runtime |
| NumPy | 1.24+ | Numerical computing |
| SciPy | 1.10+ | Statistical functions, nested sampling |
| Pandas | 2.0+ | Data manipulation |
| Matplotlib | 3.7+ | Visualization |
| emcee | 3.1+ | Ensemble MCMC sampling |
| dynesty | 2.1+ | Nested sampling for Bayesian evidence |
| Cobaya | 3.6+ | Joint MCMC with Planck likelihoods |
| classy (CLASS) | 3.2+ | CMB Boltzmann solver (patched for TEP) |

All dependencies are specified in `requirements.txt`. External dependencies (patched CLASS, AlterBBN) are included in the `external/` directory.

### Appendix Figures

![Joint SNe+CMB Background/Acoustic MCMC Diagnostic](results/figures/step_03_05_analyze_cobaya_triangle.png)

Figure A1: Joint SNe+CMB Background/Acoustic MCMC Diagnostic. This triangle plot shows the joint posterior from the Cobaya MCMC, including the homogeneous acoustic-sector amplitude $\epsilon_T^{\rm CMB}$. This is a diagnostic figure, not the SNe-only M1 evidence result. The $H_0$ boundary behaviour is separately stress-tested (Section 4.4.3). The $\epsilon_T$ shown here is the homogeneous acoustic-sector amplitude, distinct from the line-of-sight $\epsilon_{\rm shear}^{\rm los}$ fitted to supernovae.

![Minimal Conformal Perturbations vs LCDM](results/figures/step_05_09_minimal_perturbations_perturbation_spectra.png)

Figure A2: TEP Minimal Conformal Perturbations vs. $\Lambda$CDM. **Top panel:** TT power spectrum $D_\ell^{TT}$ for $\Lambda$CDM and TEP minimal conformal perturbations. **Bottom panel:** fractional residuals with quantitative gate outputs (max residual, acoustic peak shift, proxy $\chi^2$). This figure is a diagnostic cross-check of minimal conformal perturbations. TEP-HC (Paper 18) carries the full active perturbation safety claim; C0 cross-checks and imports these outputs.
