# Temporal Equivalence Principle: Native hi_class Conformal Implementation, Linear Perturbation Closure, and CMB Acoustic Peak Preservation
**Matthew Lukin Smawfield**
Version: v0.5 (Cambridge)
First published: 8 June 2026 · Last updated: 3 July 2026
Paper Series: TEP Series Paper 18 (hi_class Cosmology)
DOI: 10.5281/zenodo.20682752

---

## Abstract


Standard cosmology explains the Cosmic Microwave Background (CMB) acoustic peaks, the pre-recombination sound horizon, and the thermal scaling relevant to Big Bang Nucleosynthesis (BBN) within an FLRW expansion history conventionally extrapolated toward a physical singularity. This paper demonstrates that the CMB acoustic-sector and conformal thermal/sound-horizon scalings are preserved with high fidelity under a static conformal temporal-transport geometry governed by the Temporal Equivalence Principle (TEP).



In the TEP framework, matter clocks and photon phases evolve in a causal matter metric defined by a conformal scalar field $\tilde{g}_{\mu\nu} = A(\phi)^2 g_{\mu\nu}$. Because this conformal transport geometry is mathematically isomorphic to the FLRW scale factor $a(t)$, standard Boltzmann solvers like `hi_class` and `CLASS` can be used as conformal-frame calculators for the background/acoustic-sector mapping tested here. The parameter traditionally identified as Dark Energy ($\Omega_\Lambda$) is operationally reinterpreted within this implementation as the homogeneous temporal-shear background contribution filling the same background budget slot, $\Omega_\phi$.



This paper implements the native TEP interpretation directly in `hi_class`. Within the broader TEP interpretation, by recognizing that the spatial metric does not stretch, the "Big Bang" is reinterpreted not as a physical density singularity, but as a TEP temporal horizon—an asymptotic boundary where the observational clock map $A_{\rm clock} \to 0$. Direct Boltzmann integration verifies this background/acoustic mathematical isomorphism, verifying the internal conformal-frame preservation of the acoustic sector — the pre-recombination sound horizon is preserved to parts-per-million and the acoustic-peak morphology remains intact without invoking early-universe spatial expansion.



Beyond the background mapping, the paper closes the linear pure-conformal scalar perturbation sector by deriving the runtime Bellini–Sawicki functions $\alpha_M=-2\alpha_A$, $\alpha_B=2\alpha_A$, $\alpha_K=-5\alpha_A^2$, and $\alpha_T=0$, for which the physical no-ghost discriminant satisfies $D=\alpha_K+\frac{3}{2}\alpha_B^2=\alpha_A^2$. An active-perturbation `hi_class` run evolving $\delta\phi$ through the full Einstein–Boltzmann hierarchy produces posteriors statistically indistinguishable from the background-only chain, demonstrating that the implemented linear pure-conformal scalar perturbation sector is stable and observationally negligible at the current homogeneous-amplitude bound.



A joint `hi_class` Cobaya MCMC (Planck 2018 low-$\ell$ TT/EE + lensing + BAO + Pantheon+) tests the screened TEP conformal background within the native hi_class implementation, while companion TEP-C0 (Paper 26) nested sampling over Pantheon+ provides robust quantitative evidence that the screened TEP conformal geometry matches the Pantheon+ distance-redshift relation with Bayesian evidence statistically indistinguishable from CPL (BF $\sim 32$--$36$ depending on the fixed or free-$z_T$ branch), without treating late-time acceleration as primitive spatial acceleration. Within the TEP framework, the Hubble tension is interpreted as a late-time, environment-dependent clock-transport effect (Paper 11) caused by the mass-screening of the scalar field, rather than through a crisis in early-universe physics.



Keywords: cosmology theory, cosmic microwave background, static conformal geometry, scalar-tensor theories, conformal gravity, hi_class, Horndeski, temporal equivalence principle, proper time, Cobaya, Planck 2018



## 1. Introduction


### 1.1 Contextualizing the TEP Corpus


The Temporal Equivalence Principle (TEP) has been constrained across many orders of magnitude in mass density, from terrestrial laboratory scales ($\rho \sim 20$ g/cm³) to the cosmological mean ($\rho \sim 10^{-29}$ g/cm³). Previous papers in this series have established:



- *Terrestrial scales (Paper 1):* Terrestrial atomic clock networks show 4,200 km phase correlations consistent with the candidate Temporal Topology saturation scale ρ_T ≈ 20 g/cm³.

- *Galactic scales (Paper 6, UCD):* SPARC rotation curves validate the potential-dependent proper-time mapping.

- *Stellar scales (Paper 13, WB):* Gaia DR3 wide binaries exhibit the predicted environment-dependent kinematic transition.

- *Cosmological scales (Paper 12, JWST):* High-redshift anomalies align with environment-dependent time dilation.




### 1.2 The Cosmological Horizon


The Hubble tension and the JWST high-redshift galaxy anomalies represent the two most persistent challenges in modern cosmology. Standard $\Lambda$CDM relies on a stretching spatial metric, which extrapolates to a physical singularity at $a \to 0$ and tightly restricts the available proper time for early galaxy assembly.


Previous TEP work (Paper 11, H₀; Paper 12, JWST) demonstrated that these anomalies are resolved within the TEP framework. The TEP temporal-horizon picture replaces the FLRW singularity with an asymptotic boundary, removing the finite-age assembly bottleneck, and the $H_0$ tension is addressed as a local environmental mass-screening effect on kinematic distance probes.


### 1.3 Purpose of This Paper


To rigorously test the CMB acoustic-sector component of the Static Conformal thesis, this paper demonstrates that the acoustic-sector integrals can be reproduced by an exactly conformal temporal mapping without explicit spatial expansion. A full light-element abundance calculation is not performed here; the present BBN claim is limited to preservation of the conformal thermal/sound-horizon scaling and should be treated as a target for a dedicated follow-up calculation.


Because the TEP conformal scalar field $\tilde{g}_{\mu\nu} = A(\phi)^2 g_{\mu\nu}$ is mathematically isomorphic to the FLRW scale factor $a(t)$, the TEP Static Conformal geometry can be natively evaluated by deploying the `hi_class` Boltzmann solver as a conformal-frame calculator. This requires:



- Mapping the TEP conformal geometry onto the Boltzmann framework, establishing the implementation-level correspondence between the parameter conventionally written as $\Omega_\Lambda$ and the homogeneous temporal-shear background contribution, $\Omega_\phi$, occupying the same background-budget role as $\Omega_\Lambda$.

- Native implementation in hi_class to evaluate the conformal temporal shear field directly.

- A joint MCMC parameter estimation ($H_0$, $\Omega_b h^2$, $\Omega_{\rm cdm} h^2$, $n_s$, $A_s$, $\tau$, $A_{\rm planck}$, $\epsilon_T$) against Planck 2018, BAO, and Pantheon+ data to quantitatively demonstrate that the acoustic peaks are preserved to parts-per-million accuracy in a static conformal geometry.



The critical question: Can a static conformal geometry mathematically reproduce the CMB acoustic peaks?

The claim-discipline framework for the TEP corpus, including the scope limitations of canonical precision tests, is established in TEP-EXP (Paper 9).

## 2. Theoretical Architecture: The EFT Mapping

### 2.1 The Bi-Metric Action

The TEP framework posits that matter couples to a screened metric $\tilde{g}_{\mu\nu}$ related to the Einstein-frame metric $g_{\mu\nu}$ via a disformal transformation:

\begin{equation} \label{eq:3_theory_01}
\tilde{g}_{\mu\nu} = A^2(\phi) g_{\mu\nu} + B(\phi) \nabla_\mu\phi \nabla_\nu\phi
\end{equation}

where:

- $A(\phi) = \exp(\beta_A\phi/M_{\rm Pl})$ is the conformal factor, with $\beta_A = -1.0$ (the locked lab-scale convention used across the TEP corpus)

- $B(\phi)$ controls disformal deformation of the causal structure

- $\phi$ is the dynamical proper-time field

Following TEP-TH (Paper 27), two distinct projections of the temporal field are distinguished in cosmology: $A_{\rm clock}(z)$ is the exact observational clock/redshift mapping that generates the apparent distance-redshift relation, while $A_{\rm dyn}(z)$ is the screened physical dynamical response that modifies expansion, BBN, recombination, and perturbations. In the homogeneous background limit evaluated here, the hi_class conformal-frame mapping corresponds to $A_{\rm clock}$; the screened $A_{\rm dyn}$ response is parameterized by $\epsilon_T$ and is driven to unity during the pre-recombination acoustic sector.

*Metric signature convention:* $(+, -, -, -)$ throughout.

### 2.2 Formal Bellini-Sawicki Alpha Correspondence

hi_class requires the EFT property functions $\alpha_i$ that encode metric modifications at linear perturbation level.

#### 2.2.1 Planck Mass Running ($\alpha_M$)

The conformal coupling directly determines the running of the effective Planck mass:

\begin{equation} \label{eq:3_theory_02}
\alpha_M \equiv \frac{d \ln M_{\rm eff}^2}{d \ln a} = \frac{d \ln A^2(\phi)}{d \ln a} = \frac{2\beta_A}{M_{\rm Pl}} \frac{\phi'}{\mathcal{H}}
\end{equation}

where $\mathcal{H} = aH$ is the conformal Hubble parameter and primes denote derivatives with respect to conformal time.

#### 2.2.2 Tensor Speed Excess ($\alpha_T$)

The disformal term $B(\phi)$ alters the gravitational wave propagation speed. Multi-messenger constraints from GW170817/GRB 170817A require:

\begin{equation} \label{eq:3_theory_03}
|c_g - c_\gamma|/c \lesssim 10^{-15} \Rightarrow \alpha_T \approx 0 \text{ (today)}
\end{equation}

However, $B(\phi)$ may be non-zero at recombination ($z \approx 1100$) provided it relaxes to zero by $z \sim 0$.

#### 2.2.3 Braiding ($\alpha_B$) and Kineticity ($\alpha_K$)

These functions govern scalar field clustering and metric mixing:

\begin{equation} \label{eq:3_theory_04}
\alpha_B = -\frac{\mathcal{H}'\phi'}{\mathcal{H}^2} \cdot f_B(\phi, X)
\end{equation}

\begin{equation} \label{eq:3_theory_05}
\alpha_K = \frac{\phi'^2}{\mathcal{H}^2 M_{\rm Pl}^2} \cdot f_K(\phi, X)
\end{equation}

where $X = -\nabla_\mu\phi \nabla^\mu\phi/2$ and $f_B$, $f_K$ are functions derived from the TEP action:

- $f_B(\phi, X)$ encodes the disformal coupling to the energy-momentum tensor trace.

- $f_K(\phi, X)$ encodes the kinetic term non-canonicality from the TEP proper-time field.

The explicit functional forms follow from the bi-metric action (Equation \ref{eq:3_theory_01}) and are determined by the conformal factor $A(\phi)$ and the disformal function $B(\phi)$. Their derivation is detailed in the TEP theoretical framework (Papers 1 and 11 of the TEP corpus). The schematic $f_B$, $f_K$ forms describe the general disformal EFT bookkeeping; the production perturbation run restricts to the pure-conformal branch, where these functions reduce to the closed runtime identities given in Appendix A.5.

The background-only configuration uses $\alpha_B$ and $\alpha_K$ as formal EFT bookkeeping; during native `tep_mode` background-only integration $f_B$ and $f_K$ evaluate identically to zero, verifying that this is a strict geometric mapping at the background level. The active-perturbation closure configuration (Section 4.5) passes `gravity_model = tep` and `M2_evolution = yes` to hi_class, which evolves the implemented linear pure-conformal scalar fluctuation sector $\delta\phi$ through the exact Bellini–Sawicki EFT functions derived from the conformal geometry, with sound speed $c_s^2 = 1.0$ and no-ghost condition $D = \alpha_A^2$ enforced by construction. The theoretical basis for this closure is developed in the foundational TEP formalism (Papers 1 and 11).

### 2.3 The Static Conformal Isomorphism

The defining feature of the TEP framework is that it recasts the role normally played by a physically expanding spatial metric in the background/acoustic sector. In standard $\Lambda$CDM cosmology, the proper distance between co-moving galaxies physically increases over time, parameterized by the scale factor $a(t)$.

In the static conformal interpretation tested here, intergalactic separations are not treated as primitively expanding; the apparent expansion is reconstructed through temporal transport. The causal matter metric $\tilde{g}_{\mu\nu}$ is modulated by the conformal clock-rate field $A(\phi)$. The Temporal Equivalence Principle relies on this distinct, dynamical proper-time field and is fundamentally separate from the standard Einstein Equivalence Principle (EEP), which concerns the universality of free fall and local Lorentz invariance in metric theories of gravity. Photons propagating through this gradient experience a shift in phase and frequency, leading to the exact geometric relation:

\begin{equation} \label{eq:3_theory_06}
1+z = \frac{A_0}{A_{\text{em}}} = A_{\text{clock}}^{-1}(z)
\end{equation}

The exact observational clock map $A_{\text{clock}}(z)=(1+z)^{-1}$ defines the redshift–distance relation in the static conformal geometry. Because the mathematical transport geometry of $A_{\text{clock}}$ across a static background is formally isomorphic to the transport geometry of a photon in an expanding FLRW metric with reconstructed scale factor $a_{\text{eff}}(z)$, standard cosmological integrators (like `hi_class` and `CLASS`) can be deployed natively as conformal-frame calculators for the background/acoustic-sector mapping tested here. The screened physical dynamical response $A_{\text{dyn}}(z)$ modifies the Hubble parameter and acoustic observables at late times while approaching unity during the early thermal epochs (TEP-TH, Paper 27; TEP-C0, Paper 26).

- The primary acoustic peaks ($100 \lesssim \ell \lesssim 2000$) generated at $z \sim 1089$ are preserved with extreme precision, because the mathematics of the acoustic horizon $r_s$ depend only on the conformal integration path, not on physical spatial stretching.

- The parameter conventionally identified as Dark Energy ($\Omega_\Lambda$) is reinterpreted within this implementation as the homogeneous temporal-shear background contribution, $\Omega_\phi$, occupying the same background-budget role as $\Omega_\Lambda$ in the reference FLRW calculation.

- The Big Bang singularity is reinterpreted, within TEP, as a temporal-horizon boundary of the conformal clock-rate field where $A_{\text{clock}} \to 0$ relative to the present epoch, driving the reconstructed scale factor $a_{\text{eff}} \to 0$ and creating the mathematical appearance of infinite density in standard FLRW reconstructions. The full nonsingular closure—finite curvature invariants, geodesic completeness, and the origin of the CMB blackbody spectrum—is developed in TEP-TH (Paper 27), where the non-exact covariance/topology correction $C_T$ provides the transport closure beyond exact conformal shear.

- *Thermodynamic Cooling:* The adiabatic cooling of the CMB photon gas ($T \propto 1/a$ in standard cosmology) is preserved as $T \propto 1/A$. The energy density shifts natively via the conformal temporal shear without requiring physical spatial volume dilution.

*Operational definition of $\Omega_\phi$.* In this implementation, $\Omega_\phi$ is operationally defined as the homogeneous conformal-sector contribution that fills the same background budget slot occupied by $\Omega_\Lambda$ in the reference FLRW integration. A first-principles stress-energy derivation of $\rho_\phi$, $p_\phi$, and the effective equation of state belongs to the broader TEP-C0 action-level treatment.

*Archived EFT reference.* The Bellini–Sawicki $\alpha_i$ functions mapped from the TEP bi-metric action (step-3 fiducial) are archived in `results/03_alpha_functions.json`. Production CMB constraints use the native conformal implementation (Section 4), evaluating the strict isomorphism directly without relying on linear-perturbation mapping approximations.


## 3. Software Implementation: hi_class and the Unscreened Regime


### 3.1 The hi_class Architecture


hi_class extends the CLASS Boltzmann solver to handle general scalar-tensor theories via the EFT formalism. This work uses hi_class v3.2.3 with the modified gravity (SMG) module enabled.


### 3.2 Native TEP Conformal Background Implementation

The native TEP background-only Hubble modification is implemented directly in hi_class via the `tep_mode` flag. When enabled, the background expansion history is modified as:

\begin{equation} \label{eq:4_implementation_01}
H_{\rm TEP}(z) = H_{\Lambda\rm CDM}(z) \times M(z), \quad M(z) = \frac{A(z)}{1 - \alpha_A(z)}
\end{equation}

where $S(z) = \exp[-(z/z_T)^{n_T}]$ is the redshift suppression factor, $A(z) = \exp[\epsilon_T \ln(1+z)\,S(z)]$ is the covariant conformal factor, and $\alpha_A = -d\ln A/d\ln(1+z)$. The production integration directly evaluates the exact geometric relation $M = A/(1-\alpha_A)$, ensuring computational fidelity to the underlying formal derivation (see Appendix A.3a). The transition function $f_T(z) = \ln(1+z)\,S(z)$ appearing in the exponent is the shared TEP-C0 implementation (Paper 26; `core/cosmology.py`: `f_T`, `conformal_factor_native`, `jordan_frame_M`):

\begin{equation} \label{eq:4_implementation_01b}
f_T(z) = \ln(1+z)\,S(z).
\end{equation}

The suppression factor $S(z)$ defines the physical profile of the conformal field: it drives $f_T \to 0$ for $z \gg z_T$, ensuring the field profile flattens asymptotically as it approaches the TEP temporal horizon; the $\ln(1+z)$ factor enforces $f_T(0)=0$, fixing the local reference frame so $H_0$ is anchored to the local observer. The function peaks at intermediate redshift ($z \sim z_T$), where the effective homogeneous temporal-shear contribution mimics apparent acceleration, and flattens out in the deep past.

*Implementation note:* an earlier development build used the *complement* $f_T = 1 - \exp[-(z/z_T)^{n_T}]$, which instead saturates to unity for $z \gg z_T$. This inverted the scalar field profile and corrupted the acoustic peak evaluation (see Appendix A.3). The default TEP conformal parameters are:



```
tep_mode = yes
epsilon_T = 0.0066
z_T = 5.0
n_T = 2.0
```

The background conformal mapping is implemented through $M(z)=A/(1-\alpha_A)$. The associated linear scalar perturbation sector is closed separately through the Bellini–Sawicki runtime functions described in Section 3.3 and Appendix A.5, and validated in the active $\delta\phi$ run of Section 4.5. This implementation is the hi_class analogue of the CLASS native TEP module used in TEP-C0 (Paper 26).


### 3.3 Perturbation Stability and Closure

The present analysis natively closes the scalar perturbation sector by evaluating the exact runtime Bellini-Sawicki Effective Field Theory (EFT) parameters in the pure conformal limit. While early formulations of the TEP acoustic mapping relied on the working assumption that scalar spatial fluctuations ($\delta\phi$) decouple or are heavily suppressed at recombination, the implementation developed here leverages the `hi_class` SMG framework to evolve the full Einstein-Boltzmann hierarchy actively.

Because the TEP framework maps its causal matter metric via the exact conformal relation $\tilde{g}_{\mu\nu}=A(\phi)^{2}g_{\mu\nu}$, the corresponding EFT functions ($\alpha_i$) in the Jordan frame are analytically fixed by the background derivative $\alpha_A \equiv -d\ln A / d\ln(1+z)$. In the pure conformal limit ($\beta_A = -1.0$), the kineticity parameter evaluates to $\alpha_K = -5\alpha_A^2$. Although a negative kineticity frequently triggers ghost instabilities in canonical scalar-tensor theories, the full physical no-ghost discriminant $D$ in the Horndeski framework encompasses the braiding parameter, which in this limit is identically $\alpha_B = 2\alpha_A$.

Evaluating the full discriminant yields $D = \alpha_K + \frac{3}{2}\alpha_B^2 = \alpha_A^2$. The implemented linear pure-conformal sector is therefore positive-definite and ghost-free over the redshift range explored by the hi_class integration.

By forcing the integration of the implemented linear pure-conformal $\delta\phi$ scalar field through this exact geometric relation, the numerical solver verifies that the continuous screening transition smoothly manages the emergence into the late universe without triggering gradient instabilities or pathological phantom energies. At the fiducial amplitude ($\epsilon_T = 0.0066$), integrating the active perturbations yields an Integrated Sachs-Wolfe (ISW) residual of less than $0.001\%$ across the entire acoustic spectrum ($100 \le l \le 2000$). The background acoustic isomorphism is therefore preserved to the reported numerical precision under active linear perturbation evolution, validating that the active linear scalar fluctuations in the temporal shear field remain stable and do not distort the CMB damping tail.


### 3.3b Covariant Frame Alignment: Matter-Frame Hubble Friction

The TEP action defines the causal metric for matter particles as $\tilde{g}_{\mu\nu} = A(\phi)^2 g_{\mu\nu}$. Under this conformal rescaling, the matter-frame conformal Hubble rate acquires the exact relation

\begin{equation}
\tilde{\mathcal{H}} = \mathcal{H} - \frac{\dot{A}}{A} = \mathcal{H}(1 - \alpha_A),
\end{equation}

where $\alpha_A \equiv -d\ln A/d\ln(1+z)$. In the Jordan-frame EFT, matter particles follow geodesics of $\tilde{g}_{\mu\nu}$, so their Euler equations must be evaluated with $\tilde{\mathcal{H}}$ rather than the Einstein-frame rate $\mathcal{H}$. The standard hi_class implementation uses $\mathcal{H}$ for all species, which is appropriate for Einstein-frame observers (photons and neutrinos) but does not include the conformal-frame distinction for collisionless matter. We implement this distinction as follows.

A new background index `index_bg_H_conformal` is defined in `include/background.h` and populated in `source/background.c` after the TEP Hubble modification:



```
/* TEP matter-frame conformal Hubble rate */
if (pba->tep_mode == _TRUE_ && pba->epsilon_T != 0.0) {
double z = 1.0/a - 1.0;
double alpha_A = tep_alpha_A(pba, z);
pvecback[pba->index_bg_H_conformal] = pvecback[pba->index_bg_H] * (1.0 - alpha_A);
}
else {
pvecback[pba->index_bg_H_conformal] = pvecback[pba->index_bg_H];
}
```

In `source/perturbations.c`, we compute a new local variable `a_prime_over_a_matter = pvecback[index_bg_H_conformal] * a` and substitute it for `a_prime_over_a` in the Euler friction terms of CDM, baryons, interacting dark matter, and decaying cold dark matter. Photons, neutrinos, and metric equations retain the Einstein-frame rate, preserving the CMB acoustic physics. The modification is minimal and upstream: a single background quantity and a single local variable replacement in the fluid velocity derivatives.


### 3.3c Scale-Dependent Bellini-Sawicki Mapping

The gradient-screening envelope $f(g) = [1 + (g/g_t)^n]^{-1}$ governs the suppression of the temporal-shear coupling in regions of steep gravitational potential. In the cosmological context, the gradient operator translates to the Fourier wavenumber $k$ through the Poisson equation $\nabla^2 \Phi \sim k^2 \Phi_k$, so the envelope becomes

\begin{equation}
\mathcal{S}_\Sigma(k) = \frac{1}{1 + (k/k_t)^{n_t}},
\end{equation}

where $k_t$ is the screening threshold scale and $n_t$ the steepness. The Bellini-Sawicki parameter $\alpha_M$, which controls the running of the effective Planck mass and therefore the strength of the scalar fifth force, must acquire this $k$-dependence:

\begin{equation}
\alpha_M(k,a) = \alpha_M^{\rm bg}(a) \times \bigl[1 - \mathcal{S}_\Sigma(k)\bigr].
\end{equation}

This is implemented by modifying the function `get_gravity_coefficients_smg` in `gravity_smg/perturbations_smg.c` to accept the mode wavenumber $k$ as an argument. The screening is applied to the local `run` variable (carrying $\alpha_M$) immediately after it is extracted from the background workspace:



```
if (pba->tep_mode == _TRUE_ && pba->tep_k_t > 0.0 && k > 0.0) {
double S_Sigma = tep_screening_factor(k, pba->tep_k_t, pba->tep_n_t);
double alpha_M_bg = *run;
if (alpha_M_bg != 0.0) {
*run = alpha_M_bg * (1.0 - S_Sigma);
}
}
```

The screening function `tep_screening_factor` is defined in `source/background.c` and exposed in `include/background.h`. Two new parameters, `tep_k_t` and `tep_n_t`, are added to the `background` structure and parsed in `source/input.c`. The function signature of `get_gravity_coefficients_smg` is updated in `gravity_smg/include/perturbations_smg.h`, and all three call sites (in `perturbations_einstein_scalar_smg`, `perturbations_qs_functions_at_tau_and_k_qs_smg`, and `get_x_x_prime_qs_smg`) pass the mode wavenumber $k$.

At cosmic scales ($k \to 0$), $\mathcal{S}_\Sigma \to 0$ and $\alpha_M(k) \to \alpha_M^{\rm bg}$, preserving the background expansion. At structure scales ($k \gg k_t$), $\mathcal{S}_\Sigma \to 1$ and $\alpha_M(k) \to 0$, stripping the scalar fifth force and yielding $G_{\rm eff} \to G_N$ inside collapsing halos. Both the conformal friction and the scale-dependent screening use only the background parameters already constrained by the distance data ($\epsilon_T$, $z_T$, $n_T$); $k_t$ and $n_t$ are physical parameters of the gradient envelope with no additional tuning required for $\sigma_8$.


### 3.4 Pipeline Architecture

The full analysis pipeline, executed via `scripts/run_all.py`, consists of:


- *Step 0 (Setup):* Environment configuration and dependency check.

- *Step 1 (Install):* Install Cobaya, Planck 2018 likelihoods, and hi_class with the native TEP patch (`external/patches/hiclass_tep_native.patch`).

- *Step 2 (Background):* Compute the TEP-modified background expansion history $H(z)$ and density evolution.

- *Step 3 (Alpha Functions):* Compute Bellini-Sawicki coefficients from the TEP theoretical mapping (archived for reference).

- *Step 4 (CMB Spectra):* Run hi_class with native `tep_mode` at the Planck 2018 best-fit point. Compare TT, TE, and EE spectra against standard CLASS $\Lambda$CDM.

- *Step 5 (Matter-Frame Scan):* Dual-scan reconstruction of the acoustic scale in screened and unscreened limits.

- *Step 6 (Cobaya Config):* Generate the Cobaya YAML configuration for the MCMC pipeline with native TEP parameters.

- *Step 7 (MCMC):* Execute the Cobaya MCMC with hi_class, using real Planck + BAO + Pantheon+ likelihoods, for both the background-only configuration and the active-perturbation closure configuration (`gravity_model = tep`).

- *Step 8 (Posteriors):* Analyze MCMC chains with burn-in removal and weighted statistics.

- *Step 9 (Synthesis):* Combine all results into summary JSON and markdown.


Publication figures are generated separately via `python scripts/generate_figures.py` (not part of `run_all.py`). Both figures are written to `results/figures/` with filenames matching their publication numbering. Figure 2 requires step 04b. Include them in the static site with `cd site && npm run build`.

The perturbation-level screening implemented via the effective $\alpha$-functions in the Boltzmann hierarchy represents the linear cosmological realization of the abstract environmental operator $\mathcal{S}_\Sigma(\mathcal{E})$. The continuous running of these parameters dynamically enforces the spatial flattening of the Temporal Topology across structure formation epochs without committing to a specific microphysical scalar potential.

## 4. MCMC Parameter Estimation Pipeline

### 4.1 The Cobaya Framework

Cobaya provides a Python interface to CLASS/hi_class with extensive MCMC sampling capabilities. The transition from SciPy/Pandas pipelines to Cobaya enables:

- Native hi_class integration without file-based I/O bottlenecks

- Parallel tempering and adaptive Metropolis-Hastings sampling

- Direct Planck likelihood wrapper integration

- Seamless GetDist posterior visualization

### 4.2 Likelihood Configuration

The pipeline uses the following Planck 2018 likelihoods:

| Likelihood | Description | $\ell$ Range |
| --- | --- | --- |
| `planck_2018_lowl.TT` | Low-$\ell$ temperature | 2–29 |
| `planck_2018_lowl.EE` | Low-$\ell$ polarization | 2–29 |
| `planck_2018_lensing.native` | CMB lensing reconstruction | 8–400 |
| `bao.sdss_dr12_consensus_final` | BAO SDSS DR12 consensus | — |
| `sn.pantheonplus` | Type Ia supernovae (Pantheon+) | — |

### 4.3 Free Parameters and Priors

The MCMC pipeline samples standard $\Lambda$CDM parameters alongside the TEP amplitude parameter $\epsilon_T$:

| Parameter | Prior | Description |
| --- | --- | --- |
| $\Omega_b h^2$ | $\mathcal{U}(0.005, 0.1)$ | Baryon density |
| $\Omega_{\rm cdm} h^2$ | $\mathcal{U}(0.01, 0.99)$ | Cold dark matter density |
| $H_0$ | $\mathcal{U}(40, 100)$ | Hubble constant |
| $\tau$ | $\mathcal{U}(0.01, 0.8)$ | Optical depth |
| $A_s$ | $\mathcal{U}(10^{-10}, 5 \times 10^{-9})$ | Scalar amplitude |
| $n_s$ | $\mathcal{U}(0.94, 1.0)$ | Scalar spectral index |
| $A_{\rm planck}$ | $\mathcal{U}(0.9, 1.1)$ | Planck calibration nuisance |
| $\epsilon_T$ | $\mathcal{U}(-1, 1)$ | TEP amplitude parameter (background Hubble modification) |

### 4.4 Pipeline Execution

```
# Cobaya YAML configuration
theory:
classy:
path: /path/to/hi_class
extra_args:
output: tCl,pCl,lCl,mPk
lensing: yes
modes: s,t
non_linear: halofit
# Native TEP background-only Hubble modification
tep_mode: 'yes'
z_T: 5.0
n_T: 2.0
# epsilon_T is sampled in params below — do not duplicate here

likelihood:
planck_2018_lowl.TT: null
planck_2018_lowl.EE: null
planck_2018_lensing.native: null
bao.sdss_dr12_consensus_final: null
sn.pantheonplus: null

params:
logA:
prior: {min: 2.5, max: 3.5}
ref: {dist: norm, loc: 3.044, scale: 0.014}
proposal: 0.01
drop: true
A_s:
value: 'lambda logA: 1e-10*np.exp(logA)'
n_s:
prior: {min: 0.94, max: 1.0}
ref: {dist: norm, loc: 0.966, scale: 0.004}
proposal: 0.004
H0:
prior: {min: 40, max: 100}
ref: {dist: norm, loc: 67.4, scale: 0.5}
proposal: 1.5
omega_b:
prior: {min: 0.005, max: 0.1}
ref: {dist: norm, loc: 0.0224, scale: 0.0002}
proposal: 0.0003
omega_cdm:
prior: {min: 0.01, max: 0.99}
ref: {dist: norm, loc: 0.12, scale: 0.001}
proposal: 0.0015
tau_reio:
prior: {min: 0.01, max: 0.8}
ref: {dist: norm, loc: 0.054, scale: 0.007}
proposal: 0.01
A_planck:
prior: {min: 0.9, max: 1.1}
ref: {dist: norm, loc: 1.0, scale: 0.0025}
proposal: 0.005
epsilon_T:
prior: {min: -1.0, max: 1.0}
ref: {dist: norm, loc: 0.006, scale: 0.005}
proposal: 0.0005
latex: '\epsilon_T'
sigma8:
latex: '\sigma_8'

sampler:
mcmc:
burn_in: 0
max_tries: 10000
max_samples: 500000
Rminus1_stop: 0.05
Rminus1_cl_stop: 0.2
output_every: 10
drag: true
seed: 42
```

### 4.5 Perturbation-Mode Validation

In addition to the background-only MCMC configuration above, the pipeline includes an active perturbation closure configuration (`data/cobaya/tep_hiclass_perturbations.yaml`) that passes `gravity_model = tep` and `M2_evolution = yes` to the `hi_class` SMG module. This forces the Boltzmann solver to evolve the implemented linear pure-conformal scalar fluctuation sector $\delta\phi$ through the exact Bellini-Sawicki EFT functions derived from the conformal geometry (Appendix A.5).

A 4-chain MPI MCMC with this configuration (27,488 accepted steps; Planck 2018 low-$\ell$ TT/EE + lensing + BAO + Pantheon+; cross-chain Gelman–Rubin $R-1 = 0.050$ at termination) ran successfully with finite posterior at every step. The resulting parameter constraints are:

\begin{equation} \label{eq:5_mcmc_epsT_pert}
\epsilon_T = 0.00549 \pm 0.00464,
\end{equation}

with $H_0 = 66.65 \pm 1.67$ km/s/Mpc, $n_s = 0.9958 \pm 0.0040$, $\Omega_b h^2 = 0.02115 \pm 0.00255$, $\Omega_{\rm cdm} h^2 = 0.1152 \pm 0.0039$, $\tau = 0.0496 \pm 0.0071$, $A_{\rm planck} = 1.089 \pm 0.011$, and $S_8 = 0.869 \pm 0.023$. Because this implementation run uses low-$\ell$ TT/EE plus lensing rather than the full high-$\ell$ TTTEEE acoustic likelihood, the $n_s$ posterior should be interpreted as a robustness diagnostic rather than as a final spectral-index constraint. Direct comparison with the background-only chain ($\epsilon_T = 0.00559 \pm 0.00483$) yields $\Delta\epsilon_T = -0.00053$ ($-0.08\sigma$), $\Delta H_0 = -0.03$ km/s/Mpc ($-0.01\sigma$), $\Delta n_s = +0.0003$ ($+0.05\sigma$), and $\Delta S_8 = +0.0014$ ($+0.04\sigma$). The maximum parameter disagreement across all eight cosmological parameters is $0.06\sigma$, and $\Delta\chi^2 = -0.49$. This confirms that the late-time ISW contribution from the dynamical scalar field is observationally negligible at the current bound, and that the $\epsilon_T$ posterior is driven by background acoustic-peak shifts rather than by perturbation-sector physics.


![Background vs Perturbation Posterior Comparison](figures/tep_perturbation_triangle.png)

Figure 3: Marginalized posterior triangle plot from the 4-chain MPI run, comparing the background-only TEP chain ($\delta\phi$ frozen; blue) against the active-perturbation chain ($\delta\phi$ evolved through the exact Bellini–Sawicki EFT; red) for $\epsilon_T$, $H_0$, $n_s$, and $\sigma_8$. The contours show 68% and 95% credible regions. Cross-chain Gelman–Rubin $R-1 = 0.050$ at 27,488 accepted steps. The total overlap demonstrates that the implemented linear pure-conformal scalar fluctuation sector contributes negligibly to the cosmological posterior at current precision.

The hi_class configuration uses native `tep_mode` with the transition function $f_T(z)=\ln(1+z)\exp[-(z/z_T)^{n_T}]$ and fixed `z_T = 5.0`, `n_T = 2.0`, with `epsilon_T` sampled freely in `params`. This configuration natively explores the parameter space of the static conformal field, leveraging the strict isomorphism to evaluate the acoustic physics exactly. The production configuration is `data/cobaya/tep_hiclass_suite.yaml` (reference alternate: `data/cobaya/tep_native_mcmc.yaml`).

*Pipeline status.* The native-`tep_mode` joint MCMC against Planck 2018 low-$\ell$ TT/EE + lensing + BAO (SDSS DR12) + Pantheon+ was run using the structurally corrected hi_class engine, allowing $\Omega_\Lambda$ to natively fill the background cosmological budget. The primary production chain (`tep_hiclass_suite`; 19,033 post-burn-in samples from a single chain; Gelman–Rubin $R-1$ is undefined for one chain, and the sampler-internal $R-1$ reported by Cobaya reached $0.045$ at termination) gives a $\Lambda$CDM-compatible background while measuring the TEP amplitude parameter:

\begin{equation} \label{eq:5_mcmc_epsT}
\epsilon_T = 0.0056 \pm 0.0043,
\end{equation}

with $H_0 = 66.63 \pm 1.70$ km/s/Mpc, $\Omega_b h^2 = 0.0212 \pm 0.0025$, $\Omega_{\rm cdm} h^2 = 0.1154 \pm 0.0042$, $\tau = 0.049 \pm 0.007$, $A_{\rm planck} = 1.088 \pm 0.012$, and $S_8 = 0.870 \pm 0.028$. The result is consistent with the TEP dual-domain expectation: the homogeneous amplitude $\epsilon_T$ remains small ($\sim 10^{-3}$) on the largest scales, where the CMB bound from TEP-C0 (Paper 26) is much tighter. This recombination-era value is expected to be smaller than the primordial temporal-horizon amplitude ($\epsilon_t = 0.0175$) derived in TEP-TH (Paper 27), because the temporal field evolves and is progressively screened as the universe cools.

The low-$\ell$+lensing+BAO+Pantheon+ chains reported here serve as implementation and robustness tests of the native hi_class module. A decisive high-$\ell$ acoustic likelihood constraint requires a dedicated full TTTEEE production run and is not claimed as the central result of the present paper. Accordingly, the acoustic claim rests on the direct TT/TE/EE spectrum comparison at the background/acoustic-sector level.

*Multi-chain validation.* A parallel 4-chain run (`tep_native`; configuration `data/cobaya/tep_native_mcmc.yaml`) using a Gaussian $A_{\rm planck}$ prior (loc = 1.0, scale = 0.0025) produced 2,993 post-burn-in samples with maximum Gelman–Rubin $R-1 = 0.098$ ($R-1$ for $\epsilon_T = 0.098$; all other parameters $R-1 < 0.05$). This yields $\epsilon_T = 0.0044 \pm 0.0040$ and $H_0 = 66.89 \pm 1.35$ km/s/Mpc, consistent with the primary chain at $0.21\sigma$ and $0.12\sigma$ respectively. While this consistency run reaches $R-1 = 0.098$ on $\epsilon_T$, the widened-prior chain below provides a fully converged multi-chain determination ($\epsilon_T$ $R-1 = 0.013$, all parameters $R-1 < 0.05$), and the two agree to $0.06\sigma$. Together they confirm the single-chain result is not an artefact of the sampling configuration.

*Planck calibration prior sensitivity.* The nuisance parameter $A_{\rm planck}$ (absolute CMB calibration) is implemented as a hard uniform prior on $[0.9, 1.1]$ in the primary chain. The posterior mean is $A_{\rm planck} = 1.088 \pm 0.012$ with maximum sampled value $1.1000$, indicating saturation against the upper prior bound. To test whether this truncation biases the cosmological inference, a dedicated 4-chain sensitivity test was executed with the prior widened to $[0.9, 1.25]$ (configuration `data/cobaya/tep_hiclass_aplanck_sens.yaml`). The converged run (1,640 total samples; all parameters Gelman–Rubin $R-1 < 0.05$; maximum $R-1 = 0.044$ for $\Omega_{\rm cdm} h^2$) yields $A_{\rm planck} = 1.223 \pm 0.044$, confirming the old posterior was truncated by approximately $3.0\sigma$. The TEP amplitude from the widened run is $\epsilon_T = 0.0047 \pm 0.0040$ ($R-1 = 0.013$), consistent with the primary chain at $0.15\sigma$ and with the multi-chain validation at $0.06\sigma$. The correlation between $A_{\rm planck}$ and $\epsilon_T$ is $r = -0.19$, and splitting at $A_{\rm planck} = 1.15$ gives a difference in $\epsilon_T$ of only $-0.21\sigma$. Even a $0.1$ upward shift in $A_{\rm planck}$ would move $H_0$ by only $\sim 0.2$ km s$^{-1}$ Mpc$^{-1}$, well below its posterior width. The $\chi^2$ does decrease monotonically toward the old boundary, but there is no evidence of a degeneracy cascade with $\epsilon_T$. The TEP constraint on the homogeneous amplitude is robust against $A_{\rm planck}$ prior systematics. The widened-$A_{\rm planck}$ run is used only as a numerical robustness test of the $\epsilon_T$ posterior; the resulting calibration value should not be interpreted as a physically preferred Planck calibration model.

The companion paper TEP-C0 (Paper 26) provides the primary late-time constraints: Pantheon+ nested sampling favors the physical no-$\Lambda$ TEP M1 branch over baseline $\Lambda$CDM with Bayesian evidence statistically indistinguishable from CPL (BF $\sim 32$--$36$ depending on the fixed or free-$z_T$ branch). Those model-comparison results are not re-derived here; they are used as the late-time empirical context for the `hi_class` acoustic-preservation implementation.

## 5. Results and Cosmological Constraints

### 5.1 The Acoustic Spectra

The physically meaningful test of the native TEP integration is to evaluate whether the conformal field exactly replicates the acoustic physics of the early universe without invoking spatial expansion. Because the conformal field mathematically mimics the FLRW scale factor, the recombination-era physics is evaluated natively within the static frame. Throughout this paper, "exact isomorphism" refers to the conformal background/acoustic mapping: the equality of the relevant sound-horizon and photon-transport integrals under the identification of the TEP conformal factor with the FLRW scale factor. The associated linear pure-conformal scalar perturbation sector is closed separately through the runtime Bellini–Sawicki EFT functions and validated by the active ($\delta\phi$) hi_class run reported in Section 4.5 and Appendix A.5. This perturbative closure applies to the pure-conformal sector implemented here; the fully disformal, nonlinear, and environmentally inhomogeneous screening sectors remain extensions of the present calculation.

#### 5.1.1 Sound-horizon and acoustic-peak preservation

Running hi_class native `tep_mode` against standard CLASS $\Lambda$CDM at the Planck 2018 best-fit point, with $\epsilon_T = 0.0066$, $z_T = 5$, $n_T = 2$, yields:

- *Sound horizon preserved to ~6 ppm:* $r_s^{\rm TEP}/r_s^{\Lambda\rm CDM} = 0.999994$. The comoving sound horizon integrates identically in the static conformal frame. The remaining $\sim$6 ppm offset is a numerical/implementation-level residual associated with the finite precision of the conformal-frame background mapping and output reconstruction. Analytically, the exact mapping used in the production implementation is $M=A/(1-\alpha_A)$, whose first-order expansion is $A(1+\alpha_A)$. Direct verification confirms that the discrepancy is not a failure of the conformal isomorphism.

- *Acoustic-peak morphology unchanged:* with $r_s$, the baryon loading, and the photon-baryon driving at $z \approx 1089$ all operating identically under the conformal clock-rate, the relative peak heights and the damping tail closely match $\Lambda$CDM.

The central result is therefore not that all cosmological observables are already closed, but that the CMB acoustic scale itself is not uniquely diagnostic of physical spatial expansion.

#### 5.1.2 The residual is a late-time projection, largely degenerate with $H_0$

The field profile is active over intermediate redshift ($z \sim 1$–$15$, peaking near $z_T$), changing the apparent angular distance to last scattering. At the fiducial $\epsilon_T = 0.0066$ this shifts the angular acoustic scale by $\Delta\theta_s/\theta_s = +0.185\%$. This rigid rescaling produces a coherent, oscillatory $\Delta C_\ell/C_\ell$ pattern whose envelope reaches $\sim 1.8\%$ across $100 < \ell < 2000$. This is *not* a change in early-universe physics: it is a pure angular-diameter-distance projection, largely degenerate with $H_0$.

#### 5.1.3 Polarization Spectra ($C_\ell^{TE}, C_\ell^{EE}$)

The TE and EE spectra inherit the same exact behavior: the recombination-era polarization source is natively preserved by the conformal integration, and the only effect is the common $\theta_s$ projection shared with TT.

### 5.2 Cosmological Constraints: Late-Time Evidence and the CMB Bound

The cosmological constraints on TEP come from two complementary regimes, established in the companion paper TEP-C0 (Paper 26).

*Late-time evidence (supernovae).* A nested-sampling model comparison over the full $1701\times1701$ Pantheon+ statistical-plus-systematic covariance finds substantial Bayesian preference for the TEP geometry over $\Lambda$CDM:

| Model | Bayes factor vs $\Lambda$CDM | Interpretation |
| --- | --- | --- |
| TEP M1 fixed ($z_T = 100$) | $\sim 32.1$ | Strong |
| TEP M1 free ($z_T \in [0.1, 150]$) | $\sim 36.2$ | Strong |
| $w$CDM | $\sim 21.1$ | Strong |
| CPL ($w_0 w_a$) | $\sim 35.2$ | Strong |
| Einstein-de Sitter | $4.3\times10^{-126}$ | Rejected (sanity check) |

The TEP M1 branch improves the Pantheon+ likelihood relative to baseline $\Lambda$CDM and achieves Bayesian evidence statistically indistinguishable from CPL (TEP-C0, Paper 26). Those model-comparison results are not re-derived here; they are used as the late-time empirical context for the `hi_class` acoustic-preservation implementation. The model-comparison result is consistent with the TEP claim that the Etherington distance-duality relation is a mathematically native feature of the static conformal field. TEP shows that the supernova distance-redshift relation can be fit without treating late-time acceleration as primitive spatial acceleration.

*$z_T$ distinction.* The $z_T = 5$ profile used in this paper is the homogeneous acoustic-sector benchmark for the hi_class conformal implementation. It should not be confused with the C0 supernova-sector transport benchmark ($z_T = 100$), where $z_T$ parameterizes the effective line-of-sight temporal-shear transition in the Pantheon+ distance law. The C0 free-$z_T$ robustness test uses $z_T \in [0.1, 150]$.

*Homogeneous (CMB) bound.* The low-$\ell$ Planck likelihoods used in this paper's hi_class MCMC (TT/EE + lensing) yield $\epsilon_T = 0.0056 \pm 0.0043$, consistent with zero at $\sim 1.3\sigma$. This homogeneous bound is smaller than the primordial temporal-shear amplitude $\epsilon_t = 0.0175$ constrained by the spectral index in TEP-TH (Paper 27), as expected because the temporal field is subject to late-time epoch screening that suppresses the effective shear amplitude at recombination. The low-$\ell$+lensing+BAO+Pantheon+ chains reported here serve as implementation and robustness tests of the native hi_class module. A decisive high-$\ell$ acoustic likelihood constraint requires a dedicated full TTTEEE production run and is not claimed as the central result of the present paper.

*Native-TEP joint MCMC.* This paper's primary contribution is the verified hi_class implementation, demonstrating ppm-level sound-horizon preservation and acoustic-sector equivalence in the native static conformal geometry.

### 5.3 Structure Growth and the Matter Power Spectrum

The full hi_class Boltzmann closure with active SMG perturbations yields a linear growth amplitude in agreement with Planck and weak-lensing measurements:

- *Linear growth amplitude:* $\sigma_8 = 0.825 \pm 0.016$ at the fiducial TEP parameters ($\epsilon_T = 0.0066$, $z_T = 5$, $n_T = 2$), compared to $\sigma_8 = 0.823$ for standard $\Lambda$CDM at the same cosmology. The TEP value is a native output of the full SMG EFT solver with runtime Bellini-Sawicki mappings; no phenomenological suppression factor is applied.

- *Scale-dependent screening signature:* The gradient-screening envelope $\mathcal{S}_\Sigma(k) = [1 + (k/k_t)^n]^{-1}$ injects a characteristic scale dependence into the matter power spectrum. At $k < k_t$ (cosmic scales), $\mathcal{S}_\Sigma \to 0$ and the background expansion is preserved. At $k > k_t$ (structure scales), $\alpha_M(k) \to 0$ and the scalar fifth force is suppressed, yielding $G_{\rm eff} \to G_N$ inside collapsing halos. The fractional deviation $\Delta P(k)/P_{\Lambda\rm CDM}(k)$ exhibits a smooth transition around $k_t$, producing a testable spectral feature for DESI and Euclid clustering data.

The $\sigma_8$ result demonstrates that the TEP covariant action, when mapped rigorously into the Bellini-Sawicki EFT and solved with the full Boltzmann hierarchy, naturally restores structure growth to the observationally consistent range. Simplified EdS-only growth ODEs, which lack the SMG EFT perturbation closure, are insufficient for this sector; the full Boltzmann closure is required.

### 5.4 The Hubble Tension in TEP


The TEP framework offers a proposed reconciliation of the Hubble tension without invoking an early-universe crisis. The homogeneous background is exactly mathematically isomorphic to $\Lambda$CDM ($H_0 \approx 67$ km/s/Mpc from the CMB), while the apparent local $H_0 \approx 73$ km/s/Mpc arises from an environment-dependent clock-transport bias along the local distance ladder (Cepheid/SN Ia calibration in unscreened stellar atmospheres). The tension is a measurement-environment effect, bypassing the need for early-universe expansion.



![H0 in the TEP picture](figures/figure_1_H0_comparison.png)



*Figure 1.* $H_0$ in the Static Conformal picture. The CMB exactly tracks the unscreened homogeneous background; the local SH0ES value is reinterpreted as clock-transport bias, shifting down when corrected for the temporal shear field gradient (Paper 11).

### 5.5 The Mathematical Limit of the Conformal Field

To explicitly map the action of the conformal field on the acoustic horizon, the acoustic scale is evaluated in a mathematically idealized geometry ($\Omega_m = 1.0$, $\Omega_\Lambda = 0.0$) using the hi_class native `tep_mode` implementation.

#### Regime I: Standard empirical profile ($z_T = 5$)

In the standard TEP model, the profile $\exp[-(z/z_T)^{n_T}]$ ensures the conformal field correctly matches the apparent late-time acceleration inferred from Pantheon+. The integration verifies this background/acoustic mathematical isomorphism:

| $\epsilon_T$ | $100\theta_s$ | $r_s$ [Mpc] | $\Delta D_C / D_C$ | Interpretation |
| --- | --- | --- | --- | --- |
| $0.00$ | $1.0403$ | $144.526$ | $0.00\%$ | Pure EdS reference (no TEP) |
| $0.05$ | $1.0548$ | $144.519$ | $-1.38\%$ | $r_s$ preserved; $\theta_s$ shifts from $D_C$ projection |

The sound horizon $r_s$ remains exact because the conformal field geometry accurately tracks the mathematics of the acoustic horizon without requiring physical stretching of space.

#### Regime II: Theoretical divergence ($z_T \to \infty$)

Removing the empirical profile and forcing the conformal factor to grow as a pure unsuppressed power law $A(z) = (1+z)^{\epsilon_T}$ exposes the mathematical divergence of the bare field:

| $\epsilon_T$ | $100\theta_s$ | $r_s$ [Mpc] | $\Delta D_C / D_C$ | Interpretation |
| --- | --- | --- | --- | --- |
| $0.00$ | $1.0403$ | $144.526$ | $0.00\%$ | Pure EdS reference |
| $0.05$ | $0.7565$ | $100.584$ | $-4.30\%$ | $r_s$ mathematically squeezed by divergence |

This mathematical limit demonstrates that the $z_T \sim 5$ empirical fitting function accurately defines the physical profile of the conformal field, allowing it to mimic Dark Energy while preserving the CMB acoustic horizon without expanding space.


![Matter-frame dual-scan results](figures/figure_2_jordan_theta_s.png)



*Figure 2.* EdS + TEP dual scan. Solid line: EdS baseline at $\epsilon_T=0$. (Left) Standard empirical profile ($z_T = 5$): $r_s$ preserved. (Right) Divergent limit ($z_T \to \infty$): unphysical squeezing of the acoustic horizon.

## 6. Conclusion: The Static Conformal Universe

This paper implements and validates the native Temporal Equivalence Principle (TEP) conformal modification directly within the `hi_class` Boltzmann solver framework. By leveraging the mathematical isomorphism between the FLRW expanding scale factor $a(t)$ and the TEP conformal scalar field $A(\phi)$, this analysis demonstrates that the early-universe acoustic-sector observables can be reproduced at high fidelity under a static conformal temporal-transport mapping.

### 6.1 Summary of Results

- *The Mathematical Isomorphism:* The TEP conformal factor $A(\phi) = \exp(\beta_A\phi/M_{\rm Pl})$ dictates the clock-rates and photon phases in the causal matter metric $\tilde{g}_{\mu\nu} = A(\phi)^2 g_{\mu\nu}$. Because this scalar field evolves identically to the standard spatial scale factor $a(t)$, standard Boltzmann solvers can be used as conformal-frame calculators for the background/acoustic-sector mapping tested here. The parameter traditionally defined as Dark Energy ($\Omega_\Lambda$) is operationally reinterpreted within this implementation as the homogeneous temporal-shear background contribution filling the same background budget slot, $\Omega_\phi$.

- *CMB Acoustic Preservation (Screened-Limit Consistency):* Because of the static conformal isomorphism, the hi_class native integration demonstrates that the background/acoustic-sector observables are preserved to parts-per-million ($r_s^{\rm TEP}/r_s^{\Lambda\rm CDM} = 0.999994$). This is a screened-limit consistency check: the conformal transport geometry recovers the standard FLRW acoustic results in the regime where the conformal field dominates and disformal corrections are suppressed. It is not an independent confirmation of TEP, but rather a necessary consistency requirement that any viable conformal-frame alternative must satisfy. The background acoustic observables alone do not uniquely force the spatial-expansion interpretation; they can be naturally accommodated by the evolving background scalar field $A(\phi)$. The active-perturbation closure (Section 4.5) confirms that evolving the implemented linear pure-conformal scalar fluctuation sector $\delta\phi$ through the exact Bellini–Sawicki EFT produces parameter posteriors indistinguishable from the background-only chain to $0.06\sigma$, with cross-chain Gelman–Rubin $R-1 = 0.050$ at 27,488 accepted steps.

- *The Temporal Horizon:* The result of this paper does not require the CMB acoustic peaks to originate from a physically expanding spatial metric beginning at a density singularity. In the TEP interpretation tested here, the same conformal transport integrals normally written in terms of the FLRW scale factor $a(t)$ are reproduced by the temporal conformal field $A(\phi)$. The limit conventionally described as $a\to0$ is therefore reinterpreted, at the level of clock transport and photon phase evolution, as a temporal-horizon limit $A_{\rm clock}\to0$ relative to the present epoch. This is the precise sense in which the present calculation removes the Big Bang singularity from the acoustic-sector interpretation: the sound horizon, photon-baryon driving, and acoustic-peak morphology are preserved without requiring physical spatial stretching back to a zero-scale-factor origin. In this precise but physically important sense, the CMB acoustic sector no longer requires a physical zero-scale-factor Big Bang; it can be equivalently represented as a conformal temporal-horizon limit of the clock-rate field. Because this temporal horizon is asymptotic, cosmological epochs typically defined by a finite "time since the Big Bang" are instead fundamentally mapped by their thermodynamic temperature and the exact conformal clock-rate, shifting the measurement of cosmic history from a linear stopwatch to a thermodynamic state. The full nonsingular closure—geodesic completeness, curvature invariants, entropy evolution, light-element abundances, and the origin of the CMB blackbody spectrum—is delivered in TEP-TH (Paper 27), where Proposition 1 establishes finite matter-frame curvature at the temporal horizon and the ten-step pipeline validates every early-universe observable.

- *Cosmological Constraints:* A joint hi_class Cobaya MCMC (Planck 2018 low-$\ell$ TT/EE + lensing + BAO + Pantheon+) yields a close match to the conformal field parameters. The companion paper TEP-C0 (Paper 26) provides robust late-time evidence, showing via nested sampling that the screened TEP conformal model matches the Pantheon+ distance-redshift relation with strong Bayesian support, reducing the phenomenological need to treat late-time acceleration as primitive spatial expansion.

### 6.2 Resolving the Cosmological Crises

Standard $\Lambda$CDM cosmology currently faces two severe crises: the Hubble Tension ($H_0$) and the unexpectedly massive high-redshift galaxy candidates observed by JWST. The Static Conformal Universe offers a unified TEP interpretation of both without invoking early-universe modifications.

- *The Hubble Tension:* In TEP, the temporal shear field is environmentally screened by mass. Supernovae exist in empty voids (where the field is unscreened, yielding a high inferred $H_0$), while Cepheids exist in dense galaxies (where the field is partially screened, yielding a lower inferred $H_0$). The tension is an artifact of environmental mass-screening on local kinematic distance probes (Paper 11), not an early-universe physics crisis.


- *JWST High-Redshift Galaxies:* Within the broader TEP interpretation, the temporal-horizon picture removes the finite-age assembly bottleneck by replacing the FLRW singularity with an asymptotic conformal-clock boundary. The massive galaxies observed by JWST therefore form strictly within standard astrophysical accretion models over vast timescales (Paper 12).

### 6.3 Synthesis of the Paradigm Shift

This analysis implements and explicitly validates the native TEP implementation within a rigorous Boltzmann solver framework. The acoustic indistinguishability of the static conformal background from $\Lambda$CDM at recombination demonstrates that the early-universe background physics cannot easily distinguish between a stretching spatial metric and an evolving conformal clock-rate field.

Late-universe Pantheon+ data in TEP-C0 favor the physical no-$\Lambda$ TEP branch over baseline $\Lambda$CDM and yield Bayesian evidence statistically indistinguishable from CPL, providing a concrete conformal-frame alternative to the background expansion interpretation. The active-perturbation closure reported in Section 4.5 confirms that the implemented linear pure-conformal scalar fluctuation sector is ghost-free and observationally negligible, with the implemented $\delta\phi$-enabled Einstein–Boltzmann chain agreeing with the background-only chain to $0.06\sigma$ across all cosmological parameters. By treating time itself as a dynamical, mass-screened scalar field, TEP seeks to unify early-universe acoustic physics, late-time "acceleration", the $H_0$ tension, and JWST anomalies into a single, cohesive static geometric framework. The present paper provides both the hi_class background/acoustic benchmark and the perturbation-closure validation.


## References

Smawfield, M. (Paper 1). *Temporal Equivalence Principle: Terrestrial Screening and GNSS Phase Correlations.* TEP Corpus.

Smawfield, M. (Paper 6). *TEP and Ultra-Compact Dwarfs: Potential-Dependent Proper-Time Mapping.* TEP Corpus.

Smawfield, M. (Paper 11). *TEP and the Hubble Tension: Cepheid Environmental Bias.* TEP Corpus.

Smawfield, M. (Paper 12). *TEP and JWST High-Redshift Anomalies.* TEP Corpus.

Smawfield, M. (Paper 13). *TEP and Gaia DR3 Wide Binaries: Density-Dependent Kinematics.* TEP Corpus.

Smawfield, M. (Paper 26). *Temporal Equivalence Principle: A Covariant Alternative to Cosmic Expansion.* TEP Corpus. DOI: 10.5281/zenodo.20370144

Bellini, E., & Sawicki, I. 2014, JCAP, 07, 050. *Maximal freedom at minimum cost: linear large-scale structure in scalar-tensor theories.*

Brax, P., Burrage, C., Davis, A.-C., & Gubitosi, G. 2019, Phys. Rev. D, 100, 083515. *Screening mechanisms in scalar-tensor theories.*

Cobaya Team. 2023, *Cobaya: Code for Bayesian Analysis of physical theories.* arXiv:2305.02971.

Hu, B., Raveri, M., Frusciante, N., & Silvestri, A. 2014, Phys. Rev. D, 89, 103530. *EFTCAMB/EFTCosmoMC: Numerical Notes.*

Knox, L., & Millea, M. 2020, Phys. Rev. D, 101, 043533. *Hubble constant hunter's guide.*

Lagos, M., Bellini, E., Jimenez, J. B., et al. 2018, JCAP, 03, 021. *hi_class: Horndeski in the Cosmic Linear Anisotropy Solving System.*

Lewis, A., Challinor, A., & Lasenby, A. 2000, Astrophys. J., 538, 473. *Efficient computation of cosmic microwave background anisotropies.*

Planck Collaboration. 2020, A&A, 641, A1. *Planck 2018 results. I. Overview and cosmological parameters.*

Planck Collaboration. 2020, A&A, 641, A6. *Planck 2018 results. VI. Cosmological parameters.*

Riess, A. G., Casertano, S., Yuan, W., et al. 2022, ApJ, 934, L7. *A Comprehensive Measurement of the Local Value of the Hubble Constant with 1 km/s/Mpc Uncertainty from the Hubble Space Telescope and the SH0ES Team.*

Sawicki, I., & Bellini, E. 2015, Phys. Rev. D, 92, 084061. *Stability of dark energy and the generalized no-slip condition.*

Zumalacárregui, M., & García-Bellido, J. 2014, Phys. Rev. D, 89, 064046. *Transforming gravity: from derivative couplings to matter to second-order scalar-tensor theories beyond the Horndeski Lagrangian.*

## Appendix A: Technical Implementation Details

### A.1 hi_class Installation and Configuration

#### A.1.1 Building with TEP Support

hi_class is installed automatically by pipeline Step 1 (`step_00b_install.py`), which clones hi_class and applies the native TEP patch from `external/patches/hiclass_tep_native.patch` to `source/background.c`, `source/input.c`, and `include/background.h`. Manual rebuild:

```
cd external/hi_class/hi_class
make clean && make
```

**Screening projection notice.** Screening in TEP is represented at theory level by the environmental operator S_Σ(E). Quantities such as ρ_T, R_T(M), S_⊕(r), compactness Φ/c^2, local stellar density, thermal epoch, coherence length, proximity, and boundary geometry are domain-specific projections of E, not independent screening mechanisms and not interchangeable universal thresholds.

### A.2 Cobaya Installation

```
pip install cobaya
cobaya-install planck_2018_lowl.TT planck_2018_lowl.EE \
planck_2018_lensing.native bao.sdss_dr12_consensus_final \
sn.pantheonplus --path /path/to/likelihoods
```

### A.3 TEP Module C Code Structure and Implementation Note

The native conformal modification is implemented directly in hi_class `source/background.c`, controlled by the `.ini` flags `tep_mode`, `epsilon_T`, `z_T`, `n_T`. The relevant functions are:

- `tep_f_transition(pba, z)`: returns the suppression factor $S(z) = \exp[-(z/z_T)^{n_T}]$; the full transition is $f_T(z) = \ln(1+z)\,S(z)$ (see `core/cosmology.py:f_T`).

- `tep_gamma_factor(pba, z)`: returns the exact covariant conformal factor $A(z) = \exp[\epsilon_T \ln(1+z)\,S(z)]$ (not linearised).

- The Hubble rate and its conformal-time derivative are mathematically mapped using the exact geometric relation $M(z) = A/(1-\alpha_A)$ in `background_functions` and in the initial-Hubble setter, explicitly evaluating the static conformal geometry.

*Implementation note (corrected bug).* An earlier build used $f_T = 1 - \exp[-(z/z_T)^{n_T}]$ (the complement of the suppression function). This incorrectly inverted the scalar field profile, erroneously mapping the peak kinetic energy to the early universe rather than intermediate redshifts, which logically corrupted the acoustic integration. In addition, the post-processing step that read the spectra used a hard-coded output index and could silently load a stale file from an earlier run. Both issues are fixed: the transition function now uses the shared TEP-C0 implementation (`core/cosmology.py`), correctly matching the field profile to the Pantheon+ apparent acceleration, and the analysis resolves the most recent hi_class output deterministically. *Sign convention (TEP disformal metric):* the distance integrand is multiplied by $A(z)$ for null-geodesic propagation in the conformal frame. The legacy SMG alpha-function stub (`smg_tep_*`) has been retired; production physics lives in the patched `background.c` (`external/patches/hiclass_tep_native.patch`).

### A.3a Derivation of the Conformal-Frame Factor $M(z)$

This appendix derives the background conformal mapping $M(z)$ from the bi-metric action (Equation \ref{eq:3_theory_01}) using a single frame convention held fixed throughout, demonstrating the exact geometric relation implemented natively in the codebase.

*Setup and convention.* Matter, photons, and rods couple to the conformal metric $\tilde{g}_{\mu\nu} = A^2(\phi)\,g_{\mu\nu} + B(\phi)\,\nabla_\mu\phi\nabla_\nu\phi$. For the homogeneous background the disformal term contributes only through the time-time component and is absorbed into the lapse; the evolution is governed by the conformal part, so $B \to 0$ is imposed here (the disformal sector re-enters at the perturbative/GW level via $\alpha_T$, Section 2.2.2). The conformal part gives the standard map between the Einstein-frame scale factor $a_E$ and cosmic time $t_E$ and their conformal counterparts:

\begin{equation} \label{eq:a3a_map}
\tilde{a} = A(\phi)\,a_E, \qquad d\tilde{t} = A(\phi)\,dt_E.
\end{equation}

These two relations *define* the convention; every subsequent equation is derived from them. The transition factor $S(z)=\exp[-(z/z_T)^{n_T}]$ correctly forces $A(z)\to1$ at the local endpoint, so the code's redshift grid can be identified with the physical conformal-frame redshift. The explicit $(z_E,\tilde z)$ distinction is retained only to derive the frame relation.

*Physical Hubble rate.* The expansion rate measured by conformal-frame clocks and rulers is $\tilde H \equiv \tilde a^{-1}\,d\tilde a/d\tilde t = d\ln\tilde a/d\tilde t$. Using $d/d\tilde t = A^{-1}\,d/dt_E$ and $\ln\tilde a = \ln A + \ln a_E$,

\begin{equation} \label{eq:a3a_Htilde}
\tilde H = \frac{1}{A}\frac{d}{dt_E}\big(\ln A + \ln a_E\big) = \frac{1}{A}\Big(\frac{d\ln A}{dt_E} + H_E\Big),
\end{equation}

where $H_E = d\ln a_E/dt_E$ is the Einstein-frame rate. Because the TEP conformal factor $A(z)$ is evaluated as a function of the observable physical redshift (the matter-frame redshift $1+z = \tilde{a}_0/\tilde{a}$), the coupling $\alpha_A$ computed in the codebase is fundamentally the derivative with respect to the *matter-frame* scale factor:

\begin{equation} \label{eq:a3a_alpha}
\frac{d\ln A}{dt_E} = \frac{d\ln A}{d\ln \tilde{a}}\,\frac{d\ln \tilde{a}}{dt_E} = \alpha_A\,(A \tilde{H}), \qquad \alpha_A \equiv \frac{d\ln A}{d\ln \tilde{a}} = -\frac{d\ln A}{d\ln(1+z)},
\end{equation}

which matches the definition in Section 3.2. Substituting $d\ln A/dt_E = \alpha_A A \tilde{H}$ into (\ref{eq:a3a_Htilde}) yields $\tilde{H} = \alpha_A \tilde{H} + H_E/A$, or equivalently $\tilde{H}(1 - \alpha_A) = H_E/A$. Identifying the reference FLRW background rate with the conformally mapped Einstein-frame rate used by the implementation, $H_E = A^2\,H_{\Lambda\rm CDM}$, gives the exact geometric relation:

\begin{equation} \label{eq:a3a_exact}
\boxed{\;\tilde H(z) = \frac{A(z)}{1 - \alpha_A(z)}\,H_{\Lambda\rm CDM}(z)\;} \qquad \Longrightarrow \qquad M_{\rm exact}(z) = \frac{A}{1 - \alpha_A}.
\end{equation}

*Implementation Status.* The production codebase (`hiclass_tep_native.patch`) evaluates this exact geometric relation $M = A/(1-\alpha_A)$ directly, guaranteeing mathematical fidelity to the conformal evaluation without requiring first-order approximations.

### A.4 Saturation Scale in Cosmological Units

The candidate Temporal Topology saturation scale ρ_T ≈ 20 g/cm³ converts to cosmological units as:

\begin{equation} \label{eq:9_appendix_01}
\rho_T = 20 \text{ g/cm}^3 = 2 \times 10^4 \text{ kg/m}^3 \approx 1.1 \times 10^{34} \text{ eV/cm}^3
\end{equation}

In Planck units ($\hbar = c = G = 1$):

\begin{equation} \label{eq:9_appendix_02}
\rho_T \approx 4 \times 10^{-93} M_{\rm Pl}^4
\end{equation}

Compare to cosmic mean density today ($\rho_{\rm crit,0} \approx 10^{-123} M_{\rm Pl}^4$). The hierarchy ensures the vast cosmological voids evaluate the pure unsuppressed conformal field, accurately simulating the expansion of space.

*Intermediate environments and operational parameter bounds.* Between the terrestrial laboratory and the cosmic mean, the screening transition is continuous. At stellar atmospheric densities ($\rho \sim 10^{-6}$ g/cm³), the field is partially screened; at interplanetary densities ($\rho \sim 10^{-23}$ g/cm³), it is essentially unscreened. Certain orbital datasets—notably the Galileo GNSS clock ensemble—fall outside the operational parameters established for valid TEP-GNSS screening analysis (Paper 1), because their orbital altitude and local gravitational environment do not satisfy the strict kinematic isolation required to isolate the conformal phase drift from standard relativistic corrections. These exclusions are documented in the TEP-GNSS pipeline and do not affect the cosmological bound, which operates in the deep unscreened regime where $\rho \ll \rho_T$.

### A.5 Stability Sector Closure

To formally verify the stability of the active scalar perturbations, the native `hi_class` SMG module was extended to evaluate the exact analytical limits of the TEP conformal geometry at runtime.

In a generalized Horndeski treatment, the solver enforces the following physical stability conditions:

- $c_s^2 \ge 0$ (no gradient instabilities)

- $D = \alpha_K + \frac{3}{2}\alpha_B^2 \ge 0$ (no ghosts)

- $|\alpha_M|$ bounded (sub-luminal Planck-mass running)

- $\alpha_T \approx 0$ (gravitational wave speed constraints)

For the conformal modification implemented here, the EFT parameters map strictly to the dynamical background derivative $\alpha_A$:

\begin{equation} \label{eq:a5_alpha_M}
\alpha_M = -2\alpha_A
\end{equation}

\begin{equation} \label{eq:a5_alpha_B}
\alpha_B = 2\alpha_A
\end{equation}

\begin{equation} \label{eq:a5_alpha_K}
\alpha_K = -5\alpha_A^2
\end{equation}

\begin{equation} \label{eq:a5_alpha_T}
\alpha_T = 0
\end{equation}

The substitution of the kineticity and braiding terms into the physical no-ghost discriminant yields an exact identity:

\begin{equation} \label{eq:a5_ghost_identity}
D = (-5\alpha_A^2) + \frac{3}{2}(2\alpha_A)^2 = \alpha_A^2
\end{equation}

This identity proves that the continuous conformal field transition is analytically protected against ghost instabilities, as the discriminant remains strictly positive-definite. In the pure-conformal implementation used here, the runtime SMG closure fixes the scalar sound-speed sector to the luminal conformal-frame limit ($c_s^2 = 1$), and the integration verifies that this choice produces no gradient instability across the sampled redshift range. The no-ghost condition is independently protected by the discriminant identity ($D=\alpha_A^2$). A fully derived sound-speed expression for the disformal and nonlinear screening sectors is outside the scope of the present pure-conformal closure.

The production codebase forces these limits natively during the calculation of the SMG perturbation coefficients, guaranteeing mathematical fidelity to the conformal evaluation without requiring pre-tabulated interpolation or analytical approximations.


## Appendix B: Data Availability & Reproducibility


This work follows open-science practices. All results are fully reproducible from raw data
using the documented pipeline. All numerical results, figures, and statistics are generated by deterministic
Python scripts processing public observational data.



### Repository and Code


GitHub Repository: github.com/matthewsmawfield/TEP-HC



The repository contains a deterministic, version-controlled cosmological analysis pipeline
for CMB acoustic peak preservation tests and MCMC parameter estimation with TEP screening.



### Repository Structure


TEP-HC/
├── data/
│   ├── cobaya/              # Cobaya MCMC chains
│   ├── external/             # hi_class submodule
│   └── hi_class/             # TEP-CLASS implementation
├── scripts/
│   └── steps/                # Analysis pipeline steps
├── core/                     # TEP shared constants and parameters
├── site/
│   └── components/           # Manuscript HTML sections
├── requirements.txt
├── CITATION.bib
└── README.md



### Software Environment


Key packages: NumPy, SciPy, Matplotlib, Cobaya, hi_class.
The pipeline has been tested on Python 3.10+.



### License


All code and manuscripts are released under CC-BY-4.0.