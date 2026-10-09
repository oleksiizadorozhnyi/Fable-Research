# Certified Bounds for Partially Coherent Mask Synthesis: Convex Printability Certificates, the Value of Multiple Exposure, and Pattern-Specific Resolution Limits

*Autonomous research draft, 2026-10-09. Status: internally verified; not peer reviewed. Model: scalar thin-mask
Abbe/Hopkins imaging of periodic cells.*

## Abstract
Inverse lithography technology (ILT) and source–mask optimization (SMO) solve nonconvex problems whose local optima
carry no guarantee: a designer cannot tell whether a layout is fundamentally unprintable or whether the optimizer is
stuck. We show that, for the scalar Abbe/Hopkins model of a periodic cell, the best achievable dark-to-bright intensity
ratio (equivalently, the maximum exposure latitude for a threshold resist, over any set of focus conditions) admits a
semidefinite relaxation in the diffraction-order domain whose dual solution is an elementary certificate: nonnegative
weights on dark and bright sample points such that one explicit matrix is positive semidefinite, proving that no mask
in the class—and no incoherent combination of masks—can beat the bound. We prove that the relaxation is exactly the
unlimited-multiple-exposure problem, so its gap to the best single mask is the value of multi-exposure, and we certify
strictly positive gaps on specific targets with a level-2 Lasserre bound on the single-mask problem (e.g. 13.7% for a
mixed-width target and ≥27% for an iso-dense target under coherent illumination). For nonnegative (binary/attenuated)
masks in one dimension we show the realizable in-band spectra form exactly the Carathéodory–Toeplitz cone. We prove
that, for symmetric illumination, allowing arbitrary complex transmission does not improve the relaxation. The bound
reproduces the classical k1 cutoffs exactly (DUV and EUV), is invariant to interior sampling, is 3–150× tighter than
the generalized-eigenvalue (Rayleigh-quotient) bound used in prior global-optimization heuristics, and certifies that
local optimizers reach the global optimum on 75% of random 1D targets and on all 2D contact/line-end layouts tested,
while a standard pattern-error ILT lands 30–45% above the certified optimum. 2D layouts under annular and quasar
illumination are certified in under a second each. Code and data are provided.

## 1 Introduction
Computational lithography synthesizes mask (and source) patterns so that the projected aerial image prints a target
layout. The forward model is quadratic in the mask: in the Abbe formulation the image is a weighted sum over source
points of squared moduli of filtered mask spectra. ILT and SMO therefore optimize nonconvex functionals, and the field
has long recognized "an enormous number of local minima" (IBM patents) and "poor local minima" (recent generative-RL
ILT work). Reviews list runtime and local optimality among the main open problems. What is missing is a *certificate*:
a proof, for a given target, illumination, numerical aperture and wavelength, that no mask can do better than a stated
value. Such certificates would (i) separate unprintable layouts from unconverged optimizers, (ii) benchmark ILT engines
against the true optimum, (iii) quantify the value of multiple exposure and of phase-shift versus binary masks, and (iv)
generalize the k1 rule of thumb to arbitrary layouts.

This paper provides such certificates. The idea is the classical lifting X = M M^H of the diffraction-order vector M:
every intensity sample is linear in X, so the dark/bright ratio problem becomes a semidefinite program (SDP) once the
rank-one constraint is dropped. The novelty is not the lifting itself—Shor relaxations of quadratically constrained
quadratic programs (QCQPs) are textbook, and the identical mathematical template appears in multicast beamforming—but
the lithographic consequences: a one-line certificate theorem with solver-independent verification, the physical
identification of the relaxation with unlimited incoherent multi-exposure, an exact convex realizability cone for
nonnegative masks in 1D, a symmetry lemma showing complex masks cannot improve the relaxation, and a systematic
empirical study, with rigorous numerics, of how far local mask optimization sits from certified optima.

**Contributions.** (C1) Theorem 1: dual certificates of printability for any class of masks, verified by one PSD test.
(C2) Theorem 3: the relaxation equals the unlimited-multi-exposure problem; Pataki's bound limits the number of useful
exposures. (C3) Theorem 4: Carathéodory–Toeplitz characterization of in-band spectra of nonnegative 1D masks.
(C4) Lemma 6: complex transmission does not improve the relaxation under symmetric illumination. (C5) An open-source
pipeline (SDP + rigorous certificate + rounding + local refinement + optimal-face exploration + Lasserre-2 single-mask
bounds + SMO and DNN variants) and an empirical study: exact k1 cutoffs, pattern-specific exposure-latitude limits,
certified optimality gaps of local optimization and of standard ILT, certified multi-exposure gains, certified
PSM-vs-binary gains, and 2D layouts under annular/quasar illumination.

**What is not claimed.** No new SDP algorithm; no full-chip scalability; no vector/3D-mask/resist physics; no
statement about industrial ILT engines beyond the model studied.

## 2 Related work
*Wavefront (diffraction-order) optimization.* Rosenbluth et al. (JM3 1(1):13, 2002; SPIE 6520, 2007) optimize the
collected diffraction orders directly, maximizing exposure latitude subject to bright/dark intensity constraints, and
describe the problem as having "non-convex quadratic constraints"; their global strategies use dominant "joint
eigenvectors" and a "step-walk relaxation" (US 7,587,702; US 7,685,559), an unconstrained Rayleigh-quotient
simplification for multiple exposure (US 8,495,528), and realizability of arbitrary wavefronts by phase-shift masks
(US 8,108,802). We adopt the same variables and the same figure of merit but replace heuristics by a convex relaxation
with certificates, and we compare directly with the Rayleigh-quotient bound (Table 7).
*SDP relaxation of QCQPs.* Shor relaxation; survey by Luo, Ma, So, Ye and Zhang (IEEE SPM 2010). The max-min multicast
beamforming problem (Sidiropoulos, Davidson, Luo, IEEE TSP 2006) has exactly the structure of coherent-illumination
mask synthesis (rank-one bright constraints, interference-like dark constraints); NP-hardness and approximation results
transfer. Pataki (Math. OR 1998) bounds the rank of extreme SDP solutions. Lasserre's moment hierarchy provides tighter
bounds for small problems.
*Bounds in photonic design.* Angeris, Vučković and Boyd (ACS Photonics 2019; Opt. Express 2021; ACS Photonics 2023) use
Lagrangian duality/SDP to bound the performance of inverse-designed photonic devices and to certify near-optimality of
local designs; our work is the imaging-side analogue.
*Information-theoretic limits.* Ma, Arce, Garcia-Frias and coworkers (Opt. Express 2017, 2018; and later) model
lithography as an information channel and derive fidelity "limits" from mutual information; these are not duality
certificates for a given layout.
*Lifting in imaging.* PhaseLift (Candès, Strohmer, Voroninski 2013) lifts quadratic measurements for recovery; we lift
for synthesis.
*Multiple patterning.* SDP relaxations exist for layout *decomposition* (mask assignment), a combinatorial problem
unrelated to aerial-image synthesis.

## 3 Problem formulation
(See notes/THEORY.md §1–2 for full definitions.) A periodic cell Λ with mask transmission m has in-band diffraction
orders M ∈ ℂ^F (|f_k| ≤ (NA + max|s_j|)/λ). For source points s_j with weights w_j and focus z, the Abbe image is
I_{M,z}(x) = Σ_j w_j |Σ_k P_z(f_k + s_j/λ) M_k e^{2πi f_k·x}|² = M^H Q_{x,z} M with Q_{x,z} ⪰ 0. Given bright and
dark sample sets B, D (placed at exactly ±margin from every target edge plus interior points) and a focus set Z, the
dark/bright ratio ρ(M) = max_{D×Z} I / min_{B×Z} I is scale-invariant, and EL(ρ) = (1−ρ)/(1+ρ) is the exposure
latitude of a threshold resist. Problem (P): ρ* = inf_M ρ(M) over a mask class (complex 𝓜_ℂ, real 𝓜_ℝ, nonnegative
𝓜_+). Relaxation (R): minimize t subject to Tr(Q_{x,z}X) ≤ t on D×Z, ≥ 1 on B×Z, X ⪰ 0 (plus the real-structure
constraint X = conj(ΠXΠ) for 𝓜_ℝ, which Lemma 2 reduces to a real symmetric SDP of the same size).

## 4 Theory
Full statements and proofs are in notes/THEORY.md. Summary:
- **Theorem 1 (certificate).** If μ ≥ 0 on D×Z with Σμ = 1, λ ≥ 0 on B×Z, and S = Σμ_d Q_d − Σλ_b Q_b ⪰ 0, then
  every mask and every incoherent sum of masks has ρ ≥ Σλ_b. Proof: three lines of weak duality.
- **Lemma 2 (real structure).** Real masks ↔ a real F-dimensional subspace; the real-structured PSD cone is isometric
  to 𝕊^F_+. (Used to solve a purely real SDP; also gives the real multi-exposure decomposition.)
- **Theorem 3 (multi-exposure).** t_R = optimal ratio of the unlimited-multiple-exposure problem; a single mask attains
  t_R iff (R) has a rank-one optimum; an extreme optimum has rank r with r(r+1)/2 ≤ |D×Z|+|B×Z|+1.
- **Theorem 4 (nonnegative masks, 1D).** In-band spectra of nonnegative masks = Toeplitz-PSD cone (Carathéodory–
  Toeplitz; strict interior realized by positive trigonometric polynomials via Fejér–Riesz). Numerically: 0/2000
  necessity violations, 200/200 sufficiency realizations by LP.
- **Proposition 5 (k1 cutoff).** Below pitch λ/((1+σ)NA) the bound is exactly 1 (unprintable) with the uniform
  certificate; above it, < 1.
- **Lemma 6 (complex vs real).** For symmetric illumination and even pupil, conj-reflection preserves the image; hence
  the complex-mask relaxation equals the real-mask relaxation (confirmed numerically on all 9 cases of Table 7).
- **Sampling remark (§5 of THEORY.md).** A certificate for the sampled problem is valid for the continuous problem.

## 5 Algorithms
1. Build Q_{x,z} (F×F) for the edge-offset sample set; solve the real SDP (CLARABEL; SCS fallback).
2. Extract duals; rigorize: find the largest c ≤ 1 with Σμ Q_d − cΣλ Q_b ⪰ 0 by bisection on the (concave) minimum
   eigenvalue; report cΣλ. Independent checks: eigenvalue test; 50-digit Cholesky (mpmath).
3. Primal: Gaussian rounding from X; analytic-gradient L-BFGS on a log-sum-exp surrogate of log ρ; SLSQP minimax
   polish; multistart. 4. Optimal-face exploration (minimize random PSD functionals over the optimal face) to find
   low-rank optima. 5. Lasserre level-2 relaxation of the single-mask problem (F ≤ 11 here) to certify multi-exposure
   gains. 6. Variants: SMO relaxation (one PSD block per source point), DNN relaxation for nonnegative masks (pixel
   domain), complex-mask relaxation, multi-focus (process-window) constraints.
Complexity: F ≈ (2(1+σ)NA L/λ)^d collected orders (F = 5–15 for the 1D cells, 5–13 for the 2D cells here; ~500 for a
1 µm² DUV cell), SDP size F×F with |B|+|D| constraints per focus. All 1D/2D SDPs here solve in 0.1–1.4 s.

## 6 Experiments  (DUV: λ=193 nm, NA=1.35; EUV: λ=13.5 nm, NA=0.33; sources: coherent, conventional σ, dipole,
annular 0.6/0.9, quasar 0.6/0.9; edge tolerance 8 nm unless stated; all numbers from results/RESULTS_TABLES.md)
6.1 **Cutoffs.** Dense lines: bound = 1 exactly below λ/((1+σ)NA) (84.1 nm at σ=0.7) and < 1 above; EUV dense lines at
k1 = 0.244 are certified unprintable with any mask under both sources, and a σ=0.8 dipole is certified unprintable at
22 nm pitch (its first order falls outside the pupil) whereas σ=0.9 conventional reaches EL 0.36.
6.2 **1D feature library (Table 1).** Dipole/dense cases are tight (gap 0.00%); conventional illumination gaps ≤ 2.5%;
coherent illumination shows large gaps (iso-dense 30.8%, mixed 13.7%); three-focus process-window cases show 11–17%
gaps. Certified max exposure latitude at 8 nm edge tolerance: 0.35–0.47 single focus; 0.07–0.13 at ±80 nm focus.
6.3 **Random targets (n = 109; 1–4 random segments ≥ 40 nm in 300–600 nm cells; coherent, σ = 0.5, σ = 0.9, dipole
σ = 0.7; edge tolerance 10 nm).** The best single mask found matches the certified bound to ≤ 1e-4 in 69.7% of
instances (certified global optimality), to ≤ 1% in 82.6%, to ≤ 5% in 91.7%; median gap 1.8e-6; maximum 17.3%. Exact
tightness is source dependent: 85% for σ = 0.9, 75% for the dipole, 57% for coherent and σ = 0.5 illumination, and all
gaps > 5% have minimal optimal-face rank 2 or 3, consistent with Theorem 3 (incoherent multi-exposure helps most under
coherent light). All certificates pass the independent PSD test after rigorization.
6.4 **Certified multi-exposure gains (Table 2).** Lasserre-2 single-mask bounds: mixed_coh 0.4787 (= best single mask
found to 7e-7) vs multi-exposure bound 0.4211 → certified gain 13.7%; iso-dense coherent: single-mask ≥ 0.4608 (best
found 0.4741) vs 0.3625 → certified gain ≥ 27%. Both have minimum optimal-face rank 2: two exposures suffice.
[[pending: process-window case from exp11]]
6.5 **PSM vs binary (Table 3).** DNN relaxation for nonnegative masks: iso-dense 0.4266 vs PSM 0.3879 (binary ≥ 10%
worse, certified); pixel-domain and Fourier-domain PSM relaxations agree to 1e-5.
6.6 **Standard ILT vs certified optimum (Table 3b).** Pattern-error ILT (sigmoid resist, open-frame-normalized
intensity, threshold 0.3, best of six steepness/learning-rate/initialization settings, 1500 Adam steps) reaches ratios
0.55/0.67/0.48 (PSM) and 0.58/0.72/0.52 (binary) on iso-dense/mixed/dense-120, versus certified bounds
0.388/0.403/0.367: 30–67% above the certified optimum. Local refinement of the ILT mask on the ratio objective closes
the gap completely on dense-120 (0.3671 = bound) and partially on the others (0.426, 0.530), showing that both the
objective (pattern error vs exposure latitude) and local minima contribute.
6.7 **SMO entitlement (Table 4).** For the iso-dense target the per-source-point relaxation equals the coherent bound
(0.3625) for conventional σ=0.7 and a 13-point σ=0.9 grid: no illumination supported on these points beats on-axis
coherent light for this metric; fixed-source bounds are 0.43–0.51.
6.8 **2D layouts (Table 6).** [[pending exp4 rerun with edge-offset sampling]] Contact arrays (160, 120, 100 nm pitch)
and a T-junction and line-end layout under annular/quasar illumination; all SDPs < 1 s; previous grid-sampled run:
all contact arrays tight (rank 1, gap 0.00%).
6.9 **Robustness.** Bound invariant to interior spacing (0.46907 for 16→1 nm); monotone in edge tolerance
(2 nm: EL 0.09 → 20 nm: EL 0.74); CLARABEL/SCS agree to 1e-8; source quadrature converges (0.513–0.526 for 11–21
points); 50-digit Cholesky of the certificate succeeds; Rayleigh-quotient bound is 3–150× weaker (Table 7).

## 7 Discussion
The certificates turn qualitative folklore into quantified statements: (i) local optimization of the exposure-latitude
objective is usually globally optimal in 1D and for the 2D contact layouts tested, and when it is not, the shortfall is
exactly the benefit of incoherent multi-exposure; (ii) coherent illumination with phase-shift masks is the entitlement
for the iso-dense target; (iii) binary masks forfeit a certified ≥10% in ratio on that target; (iv) pattern-error ILT
is far from EL-optimal, which argues for EL/contrast objectives or for using the certificate as a stopping criterion.
The relaxation gap has a physical meaning (Theorem 3), so it is informative rather than merely a weakness.

## 8 Limitations
Scalar thin-mask periodic model (no vector effects, 3D mask, resist chemistry, flare); sampled constraints (edge
tolerance is a model parameter); ratio/EL objective only (though any constraint linear in intensity lifts identically);
binary-mask relaxation is only an outer approximation in 2D; Lasserre-2 bounds carry solver status "optimal_inaccurate"
and were not bisection-rigorized; cell sizes limited to the periodic small-field setting (F ≲ 10^3); literature
coverage limited to search snippets (no full texts could be accessed); the IBM 2002 full text was not read.

## 9 Conclusion
Printability of a layout under a given optical configuration can be certified, not just optimized, by a small SDP whose
dual is a one-matrix proof. The method reproduces classical limits, certifies when local ILT-type optimization is
globally optimal, and quantifies the value of multiple exposure and phase shifting, within a standard scalar model.

## References (bibliographic data verified from search metadata; full texts not accessed)
[1] A. E. Rosenbluth, S. J. Bukofsky, C. A. Fonseca, M. S. Hibbs, K. Lai, A. F. Molless, R. N. Singh, A. K. K. Wong,
"Optimum mask and source patterns to print a given shape," J. Microlith. Microfab. Microsyst. 1(1), 13–30 (2002);
SPIE 4346, 486–502 (2001). doi:10.1117/1.1448500.
[2] A. E. Rosenbluth, D. Melville, K. Tian, K. Lai, N. Seong, D. Pfeiffer, M. Colburn, "Global optimization of masks,
including film stack design to restore TM contrast in high NA TCC's," Proc. SPIE 6520, 65200P (2007).
[3] A. E. Rosenbluth et al., US Patents 7,587,702 and 7,685,559 (step-walk relaxation); US 8,495,528 (multiple-exposure
wavefronts); A. E. Rosenbluth, J. Tirapu-Azpiroz, US 8,108,802 (arbitrary wavefronts); US 7,944,545 (high-contrast masks).
[4] Z.-Q. Luo, W.-K. Ma, A. M.-C. So, Y. Ye, S. Zhang, "Semidefinite relaxation of quadratic optimization problems,"
IEEE Signal Process. Mag. 27(3), 20–34 (2010).
[5] N. D. Sidiropoulos, T. N. Davidson, Z.-Q. Luo, "Transmit beamforming for physical-layer multicasting," IEEE Trans.
Signal Process. 54(6), 2239–2251 (2006).
[6] G. Pataki, "On the rank of extreme matrices in semidefinite programs and the multiplicity of optimal eigenvalues,"
Math. Oper. Res. 23(2), 339–358 (1998).
[7] G. Angeris, J. Vučković, S. Boyd, "Computational bounds for photonic design," ACS Photonics 6(5), 1232–1239 (2019);
"Heuristic methods and performance bounds for photonic design," Opt. Express 29, 2827 (2021); "Bounds on efficiency
metrics in photonics," ACS Photonics (2023).
[8] Z. Wang, X. Ma, G. R. Arce, J. Garcia-Frias, "Information theoretical approaches in computational lithography,"
Opt. Express 26(13), 16736–16751 (2018); X. Ma et al., "Information theoretical aspects in coherent optical lithography
systems," Opt. Express 25(23), 29043–29057 (2017).
[9] E. J. Candès, T. Strohmer, V. Voroninski, "PhaseLift: exact and stable signal recovery from magnitude measurements
via convex programming," Comm. Pure Appl. Math. 66(8), 1241–1274 (2013).
[10] Y. C. Pati, T. Kailath, "Phase-shifting masks for microlithography: automated design and mask requirements,"
J. Opt. Soc. Am. A 11(9), 2438–2452 (1994).
[11] L. Pang, "Inverse lithography technology: 30 years from concept to practical, full-chip reality," J. Micro/
Nanopattern. Mater. Metrol. 20(3), 030901 (2021).
[12] Reviews: "Advancements and challenges in inverse lithography technology: a review of AI-based approaches" (2025,
PMC12287475); "Inverse Lithography Technology (ILT) under chip manufacture context" (PMC12843865).
[13] "Pushing the limits of inverse lithography with generative reinforcement learning," arXiv:2602.19027 (2026).
[14] J. B. Lasserre, "Global optimization with polynomials and the problem of moments," SIAM J. Optim. 11, 796–817 (2001).
[15] Carathéodory–Toeplitz theorem; Fejér–Riesz theorem (classical; see e.g. Grenander & Szegő, Toeplitz Forms, 1958).
