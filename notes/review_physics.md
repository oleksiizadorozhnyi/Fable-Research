# Physics / practical-relevance referee report — paper/main.tex (v2, 2026-10-12)

Reviewer scope: Abbe/Hopkins implementation, source sampling, vector model, process-window metric, entitlement
claims, EUV, multi-exposure semantics, adversarial cases, prior art. Files inspected: paper/main.tex, notes/THEORY.md,
code/certlitho.py (Model, ModelX, source helpers), exp5/8/16/17, results/RESULTS_TABLES.md. Numerical checks were run
with the authors' own solver stack (script reproduced in Appendix A; no repository files modified).

Severity key: CRITICAL = a headline number or claim is wrong as stated; MAJOR = must be fixed/reworded before
publication; MINOR = clarify; OK = checked, no action.

---------------------------------------------------------------------------------------------------------------------
## 1. Abbe implementation and source sampling

**F1 (OK) Coherent-channel construction.** `Model.__init__`: g = f_k + s/λ, pupil 1{|g|λ ≤ NA}, order set
|f_k| ≤ (NA+max|s|)/λ, image = Σ_j w_j |Σ_k P(f_k+s_j/λ) M_k e^{2πi f_k x}|². This is the standard Abbe (stacked
pupil-shift) model; orders in K that no *discrete* source point collects produce zero columns — this is the
"common null space of dimension 4" the paper verifies on the complement. Say so explicitly (it is a sampling artefact,
not a property of the layout): "orders inside the continuous collection radius that none of the sampled source points
collects".

**F2 (CRITICAL) The "conventional σ" source is a slit, not a disk, and the headline bounds are not converged in the
quadrature.** `conventional_source_1d(σ,NA,n)` places n uniformly weighted points on a *line* |s_x| ≤ σNA, s_y = 0.
A disk illuminator seen by a 1D (y-invariant) mask has (i) weights ∝ chord length √(σ²NA² − s_x²) (zero at the rim,
whereas the paper's quadrature gives the rim points weight 1/n = 14 %), and (ii) a y-component s_y/λ that enters the
pupil cutoff and the defocus phase. The paper's own sweep (RESULTS_TABLES §Robustness: 3/5/7/11/15/21 pts →
0.413/0.473/0.469/0.516/0.515/0.526) already shows the 7-point value used in every table is 11 % below the 21-point
value and that the sequence has not converged; the text ("source quadrature converges (0.513–0.526 for 11–21 points)")
hides that all reported σ-source numbers use n = 7. Disk-faithful recomputation (Appendix A, same code path, 2D source
points on a 1D mask, uniform area weights):

| case | paper (7-pt line) | 21-pt line | disk 9×9 | disk 15×15 | disk 21×21/31×31 |
|---|---|---|---|---|---|
| DUV iso-dense σ=0.7, z=0, bound ρ | 0.4691 (EL 0.361) | 0.5261 | 0.5776 | 0.5741 | 0.5889 (EL 0.259) |
| EUV p22 conv0.9 ±50 nm, bound ρ | 0.4685 (EL 0.362) | 0.6543 | — | 0.8545 | 0.9082 (EL 0.048) |
| EUV p24 conv0.9 | 0.4733 (EL 0.357) | | | | 0.6827 (EL 0.189) |
| EUV p26 conv0.9 | 0.4123 (EL 0.416) | | | | 0.5758 (EL 0.269) |
| EUV p28 conv0.9 | 0.3951 (EL 0.434) | | | | 0.5045 (EL 0.329) |
| EUV p24 dipole 0.8: 2 points vs finite poles (r 0.7–0.9, ±35°) | 0.2007 (EL 0.666) | | | | 0.3590 (EL 0.472) |

Consequences: the abstract's "σ=0.9 conventional source reaches EL=0.36 at 22 nm pitch" is an artefact (a real disk
gives ≈0.05 because only a sliver of the disk near |s_x| ≈ 0.86 NA produces two-beam interference at k1 = 0.27); the
Fig. 2 EUV curve, Table 1 σ-rows, Table 4 (SMO), the DNN/ILT comparisons, the aberration/vector sweeps and the
random-target tightness statistics for σ-sources all inherit a 10–25 % model error in ρ. The *certificates* remain valid
for the slit illuminator actually modelled; the *physics labels* ("conventional σ = 0.7", "σ = 0.9") are wrong.
Required: (a) implement 2D source points for 1D masks (two-line change: g = (f_k + s_x/λ, s_y/λ)); (b) use a
quadrature with ≥ 150–300 points (or chord-weighted Gauss rule) and report convergence of ρ to < 1 %; (c) re-run all
σ-source tables; (d) replace "converges" by the actual convergence statement. The 2D helpers have the same problem in
milder form: `annular_source_2d(0.6,0.9,n=9)` samples a 0.3 NA-wide ring with a 0.225 NA grid (1–2 rings of points);
`quasar_source_2d` defines the poles by |‖a|−|b‖ < 0.35 σ_out NA (a band about the diagonals) rather than an opening
angle; `dipole_source_1d` is two δ-points (unrealisable; pole size matters, see EUV row above). State pole geometry
(σ_in/σ_out/opening angle) and point counts in every table caption.

**F3 (MAJOR) Defocus phase is computed in air (n = 1) and clipped at |g|λ = 1, which is unphysical for immersion
NA = 1.35.** `Model`: `gl = clip(|g|λ, 0, 1); ph = 2π z (√(1−gl²) − 1)/λ`. For NA = 1.35 all orders with
1 < |g|λ ≤ 1.35 (a large part of the pupil, e.g. k = 4 of the iso-dense cell on axis) receive a constant phase as
if propagating at 90°, and orders with |g|λ < 1 receive the vacuum defocus, ≈ n× too strong. Correct form:
φ = (2π z/λ)(√(n² − |g|²λ²) − n) with n the index of the medium in which z is measured (water 1.44, or resist 1.7 if
z is depth in resist). Recomputed PW80 iso-dense σ=0.7 bound: 0.8131 (paper) → 0.7352 (n = 1.44) → 0.6138 (n = 1.7)
with the paper's source; 0.7243 / 0.6677 with a disk source. All three-focus DUV rows of Table 1, exp9/exp11 PW cases,
and the "three-focus cases 11–16 % gap" statement must be recomputed; ranks/gaps may change. `ModelX` has the same
clip while simultaneously using n_img = 1.7 for the direction cosines — internally inconsistent. EUV (NA 0.33, vacuum,
n_resist ≈ 0.97) is unaffected.

**F4 (OK/MINOR) Cutoff recovery is a consistency check, not a prediction.** Below pitch λ/((1+σ)NA) the order set K
is {0} by construction, so t_R = 1 trivially. Reword "The bound reproduces the classical k1 cutoffs exactly" → "The
bound is consistent with the k1 cutoff built into the collected-order set (Prop. 5)". Add the closed-form two-beam
check which the code does reproduce: for dense lines with two collected beams and edge tolerance τ,
ρ_2beam = (1 − sin(2πτ/p))/(1 + sin(2πτ/p)), EL = sin(2πτ/p); dense-90, τ = 6 nm gives 0.42173 = the reported
0.4217 (coherent and dipole). This is the classic entitlement of dipole/alt-PSM imaging and is a stronger sanity
check than the cutoff.

---------------------------------------------------------------------------------------------------------------------
## 2. Vector model (ModelX)

**F5 (OK) Rodrigues rotation.** With k̂ = (α,β,γ) and v = ẑ×k̂, R x̂ = (1 − α²/(1+γ), −αβ/(1+γ), −α) is correct
(Mansuripur 1986; Flagello–Milster–Rosenbluth 1996). TE/TM assignment (azimuthal/radial per source point) and
|E_x|²+|E_y|²+|E_z|² in the resist are standard. The claim that TM loses contrast at NA 1.35 is physically right: for
the dense-90/σ=0.7 case the interfering beams are at ≈ +45° and −34° in the resist, so the TM interference term carries
cos(79°) ≈ 0.19 — ratios 70–120 % above scalar are the expected magnitude.

**F6 (MAJOR, wording) What the vector model omits, and what it does to the claims.**
(i) Fresnel transmission at the fluid/resist interface (TE ≠ TM; Brewster at tan θ = 1.7/1.44 → 49.7° in water, i.e.
inside the pupil) and the energy/obliquity factor across the interface are omitted; the 1/√γ aplanatic factor is applied
with γ in the *resist* (n_img = 1.7) whereas it belongs to the lens exit medium (fluid, 1.44), followed by transmission.
(ii) No thin-film/standing-wave/BARC treatment (homogeneous resist, no substrate reflection). (iii) `n_img` is a single
number; immersion fluid vs resist are not distinguished. (iv) Defocus clip (F3). None of these change the qualitative
statements (TM contrast loss; certificate validity; tightness statistics), because every omitted factor is a per-channel
complex amplitude that keeps the image a sum of rank-one Hermitian forms (Lemma 7 still holds). They do change TM/unpol
numbers by O(10 %). Recommended wording in §Robustness: "The vector model includes the polarization rotation and
longitudinal field in a homogeneous image medium of index 1.7 and the aplanatic obliquity factor; it omits Fresnel
transmission into the resist, thin-film interference and substrate reflection. These are per-channel amplitude factors
that preserve the quadratic structure (Lemma 7) and alter TM/unpolarized bounds quantitatively, not the certificate
or its tightness." Also rename bib key `mack_vector` (it is Flagello et al.) and add Brunner et al., "High-NA lithographic
imagery at Brewster's angle", JM3 1(3) 2002 for TM contrast loss.

**F7 (MINOR) On-axis TE/TM is undefined.** For s = 0 the code silently chooses the radial direction x̂; for 1D
lines along y this makes "TM" on axis = x-polarized, which is fine for this geometry but should be stated.

---------------------------------------------------------------------------------------------------------------------
## 3. The dark/bright ratio and EL = (1−ρ)/(1+ρ)

**F8 (MAJOR) EL is misdefined by a factor 2 and quoted at a non-standard CD tolerance.** For a threshold T and dose E:
print window E ∈ [T/I_b, T/I_d), relative width about its midpoint = 2(I_b − I_d)/(I_b + I_d) = 2(1−ρ)/(1+ρ). The paper's
EL = (1−ρ)/(1+ρ) = (I_b − I_d)/(I_b + I_d) is the Michelson contrast of the image sampled at the ±τ offsets (half-range,
"±EL"). Moreover τ = 8 nm on a 50 nm feature is a ±16 nm (±32 %) CD tolerance; the industry convention is EL at ±10 % CD,
i.e. τ = 0.05·CD = 2.5 nm, where the paper's own sweep gives EL_paper ≈ 0.11 (total ≈ 22 %). Headline EL values
(0.36–0.47, "EL = 0.36 at 22 nm pitch") are therefore not comparable to any tool or OPC-report EL. Required: define
EL_total = 2(1−ρ)/(1+ρ) (or call the quantity "image contrast at tolerance τ"), state "τ is the per-edge placement
tolerance, so the CD tolerance is ±2τ for a symmetric feature", and report headline numbers at τ = 5 % of the minimum CD
in addition to τ = 8 nm. Relation to NILS worth stating: if the image is linear across ±τ, slope = (I_b−I_d)/(2τ) and
NILS = CD·slope/I_mid = (CD/τ)·(1−ρ)/(1+ρ); the certificate is thus a lower bound on a finite-difference NILS at scale τ.

**F9 (OK) Is the metric standard?** Yes: maximizing min bright / max dark intensity with sampled bright/dark points is
exactly the objective of Rosenbluth et al. 2002/2007 and of IBM's "high-contrast masks" patents; say "we adopt Rosenbluth's
bright/dark constraint formulation" rather than presenting the ratio as new. It is *not* the PV-band/EPE metric fabs
report; the paper says so (good).

**F10 (MINOR, correct but needs care) "Any constraint linear in intensity lifts identically."** Correct for: multiple
foci, multiple dose corners, image slope at an edge with known orientation (∂_x Tr(Q_x X) = Tr(∂_x Q_x X) is linear),
dose-to-size equalities Tr(Q_edge X) = T, and ILS/NILS ≥ c as the homogeneous inequality Tr(Q'_x X) − (c/CD)·Tr(Q_x X) ≥ 0.
Not linear: two-sided |slope| constraints, CD measured by threshold crossing (nonconvex), MEEF (needs a derivative with
respect to the mask, i.e. the X variable itself), and any Boolean combination of exposures (F15). Reword: "constraints
affine in the intensity (multi-focus, multi-dose, signed slope/ILS at an edge, dose-to-size) lift to linear constraints in
X and admit certificates of the same form; EPE enters only through τ; MEEF and Boolean multi-patterning do not lift."

**F11 (MAJOR) Scale invariance hides dose/throughput and is the root of the 20 nm-space degeneracy (F17).** ρ is
invariant to M → αM, so the bound never penalises masks whose bright intensity is a tiny fraction of open frame (extreme
dose, extreme MEEF, phase-edge solutions). A dose-aware certificate is one line away and still convex: with |m| ≤ 1,
Parseval gives Σ|M_k|² ≤ 1, i.e. Tr X ≤ 1; adding min_B Tr(Q_b X) ≥ κ (κ = required bright intensity relative to open
frame, e.g. 0.3) keeps (R) an SDP and Theorem 1 acquires one extra dual term. Recommend adding this as a variant and
reporting at least the iso-dense and 20 nm-space cases with κ ∈ {0.1, 0.3}.

---------------------------------------------------------------------------------------------------------------------
## 4. "Coherent + strong PSM is the entitlement" (iso-dense)

**F12 (MAJOR, misleading as stated).** Table 4 shows the *per-source-point unlimited-exposure* relaxation equals the
on-axis coherent relaxation (0.3625). But Table 1 shows the best *single* coherent mask is 0.4741, worse than the single
dipole mask (0.4313) and essentially equal to σ=0.7 (0.4808). So for one exposure the dipole is the entitlement, and the
coherent value is only reachable with ≥ 2 incoherently summed exposures (rank-2 face). Reword Discussion/§Source–mask:
"For the iso-dense target, on-axis coherent illumination with a real (phase-shifting) mask achieves the lowest
*multi-exposure* relaxation value (0.3625), but no single coherent mask does better than 0.474; for a single exposure
the 0.8 dipole (0.431, tight) is the best illumination tested." Physical plausibility: yes — on-axis coherent light with
a real band-limited mask is the two/three-beam alt-PSM/chromeless regime, and partial coherence washes out phase-edge
contrast; this is textbook (Levenson 1982; alt-PSM wants low σ). Caveats to state: σ → 0 is unrealisable (practical
σ ≥ 0.2–0.3, throughput, laser coherence/speckle); the optimal mask is a *gray-tone* real transmission in [−1,1] (not a
0/π chrome alt-PSM), which raises manufacturability, phase-conflict (2D) and MEEF questions; alt-PSM intensity imbalance
and best-focus shift (mask 3D) are absent in the thin-mask model; the ratio ignores dose (F11). Also note Lemma 6 only
says complex transmission does not help the *relaxation*; a single complex mask can still beat a single real mask, so the
"PSM" in the headline means "real transmission of either sign".

**F13 (MINOR) Att-PSM is misclassified.** §3 writes "nonnegative masks 𝓜₊ ⊂ 𝓜_ℝ (binary/attenuated)". A 6 % att-PSM
has transmission ≈ −0.245 (negative real), so it is *not* in 𝓜₊; it lies in 𝓜_ℝ with a lower bound on m. Replace
"(binary/attenuated)" by "(binary, gray-tone)" and say the PSM class is "arbitrary real transmission in [−1,1] (strong
PSM / chromeless); attenuated PSM is intermediate and not treated". The "binary masks forfeit 3–14 %" statement should
read "relative to unconstrained real-transmission masks".

---------------------------------------------------------------------------------------------------------------------
## 5. EUV (λ 13.5 nm, NA 0.33, σ 0.9 / dipole 0.8, ±50 nm)

**F14 (CRITICAL for the EUV numbers; MAJOR for wording).** Beyond F2 (the 22 nm "EL = 0.36" collapses to ≈ 0.05 with a
disk; the dipole numbers drop ~45 % with finite poles), the following make every certified EUV EL optimistic and must be
listed where the numbers appear, not only in Limitations: (i) mask shadowing from the 6° chief-ray angle (reflective
mask, ~60–70 nm absorber): H/V bias of several nm and orientation-dependent best focus; (ii) mask-3D phase effects:
pitch-dependent best-focus shift, pattern shift through focus (non-telecentricity), contrast loss — the thin-mask
spectrum is not incidence-angle independent; (iii) flare (several % DC background): ρ_eff = (I_d + φ)/(I_b + φ) > ρ; the
bound remains valid as a lower bound, but the EL numbers do not; (iv) resist blur (~3–5 nm) and image log-slope loss at
10–11 nm half-pitch; (v) stochastics: at 20–24 nm pitch photon shot noise and resist stochastics, not aerial-image EL,
set the failure rate, so a threshold-resist EL is not a process-window statement; (vi) realisable sources: NXE
illuminators require a minimum pupil-fill ratio (~20 %), so a two-point dipole and the sliver-of-disk two-beam regime are
both unrealisable; (vii) 4× anamorphic/high-NA (0.55) is the relevant regime for 20–24 nm pitch, not 0.33; (viii) EUV
"masks" are reflective att-PSM-like structures (absorber reflectance a few %, phase ≈ 180–200°), so the 𝓜_ℝ class is an
upper envelope. Note a genuine strength the authors miss: the per-source-point relaxation `sdp_bound_smo` is valid even
for angle-dependent (thick-mask) diffraction because it gives every incidence angle its own order vector; the fixed-source
relaxation is not. Recommended sentence for §Cutoff recovery / Table 5: "EUV values are aerial-image entitlements of the
thin-mask scalar model with an idealised illuminator; shadowing, mask-3D best-focus shifts, flare, resist blur and
stochastic failure all reduce the achievable window, so the certified EL is a necessary, not sufficient, condition for
printability." Remove "reaches EL = 0.36" from the abstract until F2 is fixed.

---------------------------------------------------------------------------------------------------------------------
## 6. Multi-exposure semantics

**F15 (MAJOR, wording) Theorem 3 models one specific process.** "Incoherent superposition of aerial images with a common
threshold" is the correct model of *double exposure of a single resist layer without intermediate processing* (latent
image ∝ summed dose; valid for a linear-absorbing CAR to first order, same resist, dose split free, no inter-exposure
alignment error). It does *not* model litho-etch-litho-etch / litho-freeze-litho-etch (each exposure developed
separately; the final pattern is a Boolean union/intersection, which can beat the incoherent sum — e.g. pitch splitting
below the cutoff), nor nonlinear/threshold-per-exposure processes (two-photon, CEL, dual-tone). Add after Theorem 3:
"Here 'multiple exposure' means incoherent dose summation in one resist (double exposure, or equivalently a source
decomposed into coherent exposures, each with its own mask); multi-patterning with intermediate development or etch
(LELE, LFLE) composes images nonlinearly and is not bounded by t_R." Also state that the SMO relaxation (one PSD block per
source point) is the same object Rosenbluth et al. 2002 call per-direction optimal wavefronts and that its value is a
bound on *any* source supported on the sampled points *with unlimited exposures*, not on single-exposure SMO.

---------------------------------------------------------------------------------------------------------------------
## 7. Isolated 20 nm space with bound ≈ 0

**F16 (OK as mathematics; MAJOR as interpretation).** `extreme_duty_95pct` = 380 nm bright, 20 nm dark line, τ = 6 nm:
dark samples sit at ±4 nm about the line centre, bright samples at ±16 nm. A real mask with a sign change (phase edge)
produces a field zero at the edge, so I(0) = 0 and I(±4 nm) = O(slope²·16 nm²) while I(±16 nm) is large; ρ → 0 is correct
*for this constraint set*. The physics it reports is the chromeless/phase-edge line, whose printed width is set by dose,
not by the mask — exactly what a metric without dose-to-size or CD constraint cannot see: the dark line is allowed to
print anywhere between 8 nm and 32 nm wide (τ = 6 nm is a ±60 % CD tolerance here). Reword "a phase-edge mask prints an
arbitrarily dark narrow line, which the relaxation recognises" → "a phase-edge mask yields a field zero, so the sampled
constraints (dark within ±4 nm, bright beyond ±16 nm) are satisfiable with arbitrarily small dark intensity; the ratio
metric does not fix the printed CD of such a line (its width is dose-determined), so this bound is a statement about the
constraint set, not a printability entitlement. With a dose-to-size constraint (F11) or τ ≤ 10 % of CD the degeneracy
disappears." Also clean Table 12: gaps "−100.00 %" and "1592.09 %" and ranks 3/7 on bounds of 0.0000 are 0/0 artefacts;
report "degenerate (bound 0)" and exclude from tightness statistics.

**F17 (MINOR) Same mechanism in `sraf_like_narrow_space`** (20 nm dark lines flanking a 20 nm bright space; ρ = 0.2426
coherent/dipole, tight): say that the solution is a phase-edge pair and give the printed-CD caveat.

---------------------------------------------------------------------------------------------------------------------
## 8. Prior art the manuscript should position against (verify full texts; items marked † are from the reviewer's
recollection and must be checked against the originals)

**F18 (MAJOR).**
- Rosenbluth et al., JM3 1(1) 13–30 (2002) "Optimum mask and source patterns to print a given shape": sampled bright/dark
  constraints, minimum-bright/maximum-dark objective, per-source-direction optimal wavefronts, and — † — the use of
  pairwise products of diffraction orders ("interference terms") as the optimization variables with the rank-one
  condition relaxed; if confirmed, the lifted variable X is already there and the paper's novelty is the *dual
  certificate, rigorization and the exact multi-exposure identification*, which should be said precisely. Obtaining the
  full text is mandatory; "full texts were not accessible" is not acceptable for the central prior-art item.
- Rosenbluth & Seong, Proc. SPIE 6154 (2006) "Global optimization of the illumination distribution to maximize integrated
  process window": LP over the source for fixed mask (the source step is convex); and the IBM patents already cited.
- Granik, JM3 3(4) 2004 "Source optimization for image fidelity and throughput" and JM3 5(4) 2006 "Fast pixel-based mask
  optimization for inverse lithography"; Poonawala & Milanfar, IEEE TIP 16(3) 2007 — standard ILT baselines; the "pattern
  error ILT" strawman in §PSM vs binary should cite one of them and use their threshold/normalization conventions.
- Socha et al., Proc. SPIE 5377 (2004) "Contact hole reticle optimization by using interference mapping lithography" and
  Socha, Proc. SPIE 5853 (2005) "Simultaneous source mask optimization": coherent "ideal image"/interference-map
  entitlement, the same object as the coherent relaxation here.
- Yamazoe, JOSA A 25(12) 2008 "Computation theory of partially coherent imaging by stacked pupil shift matrix": the Abbe
  channel matrix C (field = C·M) and its relation to the Hopkins/TCC eigen-decomposition (Pati–Kailath 1994 cited; add
  Cobb 1998 SOCS). The paper's "rank ≤ #source points" remark is Yamazoe's.
- Levenson, Viswanathan, Simpson, IEEE TED 29 (1982) alt-PSM; Brunner et al., JM3 2002 (Brewster/TM); Mansuripur, JOSA A 3
  (1986) (vector rotation).
- Two-beam/NILS entitlement: Mack, *Fundamental Principles of Optical Lithography* (2007) ch. 2/8 (two-beam image, NILS
  ≤ π·contrast, ILS/EL relations) — supports F4/F8; Brunner, JVST B 21 (2003) "Why optical lithography will live forever"
  for k1 limits.
- Lifting/SDP in photonics is cited; add Angeris et al. is fine. Also note Liu–Zakhor 1992 already cited.
- EUV: Krautschik et al. (SPIE 2001, shadowing), Erdmann et al. JM3 (mask-induced best-focus shifts), De Bisschop JM3
  16(4) 2017 (stochastics), for F14.

---------------------------------------------------------------------------------------------------------------------
## 9. Other physics-adjacent items

**F19 (MINOR) Direction-cosine bookkeeping.** "|s_j| ≤ σNA (direction cosines)" exceeds 1 for σ = 0.9, NA = 1.35; call
them normalized spatial frequencies ×λ (sine of the angle in the image medium, NA = n sin θ) and state the medium.

**F20 (MINOR) Vector-model polarization ensemble weights.** `ws += [wsrc/len(pols)]*nch` — unpolarized = ½(x + y), fine;
but TE/TM for an *annular/quasar* source changes direction per point, so Lemma 6's symmetry argument needs w(−s) = w(s)
*and* the polarization map to be odd/even consistently — true for azimuthal/radial fields; one sentence suffices.

**F21 (OK) Rayleigh-quotient comparison and complex=real check**: fine; Lemma 6 requires an even pupil, which excludes
odd aberrations (coma) — say "for aberration-free or even-aberrated pupils".

---------------------------------------------------------------------------------------------------------------------
## Summary verdict

Mathematics (certificate, Theorem 3, Theorem 9) is sound and the code implements the stated model faithfully. The
*physics labels* are not: the "conventional σ" illuminator is a 7-point slit (F2), immersion defocus is computed in air
with a clip (F3), and EL is half the standard quantity at a non-standard CD tolerance (F8). These three items change every
σ-source and every through-focus number in the paper by 10–25 % (and the EUV 22 nm headline from EL 0.36 to ≈ 0.05), so
the empirical sections must be re-run before any claim about "conventional illumination", EUV process windows or
three-focus gaps can stand. The qualitative conclusions most likely to survive: tightness for dipole/dense lines (the
two-beam closed form), the certified multi-exposure gains under coherent light (source-independent), the vector-model
TM/tightness statements, and the structural theorems.

---------------------------------------------------------------------------------------------------------------------
## Appendix A — reviewer check script (run with the repository's certlitho.py; CLARABEL)

```python
import sys, numpy as np; sys.path.insert(0, '/home/user/Fable-Research/code')
from certlitho import *
class Model1D2S(Model):   # 1D mask, 2D source points (sx, sy); defocus in medium of index n_med
    def __init__(self, L, lam, NA, src2, w, xs, focus=(0.0,), n_med=1.0):
        self.dim = 1; self.L = np.array([float(L)]); self.lam, self.NA = lam, NA
        self.src = np.asarray(src2, float); self.w = np.asarray(w, float)/np.sum(w); self.focus = tuple(focus)
        smax = np.max(np.linalg.norm(self.src, axis=1)); R = (NA + smax)/lam
        K = int(np.floor(R*L + 1e-9)); self.k = np.array([(k,) for k in range(-K, K+1)]); self.F = len(self.k)
        self.f = self.k/self.L; self.kindex = {tuple(kk): i for i, kk in enumerate(self.k)}
        self.x = np.asarray(xs, float).reshape(-1, 1); self.nx = len(self.x); self.C = {}
        for z in self.focus:
            Cz = []
            for s in self.src:
                gx = self.f[:, 0] + s[0]/lam; gy = np.full(self.F, s[1]/lam); gn = np.hypot(gx, gy)
                inside = gn*lam <= NA + 1e-12
                ph = 2*np.pi/lam*z*(np.sqrt(np.clip(n_med**2 - (gn*lam)**2, 0, None)) - n_med)
                Cz.append(np.exp(2j*np.pi*(self.x @ self.f.T))*(inside*np.exp(1j*ph))[None, :])
            self.C[z] = np.stack(Cz)
        self.flip = np.array([self.kindex[tuple(-kk)] for kk in self.k])
def disk_pts(sigma, NA, n):
    g = np.linspace(-sigma*NA, sigma*NA, n); pts = [(a, b) for a in g for b in g if np.hypot(a, b) <= sigma*NA + 1e-12]
    return np.array(pts), np.ones(len(pts))
def dipole_pts(sc, NA, n, halfw=0.1, opening_deg=35):
    g = np.linspace(-NA, NA, n); pts = [(a, b) for a in g for b in g
          if sc-halfw <= np.hypot(a, b)/NA <= sc+halfw and np.degrees(np.arctan2(abs(b), abs(a))) <= opening_deg]
    return np.array(pts), np.ones(len(pts))
def bound(L, lam, NA, pts, w, segs, margin, focus=(0.0,), n_med=1.0, spacing=4.0):
    x, B, D, _ = samples_1d(L, segs, margin, spacing); rm = RealModel(Model1D2S(L, lam, NA, pts, w, x, focus, n_med))
    r = sdp_bound_real(rm, B, D); return rigorize_certificate(rm, r['mu'], r['lam'])[0]
# DUV iso-dense sigma 0.7: paper line7 -> 0.4691; disk 21x21 -> 0.5889.  PW80: n_med 1.0/1.44/1.7 -> 0.8131/0.7352/0.6138 (line7)
# EUV p22 conv0.9 +-50: line7 -> 0.4685; disk 31x31 -> 0.9082.  EUV p24 dipole: 2pt -> 0.2007; finite poles -> 0.3590
# dense-90 two-beam closed form: (1-sin(2*pi*6/90))/(1+sin(2*pi*6/90)) = 0.42173 (paper: 0.4217)
```
