# Mathematical review of notes/THEORY.md and paper/main.tex §3–§4

Scope: Theorems 1, 3, 4, 9; Lemmas 2, 6, 7, 8; Proposition 5; sampling remark; §3 definitions (τ, EL, pupil).
Code consulted where a proof leans on it: `code/certlitho.py` (Model, RealModel, ModelX, samples_1d/2d, sdp_bound_pixel*), `code/bnb_cert.py`.
Labels: CRITICAL = claim false as used by a headline result; MAJOR = statement/proof wrong or materially incomplete; MINOR = imprecision; OK = checked, no change needed.

---

## Theorem 1 (certificate, weak duality)

1. **OK** — Main chain. δ = δΣμ_d ≥ Σμ_d I_d (needs μ ≥ 0, I_d ≤ δ) = M^H S M + Σλ_b M^H Q_b M ≥ βΣλ_b (needs S ⪰ 0, λ ≥ 0, I_b ≥ β). Division by β needs β > 0. Correct. Note the proof never uses Q ⪰ 0 or Hermitian-ness beyond M^H Q_b M = I_b real.

2. **OK** — Real-mask variant. For M = Ur, r ∈ ℝ^F: M^H S M = r^T(U^H S U)r = r^T Re(U^H S U) r because r^T(Im U^H S U)r = 0 (antisymmetric). S̃ ⪰ 0 suffices and bounds only M ∈ 𝓜_ℝ; the paper states this correctly.

3. **OK** — Multi-exposure variant. With J = Σ_i I_{M_i}: max_D J ≥ Σ_d μ_d J_d = Σ_i[M_i^H S M_i + Σ_b λ_b M_i^H Q_b M_i] ≥ Σ_b λ_b J_b ≥ (min_B J)Σλ_b. Doses are absorbed in ‖M_i‖. Hidden assumption worth one sentence: every exposure uses the same source, pupil and focus set Z (per-exposure sources are the separate SMO variant).

4. **MINOR** — ρ is 0/0 on the common null space. The paper proves ρ(M) ≥ Σλ_b "for every M ∈ ℂ^F∖{0}", but masks supported on orders invisible to every constraint exist (the paper itself reports a 4-dimensional common null space at F = 61), and for them min_B I = max_D I = 0.
   Fix (§3, after the definition of ρ): "with the convention ρ(M) = +∞ when min_{B×Z} I_{M,z}(x) = 0".

5. **MINOR** — Rigorization (§5) uses a *near*-null space in floating point. The argument "on N every form vanishes identically, hence PSD on N^⊥ suffices" is exact only when S_D v = S_B v = 0 exactly; with eigenvalues ~1e-12 discarded, the cross terms v_⊥^H S v_N are not exactly zero and are not covered by the 50-digit Cholesky of the projected matrix. State the residual bound (‖S v_N‖ ≤ null_tol·‖S‖) and that it is absorbed by the 1e-10 relative margin, or project with the exact null space of the floating-point matrices.

## Lemma 2 (real structure)

6. **OK** — J antilinear involution (Π real, Π² = I); Fix(J) has real dimension F; u^H v is real on Fix(J) (conj(u^H v) = u^TΠv = (Πu)^Tv = conj(u)^Tv = u^Hv); XJ = JX ⇔ X = conj(ΠXΠ) (JX M = conj(ΠXΠ)·conj(ΠM)); Y = U^H X U is real symmetric PSD; Tr(QX) = Tr(Q̃Y) because Tr(Im(U^HQU)·Y) = 0 (antisymmetric × symmetric). The code's U (columns (e_i+e_j)/√2, i(e_i−e_j)/√2, e_i for self-conjugate orders) is unitary and spans Fix(J).

## Theorem 3 (relaxation = unlimited multiple exposure; Pataki)

7. **OK** — Equivalence t_R = inf_r ρ*(ME_r) = ρ*(ME_F) as infima: X ⪰ 0 ⇔ X = Σ_{i≤F} M_iM_i^H, Tr(QX) = Σ_i I_{M_i}; Lemma 2 gives real M_i in the real-structured case; rank-one ⇔ single mask both ways.

8. **MAJOR** — Pataki bound is stated for 𝓜_ℂ with the real-symmetric count r(r+1)/2. For a Hermitian PSD block the face dimension is r², so the extreme-point bound is **r² ≤ m** (Pataki; Huang–Palomar/Huang–Zhang for the complex case). r(r+1)/2 ≤ m is correct only for 𝓜_ℝ (via Lemma 2, a real SDP).
   Fix (statement): "An extreme optimal solution has rank r with r(r+1)/2 ≤ |D×Z|+|B×Z| for 𝓜_ℝ and r² ≤ |D×Z|+|B×Z| for 𝓜_ℂ."

9. **MINOR** — Slack and free-variable bookkeeping. With equalities Tr(Q_dX)+s_d = t⁺−t⁻ and Tr(Q_bX)−u_b = 1, Pataki gives r(r+1)/2 + #{s_d>0} + #{u_b>0} + #{t^±>0} ≤ m. The stated bound drops nonnegative terms, so it is valid but loose: since t_R > 0 forces one of t^± > 0, r(r+1)/2 ≤ m − 1 − #(inactive constraints). Splitting t = t⁺−t⁻ is legitimate (the split feasible set is line-free, so extreme points exist). Suggest stating the sharper count or at least "≤ m − 1".

10. **MAJOR (gap)** — Pataki requires an optimal solution to exist, and neither text argues that (R) attains its minimum; the adjacent Remark on non-attainment of (P) invites the reader to doubt it. It does hold, and the argument is short: let N = ∩_d ker Q_d. Since each Q_d ⪰ 0, Q_d P_N = 0, so dark constraints and any bright constraint with Q_b P_N = 0 see only the N^⊥-block of X, on which Σ_d Q_d ≻ 0 bounds Tr(X_⊥⊥) ≤ C·t; bright constraints with Q_b P_N ≠ 0 are met by adding a large multiple of P_N at zero dark cost. Hence the level sets of the reduced problem are compact and the minimum is attained.
   Fix: add one sentence "(R) always attains its optimum: [argument above]; the attainment caveat concerns only the rank-constrained problem (P)."

11. **MINOR** — "(P_𝓜) ... = min{t : ...}" in §3 should read inf (the Remark itself says attainment may fail).

## Theorem 4 (nonnegative masks in 1D)

12. **MAJOR (THEORY.md statement; paper is fine)** — "If T(M) ≻ 0 the measure may be taken to be ... a strictly positive trigonometric polynomial of degree Kmax, by the Fejér–Riesz theorem" is false. T(M) ≻ 0 does not make the in-band polynomial p(x) = Σ_{|k|≤Kmax} M_k e^{2πikx/L} nonnegative: M_0 = 1, M_{±1} = 0.6 gives T = [[1,.6],[.6,1]] ≻ 0 but p = 1+1.2cos(·) < 0 somewhere. Fejér–Riesz is irrelevant here. The proof's actual construction (maximum-entropy AR density σ²/|a(e^{iθ})|²) is a rational, not polynomial, density with infinitely many nonzero Fourier coefficients, of which only M_0..M_Kmax are prescribed.
   Fix: replace the parenthetical by "(a bounded, strictly positive continuous function, e.g. the maximum-entropy autoregressive extension)". Delete the Fejér–Riesz reference.

13. **OK** — AR extension: T ≻ 0 ⇒ Levinson–Durbin prediction-error variance σ² > 0 and a(z) zero-free on the closed disk; f = σ²/|a|² is continuous on the circle, 0 < min f ≤ f ≤ max f < ∞, and matches M_0..M_Kmax (Yule–Walker). Scaling by 1/max f gives 0 < m ≤ 1; legitimate for the ratio problem because ρ is scale-invariant.

14. **MINOR** — Closure/class bookkeeping. (a) T(M) is Hermitian only if M_{−k} = conj(M_k); the hypothesis should read M ∈ 𝓜_ℝ, not M ∈ ℂ^F. (b) If 𝓜_+ means nonnegative *functions* (not measures), then 𝓜_+ = {T ≻ 0} ∪ {0} exactly: a singular T with M ≠ 0 forces ∫ m|p_v|² = 0, i.e. m = 0 a.e. — stronger and cleaner than "int ⊂ 𝓜_+ ⊂ closure". (c) If 𝓜_+ is meant with 0 ≤ m ≤ 1 it is not a cone and "closure(𝓜_+) = {T ⪰ 0}" is false; what is true is cone(𝓜_+) closure. State which.
   Fix: "Let 𝓜_+ be the in-band coefficients of nonnegative integrable m (the bound 0 ≤ m ≤ 1 is immaterial for the scale-invariant ρ). Then 𝓜_+ = {M ∈ 𝓜_ℝ : T(M) ≻ 0} ∪ {0} and closure(𝓜_+) = {T(M) ⪰ 0}."

15. **MINOR** — inf over 𝓜_+ equals inf over {T ⪰ 0}∖{0}: true, but needs the observation that ρ is continuous where min_B I > 0 and +∞ elsewhere, so boundary values are limits of interior values. One clause.

16. **MINOR** — 2D remark: conclusion correct (block-Toeplitz PSD necessary, not sufficient), reason misattributed. The relevant fact is the failure of the 2D truncated trigonometric moment problem (Rudin, "The extension problem for positive-definite functions", 1963; Calderón–Pepinsky), not "no 2D Fejér–Riesz" (which concerns squares of polynomials). Cite Rudin.

17. **MAJOR (use in §6)** — "binary masks forfeit a certified 10%/14%/3%" is derived from `sdp_bound_pixel(..., nonneg=True)` with N = 120 *point-sampled* pixels (no sinc). That DNN value is a lower bound only for nonnegative masks whose in-band spectrum lies in the polyhedral cone generated by the 120 point-mass spectra, i.e. masks on that 5-nm grid; a binary mask with edges off the grid is not covered (the Toeplitz cone of Theorem 4 is not finitely generated). Theorem 4 provides the tool for *all* nonnegative masks (LMI T(M) ⪰ 0, or its lifted analogue) but was not used for this claim.
   Fix: either (i) qualify: "binary/attenuated masks on a 5-nm pixel grid forfeit a certified ...", or (ii) recompute with the Toeplitz-cone relaxation and keep the unqualified claim.

## Proposition 5 (cutoff)

18. **MAJOR** — "i.e./equivalently" is wrong in 2D and the converse is false. (a) "No two collected orders pass for a common source point" ⇔ no order k ≠ 0 passes for any source point ⇔ Q independent of x. The condition L(NA+max|s|)/λ < 1 (⇔ F = 1) implies this, and in 1D with source points at ±s_max it is equivalent to it; in 2D it is only sufficient (an axis order (1,0) is invisible under a quasar source until L ≥ λ/((1+σ/√2)NA); for a rectangular cell the condition must use max(L_1,L_2)). (b) "Above this pitch t_R < 1 for any target with a bright and a dark point" is false. Just above cutoff only orders 0 and ±1 interfere (for the extreme source points), so each focus image is c + A cos(2πx/L + φ_z): a single sinusoid cannot separate, e.g., B = {0, L/2} from D = {L/4, 3L/4} (a period-L/2 target in cell L), and a wide focus spread can also make separation impossible; then t_R = 1. (c) The proof's mechanism "two interfering orders ±1" is wrong near cutoff: ±1 interfere only if both pass for one source point, i.e. L ≥ λ/NA, far above cutoff; near cutoff the interfering pair is (0, ±1).
   Fix (statement): "If max_i L_i (NA + max_j|s_j|)/λ < 1 then K = {0}, Q_{x,z} = Σ_j w_j e_0e_0^H is independent of (x,z), ρ(M) = 1 for every M with M_0 ≠ 0, t_R = 1 and uniform (μ,λ) is a certificate (S = 0). Conversely, in d = 1, if L(NA+max|s|)/λ ≥ 1 then orders 0 and ±1 interfere for the extreme source points and t_R < 1 for targets whose bright and dark samples are separated by a single sinusoid of period L at every focus (in particular for dense lines of pitch L with Z = {0})." Drop "for any target".

19. **OK** — Uniform (μ,λ) gives S = Q − Q = 0 ⪰ 0 and Σλ = 1; the defocus phase drops out in |P|², so Q is also z-independent (needed for S = 0 with |Z| > 1). Requires |s_j| ≤ NA (true for σ ≤ 1); and s_max = σNA holds for the code's source quadratures (endpoints included).

## Lemma 6 (complex vs real masks)

20. **MAJOR** — Proof drops a conjugate on the pupil, and the hypothesis is wrong. Substituting k → −k: I_{M̃}(x) = Σ_s w_s |Σ_k P(−f_k + s/λ) conj(M_k) e^{−2πif_kx}|² = Σ_s w_s |Σ_k conj(P(−f_k + s/λ)) M_k e^{2πif_kx}|². The text writes conj(Σ_k P(f_k − s/λ) M_k e^{…}) — the pupil should be conj(P), not P. Hence I_{M̃} = image of M under the pupil P̂(g) := conj(P(−g)) and source −s. The identity I_{M̃,z} = I_{M,z} holds iff P̂ = P, i.e. the pupil is **Hermitian-symmetric**: |P| even and phase Φ **odd** (Φ(−g) = −Φ(g)). Evenness of P (as stated) is the right hypothesis only when P is real (z = 0, W = 0). With the defocus phase Φ_z (even in g): P̂_z = P_{−z}, so I_{M̃,z} = I_{M,−z}. The per-(x,z) claim "for all x, z" is therefore false for z ≠ 0; the consequence t_R(ℂ) = t_R(ℝ) survives iff Z = −Z (constraints are permuted; X_r = ½(X + conj(ΠXΠ)) has constraint values ½(Tr Q_{x,z}X + Tr Q_{x,−z}X)). Aberrations: odd Zernikes (tilt, coma Z7/Z8) preserve the identity; even ones (Z4 defocus, Z5/Z6 astigmatism, Z9 spherical) map to negated coefficients and break it; with even aberrations present, t_R(ℂ) < t_R(ℝ) is possible. The numerical confirmation (Table 7) covers only z = 0, W = 0.
   Fix (statement): "If w(−s) = w(s) and the pupil satisfies P(−g) = conj(P(g)) (even amplitude, odd phase; this holds for the aberration-free pupil at z = 0 and for odd aberrations), then I_{M̃,z} = I_{M,z}. For the defocus phase, I_{M̃,z} = I_{M,−z}; hence t_R(𝓜_ℂ) = t_R(𝓜_ℝ) whenever the focus set is symmetric (Z = −Z) and all aberrations are odd." Fix the proof line to carry conj(P) and use conj(P_z(f_k − s/λ)) = P_{−z}(f_k − s/λ).

21. **OK** — X = ½(MM^H + M̃M̃^H) is PSD and real-structured: conj(ΠMM^HΠ) = M̃M̃^H and conj(ΠM̃M̃^HΠ) = MM^H (Π^T = Π, J² = I). Equivalently X = UYU^H with Y real.

## Lemma 7 (vector imaging, aberrations)

22. **OK** — Rodrigues factors. With v = ẑ × k̂ = (−β, α, 0), R = I + [v]_× + [v]_×²/(1+γ): R x̂ = (1 − α²/(1+γ), −αβ/(1+γ), −α), R ŷ = (−αβ/(1+γ), 1 − β²/(1+γ), −β). Matches text and code. Under g → −g ((α,β) → (−α,−β), γ fixed): E_x, E_y even, E_z odd, obliquity 1/√γ even; the odd sign is global per channel so |·|² is invariant. TE/TM per source point: rad(−s) = −rad(s), azi(−s) = −azi(s), again a global sign. Unpolarized = ½(x̂ + ŷ) incoherent, weights w_s/|𝒫_s| as in code. The channel decomposition (s,p,c) is rank-one in M, so Theorems 1, 3 (constraint count unchanged) and Prop. 5 carry over.

23. **MINOR** — "Lemma 6 extends" inherits finding 20: requires the Hermitian pupil condition (odd phase) or Z = −Z with odd aberrations only; "any aberration" in Lemma 7's first sentence is fine for the quadratic structure, not for Lemma 6.

24. **MAJOR (model, §3 pupil)** — P_z(g) = 1{|g|λ ≤ NA} exp(i2πz(√(1−|g|²λ²) − 1)/λ) is ill-defined for immersion: at NA = 1.35 the orders with |g|λ ∈ (1, 1.35] are inside the pupil but the square root is imaginary (evanescent), and `Model` silently clips gl = min(|g|λ, 1), making the defocus phase constant on that annulus. The physical scalar defocus in an image medium of index n is (2πnz/λ)(√(1 − (|g|λ/n)²) − 1). This affects every DUV multi-focus number (Table 1 three-focus rows, ±80 nm) and the ModelX vector model (direction cosines use n_img = 1.7 while the defocus term uses n = 1 with clipping). The theorems are unaffected (any phase works); the quantitative multi-focus DUV results are for a nonphysical pupil.
   Fix: define P_z with the image-medium index and remove the clip; rerun the ±80 nm rows; or state explicitly that defocus is modelled in vacuum with the annulus clipped and treat those rows as illustrative.

## Lemma 8 (pixelation)

25. **OK/MINOR** — For pixels of width L/N *centred* at nL/N: M_k = sinc(k/N)·(1/N)Σ_n m_n e^{−2πikn/N} (code `sdp_bound_pixel_int`). For pixels starting at nL/N there is an extra factor e^{−iπk/N} (half-pixel shift), harmless but the convention should be stated. Surjectivity: N ≥ F = 2Kmax+1 ⇒ residues k mod N distinct ⇒ DFT restricted to K is onto ℂ^F; sinc(k/N) ≠ 0 for |k| < N. Real pixels map onto 𝓜_ℝ. Pixel-domain SDP = Fourier-domain SDP because G X_p G^H ranges over all of 𝕊^F_+ (X_p = G^+XG^{+H}). N < F only shrinks the class. Correct. Note: surjectivity does **not** hold for the nonnegative/DNN class (finding 17), and `sdp_bound_pixel` (used in exp3) omits the sinc factor; say "point-sampled" where it is used.

## Theorem 9 (branch-and-bound)

26. **OK** — (i) ρ(αr) = ρ(r) for α ≠ 0; pick i = argmax|r_j| and scale by 1/r_i (sign included), so ∪Φ_i with r_i = +1 covers ℝ^F∖{0} up to scale; nothing is lost. (ii) q(c+δ) = q(c) + 2(Qc)·δ + δ^TQδ, |2(Qc)·δ| ≤ 2|Qc|·h, 0 ≤ δ^TQδ ≤ h^T|Q|h; Q̃ = Re(U^HQU) is PSD (real part of a PSD Hermitian matrix). (iii) numerator ≥ max_i q_i^lo ≥ 0, denominator ≤ min_j p_j^hi; multi-focus handled as max over (d,z), min over (b,z) in code, consistent with the common-dose ρ. (iv) Splitting preserves the cover, so min of LB over *all* leaves (pruned and still-active) is a valid lower bound at any stopping time, independently of the incumbent U and of pruning decisions; on convergence L ≥ (1−tol)U ≥ (1−tol)ρ*. Code `bnb_cert.py` implements exactly this (L = min(pruned leaves, active heap top)).

27. **MINOR** — (i) says "min over the faces"; the faces are compact but ρ is not lower semicontinuous at points where min_B I = max_D I = 0 (the Remark's "vanishing together" case), so write inf. The B&B is unaffected (LB bounds the inf; boxes near such points have d_lo, b_hi → 0 together).

28. **MINOR** — Make explicit in the statement that "leaves" includes unexpanded active nodes, and that U is only needed for pruning/termination, not for validity of L.

29. **MAJOR (wording)** — "Rigorous" is claimed in the abstract and title of Theorem 9, but the floating-point safeguard (error terms × (1+1e-9), LB − 1e-12) is a heuristic, not a proof: rounding in q(c), Qc, the division, and the ≈1e-16 indefiniteness of the symmetrized Q̃ are not bounded by an argument, only plausibly dominated. The Numerical note is honest about this; the theorem statement and abstract are not.
   Fix: "Theorem 9 is rigorous in exact arithmetic; the implementation uses inflated error terms (not interval arithmetic), so the certified gains (≥ 6.4%) should be read as exceeding the floating-point perturbation (~1e-12) by many orders of magnitude rather than as formally verified." Replace "rigorous branch-and-bound certifier" in the abstract by "branch-and-bound certifier (exact-arithmetic proof, floating-point implementation)".

## Sampling remark

30. **OK** — Logic correct: sampled constraints ⊂ continuous constraints ⇒ ρ*_cont ≥ ρ*_sampled ≥ certificate, provided every sample point lies in the continuous constraint set.

31. **MINOR** — That proviso is violated at the 1e-9 nm level in 1D (`keep = d ≥ margin − 1e-9`) and at the ≈0.5–0.7 nm level in 2D: `samples_2d` takes points with EDT distance in [τ, τ+0.6) computed on a 1-nm grid against pixel *centres*, so true boundary distance can be < τ by up to half a pixel plus band. The "placed at exactly ±τ" wording in §3 is accurate for 1D only. State: 2D samples lie at distance τ ± 0.5 nm; the bound is valid for the continuous problem with tolerance τ − 0.5 nm.

32. **MINOR** — Abstract: "valid for the continuous constraint set when the sample set contains the edge-offset points" conflates validity with tightness. Any sample subset gives a valid bound; the edge-offset points are what make it (empirically) tight. Reword: "valid for the continuous constraint set for any sample subset; empirically tight once the edge-offset points are included".

## §3 definitions: τ and exposure latitude

33. **MAJOR (definition/wording)** — EL formula. Derivation: with dose E and resist threshold θ, all B print and no D prints iff E·β ≥ θ and E·δ < θ, i.e. E ∈ [θ/β, θ/δ) with β = min_{B×Z} I, δ = max_{D×Z} I (one dose across Z is exactly what makes β, δ run over the whole focus set). The window is nonempty iff ρ = δ/β < 1. Then
   - half-width relative to the window-centre dose E_c = (θ/β + θ/δ)/2: (E_max − E_min)/(E_max + E_min) = (β − δ)/(β + δ) = **(1−ρ)/(1+ρ)** — the paper's formula, exact;
   - full width relative to E_c (the conventional "EL = (E_max − E_min)/E_nom"): 2(1−ρ)/(1+ρ);
   - relative to E_min: 1/ρ − 1; relative to the geometric-mean dose: (1−ρ)/√ρ.
   So the formula is exact only for the ±EL (half-range) convention about the window centre. THEORY.md's "threshold resist placed at the midpoint" is wrong: at nominal dose E_c the effective intensity threshold is θ/E_c = 2βδ/(β+δ), the *harmonic* mean of δ and β, not the midpoint; a threshold at the arithmetic midpoint (β+δ)/2 gives an asymmetric dose window (−(1−ρ)/2, +(1−ρ)/(2ρ)).
   Fix (§3): "For a threshold resist exposed with one dose across the focus set, the doses for which every bright sample prints and no dark sample prints form the interval [θ/β, θ/δ); its half-width relative to its centre is EL(ρ) = (1−ρ)/(1+ρ), i.e. the dose may deviate by ±EL from the window-centre dose (equivalently, intensity threshold at the harmonic mean of δ and β). The conventional full-range exposure latitude is 2·EL." Replace "midpoint" in THEORY.md accordingly. Also note that ρ ↦ EL is decreasing, so max EL ⇔ min ρ (used throughout).

34. **OK** — Edge tolerance τ: samples at exactly e ± τ in 1D (code `samples_1d`) and interior points at distance ≥ τ; this is the CD-tolerance interpretation (edge within ±τ ⇔ bright point at +τ prints and dark point at −τ does not) that makes the EL above the standard CD-spec EL. Monotonicity of the bound in τ is consistent.

35. **MINOR** — 𝓜_ℝ is called "phase-shift masks" and the "best single mask" is a continuous-tone real band-limited mask (grayscale amplitude, 0/π phase), not a manufacturable {−1,0,1} alt-PSM; 𝓜_+ is grayscale nonnegative, not binary. Bounds remain valid for the manufacturable subclasses (smaller class ⇒ larger ρ*), but the reported "gaps" are to continuous-tone masks. Add one sentence in §3.

36. **MINOR** — Collected-order set K: correct that orders outside K never contribute (|f_k + s_j/λ| ≤ NA/λ ⇒ |f_k| ≤ (NA+|s_j|)/λ); but orders inside K may contribute to no constraint (finding 4), which is why the common null space in §5 arises. One clause linking the two.

---

### Summary of required changes (by severity)
- MAJOR: 8 (complex Pataki r²), 10 (attainment of (R) — add argument), 12 (Fejér–Riesz claim false), 17 (binary-mask claim only for 5-nm-grid masks), 18 (Prop. 5 "equivalently"/converse false; proof mechanism wrong), 20 (Lemma 6 proof drops conj(P); hypothesis must be Hermitian pupil / Z = −Z; even aberrations break it), 24 (pupil formula ill-defined for NA > 1; code clips), 29 ("rigorous" overclaim for floating-point B&B), 33 (EL convention: half-range about window centre; "midpoint" wrong).
- MINOR: 4, 5, 9, 11, 14, 15, 16, 23, 25, 27, 28, 31, 32, 35, 36.
- No CRITICAL: no headline certified bound is invalidated (Theorem 1 and the sampled-⊂-continuous argument are sound); the affected headline items are the complex rank count, the Lemma 6 scope, the binary-mask forfeit claim's class, the DUV multi-focus pupil, and the EL convention.
