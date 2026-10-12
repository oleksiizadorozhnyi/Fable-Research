# Formal problem construction and theory (Direction A)

## 1. Definitions
- Periodic cell Λ = [0,L_1)×...×[0,L_d), d∈{1,2}. Mask transmission m: Λ→ℂ (periodic), Fourier coefficients
  M_k = (1/|Λ|)∫ m(x) e^{-2πi f_k·x} dx, f_k = (k_1/L_1,...,k_d/L_d), k∈ℤ^d.
- Illumination: finite set of source points s_j ∈ ℝ^d (direction cosines, |s_j| ≤ σ NA), weights w_j>0, Σw_j=1.
- Pupil P_z(g) = 1{|g|λ ≤ NA}·exp(i 2π z (√(1-|g|²λ²) - 1)/λ) (scalar, defocus z).
- Collected orders: K = {k : |f_k| ≤ (NA+max_j|s_j|)/λ}; F=|K|. M ∈ ℂ^F denotes (M_k)_{k∈K}.
- Abbe image at x, focus z: I_{M,z}(x) = Σ_j w_j |Σ_{k∈K} P_z(f_k + s_j/λ) M_k e^{2πi f_k·x}|² = M^H Q_{x,z} M,
  Q_{x,z} = Σ_j w_j c_{j,x,z} c_{j,x,z}^H ⪰ 0 (Hermitian PSD, rank ≤ #source points).
  (Standard Abbe/Hopkins thin-mask scalar model; orders outside K never contribute.)
- Target: finite sample sets B (bright) and D (dark) ⊂ Λ, and focus set Z. Constraint set 𝒞 = (B∪D)×Z.
- Mask classes: 𝓜_ℂ = ℂ^F (complex transmission); 𝓜_ℝ = {M : M_{-k} = conj(M_k)} (real transmission, e.g. PSM);
  𝓜_+ ⊂ 𝓜_ℝ: in-band orders of some real m ≥ 0 (binary/attenuated masks).
- Dark/bright ratio of a mask: ρ(M) = max_{(x,z)∈D×Z} I_{M,z}(x) / min_{(x,z)∈B×Z} I_{M,z}(x). Scale-invariant.
- Exposure latitude associated with ratio ρ: for a threshold resist (threshold θ, dose E, printed where E·I ≥ θ) every
  B point prints and no D point prints iff θ/β ≤ E < θ/δ (β = min_B I, δ = max_D I). The full-range dose latitude,
  window width over its centre dose, is EL(ρ) = 2(1−ρ)/(1+ρ) with ρ = δ/β; the half-range (±) latitude about the
  centre dose (θ(β+δ)/(2βδ), i.e. threshold at the harmonic mean of δ and β) is (1−ρ)/(1+ρ). The manuscript reports
  the full-range value. The edge tolerance τ is quoted as a fraction of the minimum CD (10% per edge).

## 2. Problems
(P_𝓜)  ρ*(𝓜) = inf_{M∈𝓜\{0}} ρ(M)  (inf, not min: ρ need not be lower semicontinuous where bright and dark
       intensities vanish together; ρ := +∞ when min_B I = 0)  ⇔  min t s.t. M^H Q_{x,z} M ≤ t ∀(x,z)∈D×Z,  M^H Q_{x,z} M ≥ 1 ∀(x,z)∈B×Z.
       This is a nonconvex homogeneous QCQP (the bright constraints are reverse-convex).
(R_ℂ)  t_R(ℂ) = min t s.t. Tr(Q_{x,z} X) ≤ t (D×Z), Tr(Q_{x,z} X) ≥ 1 (B×Z), X ⪰ 0 Hermitian F×F.
(R_ℝ)  same with X restricted to the real-structured cone: X = conj(Π X Π), Π = order-flip permutation.
       Equivalently (Lemma 2) X = U Y U^H with Y real symmetric PSD, Q̃ = Re(U^H Q U).
(ME_r) multi-exposure problem: masks M_1..M_r ∈ 𝓜, summed image Σ_i I_{M_i,z}(x); ratio of summed image minimized.
(Dual) max Σ_{b} λ_b  s.t. λ ≥ 0, μ ≥ 0, Σ_d μ_d = 1, S(μ,λ) := Σ_d μ_d Q_d − Σ_b λ_b Q_b ⪰ 0  (indices run over D×Z, B×Z).

## 3. Results

**Theorem 1 (certificate; weak duality).** Let μ ≥ 0 on D×Z with Σμ = 1 and λ ≥ 0 on B×Z satisfy S(μ,λ) ⪰ 0
(for 𝓜_ℝ: S̃ := Re(U^H S U) ⪰ 0 suffices). Then every mask M ∈ 𝓜 (any class ⊂ ℂ^F; for the real version, any
M ∈ 𝓜_ℝ) satisfies ρ(M) ≥ Σ_b λ_b. In particular no single mask, and no incoherent combination of masks, achieves a
dark/bright ratio below Σλ_b.
Proof. Let M be any mask with min_B I = β > 0 and max_D I = δ. Then δ = δ Σμ_d ≥ Σ_d μ_d I_d = Σ μ_d M^H Q_d M
 = M^H S M + Σ_b λ_b M^H Q_b M ≥ 0 + β Σ_b λ_b. Hence ρ(M) = δ/β ≥ Σλ_b. For M ∈ 𝓜_ℝ, M = U r with r real and
 M^H S M = r^T Re(U^H S U) r, so S̃ ⪰ 0 suffices. For a sum of images Σ_i I_{M_i}: the same chain holds with
 M^H Q M replaced by Σ_i M_i^H Q M_i (each term ≥ 0 after using S ⪰ 0). ∎
Remark. The certificate is verified by one PSD test (Cholesky/eigenvalue) — independent of the SDP solver.

**Lemma 2 (real structure).** Let J(M) = conj(ΠM) (antilinear involution). V_ℝ := Fix(J) = 𝓜_ℝ is a real
F-dimensional subspace and the Hermitian form ⟨u,v⟩ is real on V_ℝ. A Hermitian PSD X satisfies X = conj(ΠXΠ) iff
XJ = JX, iff X maps V_ℝ to itself; then X|_{V_ℝ} is a real symmetric PSD operator and X = Σ_i θ_i v_i v_i^H with
v_i ∈ V_ℝ, θ_i ≥ 0. Writing V_ℝ = U ℝ^F with U unitary gives X = U Y U^H, Y ∈ 𝕊^F_+, and Tr(Q X) = Tr(Q̃ Y).
Proof. Standard: Πconj(·) commutes with X iff X = conj(ΠXΠ); for M ∈ V_ℝ, J(XM) = X(JM) = XM. Restricting a
Hermitian operator to an invariant real form gives a real symmetric operator; spectral theorem. ∎

**Theorem 3 (relaxation = unlimited multiple exposure).** For 𝓜 ∈ {𝓜_ℂ, 𝓜_ℝ}: t_R(𝓜) = inf_r ρ*(ME_r) = ρ*(ME_F),
i.e. the SDP relaxation value is exactly the optimal ratio achievable by incoherently superposing images of (at most
F) masks of the class. Moreover, if (R) has a rank-one optimal solution then ρ*(𝓜) = t_R(𝓜) and the infimum in (P_𝓜) is attained;
conversely, if ρ*(𝓜) = t_R(𝓜) is attained by some mask then (R) has a rank-one optimal solution. (Attainment is not
automatic: ρ is scale-invariant and continuous on the unit sphere away from {min_B I = 0}, and the infimum can be
approached by masks whose bright and dark intensities vanish together; in that case ρ* is still certified by Theorem 1.)
Writing the inequalities of (R) with slack variables, (R) is a standard-form SDP with m = |D×Z| + |B×Z| equality
constraints and block-diagonal cone (X, slacks, t = t⁺ − t⁻); (R) attains its optimum: the feasible set is closed, t is bounded below by 0, and along any minimizing sequence X_n
the constraints Tr(Q_b X_n) ≥ 1 and Tr(Q_d X_n) ≤ t_n ≤ t_1 bound Tr(X_n) on the complement of the common null space of
all Q's (on that null space X may be set to zero without changing any constraint), so a convergent subsequence exists.
Pataki's theorem then gives an optimal extreme solution whose blocks satisfy Σ_blocks r_b(r_b+1)/2 ≤ m for REAL
symmetric blocks (the 𝓜_ℝ problem in the real basis of Lemma 2), hence rank(Y) = r with r(r+1)/2 ≤ |D×Z| + |B×Z|;
for the complex Hermitian problem (𝓜_ℂ) the corresponding count is r² ≤ |D×Z| + |B×Z|. At most that many exposures
are ever needed to attain t_R.
Proof. X ⪰ 0 ⇔ X = Σ_{i=1}^{r} M_i M_i^H (r ≤ F), and Tr(Q X) = Σ_i M_i^H Q M_i = summed image. Hence feasible X of
(R) ↔ exposure ensembles (M_i) with summed-image constraints; objectives coincide. For 𝓜_ℝ use Lemma 2 to choose
M_i ∈ V_ℝ. Rank-one optimum ↔ a single mask attaining t_R; conversely any single optimal mask gives a rank-one
feasible X with the same value. Pataki's theorem on extreme points of spectrahedra bounds the rank. ∎

**Theorem 4 (exact convex description of 𝓜_+ in 1D).** d = 1, K = {-Kmax..Kmax}. M ∈ ℂ^F is the in-band Fourier
coefficient vector of some nonnegative measure on the circle iff the Toeplitz matrix T(M) = (M_{a-b})_{a,b=0..Kmax} ⪰ 0
(Carathéodory–Toeplitz). If T(M) ≻ 0 the measure may be taken to be a bounded, strictly positive continuous function (the maximum-entropy
autoregressive density σ²/|a(e^{iθ})|², a rational — not polynomial — function), hence a legitimate mask after scaling. Consequently the
closure of 𝓜_+ equals {M : T(M) ⪰ 0}, and (P_{𝓜_+}) is a QCQP over an SDP-representable cone.
Proof. (Necessity) For m ≥ 0 and any v ∈ ℂ^{Kmax+1}, v^H T(M) v = (1/|Λ|)∫ m(x) |Σ_a v_a e^{2πi a x/L}|² dx ≥ 0.
(Sufficiency, T(M) ⪰ 0) Carathéodory–Toeplitz [Grenander–Szegő, Ch. 1–2]: a finite sequence with PSD Toeplitz matrix
is the truncated moment sequence of a nonnegative measure on the circle; such measures are weak-* limits of bounded
nonnegative functions (Fejér means), so M lies in the closure of 𝓜_+. (Interior, T(M) ≻ 0) The Levinson–Durbin
recursion yields the maximum-entropy autoregressive extension: a polynomial a(z) = 1 + a_1 z + … + a_{Kmax} z^{Kmax}
with all zeros in |z| > 1 and σ² > 0 such that f(θ) = σ²/|a(e^{iθ})|² has Fourier coefficients M_0..M_{Kmax} (Burg;
Grenander–Szegő §5). Since a has no zeros on the unit circle, f is continuous, bounded and bounded away from zero;
scaled by its maximum it is an admissible mask with |m| ≤ 1. Hence int{T ⪰ 0} ⊂ 𝓜_+ ⊂ {T ⪰ 0}, and the closure of
𝓜_+ is {T(M) ⪰ 0}, an SDP-representable cone. ∎
Note. In d = 2 the block-Toeplitz PSD condition is necessary but not sufficient (the two-dimensional truncated
trigonometric moment problem is not characterised by block-Toeplitz positivity alone; see Dumitrescu, Positive
Trigonometric Polynomials and Signal Processing Applications, 2007); it still gives a valid relaxation. T(M) is
Hermitian for M ∈ 𝓜_ℝ (M_{-k} = conj M_k), which is the case of interest.

**Proposition 5 (resolution limit recovery).** If no two collected orders k ≠ k' both pass the pupil for a common
source point — in d = 1 this is exactly L·(NA+max|s|)/λ < 1, i.e. pitch below λ/((1+σ)NA); in d = 2 it holds when
max(L_1, L_2)·(NA+max|s|)/λ < 1 — then Q_{x,z} is independent of x, so ρ(M) ≥ 1 for all M, t_R = 1, and
(μ,λ) = (uniform, uniform) is a certificate.
Proof. With a single order per source point the image is Σ_j w_j |M_{k_j}|² — x-independent. The dual certificate with
μ uniform, λ uniform gives S = (1/|D|)Σ Q − (1/|B|)Σ Q = 0 ⪰ 0, Σλ = 1. ∎
Remark. The converse is not true for every target: just above the cutoff only the orders 0 and ±1 interfere for some
source points and the image is a sinusoid at the cell period, so a target whose bright and dark samples cannot be
separated by such a sinusoid (e.g. B = {0, L/2}, D = {L/4, 3L/4}) still has t_R = 1. For dense lines (one line per
period) t_R < 1 immediately above the cutoff, which is the classical statement verified in the experiments.

## 4. Hypothesis, success and failure criteria
H1 (certification): on realistic 1D/2D periodic targets the certified bound t_R and the best feasible mask found by
local optimization agree within a few percent; i.e. local ILT-type optima are near-globally optimal and this can be
PROVEN per instance.  Falsified if typical gaps are large (>20%) AND rounding cannot close them (then the SDP is a weak bound).
H2 (multi-exposure value): rank(X*) > 1 instances exist where the bound is strictly below any single mask, quantifying
the benefit of incoherent multi-exposure. Falsified if all instances are rank-one (then Theorem 3 is vacuous in practice).
H3 (sanity): k1 cutoffs recovered exactly (Prop. 5).
Success: H3 holds; H1 or H2 supported with independent certificate verification; 2D instances solvable in minutes.
Failure: certificates fail independent verification; SDP intractable for 2D; bounds trivially loose.

## 5. Remark on sampled vs continuous constraints
The certified bound is computed for finite sample sets B, D (edge-offset points at exactly ±margin from every target
edge, plus interior points). The continuous problem (constraints for all x in the bright/dark regions at distance
≥ margin from edges) has MORE constraints, hence a larger optimal ratio; therefore any certificate for the sampled
problem is also a valid lower bound for the continuous problem. Numerically the interior points are inactive (dual
weight 0): the bound is determined by the edge-offset points, and is invariant to interior spacing (verified).

**Lemma 6 (complex vs real masks under symmetric illumination).** Let the source be symmetric (w(s) = w(−s)) and the
pupil satisfy P(−g) = conj(P(g)) (even amplitude, odd phase: e.g. no aberration, or odd aberrations such as coma).
Then for every complex order vector M the conjugate-reflected vector M̃ := conj(ΠM) produces the same image at best
focus, I_{M̃,0}(x) = I_{M,0}(x). Since the defocus phase is even in g, with defocus the images are related by
I_{M̃,z}(x) = I_{M,−z}(x). Consequently, for a focus set with Z = −Z (e.g. {0, ±z_0}) and no even aberrations,
t_R(𝓜_ℂ) = t_R(𝓜_ℝ): the relaxation over complex-transmission masks equals that over real masks
(X := ½(MM^H + M̃M̃^H) is real-structured with the same constraint values). Even aberrations (astigmatism, spherical,
Z4) break the identity.
Proof. I_{M̃}(x) = Σ_s w_s |Σ_k P(f_k + s/λ) conj(M_{-k}) e^{2πi f_k x}|² ; substitute k → −k:
= Σ_s w_s |Σ_k P(−f_k + s/λ) conj(M_k) e^{−2πi f_k x}|² = Σ_s w_s |conj(Σ_k conj(P(−f_k + s/λ)) M_k e^{2πi f_k x})|²
= Σ_s w_s |Σ_k P(f_k − s/λ) M_k e^{2πi f_k x}|² using conj(P(−g)) = P(g) (for the defocus factor, conj flips the sign of z)
= Σ_s w_{s} |Σ_k P(f_k + (−s)/λ) M_k e^{2πi f_k x}|² = Σ_{s'} w_{s'} |Σ_k P(f_k + s'/λ) M_k e^{2πi f_k x}|² = I_M(x),
using w(−s) = w(s). The matrix ½(MM^H + M̃M̃^H) is PSD, satisfies X = conj(ΠXΠ) (since Π conj(M̃ M̃^H) Π = MM^H),
and has image ½(I_M + I_{M̃}) = I_M. Hence every complex-feasible X yields a real-structured feasible X with the same
constraint values; the converse inclusion is trivial. ∎  (Numerically confirmed: Table 7, all cases equal.)
Note: for SINGLE masks the complex optimum can still be better than the real optimum; the lemma concerns the relaxation
(equivalently, a two-exposure ensemble {M, M̃}).

## 6. Model extensions (session 2)
**Lemma 7 (vector imaging and aberrations preserve the quadratic structure).** With a pupil phase W(g) (any
aberration) and the high-NA vector model in which the image-side plane wave of direction cosines (α,β,γ) carries the
mask-side polarization p rotated by the minimal rotation taking ẑ to (α,β,γ) (Rodrigues; for p = x̂:
E = (1−α²/(1+γ), −αβ/(1+γ), −α)), the image I(x) = Σ_s w_s Σ_{p∈𝒫_s} Σ_{c∈{x,y,z}} |Σ_k V_{c,p,s}(k) M_k e^{2πi f_k·x}|²
is a sum of rank-one Hermitian forms in M. Hence Theorems 1, 3 and Proposition 5 and the SDP relaxation apply verbatim
with the channel set (s,p,c) in place of s. Lemma 6 extends when the polarization ensemble is reflection-symmetric
(x, y, unpolarized, or TE/TM with a symmetric source): the rotation factors satisfy E_{x,y}(−g) = E_{x,y}(g),
E_z(−g) = −E_z(g), so |·|² is invariant. Proof: direct substitution. ∎
**Lemma 8 (pixelation).** If the mask is piecewise constant on N ≥ F pixels per period (1D), the map from pixel values
to collected orders, M_k = sinc(k/N)·DFT_k(m), is surjective onto ℂ^F (sinc(k/N) ≠ 0 for |k| < N), so for the scale-free
ratio problem the pixel-domain relaxation equals the Fourier-domain relaxation; pixelation below the collected-band
count (N < F) restricts the mask class and can only raise the bound. (Verified: Table "pixelation".) ∎

**Theorem 9 (rigorous single-mask lower bound by branch and bound).** Let Q_1..Q_m ⪰ 0 (real symmetric, F×F) be the
dark forms (over D×Z) and P_1..P_n ⪰ 0 the bright forms (over B×Z), and ρ(r) = max_i r^TQ_i r / min_j r^TP_j r for r ≠ 0.
(i) ρ is invariant under r ↦ αr (α ≠ 0), so ρ*(𝓜_ℝ) = min over the union of faces Φ_i = {r : r_i = 1, |r_j| ≤ 1 ∀j}.
(ii) For a box 𝔅 = {c + δ : |δ_j| ≤ h_j} ⊂ Φ_i and q(r) = r^TQr with Q ⪰ 0:
     q(c) − 2|Qc|·h ≤ q(r) ≤ q(c) + 2|Qc|·h + h^T|Q|h  for all r ∈ 𝔅,
   because q(c+δ) = q(c) + 2(Qc)·δ + δ^TQδ, |2(Qc)·δ| ≤ 2|Qc|·h, 0 ≤ δ^TQδ ≤ h^T|Q|h.
(iii) Hence inf_{𝔅} ρ ≥ LB(𝔅) := max_i max(q_i(c) − 2|Q_ic|·h, 0) / min_j (p_j(c) + 2|P_jc|·h + h^T|P_j|h).
(iv) If a finite family of boxes covers ∪Φ_i, then ρ*(𝓜_ℝ) ≥ min_𝔅 LB(𝔅). Branch and bound refines boxes whose LB is
   below the incumbent U(1−tol) (U = ρ at a box centre, an upper bound on ρ*); at any stopping time the minimum of LB
   over all leaves is a valid lower bound, and on convergence it is within tol of ρ*.
Proof. (i) sign and scale invariance; every r ≠ 0 has a coordinate of maximal modulus, which can be scaled to +1.
(ii)–(iv) as written. ∎
Numerical note. LB is evaluated in floating point; error terms are multiplied by (1+1e-9) and 1e-12 is subtracted; the
forms Q, P are themselves floating-point images of exact complex exponentials with relative error ~1e-15. These
safeguards are heuristic, not a proof: the certificate is exact in exact arithmetic and 'certified' in the manuscript
means 'exact modulo double-precision rounding of order 1e-12 relative', many orders below the reported gains (≥ 6%).
A directed-rounding (interval) implementation is future work.

### Session-3 corrections (2026-10-12, code audit)
- Lemma 8: the pixel-to-order map carries the half-pixel phase, M_k = sinc(k/N) e^{-iπk/N} DFT_k(m), for pixel n covering
  [nL/N, (n+1)L/N); the earlier implementation modelled pixels centred at nL/N (N ≥ F bounds unaffected; N < F bounds change).
- Model: orders that no source point brings into the pupil ("dead orders") are discarded from K so that every order in K
  is constrained; F counts live orders only (bounds unchanged, reported ranks no longer inflated by the unconstrained block).
- Two-dimensional sampling: targets are unions of axis-aligned rectangles; samples lie exactly on the offset contours at
  distance τ (asserted), replacing the pixel-quantised distance transform whose rings were empty for τ not a multiple of
  the raster step and biased by up to one raster step otherwise. All 2D numbers were regenerated.
- Certified multi-exposure gain (Theorem 9 use): the reference is the ratio realized by the SDP ensemble (an upper bound
  on t_R), not the rigorized lower bound, so gain = (L − t̄)/t̄ ≤ (L − t_R)/t_R is itself certified.

### Referee-driven corrections (2026-10-12, major revision)
- Theorem 3: the attainment argument via boundedness of Tr X was wrong (bright constraints are lower bounds: X_n =
  diag(1/n, n) with Q_d = diag(1,0), Q_b = diag(0,1) is a feasible minimizing sequence with unbounded trace). Replaced by
  the Attainment Lemma: reduce to N^⊥ (N = common dark null space); bright forms vanishing on N give a reduced SDP whose
  dual is strictly feasible, hence attained (Vandenberghe–Boyd 1996); the others are satisfied by adding c(I−P).
- Single-mask attainment: ρ* IS attained (generic large component in N makes the B_N bright values large; the reduced
  ratio is lower semicontinuous on the compact sphere of N^⊥). The remark claiming possible non-attainment was wrong.
- Theorem 4: {T ≻ 0} ⊂ C_+ ⊂ {T ⪰ 0} = measures, both strict (point mass M_0 = M_±1 = a is not a bounded mask);
  cl(C_+) = {T ⪰ 0}; same optimal value for the nonnegative-mask problem and its measure relaxation.
- Proposition 5: the radial condition L(NA + max|s|) < λ is sufficient, not necessary (2D sources for 1D masks:
  exact passing condition |k/L + s_x/λ| ≤ sqrt(NA² − s_y²)/λ; exact criterion and the sufficient condition
  L max_j(|s_jx| + sqrt(NA² − s_jy²)) < λ, equal to the classical disk cutoff).
- K = −K is enforced: an order is dropped only when neither it nor its mirror is ever collected.
- Certificates: exact 50-digit rebuild from the model equations (exp21) replaces the "assembled from rounded
  coefficients" check; wording "safe against round-off" removed.
