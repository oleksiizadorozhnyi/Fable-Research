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
- Exposure latitude associated with ratio ρ: EL(ρ) = (1-ρ)/(1+ρ) (relative dose window within which every B point
  prints and no D point prints, for a threshold resist placed at the midpoint).

## 2. Problems
(P_𝓜)  ρ*(𝓜) = inf_{M∈𝓜\{0}} ρ(M)  ⇔  min t s.t. M^H Q_{x,z} M ≤ t ∀(x,z)∈D×Z,  M^H Q_{x,z} M ≥ 1 ∀(x,z)∈B×Z.
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
F) masks of the class. Moreover ρ*(𝓜) = t_R(𝓜) iff (R) has an optimal solution of rank one; an extreme optimal
solution has rank r with r(r+1)/2 ≤ |D×Z| + |B×Z| + 1 (Pataki), so at most that many exposures ever help.
Proof. X ⪰ 0 ⇔ X = Σ_{i=1}^{r} M_i M_i^H (r ≤ F), and Tr(Q X) = Σ_i M_i^H Q M_i = summed image. Hence feasible X of
(R) ↔ exposure ensembles (M_i) with summed-image constraints; objectives coincide. For 𝓜_ℝ use Lemma 2 to choose
M_i ∈ V_ℝ. Rank-one optimum ↔ a single mask attaining t_R; conversely any single optimal mask gives a rank-one
feasible X with the same value. Pataki's theorem on extreme points of spectrahedra bounds the rank. ∎

**Theorem 4 (exact convex description of 𝓜_+ in 1D).** d = 1, K = {-Kmax..Kmax}. M ∈ ℂ^F is the in-band Fourier
coefficient vector of some nonnegative measure on the circle iff the Toeplitz matrix T(M) = (M_{a-b})_{a,b=0..Kmax} ⪰ 0
(Carathéodory–Toeplitz). If T(M) ≻ 0 the measure may be taken to be a bounded nonnegative function (a strictly
positive trigonometric polynomial of degree Kmax, by the Fejér–Riesz theorem), hence a legitimate mask. Consequently the
closure of 𝓜_+ equals {M : T(M) ⪰ 0}, and (P_{𝓜_+}) is a QCQP over an SDP-representable cone.
Proof sketch. Carathéodory–Toeplitz: a finite sequence is a truncated trigonometric moment sequence of a nonneg
measure iff its Toeplitz matrix is PSD. Riesz–Fejér: a nonneg trig polynomial is |h(e^{iθ})|²; with T(M) ≻ 0 the
sequence extends to a positive-definite function with strictly positive density (take the Toeplitz Gaussian/maximum-
entropy extension, whose spectral density is positive and bounded). ∎  [Open proof obligation: boundedness statement
to be checked against a reference; the SDP-cone statement for closures is standard.]
Note. In d = 2 the block-Toeplitz PSD condition is necessary but not sufficient (no 2D Fejér–Riesz); it still gives
a valid relaxation.

**Proposition 5 (resolution limit recovery).** If no two collected orders k ≠ k' both pass the pupil for a common
source point (i.e., L·(NA+max|s|)/λ < 1 — equivalently pitch below λ/((1+σ)NA)), then Q_{x,z} is independent of x,
so ρ(M) ≥ 1 for all M, t_R = 1, and (μ,λ) = (uniform, uniform) is a certificate. Above this pitch, t_R < 1.
Proof. With a single order per source point the image is |M_{k_j}|² per source point — x-independent. The dual
certificate with μ uniform, λ uniform gives S = (1/|D|)Σ Q − (1/|B|)Σ Q = 0 ⪰ 0, Σλ = 1. For the converse, choose
M with two interfering orders ±1 to get modulation. ∎

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

**Lemma 6 (complex vs real masks under symmetric illumination).** If the source is symmetric (w(s) = w(−s)) and the
pupil is even (P_z(g) = P_z(−g)), then for every complex order vector M the conjugate-reflected vector M̃ := conj(ΠM)
produces the same image, I_{M̃,z}(x) = I_{M,z}(x) for all x, z. Consequently t_R(𝓜_ℂ) = t_R(𝓜_ℝ): the relaxation over
complex-transmission masks equals that over real masks (X := ½(MM^H + M̃M̃^H) is real-structured with the same image).
Proof. I_{M̃}(x) = Σ_s w_s |Σ_k P(f_k + s/λ) conj(M_{-k}) e^{2πi f_k x}|² ; substitute k → −k and use evenness of P:
= Σ_s w_s |Σ_k P(−f_k + s/λ) conj(M_k) e^{−2πi f_k x}|² = Σ_s w_s |conj(Σ_k P(f_k − s/λ) M_k e^{2πi f_k x})|²
= Σ_s w_{s} |Σ_k P(f_k + (−s)/λ) M_k e^{2πi f_k x}|² = Σ_{s'} w_{s'} |Σ_k P(f_k + s'/λ) M_k e^{2πi f_k x}|² = I_M(x),
using w(−s) = w(s). The matrix ½(MM^H + M̃M̃^H) is PSD, satisfies X = conj(ΠXΠ) (since Π conj(M̃ M̃^H) Π = MM^H),
and has image ½(I_M + I_{M̃}) = I_M. Hence every complex-feasible X yields a real-structured feasible X with the same
constraint values; the converse inclusion is trivial. ∎  (Numerically confirmed: Table 7, all cases equal.)
Note: for SINGLE masks the complex optimum can still be better than the real optimum; the lemma concerns the relaxation
(equivalently, a two-exposure ensemble {M, M̃}).
