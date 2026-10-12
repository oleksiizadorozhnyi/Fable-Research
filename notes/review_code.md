# Adversarial code audit — `code/certlitho.py`, `code/bnb_cert.py`, `code/bnb_sdp.py`, `code/pipeline.py`

Audited 2026-10-12 (commit 3e94069 + working tree; the session-3 `Model`/`ModelX`/`build_*` appended at the end of
`certlitho.py` are the effective definitions, as is `rigorize_certificate` v3 and the vectorised `exact_minimax_polish`).
Test scripts: `code/audit_*.py` (run from `code/` with `python3 -W ignore audit_X.py`). No existing file was modified.
A clean reproduction (`run_all_par.sh`) was running during the audit; numbers below marked "stored" refer to
`results/*.json` as found, "as coded" to a fresh run of the current code.

## Summary
| # | label | item |
|---|---|---|
| 1 | CRITICAL | `samples_2d`: edge-offset ring is built from a pixel-quantised EDT; for the new tolerances (τ = 7, 5, 4.4 nm with `fine=2`) the sample set contains **no bright edge-offset points at all** → 2-D bounds far too low (0.385 vs 0.815 at the true τ = 7) |
| 2 | MAJOR | `samples_2d`: even for τ a multiple of `fine`, samples sit 1–2 nm (one `fine`) closer to the true edge than τ → published 2-D bounds are for τ_eff ≈ τ − 2 nm (0.584 vs 0.524 at exact τ = 12) and the "valid for the continuous problem" claim fails at the stated τ |
| 3 | MAJOR | Immersion defocus: `Model` defaults to `n_img=1.0`; several DUV multi-focus experiments (exp14, exp19, exp11, exp11b, exp15 PW80 rows) omit `n_img=1.44` → 712 passing (order, source) pairs have |g|λ > 1 and get a clipped, non-physical defocus phase; exp15 then "verifies" a different model than exp1 |
| 4 | MAJOR | `exp15_verify_all.py` no longer rebuilds the reported instances (τ = 8/10 vs TAU_DUV = 5; no `n_img`) → the "every reported certificate independently verified" claim is false for the new run |
| 5 | MAJOR | Dead orders (never inside the pupil for any source point) are kept in K → unconstrained block in Y inflates the reported rank (T-junction annular: rank 6 reported, 2 on live coordinates) and forces the null-space projection in `rigorize_certificate`, whose soundness argument is incomplete |
| 6 | MAJOR | Degenerate instances: the stored exp17 record `extreme_duty_95pct/coh` has rigorous = 1.48e-5 with best_feasible = 1.3e-24 (a feasible mask below the "rigorous" bound; produced by the v2 relative-tolerance test). v3 correctly returns 0; the row is already shown as 0.0000 in Table 12, but the JSON still carries the invalid value |
| 7 | MINOR | `bnb_cert`/`exp14`: `certified_gain = (L − rig)/rig` uses the *lower* bound rig as the multi-exposure reference; the correct reference is an upper bound (ensemble ratio of the verified-feasible X*) |
| 8 | MINOR | `sdp_bound_pixel_int`: `G` omits the half-pixel phase e^{−iπk/N} → models pixels *centred* at nL/N, not edges at nL/N as the docstring says; changes the N < F bounds (N = 12: 0.5407 vs 0.5265) |
| 9 | MINOR | `bnb_sdp.node_certificate`: dual recovery is weak (LB = −inf or < 0 while the primal s* > 0) → sound but needs many nodes; `LB > 0` has no floating-point slack on Σ y c |
| 10 | MINOR | `ModelX`/`Model` accept NA > n_img silently (direction cosines clipped) |
| 11 | OK | RealModel basis U, Re(UᴴQU), DFT conventions, K/pupil cutoff, dual extraction, `rigorize_certificate` v3 bound ≤ ρ(r) on random masks (0 violations), `bnb_cert` interval bounds and final L (0 violations), common-dose metric consistency, `samples_1d` exactness, ModelX channel weights/rotation, ModelX(vector=False) ≡ Model, new sources symmetric, new defocus formula |

---

## 1. CRITICAL — `samples_2d` produces (almost) no edge-offset points for the new tolerances
`certlitho.samples_2d` labels a fine-grid pixel centre as an edge-offset sample iff its EDT distance lies in
`[margin, margin+band)` (`band=0.6`). The EDT is pixel-centre-to-pixel-centre, so along straight edges it takes only the
values `fine, 2·fine, 3·fine, …`. With `fine=2` the ring is non-empty only when τ is (close to) an even integer; for the
session-3 settings (exp4: 7/5/4.4 nm, exp15 2-D: same, exp16 2-D: 5, exp17 random 2-D: 5) the bright ring is **empty**
and the dark ring consists of a few corner pixels:
```
$ python3 -W ignore audit_2dmargin3.py      # and the ring-count snippet in the transcript
square hole P=160 hw=35, fine=2:
 tau=4.4: ring B=  0 ring D=  8 | true dist ring D: [2.]    | nB=9  nD=83
 tau=5  : ring B=  0 ring D=  0                              | nB=9  nD=75
 tau=6  : ring B=116 ring D=144 | true dist B: [6.] D: [4.]  | nB=125 nD=219
 tau=7  : ring B=  0 ring D=  8 | true dist ring D: [4.47]  | nB=9  nD=83
 tau=8  : ring B=108 ring D=148 | true dist B: [8.] D: [5.66 6.]
contacts160_annular tau=7 : nB=9 (interior grid only), min true dist B=19.0, D=4.47 ; bound as coded 0.38519 ; exact-contour tau=7: 0.81467
contacts120_quasar  tau=5 : nB=9, min true dist B=9.0, D=7.0                       ; bound as coded 0.60408 ; exact-contour tau=5: 0.84955
contacts100_quasar  tau=4.4: nB=9, min true dist B=6.0, D=3.16                     ; bound as coded 0.80513 ; exact-contour tau=4.4: 1.0 (unprintable)
EUV_contacts40_annular tau=2 (fine=1): 81% of samples at 1.5 nm (< tau)            ; as coded 0.902 ; exact-contour tau=2: 0.764
```
The bright constraint set degenerates to 9 interior grid points, so the SDP is solving a much weaker problem; the
resulting "certified bounds" are valid only for that degenerate sample set, not for the τ = 7/5/4.4 nm problem. Every
2-D number produced by the running reproduction (exp4, exp13/exp18 with τ = 6 → true dark distance 4 nm, exp15 2-D,
exp16 2-D, exp17 random 2-D and EUV 2-D) is affected.
**Fix.** Sample the offset contours exactly, as `samples_1d` does: compute the signed distance analytically for the target
primitives (union of axis-aligned rectangles / holes → distance to the boundary of the union, periodic), place points
*on* the curves {dist = ±τ} with step ≤ 2 nm (rounded corners included), label them with the exact `target_fn`, add the
interior grid with dist ≥ τ + band, and `assert d.min() >= margin - 1e-9` (as `samples_1d`). A reference implementation
for square holes is `samples_exact` in `code/audit_2dmargin2.py`. If the EDT route is kept: use a sub-pixel distance
(EDT of a 4–8× finer grid, or EDT − fine/2 for axis-aligned edges), restrict τ to multiples of `fine`, and assert that
both rings are non-empty and that every kept sample has exact distance ≥ τ.

## 2. MAJOR — 2-D samples are 1–2 nm closer to the true edge than τ (session-2 numbers)
For τ a multiple of `fine` the ring exists but the pixel-centre EDT plus the pixelised target boundary place samples at
true distance τ − fine … τ:
```
$ python3 -W ignore audit_samples.py
2D holes P=160 hw=35 margin=12 fine=2.0: true dist min B=12.000 min D=10.000 | #samples closer than margin: 144/320
2D holes P=120 hw=25 margin=10 fine=2.0: min B=10.000 min D=7.211 | 113 closer
2D holes P=100 hw=22 margin=8  fine=2.0: min B=7.000  min D=7.000 | 160 closer
2D holes P=40  hw=10 margin=4  fine=1.0: min B=3.500  min D=3.500 | 152 closer      (EUV 2-D, session 2)
random 2D (exp17 geometry, margin 12): 3806/5880 samples closer than margin; worst distance 10.01
$ python3 -W ignore audit_2dmargin2.py   # exact offset-contour sampler, old sources
contacts160_annular tau=12: as coded 0.58434 | exact tau=10: 0.64375, 10.5: 0.61282, 11: 0.58086, 12: 0.52429
contacts120_quasar  tau=10: as coded 0.57395 | exact tau=8: 0.62338, 9: 0.56235, 10: 0.50645
contacts100_quasar  tau=8 : as coded 0.87506 | exact tau=6: 0.92169, 7: 0.87564, 8: 0.83370
```
The as-coded bound equals the exact bound at τ − 1 nm (contacts100: 0.8751 ≈ 0.8756 at τ = 7). Since points at
distance < τ are *not* constraints of the continuous τ-problem, the certificate is **not** a lower bound for the
continuous problem at the stated τ (THEORY §5, PAPER §3 "placed at exactly ±τ"); it is one for τ − fine. Gaps
(bound vs best mask on the same samples) are unaffected. Fix: as in item 1.

## 3. MAJOR — immersion defocus not applied consistently
Old `Model` (session 2) used `ph = 2π z/λ (√(1−min(|g|λ,1)²) − 1)`: for NA = 1.35 (water immersion) 12/49 passing
(order, source) pairs of the conv0.7 1-D case have |g|λ > 1 and received the constant phase −2πz/λ; the paraxial defocus
was 1.44× too strong (n = 1 formula). All session-2 "PW80" rows (Table 1, Table 2 PW80, exp11b, exp14) are for this
non-physical pupil. The session-3 `Model` has the right formula (verified: phase error 0.0 vs
2π n z/λ(√(1−(|g|λ/n)²)−1)) **but** `n_img` defaults to 1.0 and is omitted in:
```
exp14_bnb.py:11  ("mixed_conv0.7_PW80_edgeonly", ..., (0.0, 80.0, -80.0))          # no n_img
exp19_bnb_sdp.py:5, exp11_lasserre_pw.py:4, exp11b_lasserre_pw.py:4              # no n_img
exp15_verify_all.py:27  build_1d(600, ..., 8, spacing=4.0, focus=(0.0,80.0,-80.0))  # no n_img, tau=8
audit_model2.py: n_img=1.0 -> passing pairs with |g|lam > n_img (clipped): 712 (conv0.7, 349 pts); n_img=1.44 -> 0
```
Fix: `Model.__init__`: `if NA > n_img: raise ValueError("NA exceeds image-medium index; pass n_img")` (or default
`n_img = max(1.0, NA)` with a warning); pass `n_img=NIMG_DUV` in every DUV call (make `build_1d/2d` take it from a
single config); store `n_img` in every result record.

## 4. MAJOR — exp15 does not verify the reported instances
`exp15_verify_all.py` builds 1-D cases with τ = 8 (lines 25, 27) and random cases with τ = 10 (line 31) while exp1/exp2
now use `TAU_DUV = 5`; the PW80 rows omit `n_img`. The 2-D rows share the τ = 7/5/4.4 degenerate sample sets of item 1.
Hence "independent verification of every reported certificate" (paper §6.9, exp15) verifies other instances. Fix: save
the instance spec (L, segs/target, src name, τ, spacing, focus, n_img, fine) in each result record and let exp15 rebuild
from `results/*.json`, or import the case lists from the experiment scripts.

## 5. MAJOR — dead orders kept in K
Orders with P(f_k + s/λ) = 0 for every source point contribute to no constraint. They exist for discrete off-axis sources:
```
$ python3 -W ignore audit_model.py / audit_model2.py / audit_misc.py
2D Tjunction320 annular (old src): F=61 dead orders=4 -> [[-3,-3],[-3,3],[3,-3],[3,3]]
2D contacts160 quasar (new src):   F=13 dead orders=4
Tj320: saved Y eigs top6=[1.42453, 0.52974, 0.11654 x4]; rank_eps=6 ; ||Y_dead||=0.233 ; rank_eps on live coords = 2
```
Consequences: (a) Table 6 / exp13 "rank 6" for Tjunction320_annular and lineend320 (0.1259 × 4) is an artefact of the
unconstrained block (true numerical rank 2); (b) SD has an exact null space, so `rigorize_certificate` v3 has to project.
The projected test is sound only if SB vanishes on N; the implemented check `‖VnᴴSB Vn‖ ≤ 1e-12‖SB‖` bounds the cross
term `VkᴴSB Vn` only to √(1e-12)·‖SB‖ = 1e-6‖SB‖ by Cauchy–Schwarz, far above the 1e-10 margin (in practice the block
is exactly zero: resid 3e-31, full-space rel. λ_min −8.7e-18, so no reported bound is wrong). Fix: in `Model.__init__`
drop k with `max over (z, channel) |P| == 0` from `ks` (after computing `inside`), then `null_dim = 0` and the full-space
test applies; recompute ranks.

## 6. MAJOR — degenerate instance with a feasible mask below the stored "rigorous" bound
```
results/exp17_dataset.json  extreme_duty_95pct/coh: t_sdp 6.20e-4, rigorous 1.48e-5 (cert_scale 0.0239), best_feasible 1.29e-24, gap_rel -1
$ python3 -W ignore audit_cert.py  (v3, same instance): t 6.20e-4, rig 0.0, c 0.0, resid_null 1.3e-7, best_local 1.86e-5
```
CLARABEL reports `optimal` at t = 6.2e-4 although masks with ratio ≈ 0 exist (true t_R ≈ 0): the SDP solution is
wrong and the dual is not a certificate. v2 (relative tolerance −1e-11 on λ_min) accepted a scaled version → 1.48e-5,
an invalid "rigorous" bound (a −ε eigenvalue is amplified by |M|²/β, unbounded when β → 0). v3 returns 0 (correct,
conservative); Table 12 already shows 0.0000 but the JSON and any downstream statistics still carry 1.48e-5 and
`gap_rel = -1`. Fix: regenerate exp17 with v3 (in progress), never report a bound whose certificate has λ_min < 0, flag
`status != optimal` and degenerate rows (`analyze_el` now does the latter), and treat `best_feasible < rigorous` as a
hard assertion failure in `analyze`.

## 7. MINOR — B&B "certified gain" reference
`exp14_bnb.py`: `certified_gain = (L − rig)/rig`. L is a lower bound on the single-mask optimum (fine), but the
multi-exposure reference must be an *upper* bound on t_R, i.e. the ratio of a verified-feasible ensemble
X* (`max_d Tr(Q_d X*)/min_b Tr(Q_b X*)`), not the rigorous lower bound. Here they agree to 1e-6 (trial82: rig 0.359679,
ensemble 0.359680), so no number changes, but the definition should be fixed: `gain = (L − t_ens)/t_ens`.

## 8. MINOR — `sdp_bound_pixel_int` sinc convention
```
$ python3 -W ignore audit_misc.py
pixel_int G: err vs exact(edges at nL/N) 2.2e-01; err vs exact(centers at nL/N) 2.3e-16
N=15: 0.46907 / 0.46907 ; N=12: as coded 0.54074, edges-at-nL/N 0.52651 ; N=9: 0.68637 vs 0.67480
```
The exact coefficient of a mask piecewise constant on [nL/N,(n+1)L/N) is `sinc(k/N) e^{−iπk/N} (1/N) Σ m_n e^{−2πikn/N}`;
the code omits `e^{−iπk/N}`, i.e. models pixels centred at nL/N. Lemma 8 (N ≥ F ⇒ equal to Fourier bound) is unaffected;
the N < F "pixelation" rows depend on grid alignment. Fix: multiply `G` by `np.exp(-1j*np.pi*k/N)[:, None]` or state
the centre convention.

## 9. MINOR — `bnb_sdp.py` (new SDP-bounded B&B) is sound but weak
Checked: RLT cut matrices and constants (all four products) are algebraically correct; weak-duality bound
`LB = Σ y c − α' − β'` is valid for *any* (α', β') with Λ ⪰ 0, so trying sign combinations is legitimate; restricting to
active (d, b) pairs only enlarges the set proved empty (valid); heap/budget logic cannot certify with open or failed
nodes; Λ negativity is compensated by `slack·(1+F)` using Tr Z ≤ 1+F (valid for lifted box points).
```
$ python3 -W ignore audit_bnbsdp2.py   (F=5 trial82, tiny boxes around the best mask r*, rho(r*) = U = 0.406553)
box hw=1e-3: T=1.01U: LB=-inf verified=False s*=6e-11   | T=0.995·min rho: LB=-1.6e-03 verified=True  s*=1.7e-03
box hw=1e-2: T=1.01U: LB=-inf verified=False s*=1e-10   | T=0.995·min rho: LB=+2.0e-04 verified=True  s*=1.7e-03
box hw=5e-2: T=1.01U: LB=-inf verified=False s*=5e-12   | T=0.995·min rho: LB=-inf     verified=False s*=1.2e-03
$ python3 -W ignore audit_bnbsdp.py   (certify_target, 400 nodes): nothing certified for T = 0.98 L_bnb, 1.02 U, 0.9 rig
```
No unsound certificate was produced, but the recovered dual is often far from the primal s* (LB ≤ 0 while s* > 0), so
nodes the SDP does certify are split instead → budget exhaustion (exp19 will likely report "not certified"). Fix: for the
extracted y solve the 2-variable problem max −α−β s.t. Σ y G + αE₀₀ + βE_ii ⪰ 0 (tiny SDP, or bisection on β after
setting α by λ_min), or solve the dual SDP explicitly; add a floating-point slack to `LB > 0`
(`LB > 1e-12·(1+Σ|y c|)`).

## 10. MINOR — missing guards
`ModelX`/`Model` accept `n_img < NA` (direction cosines a²+b² > 1, γ clipped to √1e-9) without warning;
`rank_eps` uses an IPM (maximal-rank) solution, so published ranks should come from `face_min_rank` after item 5;
`samples_2d` labels/distances coarse grid *corners* by the containing pixel *centre* (off by ≤ 1.4 nm; subsumed by item 1).

## 11. OK — verified items (evidence)
```
audit_model.py / audit_model2.py
  RealModel U: unitary err 2.2e-16 (1-D F=11, 2-D F=13, new model F=11 with 349 src pts); U r Hermitian-symmetric;
  r^T Re(U^H Q U) r == M^H Q M to 1e-14 (random r, all foci); real pixel masks lie in span(U) (round trip 1e-12)
  orders_from_pixels == (1/L)∫ m e^{-2πi k x/L} (1.4e-17); pixels_from_orders∘orders_from_pixels = id (1-D 3.9e-16, 2-D 2.0e-16)
  translation covariance of the Abbe image (shifted M vs shifted x): 8e-15
  dead orders: 0 for 1-D and old annular contacts (see item 5 for exceptions)
  session-3 defocus phase == 2π n z/λ (√(1-(|g|λ/n)²)-1) exactly; n_img=1.44 -> 0 clipped pairs
  ModelX(vector=False) == Model: C err 0.0; vector channels: Σ_c |E_c|² == |P|² (≤1.1e-15) for x,y,unpol,TE,TM
  (Rodrigues rotation unitary); channel weights w_s/|pols| per component, Σw = 3 (=components) as intended
  Lemma 6: real == complex relaxation for Z=-Z (0.519791 both); differs for Z={0,50} (0.519791 vs 0.519631) as the
  revised lemma states; all catalogue sources symmetric, weights uniform (area quadrature), smax/NA correct
audit_cert.py (old model; 8 instances: 1-D coh/conv/dipole, PW80 x2, 95%-duty, 2-D contacts quasar, contacts F2, Tj320)
  Σμ = 1 ± 4e-6, Σλ = t ± 3e-6 (dual mapping/sign correct); c ≥ 0.99999 except degenerate 95%-duty (c = 0, bound 0)
  rigorized bound ≤ ρ(r) for 3000 random real masks per instance (gaussian, sparse, Y-rounded) + local optimum: 0 violations
  Tj320: null_dim 4, resid 3.2e-31, full-space rel λ_min -8.7e-18 (projection harmless here)
audit_cert3.py (session-3 model, τ=5, n_img=1.44, 349-point sources): see appendix line below
audit_bnb.py (bnb_cert.py, F=5 coherent, incl. 3 foci)
  box bounds q(c)-2|Qc|·h ≤ q ≤ q(c)+2|Qc|·h+h^T|Q|h: 0 violations / 30000 random (box, point) pairs; fixed coordinate h=0
  L ≤ ρ(r) for 20000 random masks at budgets 10/500/20000/∞ (0 violations); L ≤ U; converged: U−L ≤ tol·U;
  no box lost (popped boxes are either leaves or replaced by both children; heap min = min active LB)
audit_misc.py
  common-dose metric: exact_minimax_polish ratio 1.194053 == max_{d,z} I / min_{b,z} I; per-focus ratios equalised;
  local_refine_real objective = log max_D I − log min_B I over all foci (common dose) — consistent with rm.ratio and the SDP
audit_samples.py: samples_1d — all 2|E| edge-offset points present at exactly ±τ, min distance = τ, labels correct,
  wrap-around correct (segments touching 0; 95% duty; near-cutoff mixed)
```
Informational: the running reproduction's `exp1_1d_features.py` failed at 07:27 (`run() got an unexpected keyword
argument 'n_img'`) and was rerun manually after a signature fix (results/logs/run_all.log).

### Appendix — session-3 model certificate check (`audit_cert3.py`, catalogue sources, n_img = 1.44, τ = 5)
```
isodense_conv      F=15 t=0.590411 rig=0.590409 c=0.999997 null_dim=2 (dead orders also appear in 1-D with the 349-pt disk source!)  rank_eps=6  worst random/local mask 0.598371  violations 0/400
mixed_coh          F=7  t=0.421083 rig=0.421081 c=0.999994 null_dim=0  worst 0.485018  violations 0/400
mixed_conv_PW80    F=11 t=0.739341 rig=0.739340 c=0.999998 null_dim=0  worst 0.759025  violations 0/400   (immersion defocus; old clipped model gave 0.8691)
duty95_coh         F=5  t=6.20e-4  rig=0.0      c=0        (invalid SDP solution correctly rejected; feasible mask 1.86e-5)
contacts120_quasar F=9  t=0.602104 rig=0.602103 c=1.0000   null_dim=0  worst 0.602115  violations 0/400   (τ=5 as coded = degenerate sample set, item 1)
```
Note for item 5: with the new area-quadrature sources, dead orders (and hence `null_dim > 0`, inflated `rank_eps`) now
occur in 1-D too (isodense_conv: null_dim 2, rank_eps 6) — dropping dead orders from K is the right fix for both.
