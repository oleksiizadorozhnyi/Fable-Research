# FINAL SCIENTIFIC DELIVERABLE — Autonomous discovery run, 2026-10-09
Classification: **VALIDATED RESEARCH CANDIDATE** (internal verification only; not peer reviewed; see Gate status §0).

## 0 Completion-gate status
| Gate | Status | Evidence |
|---|---|---|
| A Precision | PASS | notes/THEORY.md §1–2: problem, classes, relaxation, certificate defined exactly |
| B Correctness | PASS (model-internal) | Theorems 1–6 proved (elementary/classical); certificates verified by independent PSD test, bisection rigorization, 50-digit Cholesky; two solvers agree to 1e-8; 6 regression tests pass |
| C Originality | PASS with caveats | notes/NOVELTY.md: no lifted/SDP certificate for Hopkins mask synthesis found in >25 targeted searches; closest: IBM wavefront heuristics, Luo et al. SDR, Angeris–Boyd photonic bounds, Ma–Arce information limits. Full texts inaccessible (egress) ⇒ residual risk |
| D Significance | PASS (medium) | Pattern-specific printability certificates, ILT optimality-gap benchmarking, certified multi-exposure/PSM gains; scalar periodic model limits direct fab relevance |
| E Reproducibility | PASS | README.md, code/, results/*.json, logs, tests, figures |
| F Adversarial review | PASS with documented limitations | notes/ADVERSARIAL.md; bug log in RESEARCH_STATE.md (edge-detector bug found and fixed; ILT gradient bug found and fixed; all affected results rerun) |

## PART I — Executive research summary
Problem: ILT/SMO are nonconvex; nobody can tell an unprintable layout from a stuck optimizer. Contribution: for the
scalar Abbe/Hopkins model of a periodic cell, the maximum exposure latitude (min dark/bright ratio) over all masks is
bounded by a small SDP in the diffraction-order domain whose dual is a one-matrix certificate (Theorem 1). The
relaxation is exactly the unlimited-multi-exposure problem (Theorem 3), nonnegative 1D masks have an exact convex
realizability cone (Theorem 4), complex phase cannot improve the relaxation under symmetric illumination (Lemma 6), and
the classical k1 cutoffs are reproduced exactly (Proposition 5). Empirically: local optimization is certified globally
optimal on most targets (75% of random 1D targets exactly; all 2D contact layouts within 0.3%); certified
multi-exposure gains of 13.7% and ≥27% exist for coherent illumination; binary masks are certified ≥10% worse than PSM
on an iso-dense target; standard pattern-error ILT is 30–67% above the certified optimum; the SDP bound is 3–150×
tighter than the Rayleigh-quotient bound underlying prior global heuristics.

## PART II — Five research directions
See notes/PHASE2_DIRECTIONS.md (A: certified mask synthesis bounds — selected, weighted score 7.75; B: shot-noise
capacity of EUV patterning 5.05; C: TCC eigenvalue-decay theorem 6.15; D: superposition theory 5.25; E: BTSP credit
assignment theory 5.30). Rationale: A had the strongest verification pathway with the tools available and a clear gap.

## PART III — Literature and prior art
notes/NOVELTY.md and manuscript §2 and References. Limitation: all bibliographic verification via search snippets.

## PART IV — Formal contribution
notes/THEORY.md (Definitions; Theorem 1 certificate; Lemma 2 real structure; Theorem 3 multi-exposure equivalence +
Pataki; Theorem 4 Carathéodory–Toeplitz realizability; Proposition 5 cutoff; Lemma 6 complex=real; §5 sampling validity).
Open proof obligation: Theorem 4's "bounded function" refinement relies on the standard positive-definite extension
argument (Fejér–Riesz); stated as sketch. Lasserre-2 single-mask bounds are numerical (solver status optimal_inaccurate).

## PART V — Validation
results/RESULTS_TABLES.md (Tables 1–7, robustness), results/figs/*.png, results/logs/*. Summary in manuscript §6.

## PART VI — Adversarial analysis
notes/ADVERSARIAL.md. Key failed/retracted items: (1) initial 1D results (pre-fix) were invalid due to the wrap-around
edge bug and were discarded; (2) first ILT baseline results were invalid (gradient sign, normalization) and discarded;
(3) relaxation is not always tight (documented, physically interpreted); (4) grid-sampled bounds were discretization-
sensitive, replaced by exact edge-offset sampling (bound invariant to interior spacing).

## PART VII — Novelty assessment
notes/NOVELTY.md: Claims 1–4 "potentially novel; not found in the literature reviewed"; mathematics of Claims 2–4 is
classical, their lithographic use is what appears new; residual overlap risk with unseen SPIE/thesis literature.

## PART VIII — Manuscript
notes/PAPER_DRAFT.md.

## PART IX — Reproducibility package
README.md; code/certlitho.py, pipeline.py, exp*.py, test_certlitho.py, assemble.py, figs.py, fig2d.py; results/.

## PART X — Next scientific steps
1. Rigorize Lasserre-2 bounds (extract SOS dual certificate; or exact rational verification).
2. Vector imaging + 3D mask (rigorous EM) — the lifting applies unchanged to any model quadratic in mask orders; thick-
   mask models are not quadratic and need a different treatment.
3. Non-periodic clips via windowed/overlapping cells; scaling to F ~ 10^3–10^4 with Burer–Monteiro or chordal SDP.
4. Resist/NILS/MEEF constraints (linear in intensity derivatives — lift identically); EPE-based objectives.
5. Benchmark commercial/academic ILT engines (ICCAD-13 clips) against certificates; use certificates as stopping rules.
6. Theory: characterize targets for which the relaxation is tight (rank-one); connection to the multicast-beamforming
   literature's approximation bounds; 2D realizability cones for binary masks (sum-of-squares hierarchies).
7. Independent replication and submission (JM3 / Optics Express / SIAM Imaging).
