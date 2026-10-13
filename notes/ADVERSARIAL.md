# Scientific red-team review (self-administered)

## Strongest rejection argument
"This is a straightforward application of the textbook Shor SDP relaxation (Luo et al. 2010) to a QCQP that IBM
already wrote down in 2002. The physics model is a scalar thin-mask periodic toy; the 'ratio' metric is not what fabs
optimize (EPE/PV-band); the bounds are only as good as the sample-point discretization; nothing is validated on a real
imaging tool or against commercial ILT; and the 'theorems' are either weak duality or classical moment theory."

## Responses and what remains open
1. Textbook relaxation? — Yes, the relaxation itself is standard, and we say so. The contribution is (a) the
   certificate theorem in the lithography setting with solver-independent verification, (b) the multi-exposure
   identification (Theorem 3) with CERTIFIED strictly positive gains via Lasserre-2, (c) the realizability cone, (d) the
   first (to our knowledge) measurements of global-optimality gaps of local mask optimization against certified bounds,
   (e) pattern-specific exposure-latitude limits generalizing the k1 rule. We do not claim algorithmic novelty of SDP.
2. Toy physics — The scalar thin-mask Abbe model is the standard ILT model (SOCS/Hopkins). Vector/3D-mask/resist effects
   are outside scope; the bound's VALIDITY for the model is exact; its RELEVANCE to silicon is model-limited. Stated.
3. Metric — dark/bright ratio ⇔ exposure latitude for a threshold resist; Rosenbluth 2002 also maximizes exposure
   latitude. Any constraint linear in intensity (slope/NILS, multi-focus, multiple dose corners) lifts identically.
   EPE-type objectives are nonlinear in intensity only through the threshold crossing; with sampled edge-offset
   constraints they are captured by the margin parameter. Stated as a limitation.
4. Discretization — Fixed by exact edge-offset sampling; bound invariant to interior spacing; a certificate for the
   sampled problem is valid for the continuous one (Theory §5). The margin is the physical edge tolerance.
5. Numerical rigor — certificate = one PSD test; bisection rigorization; 50-digit Cholesky check; CLARABEL/SCS agree 1e-8.
   Lasserre-2 bounds are reported with solver status (optimal_inaccurate) and are NOT bisection-rigorized (open).
6. Scale — periodic cells up to ~1 µm are tractable (F ≈ 10^2–10^3); full-chip is out of scope; the method targets
   clip-level entitlement analysis and benchmarking of ILT engines. Stated.
7. Could a simpler method give the same? — Rayleigh-quotient (generalized eigenvalue) bounds (IBM's unconstrained
   simplification) are weaker than the SDP (they drop the min/max structure); checked numerically (TODO: include).
8. Hidden counterexamples — the relaxation is NOT always tight (documented rank-2 cases); the certificate is still valid.
   Binary (DNN) relaxation in 2D is only an outer approximation (no 2D Fejér–Riesz). Stated.

9. Large-2D gaps (9–14% on T-junction) are NOT certified as multi-exposure gains: Lasserre-2 is out of reach at F=61 and
   the primal search omitted SLSQP polish; they are reported as unresolved upper bounds on the gap. (Open.)
10. Rigorization initially failed on layouts whose weighted dark matrix is singular (orders invisible to every
   constraint); fixed by verifying on the complement and reporting the residual (1e-31). Documented in the bug log.

## Session-2 review items (2026-10-11)
11. Stale-data risk: two result files (exp7_lasserre.json, exp9_outlier_pw.json) were found to carry pre-fix numbers while
    the manuscript quoted rerun logs. Both regenerated from the correct logs; the final deliverable is produced by a single
    clean run (code/run_all.sh) so that tables, figures and text derive from one consistent set of outputs.
12. Solver reliability: on near-degenerate random instances (bound ≈ 1e-3) the first-order solver SCS returned values off
    by up to 0.14 (trial 103) and even negative; CLARABEL + independent certificate is used throughout; SCS only as fallback,
    and every reported bound is re-verified by the 50-digit Cholesky check (exp15). Two EUV instances with bound exactly 1
    have S = 0 (uniform certificate of Prop. 5), where Cholesky is undefined; they are reported as degenerate.
13. Lasserre-2 "optimal_inaccurate": no longer load-bearing. Multi-exposure gains are certified by Theorem 9 (B&B), which
    is rigorous up to floating-point rounding of order 1e-12. Lasserre-2 is reported as corroboration only.
14. 2D 320 nm gaps: addressed by exp13 (optimal-face exploration, 24-start strong primal, explicit two-exposure
    realization). Whatever remains is reported as an unresolved upper bound on the multi-exposure benefit.
15. Vector model: conclusions re-tested (Table 10–11): tightness statistics unchanged; TM polarization raises bounds
    substantially, as expected physically; aberrations up to 0.1λ change bounds by < 7%; pixelation below the collected-band
    count raises bounds (Lemma 8).

## Session-3 review items (2026-10-12): six independent adversarial reviews
Reviews (notes/review_*.md): mathematics, numerical provenance, lithographic physics, prior art, editing, code audit.
16. Mathematics (review_math.md). MAJOR items fixed: Pataki count for Hermitian blocks is r² (not r(r+1)/2); attainment
    of (R) now argued (null-space reduction + closed feasible set); the Fejér–Riesz "positive trigonometric polynomial"
    claim in Theorem 4 replaced by the maximum-entropy AR (rational) extension; Proposition 5 stated as a sufficient
    condition only, with a counter-example to the converse; Lemma 6 hypothesis corrected to P(−g) = conj P(g), Z = −Z
    (even aberrations break it); "rigorous" replaced by "certified up to floating-point rounding" everywhere; the
    pixel-grid nonnegative-mask claim restricted to 5-nm pixelated masks.
17. Physics (review_physics.md). CRITICAL/MAJOR fixed: the 1D "conventional" source was a 7-point slit, not a disk
    (bounds changed by up to 25%; all sources are now 2D area-weighted quadratures with a convergence study in exp8);
    defocus phase computed in air was unphysical for NA 1.35 immersion (now n = 1.44 water; EUV in vacuum); EL was the
    half-range value (now the full-range 2(1−ρ)/(1+ρ)); τ was 8 nm (16% of CD), now 10% of the minimum CD; the EUV
    22-nm "EL" was re-evaluated with the disk quadrature (0.36 → much smaller) and the EUV section carries the thin-mask
    caveats (shadowing, 3D mask, flare, stochastics); the "coherent entitlement" statement is now explicitly about the
    multi-exposure relaxation; the 20-nm-space bound ≈ 0 is explained as a sampled-constraint statement.
18. Numbers (review_numbers.md). CRITICAL: several text numbers were typed by hand from stale logs. Fix: every number
    and table body in main.tex is a macro generated by code/make_numbers.py from the JSON of one clean reproduction.
19. Prior art (review_priorart.md). Nearest neighbours added to Related Work and refs.bib (Nouralishahi–Wu–Vandenberghe
    2008 SDP calibration; Nonaka et al. 2025 0–1 QCQP; Kohira/Matsui PSD relaxations for multiple patterning; IBM
    wavefront-optimization patents; Liu–Zakhor B&B; Davidson–Luo–Sturm Toeplitz cone; Burer–Vandenbussche SDP-B&B).
    Absolute novelty claims replaced by search-scoped wording ("to our knowledge, based on titles/abstracts/patents").
20. Editor (review_editor.md). Self-contained captions, notation table, consistent symbols, booktabs; adopted in v3.
21. Code audit (review_code.md). See that file; findings fixed in this session are listed in RESEARCH_STATE.md
    (exp1 signature crash in the clean run; exp15 verifying instances with the old τ and air defocus; exp11/14/19
    aligned; make_numbers macro names and table bodies).
22. Still open after session 3: quasar T-junction gap (optimizer shortfall or multi-exposure gain, unresolved);
    interval B&B informative only to F = 7; SDP-B&B certifies targets only away from the incumbent; full texts of
    Rosenbluth 2002/2007, Liu–Zakhor 1992, Nonaka 2025 and Nouralishahi 2008 were not accessible.

## Referee report (2026-10-12, "major revision") — items 23–34
23. Rounding safety: wording corrected; exact 50-digit rebuild of every certificate from the model equations (exp21):
    143/143 positive definite with margin 1e-30; one instance (EUV cutoff, near-zero certificate) needed an exact kappa
    bisection (bound change < 1e-15). 24. Theorem 4: measure vs bounded-function distinction with the point-mass example.
25. K = -K enforced and explained. 26. Proposition 5: exact 1D criterion; radial condition only sufficient.
27. PSM penalty and ILT suboptimality now bracketed with feasible upper bounds (exp22). 28. Relaxation gap is an upper
    bound on the multi-exposure benefit; strict benefit only from B&B. 29. Attainment lemma (reduction to the complement
    of the common dark null space; dual Slater); single-mask attainment proved. 30. SDP-B&B thresholds in the
    informative range (none certified within 6,000 nodes; reported). 31. Geometry/caption/denominator inconsistencies
    fixed. 32. EL brackets reported. 33. Source quadrature converged (401^2 grid; residual ~1e-4); discretization scope
    stated. 34. Reproducibility: repository, per-instance archives (results/certs), versions, seeds; EUV cell stated;
    prior-art claim narrowed; reference [2] corrected.
