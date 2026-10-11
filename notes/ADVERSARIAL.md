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
