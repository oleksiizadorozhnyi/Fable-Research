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
