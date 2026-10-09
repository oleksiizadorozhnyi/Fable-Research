# Persistent Research State (autonomous discovery agent)
Started: 2026-10-09. Repo was empty at start.

## Tooling available
python3.13, numpy, scipy, sympy, cvxpy (SCS, CLARABEL), matplotlib; WebSearch/WebFetch. No GPU.

## Phase I literature sweep — log
- Q1 "semidefinite relaxation inverse lithography": NO hits for SDP relaxation of ILT mask synthesis. Hits only: SDP for multiple-patterning layout decomposition (Kuang/Young; Tokyo Tech), ADMM for ILT (arXiv 2209.10814), general SDP-relaxation theory.
- Q2 ILT open problems: reviews (PMC12287475 2025; PMC12843865; Pang JM3 2021 "ILT 30 years") stress nonconvexity, local minima, runtime; none mention optimality certificates/lower bounds.
- Q3 PhaseLift/Hopkins: PhaseLift (Candès, Strohmer, Voroninski 2011/2013), coded diffraction (Candès, Li, Soltanolkotabi 2014) — no Hopkins/TCC lifting found.
- Q4 EUV stochastic info-theoretic limit: nothing found framing patterning as channel capacity; shot-noise numbers (6.8 photons/nm^2 at 10 mJ/cm^2).
- Q5 EP/FF theory: Directed EP (MDPI Math 13(11) 2025), Hopfield-ResNet EP (arXiv 2509.26003), Laborieux 2021. Guarantees still partial.
- Q6 superposition toy models: arXiv 2608.27540 "Towards a mathematical theory of superposition" (Aug 2026) — must read; arXiv 2310.06301.
- Q8 SMO patents: IBM (Rosenbluth) patents 8730452/10657641: source step is LP; single-direction illumination -> nonconvex QP. MUST check Rosenbluth 2002 JM3 "Optimum mask and source patterns to print a given shape" for lifting/SDP use.

## Candidate directions (draft)
A. Lifted (SDP) relaxation of Hopkins-model mask synthesis -> provable lower bounds / infeasibility certificates for target patterns; BM low-rank solver.
B. Information-theoretic capacity of stochastic (Poisson) EUV pattern transfer.
C. Theory of superposition geometry (AI interpretability) — check 2608.27540 first.
D. EP/local learning convergence theory.
E. TBD (neuro/math).

## Next action
Prior-art deep-dive on Rosenbluth global mask optimization; QCQP/lifting in OPC; then score candidates.

## Tool limitation (recorded 2026-10-09)
Egress proxy blocks arxiv.org (and mirrors), PMC, IEEE Xplore, SPIE, Google Patents, USPTO, Semantic Scholar, OpenAlex, Stanford, UCLA, ResearchGate, publisher sites. Only pypi/github reachable by curl/WebFetch.
=> Literature verification is via WebSearch result summaries/quotes only; full texts cannot be read. All bibliographic entries must be marked "verified from search metadata" vs "full text read" (none).

## Checkpoint 2026-10-09 ~21:30 (selected direction A under investigation)
Code: code/certlitho.py (model, SDP bounds real/complex/pixel/DNN/SMO, certificates, rigorize, rounding, local solvers,
Lasserre-2), code/pipeline.py, code/exp1..7. Results in results/*.json, logs in results/logs.
Verified so far:
- Prop.5: Abbe/k1 cutoffs reproduced exactly (DUV coherent 143nm; sigma=0.7 -> 84.1nm; EUV k1=0.244 -> bound=1).
- Theorem 1 certificates verified independently (min eig >= -1e-8 before rigorization; rigorous after bisection scaling).
- Theorem 4 numerics: necessity 0/2000 violations; sufficiency 200/200 realizable by LP.
- Tightness: majority of random 1D instances rank-1 & zero gap. Genuine gaps exist (coherent iso-dense 1.63%, coherent
  mixed 2.78%, dipole random trial-7 2.3%, EUV conv0.9 pitch 26 1.6%). Lasserre-2 single-mask bound == best single mask
  found (to 4e-5) => multi-exposure gain CERTIFIED for mixed_coh (2.78%) and isodense_coh (1.63%).
- SMO entitlement bound (per-source-point PSD blocks) == coherent bound for isodense (0.5645): coherent+PSM is optimal there.
- 2D: contacts160 annular: F=13, bound 0.2324, rank 1, rounding attains bound (tight), SDP 0.3 s.
Open: exp3 (pixel vs Fourier, DNN binary, ILT baseline) rerunning after bug fix; exp4 2D battery status unclear.
Next: finish 2D battery; figures; paper draft; adversarial review (sampling, thin-mask, scalar, metric choice).

## Bug log (honesty record)
- 2026-10-09 ~22:10: Found operator-precedence bug in sets_1d (wrap-around edge placed at mid-cell). Effect: for targets
  with a feature touching x=0, samples adjacent to a true edge were included -> spuriously tight constraint pair ->
  bound strongly sampling-dependent (0.32 -> 0.90 as nx 24 -> 384). Random-target study (exp2) unaffected (no feature
  touches the cell boundary). FIX: replaced grid+margin sampling by exact edge-offset sampling (samples_1d/build_1d):
  constraints at exactly edge +- margin plus an interior grid. Verified: bound invariant to interior spacing
  (0.46907 for spacing 16..1 nm) and monotone in margin. ALL 1D experiments rerun (run_all_1d.sh). Pre-fix numbers
  in this file above are superseded.
- ILT baseline threshold changed from adaptive midpoint to fixed 0.4 of max normalized intensity (standard).
