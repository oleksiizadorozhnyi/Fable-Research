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
- 2026-10-09 ~23:20: Rigorization returned bound 0 on the three 320-nm 2D layouts (F=61). Cause: the weighted dark
  matrix SD has an exact null space (4 orders invisible to all constraints), so a strict eigenvalue margin can never be
  met; the bisection collapsed to c=0 although the SDP value (0.4387) was solver-consistent (CLARABEL/SCS agree 1e-7).
  FIX: null-space-aware rigorization: verify PSD with margin on the complement of the near-null space of SD and report
  the residual of SB on that space (1e-31 here). 1D results unaffected (null dim 0). Certificates recomputed.

## FINAL CHECKPOINT 2026-10-09 ~23:40
All experiments complete (exp1–exp12, exp3b, exp4b, exp11b); tables in results/RESULTS_TABLES.md; figures in results/figs;
manuscript notes/PAPER_DRAFT.md; report notes/FINAL_REPORT.md. Tests: 6/6 pass.
Continuation instruction (if resumed): remaining open items are (1) rigorize Lasserre-2 bounds via SOS dual extraction,
(2) resolve the 2D T-junction gaps (9–14%) with face exploration + stronger primal or a scalable Lasserre, (3) vector/
3D-mask models, (4) non-periodic clips. No experiment needs rerunning.

## Session 2 (2026-10-11) — hardening pass
- Found stale results/exp7_lasserre.json (pre-edge-fix values 0.6235/0.6408) feeding Table 2 while the manuscript used
  the rerun log values (0.4211/0.4787, 0.3625/0.4608). Regenerated the JSON from the rerun log; clean rerun scheduled.
- Found that results/exp9_outlier_pw.json had been rebuilt from the STALE pre-fix log (exp9.log, 21:30) instead of the
  rerun log (exp9_outlier_pw.log, 22:25-22:37). The manuscript's "7.8–16.2% certified gains on four random targets"
  therefore referred to pre-fix instances. Rebuilt from the correct log; numbers re-derived below and in the manuscript.
- 2026-10-12: BUG (multi-focus metric). RealModel.ratio/Model.ratio and the L-BFGS surrogate evaluated the ratio PER
  FOCUS (max over z of max_d I/min_b I), i.e. allowed a different dose at each focus, while the SDP and the SLSQP polish
  use the common-dose ratio max_{d,z} I / min_{b,z} I. Detected by "feasible masks" below certified bounds (EUV 32 nm
  contacts, −1.5/−1.9%). Certified bounds unaffected; all multi-focus "best single mask" values (PW80 rows, EUV Tables
  5/13, exp11b/exp14 PW) were optimistic and are recomputed by the clean rerun. Fixed in certlitho.py and bnb_cert.py;
  regression test test_multifocus_metric added.
- 2026-10-12: BUG (rigorization). The null-space-projected certificate check could accept an invalid certificate when
  the weighted bright matrix did not vanish on the near-null space of the weighted dark matrix (95%-duty coherent case:
  reported 6e-4, true bound ≈ 1.5e-5). Replaced by a full-space eigenvalue test with relative tolerance 1e-11 (bisection
  drives c -> 0 when needed). All certificates recomputed by the clean rerun (exp15).

## Session-2 checkpoint (2026-10-12, after clean reproduction)
Clean run (code/run_all_par.sh) completed: 20 scripts, all exit 0 (exp15 re-run separately with rigorization v3).
Deliverables: paper/main.tex + refs.bib (LaTeX), notes/PAPER_v2.md (Markdown), notes/FINAL_REPORT.md (v2),
notes/THEORY.md (Theorems 1–9), notes/NOVELTY.md, notes/ADVERSARIAL.md, results/RESULTS_TABLES.md, results/figs/fig1–6.
Open: exp15 v3 statistics to be inserted; TeX compilation (no engine); see FINAL_REPORT §11.
