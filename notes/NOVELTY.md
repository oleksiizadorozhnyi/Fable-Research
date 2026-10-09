# Novelty assessment (Direction A)
Literature access: WebSearch snippets only (publisher sites, arXiv, patents blocked by egress proxy). Coverage therefore
incomplete; full texts of the closest prior work were NOT read. Calibrated language used throughout.

## Claim 1 — Lifted (SDP) relaxation of Hopkins/Abbe mask synthesis with independently checkable dual certificates
Closest prior work:
- IBM/Rosenbluth et al. (JM3 1(1) 2002; SPIE 6520 2007; patents US7587702, US7685559, US8495528, US7944545): mask
  synthesis posed on diffraction orders ("wavefront") with quadratic intensity constraints; recognized as a problem
  "with non-convex quadratic constraints"; attacked by eigen-subspace ("joint eigenvectors", "step-walk relaxation")
  heuristics and local optimizers; multi-exposure version via unconstrained Rayleigh-quotient simplification.
- Luo, Ma, So, Ye, Zhang (IEEE SPM 2010); Sidiropoulos, Davidson, Luo (IEEE TSP 2006): SDP relaxation of QCQPs /
  multicast beamforming (max-min with rank-one quadratic forms) — same mathematical template, different field.
- Angeris, Vučković, Boyd (ACS Photonics 2019; Opt. Express 2021; ACS Photonics 2023): Lagrangian-dual / SDP bounds
  for photonic inverse design — same philosophy, different physics (Helmholtz constraints vs imaging quadratic forms).
- Ma, Arce, Garcia-Frias et al. (Opt. Express 2017, 2018; 2021; 2022; 2025): information-theoretic "limits" on OPC
  fidelity (mutual information based). Not duality certificates; not pattern-specific feasibility bounds.
Technical difference: (i) explicit PSD lifting in the diffraction-order domain, (ii) dual certificate theorem with a
one-PSD-test verification independent of the solver, (iii) rigorization by bisection scaling, (iv) real-mask structure
handled exactly by a real subspace (Lemma 2), (v) process-window (multi-focus) and source-optimization variants.
Potential overlap: an SDP/"lifting" treatment may exist in SPIE proceedings or theses not surfaced by search (e.g.,
work citing Rosenbluth with "convex relaxation" wording). Not ruled out.
Verification coverage: >25 targeted WebSearch queries (terms: semidefinite/SDP/convex relaxation/lifting/PhaseLift/
QCQP/rank-one/Toeplitz/multiple exposure/incoherent/certificate/lower bound × ILT/OPC/mask/wavefront/Hopkins/SMO).
Status: potentially novel; not found in the literature reviewed.

## Claim 2 — Relaxation = unlimited incoherent multi-exposure (Theorem 3); Pataki bound on useful exposure count
Closest: US8495528 (multi-exposure wavefront optimization, nonconvex); general SDP folklore (X PSD = Gram of vectors).
Difference: the identification of the relaxation gap with the physical benefit of incoherent multiple exposure, and the
use of Lasserre-2 single-mask bounds to CERTIFY strictly positive multi-exposure gains on specific targets.
Status: elementary mathematically; interpretation appears new in lithography; not found in the literature reviewed.

## Claim 3 — Exact convex (Toeplitz-PSD) description of in-band spectra of nonnegative 1D masks (Theorem 4)
Closest: Carathéodory–Toeplitz / Fejér–Riesz theorems (classical); Rosenbluth US8108802 (arbitrary wavefront
realizability by PSM); patents selecting quantized transmissions to match target orders. Status: the mathematics is
classical; its application as a realizability cone for binary/attenuated masks was not found in the literature reviewed.

## Claim 4 — Empirical findings (certified near-optimality of local optima; certified PSM-vs-binary and multi-exposure gains;
pattern-specific k1 curves)
Closest: NVIDIA arXiv 2602.19027 and ILT reviews discuss local minima qualitatively; no certified gap measurements found.
Status: new measurements within a simplified (scalar, thin-mask, periodic) model.
