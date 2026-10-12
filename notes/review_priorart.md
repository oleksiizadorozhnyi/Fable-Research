# Adversarial prior-art review — Direction A novelty claims

Date: 2026-10-12. Reviewer role: adversarial prior-art investigator (independent of the authoring session).
Method: 47 distinct WebSearch queries (extended mode), English + Chinese + Japanese + Korean, across lithography terms
(mask optimization, ILT, OPC, SMO, wavefront optimization, aerial image, Hopkins/TCC, diffraction orders, exposure
latitude/NILS, multiple exposure, printability), optimization terms (semidefinite/SDP, convex relaxation, Shor,
lifting/PhaseLift, rank-one, QCQP, Lagrangian dual, certificate, lower bound, global optimality, branch and bound,
moment/Lasserre/SOS, Toeplitz/trigonometric moment/Fejér–Riesz, low-rank/nuclear norm), neighbouring fields
(beamforming/MLS, holography/CGH, optical tweezers, DOE/metasurface, coherent-mode/mutual-intensity synthesis, pupil
mask design, photonic inverse-design bounds), patents (USPTO full-text mirror, Google Patents) and venue-restricted
queries (SPIE Optical Microlithography / Photomask Technology, JM3, Acta Optica Sinica, Chinese Optics Letters, JJAP,
arXiv 2024–2026).

**Access caveat.** Only search snippets were available; no publisher full text, no SPIE Digital Library full-text
search, no direct Google Patents full-text query. Every item below is recorded only as the snippet described it. A
"no overlap" verdict therefore means "not surfaced by snippet-level search", not "proven absent".

---

## 1. Table of hits

Overlap codes refer to the four claims: **C1** SDP/Shor lifting + dual certificates for Hopkins/Abbe mask synthesis;
**C2** relaxation = unlimited incoherent multiple exposure; **C3** Carathéodory–Toeplitz realizability cone for
nonnegative masks in lithography; **C4** rigorous branch-and-bound certificate of the continuous single-mask optimum.

| # | Title | Authors | Venue / year | URL | What it does (snippet basis) | Overlap | Reasoning |
|---|---|---|---|---|---|---|---|
| 1 | Positive Semidefinite Relaxation and Approximation Algorithm for Triple Patterning Lithography | Matsui, Kohira, Kodama, Takahashi | ISAAC 2014 (LNCS), pp. 365–375 | https://link.springer.com/content/pdf/10.1007/978-3-319-13075-0_29.pdf | "polynomial time (7/9)-approximation algorithm based on positive semidefinite relaxation and randomized rounding" for TPL layout decomposition (mask *assignment*, a colouring problem) | C1: **no** (partial on wording only) | Same tool (PSD relaxation + Goemans–Williamson-type rounding), different problem: discrete colour assignment, no imaging model, no intensity quadratic forms, no dual certificate. Should still be cited as "the only prior SDP use in lithography". |
| 2 | Fast mask assignment using positive semidefinite relaxation in LELECUT triple patterning lithography | Kohira et al. | IEEE (ASP-DAC 2015), doc 7059084 | https://ieeexplore.ieee.org/document/7059084/ | PSD relaxation + randomized rounding for LELECUT mask assignment; "builds on the Goemans–Williamson ±1 formulation for max-cut" | C1: **no** | As #1. |
| 3 | Yield-aware Mask Assignment using Positive Semidefinite Relaxation in LELECUT TPL | Kohira et al. | SPIE 9427 (2015) / JM3 2016 (per T2R2 listing) | https://t2r2.star.titech.ac.jp/rrws/file/CTT100684777/ATD100000413/ | Extension of #2 to overlay-error tolerance | C1: **no** | As #1. |
| 4 | LRSDP: Low-Rank SDP for Triple Patterning Lithography Layout Decomposition | (CUHK, B. Yu group) | DAC 2023 | https://henryhxu.github.io/share/yu-dac23.pdf | Low-rank SDP via augmented Lagrangian + Riemannian optimization for MPL decomposition; "speedups of 186×, 25×, 12×" | C1: **no** | Decomposition/colouring, not image synthesis. |
| 5 | Model calibration for optical lithography via semidefinite programming | Nouralishahi, Wu, Vandenberghe | Optimization and Engineering 9:19–35 (2008) | https://link.springer.com/article/10.1007/s11081-007-9026-y ; http://seas.ucla.edu/~vandenbe/publications/calib.pdf | SDP for *model calibration* (fit an imaging kernel to measurements); imposes PSD on the kernel matrix "so that I(x) ≥ 0"; custom interior-point exploiting "rank-one structure" | C1: **partial (weak)** | Uses PSD-constrained kernel (Hopkins-type) matrices in lithography, but for estimation, not mask design; no relaxation of a mask QCQP, no dual bound on printability. Must be cited and distinguished explicitly, because a reviewer searching "semidefinite lithography" will find it first. |
| 6 | Step-walk relaxation method for global optimization of masks (US 7,587,702; US 7,685,559) | Rosenbluth et al. (IBM) | USPTO | https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/7685559 | "method for narrowing the search over local intensity minima when globally optimizing lithographic masks"; problem has "non-convex constraints in a high-dimensioned solution space", "scaling of the problem difficulty being exponential in the number of dimensions" | C1: **no**; C4: **partial (goal)** | Global *search* over the nonconvex QCQP by partitioning/enumeration; no convex relaxation, no bound independent of the search. Shares the goal of C4 (global optimum of the single-mask wavefront problem) but not a rigorous certificate; snippets give no optimality proof. |
| 7 | High contrast lithographic masks (US 7,944,545) | Rosenbluth et al. (IBM) | USPTO | https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/7944545 | "optimal mask as a minimization over diffraction orders subject to quadratic constraints"; changes variables into "joint eigenvector space"; optimum is "near-closest point to an origin lying outside hyper-ellipsoids defined by the constraint matrices" | C1: **partial (formulation)** | Confirms the QCQP-in-diffraction-orders *formulation* is prior art (manuscript already says so). No lifting, no SDP, no dual certificate in snippet. |
| 8 | Method for generating a plurality of optimized wavefronts for a multiple exposure lithographic process (US 8,495,528) | Rosenbluth, Sakamoto, Tian et al. (IBM) | USPTO, filed 2010 | https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/8495528 | Multi-exposure wavefront optimization; "dose adjustment applied externally to the optimization during the e-th exposure" | C2: **partial (weak)** | Treats multiple exposures as a design variable in the same quadratic framework, but snippets show no statement that the relaxation (rank-unconstrained PSD) *equals* unlimited incoherent exposure, nor any rank bound. Closest prior work for C2; cite and distinguish. |
| 9 | Source-mask optimization for a lithography process (US 2014/0282291; US 8,954,898) | IBM | Google Patents / USPTO | https://patents.google.com/patent/US20140282291A1/en | Mask diffraction orders as vector variable; "moves the problem into a joint eigenvector space"; wavefront step "quadratic, pseudo-global", source step "linear, global"; "if the diffracted wave front orders are fixed, it is possible to find the globally optimum solution for the source intensities" | C1: **partial (formulation)** | Same as #7; explicitly calls the wavefront step only "pseudo-global", which supports the manuscript's framing that IBM yields no certified bound. |
| 10 | Synthesizing low MEEF lithography solutions (US 9,395,622) | IBM | Google Patents | https://patents.google.com/patent/US9395622B2/en | "constraint is inherently quadratic, because the optical intensity is the square of diffraction order amplitude"; joint source+wavefront via "local optimizer" | C1: **no** | Confirms quadratic framing only. |
| 11 | Method for forming arbitrary lithographic wavefronts using standard mask technology (US 8,108,802) | Rosenbluth (IBM) | USPTO | https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/8108802 | Realize an arbitrary set of diffraction orders with standard (PSM) mask technology | C3: **partial (weak)** | Addresses realizability of a target wavefront by a *phase-shift* mask (complex transmission), so no Toeplitz/nonnegativity cone arises. Relevant as the statement that complex masks realize any in-band spectrum (which makes C3 interesting only for nonnegative masks). |
| 12 | Optimum mask and source patterns to print a given shape | Rosenbluth et al. | JM3 1(1), 2002 (cited via patents) | (not directly indexed) | "first determining the optimum set of image forming waves ... then designing manufacturable mask shapes to produce these diffracted waves"; "does not require initial conditions" | C1: **partial (formulation)**, C4: **no** | Wavefront-first decomposition; snippets describe no relaxation or bound. |
| 13 | Global optimization of masks, including film stack design... | Rosenbluth et al. | SPIE 6520 (2007) (seen earlier in NOVELTY.md) | — | "non-convex quadratic constraints"; "solution space can be split into regions with at most one local minimum" | C4: **partial (goal)** | Partition-based global method; not a convex relaxation; no bound for externally found solutions. (Not re-surfaced in this session; recorded from NOVELTY.md.) |
| 14 | Binary and phase shifting mask design for optical lithography | Y. Liu, A. Zakhor | IEEE Trans. Semicond. Manuf. 5(2):138–152, 1992 | http://www-video.eecs.berkeley.edu/papers/YLiu/Binary%20IEEE%20SM%2092.pdf | Hopkins model; binary LP via "branch and bound and simulated annealing"; "annealing is used to avoid local minima" | C4: **partial** | Branch-and-bound for *combinatorial* binary mask (pixel MSE); not a certificate of the continuous ratio optimum; no interval bounds on PSD forms. Already cited in manuscript; wording there is accurate. |
| 15 | Optimal binary image design based on the branch and bound algorithm / Optimal binary image design for optical lithography | Liu, Zakhor | IEEE (doc 115862) 1990; SPIE 1264 (1990) | https://ieeexplore.ieee.org/document/115862/ ; https://www.spiedigitallibrary.org/conference-proceedings-of-spie/1264/1/... | "formulated as a binary linear programming problem and solved with the branch and bound algorithm"; "MSE criterion ... simplex" | C4: **partial** | As #14; earlier conference versions — add to the liu1992 citation. |
| 16 | Fast Mask Optimization Under Process Variation Using Guided Local Search on Quadratic Programming | Nonaka, Kuramochi, Kohira, Azuma, Matsui, Takahashi, Kodama | IEEE TCAD 44:4355–4366 (2025), doi 10.1109/TCAD.2025.3558148 | https://ieeexplore.ieee.org/document/10949618/ | "formulates mask optimization as a 0-1 convex QCQP ... maximize tolerance under fidelity constraints"; "relaxes the 0-1 QCQP into a 0-1 convex QP by introducing Lagrange multipliers"; local search | C1: **partial (weak)** | Explicit QCQP framing of pixel mask optimization plus a Lagrangian (not semidefinite) relaxation used as a *solver*, not a bound; no PSD lifting; no certificate. Not cited in manuscript — should be. |
| 17 | A formulation of mask optimization into QUBO model for Ising machines | (same group) | SPIE 12751 (2023) | https://www.spiedigitallibrary.org/conference-proceedings-of-spie/12751/127511D/... | Mask optimization via repeated QUBO solves on an Ising machine | C1/C4: **no** | Heuristic; no bound reported in snippet. |
| 18 | A fast process variation and pattern fidelity aware mask optimization algorithm | Awad, Takahashi, Tanaka, Kodama | ICCAD 2014 | https://www.semanticscholar.org/paper/dd915898ac6308fef830f3e0c0d1f651b4af6a7b | Earlier intensity-modeling mask optimization in the same line | **no** | Snippet has no method detail suggesting relaxation/bounds. |
| 19 | Mask design for single and double exposure (PhD) / Double-exposure mask synthesis using inverse lithography | A. Poonawala (adv. Milanfar) | UCSC thesis; JM3 | https://users.soe.ucsc.edu/~milanfar/publications/dissertations/PoonawalaThesis.pdf ; https://users.soe.ucsc.edu/~milanfar/publications/journal/JM3_Final.pdf | "combined image is the sum of individual exposure images", gradient descent; notes identical initialisation of both masks "defeats the purpose of double exposure" | C2: **no** | Models incoherent sum but gives no relaxation interpretation or bound on exposure count. |
| 20 | Diffraction-limit-breaking digital projection lithography via multi-exposure strategies | — | Microsystems & Nanoengineering (Nature), 2025 | https://www.nature.com/articles/s41378-025-01131-x | Multi-exposure DMD lithography; "coherence between sub-layouts is eliminated through ME"; gradient-descent sub-layout design | C2: **no** | Empirical; no relaxation equivalence. |
| 21 | Computational Bounds for Photonic Design / Heuristic Methods and Performance Bounds for Photonic Design / Bounds on Efficiency Metrics in Photonics | Angeris, Vučković, Boyd | ACS Photonics 6:1232 (2019); Opt. Express 2021; arXiv 2204.05243 | https://stanford.edu/~boyd/papers/pdf/comp_imposs_res.pdf | Lagrange-duality lower bounds on physical (Helmholtz) design; "guaranteed lower bound on the best possible objective value"; ADMM design "within about 8.7% of the dual bound" | C1: **partial (philosophy)** | Same certificate philosophy, different physics (PDE constraint vs. quadratic intensity forms). Already cited. No lithography/imaging application surfaced ("found no paper by these authors on lithography or imaging bounds"). |
| 22 | Fundamental Limits of Nanophotonic Design (thesis) | (Stanford; Angeris-line) | arXiv 2308.08698 | https://arxiv.org/pdf/2308.08698 | Shows a quadratic-constrained design problem can be turned into an SDP; notes "the semidefinite relaxation does not introduce any actual relaxation when there is only one global constraint" | C1: **partial (philosophy)** | General SDP-bound methodology for photonic design; not imaging/lithography. |
| 23 | Shaped beam synthesis via semidefinite relaxation (SDR) based on magnitude least squares | Hu et al. | IET Microwaves, Antennas & Propagation 2022 | https://ietresearch.onlinelibrary.wiley.com/doi/full/10.1049/mia2.12231 | SDR for antenna "beam pattern mask constraints" (upper/lower bounds on |pattern|²) | C1: **partial (math template)** | Same quadratic-form max/min template as coherent-illumination mask synthesis; antenna domain. Analogous to sidiropoulos2006 already cited. |
| 24 | Semidefinite Relaxation of Quadratic Optimization Problems | Luo, Ma, So, Ye, Zhang | IEEE SPM 2010 | https://www.researchgate.net/publication/224131024 | Tutorial on Shor relaxation | C1: textbook (already cited) | — |
| 25 | A finite branch-and-bound algorithm for nonconvex quadratic programming via semidefinite relaxations | Burer, Vandenbussche | Math. Programming 2008 | https://link.springer.com/article/10.1007/s10107-006-0080-6 | SDP-bounded B&B for nonconvex QP | C4: **partial (method family)** | Generic SDP-based B&B; not rigorous (floating-point) and not lithography. Shows C4's *method* is standard; novelty is the lithography application + interval/rigorous PSD-form bounds. |
| 26 | An Enhanced SDR based Global Algorithm for Nonconvex Complex Quadratic Programs | (Lu, Liu, Dai?) | arXiv 1902.04287 | https://arxiv.org/pdf/1902.04287 | B&B with SDR bounds for complex QCQP (virtual beamforming, MIMO); "guaranteed to find the global solution ... within any given error tolerance" | C4: **partial (method family)** | Same mathematical object (complex QCQP global solve); beamforming domain; tolerance-based, not rigorous. |
| 27 | New semidefinite relaxations for a class of complex quadratic programming problems | — | arXiv 2305.09934 | https://arxiv.org/pdf/2305.09934 | B&B for complex QP with modulus/phase constraints, tested on beamformer design | C4: **partial (method family)** | As #26. |
| 28 | Designing optimal masks for a multi-object spectrometer | — | arXiv 2101.03378 | https://arxiv.org/pdf/2101.03378 | B&B for *astronomical slit-mask* assignment; reports "% deviation between LB and UB" | C4: **no** | "Mask" is a slit-allocation problem; combinatorial. False-positive keyword hit. |
| 29 | Phase-shifting masks for microlithography: automated design and mask requirements | Pati, Kailath | JOSA A 11(9):2438 (1994); SPIE 2197 | https://opg.optica.org/josaa/abstract.cfm?uri=josaa-11-9-2438 | Optimal coherent approximation of Hopkins model → phase retrieval; Gerchberg–Saxton-type "(suboptimal)" masks; later cited as "sub-optimal projections onto convex sets" | C1: **no** | POCS/phase retrieval, no lifting, explicitly suboptimal. Already cited (pati1994). |
| 30 | Phase Retrieval with Masks using Convex Optimization | Jaganathan, Eldar, Hassibi | Caltech/Technion (ICASSP-era) | https://authors.library.caltech.edu/records/6n5mq-e1e21 | SDP phase retrieval using measurement masks | **no** | "Masks" are measurement modulations for recovery; not synthesis. False-positive. |
| 31 | Convex phase retrieval with multiple structured illuminations design | — | J. Optics (Springer) 2020 | https://link.springer.com/article/10.1007/s12596-020-00593-8 | PhaseLift-style recovery with physically realistic illumination design | **no** | Recovery, not synthesis. |
| 32 | Spectral mask constraints on PSD autocorrelation (filter design) | Davidson, Luo, Sturm | IEEE TSP 50(11):2702, 2002 | https://www.ece.mcmaster.ca/~davidson/pubs/Davidson_etal_spec_mask_constrs.pdf | "finite convex representation for spectral mask constraints on a positive semidefinite autocorrelation" | C3: **partial (math template)** | Autocorrelation/Toeplitz-PSD convex cone used for filter design (signal processing). Same classical mathematics as Theorem 4; "mask" here is a spectral template, not a photomask. Good citation for "the Toeplitz-PSD cone is standard in filter/sequence design". |
| 33 | McLean & Woerdeman criterion (as cited in "Solution to the Inverse Wulff Problem by Means of the Enhanced Semidefinite Relaxation Method") | McLean, Woerdeman; cited by Ševčovič et al. | arXiv 1402.5668 | https://arxiv.org/pdf/1402.5668 | "semidefinite representable criterion for non-negativity of a partial finite Fourier series" via Riesz–Fejér | C3: **partial (math template)** | Classical SDP characterisation of nonnegative trigonometric polynomials, applied to crystal-shape (Wulff) problems — not lithography. |
| 34 | The Optimization Landscape of Carathéodory Decomposition of Toeplitz Covariances | — | arXiv 2511.01605 | https://arxiv.org/html/2511.01605 | Carathéodory decomposition of Toeplitz covariances (array processing) | C3: **no** (math background) | — |
| 35 | Band-limited coronagraph masks | Kuchner, Traub | ApJ 2002 (astro-ph/0203455) | https://arxiv.org/pdf/astro-ph/0203455 | Mask "transfer function nonzero only inside a finite band and lies between 0 and 1" | C3: **no** | Band-limited nonnegative amplitude masks in astronomy; no Toeplitz cone or convex realizability statement in snippet. |
| 36 | Nonlinear compressive inverse lithography aided by low-rank regularization | Ma, Arce et al. | Opt. Express 27(21):29992 (2019?) | https://opg.optica.org/oe/fulltext.cfm?uri=oe-27-21-29992 | "sparsity and low-rank regularizations to constrain the solution space", split Bregman | C1: **no** | Low rank of the *mask image*, not of a lifted Gram matrix. |
| 37 | Efficient source and mask optimization with augmented Lagrangian methods / ADMM for inverse lithography | (Opt. Express 21:8076, 2013); (NMTMA / arXiv 2209.10814) | — | https://opg.optica.org/oe/fulltext.cfm?uri=oe-21-7-8076 ; https://arxiv.org/pdf/2209.10814 | Lagrangian machinery as *solver* ("alternating between minimizing the primal given its dual, and maximizing the dual") | C1: **no** | No bound/certificate use of duality. |
| 38 | Level-set-based inverse lithography (HKU) | Shen et al. | JVST B 2013 | https://www.eee.hku.hk/optima/pub/journal/1307_JVSTB.pdf | States "finding a rigorous global minimum is usually too expensive and is not necessary" | C1/C4: **no** | Useful quote showing the field's stance. |
| 39 | Advancements and challenges in ILT: a review of AI-based approaches | — | Light: Sci. & Appl. 2025 | https://www.nature.com/articles/s41377-025-01923-w | ILT "near-optimal", local minima; recounts Liu & Zakhor B&B | **no** | Review; no certified bounds. |
| 40 | Pushing the Limits of Inverse Lithography with Generative Reinforcement Learning | (NVIDIA) | arXiv 2602.19027 (2026) | https://arxiv.org/pdf/2602.19027 | Learned ILT; only reference list seen | **no** | Already cited (nvidia2026). |
| 41 | Lithographic exposure latitude aware source and mask optimization | — | ResearchGate 2021 | https://www.researchgate.net/publication/357013097 | Penalty-based SMO raising EL "from 5% to 11%" | **no** | Heuristic; no bound. |
| 42 | Printing a mask with maximum possible process window through adjustment of the source distribution (US 7,057,709 family) | IBM | USPTO | https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/7057709 | Max latitude = "minimum of the upper exposure limits across edges" — a worst-case feature bound; LP for source with fixed mask | **no** | Source-only global optimization (linear in source), already a known fact; not a mask bound. |
| 43 | Fühner/Erdmann GA-based source–mask optimization | Fühner, Erdmann et al. | SPIE 5182 (2003), 5377 (2004), 5754 (2005), 7122 (2008); JM3 2006 | https://www.spiedigitallibrary.org/conference-proceedings-of-spie/7122/1/... | "global optimization techniques" = genetic algorithms | **no** | "Global" means metaheuristic; no certificate. Useful to pre-empt reviewer confusion over the word "global". |
| 44 | Band Gap Optimization of 2D Photonic Crystals Using Semidefinite Programming and Subspace Methods | Men, Nguyen, Freund, Parrilo, Peraire | arXiv 0907.2267 / J. Comput. Phys. 2010 | https://arxiv.org/pdf/0907.2267 | SDP for photonic crystal band gaps | **no** | Keyword false positive. |
| 45 | Semidefinite relaxations in optimal experiment design ... hyperpolarized MRI | — | arXiv 1510.00455 | https://arxiv.org/pdf/1510.00455 | SDP exactness conditions for squared-variable QPs | **no** | Not optics. |

Not surfaced at all (negative results, each from ≥2 queries): any SPIE/JM3/IEEE/arXiv paper or patent combining
{semidefinite, SDP, Shor, lifting, PhaseLift, rank-one, moment, Lasserre, SOS} × {mask optimization, ILT, OPC, SMO,
wavefront optimization, aerial image, Hopkins, TCC, diffraction orders}; any "convex hull of achievable aerial images
under multiple exposure"; any Toeplitz/Carathéodory/Fejér–Riesz statement about photomask spectra; any
"theoretical upper bound on NILS/exposure latitude over all masks"; any ILT benchmark reporting an optimality gap
against a certified bound; any Chinese/Japanese/Korean-language hit (queries returned only generic SDP tutorials and
unrelated ILT/dataset papers).

---

## 2. Verdict per claim

### C1 — Shor/lifted SDP relaxation with dual certificates for Hopkins/Abbe mask synthesis
**Overlap found: no direct overlap. Confidence: medium (≈70%).**
- Nearest prior art: (a) IBM wavefront QCQP formulation (#7, #9, #12, #13) — same primal problem, no lifting/dual;
  (b) Nouralishahi–Wu–Vandenberghe 2008 (#5) — SDP *in lithography* but for model calibration; (c) Kohira/Kodama
  group — PSD relaxation for mask *assignment* (#1–#3) and a Lagrangian-relaxed 0-1 QCQP for pixel mask
  optimization (#16); (d) Angeris–Vučković–Boyd dual bounds in photonics (#21, #22); (e) SDR for antenna "mask"
  constraints (#23).
- Why confidence is not higher: the SPIE Digital Library, IEEE Xplore and Google Patents full texts cannot be
  searched; an SDP treatment buried in an SPIE paper body, a thesis chapter, or a patent claim set using words like
  "convex relaxation of the wavefront problem" would not appear in snippets. The IBM group in particular wrote
  dozens of patents on the wavefront QCQP and it is plausible (not shown) that one of them mentions an SDP bound.
- Action: cite #5 and #16 explicitly (both currently absent from the manuscript) and distinguish them.

### C2 — Relaxation ≡ unlimited incoherent multiple exposure; Pataki rank bound on useful exposures
**Overlap found: no. Confidence: medium-high (≈75%) for the lithography *interpretation*; the mathematics (PSD = Gram
matrix of vectors; rank bound) is folklore and must be stated as such.**
- Nearest prior art: US 8,495,528 (#8) optimizes several exposures in the same quadratic framework; Poonawala (#19)
  models incoherent sums. Neither surfaces the equivalence "SDP value = infinitely-many-exposure optimum" or a rank
  bound on the number of useful exposures.
- Risk: the equivalence is a one-line observation; a reviewer may consider it obvious rather than novel. The
  defensible novelty is the *use*: measuring the multi-exposure value as the SDP gap and certifying strictly positive
  gains via C4.

### C3 — Carathéodory–Toeplitz cone as realizability constraint for nonnegative (binary/attenuated) masks
**Overlap found: no in lithography; the mathematics and its use as a design cone exist in filter/sequence design
(#32, #33) and array processing (#34). Confidence: medium-high (≈75%) for "not previously used in lithography";
the manuscript must not claim mathematical novelty.**
- Nearest lithography item: Rosenbluth US 8,108,802 (#11) on realizing arbitrary wavefronts with PSM — complementary
  (complex masks have no such cone).

### C4 — Rigorous branch-and-bound certificate of the continuous single-mask optimum
**Overlap found: partial on method family, none on application. Confidence: medium (≈65%).**
- Branch-and-bound in lithography: Liu & Zakhor 1990/1992 (#14, #15) — combinatorial binary LP; IBM step-walk /
  partition global search (#6, #13) — non-convex global search without an independent bound.
- Branch-and-bound with SDP bounds for complex QCQPs is established in beamforming/MIMO (#25–#27), with tolerance-based
  (floating-point) guarantees rather than rigorous interval/PSD-form bounds.
- The specific combination (interval bounds on PSD quadratic forms over sphere/box cells, applied to the ratio
  objective of the single-mask wavefront problem, used to separate single- from multi-exposure optima) was not found.
- Risk: "rigorous" must be defined precisely in the paper (what is verified in exact/high-precision arithmetic and
  what relies on floating-point SDP values), because the B&B-for-QCQP community will ask.

---

## 3. Recommended wording for the novelty statements

Replace absolute claims ("no prior work", "new", "first") with calibrated, search-scoped language and name the
nearest neighbours in the same sentence. Suggested text:

**Related work / C1.**
> Semidefinite programming has entered lithography in two places: model calibration, where Nouralishahi, Wu and
> Vandenberghe fit a positive-semidefinite imaging kernel to measurements [Optim. Eng. 2008], and layout
> decomposition for multiple patterning, where positive-semidefinite relaxations with randomized rounding solve the
> combinatorial mask-assignment problem [Matsui et al. ISAAC 2014; Kohira et al. 2015/2016; LRSDP, DAC 2023]. Pixel
> mask optimization has been written as a 0-1 QCQP and attacked with a Lagrangian relaxation used as a local-search
> heuristic [Nonaka et al., IEEE TCAD 2025]. IBM's wavefront optimization poses the same diffraction-order QCQP that
> we lift [Rosenbluth et al. 2002, 2007; US 7,944,545; US 2014/0282291] and solves it by eigen-subspace and
> partition-based global search described in the patents as "pseudo-global". To our knowledge — based on a search of
> titles, abstracts and patent texts, not of full SPIE proceedings — the Shor lifting of this problem and its dual
> certificate of printability have not been reported; the mathematics is textbook [Luo et al. 2010; Boyd &
> Vandenberghe 2004] and our contribution is its application, the certificate verification protocol and the
> structural results below.

**C2.**
> That a positive semidefinite matrix is a Gram matrix, and hence that the relaxation (R) is the mask-synthesis
> problem with an unlimited number of incoherently added exposures, is elementary; Pataki's bound on the rank of
> extreme points is likewise classical. We are not aware of this identification, or of the resulting rank bound on
> the number of useful exposures, in the lithography literature, where multiple exposure has been optimized directly
> [US 8,495,528; Poonawala & Milanfar JM3].

**C3.**
> The characterization of in-band spectra of nonnegative functions by a positive semidefinite Toeplitz matrix is the
> Carathéodory–Toeplitz theorem and is routinely used as a convex design cone in filter and sequence design [Davidson,
> Luo & Sturm, IEEE TSP 2002; McLean & Woerdeman]. We use it as an exact realizability constraint for binary and
> attenuated photomasks, which we have not found in prior lithography work; for phase-shift masks no such constraint
> exists since any in-band wavefront is realizable [Rosenbluth, US 8,108,802].

**C4.**
> Branch and bound has been applied to lithographic mask design as a combinatorial method for binary pixel masks
> [Liu & Zakhor 1990, 1992], and semidefinite-bounded branch and bound for complex QCQPs is established in
> beamforming [Burer & Vandenbussche 2008; Lu et al. 2019]. Theorem 3 differs in object and in guarantee: it bounds the
> continuous single-mask ratio problem using interval enclosures of positive semidefinite quadratic forms, so that
> the resulting lower bound — and therefore the certified multi-exposure gain — holds in verified arithmetic rather
> than to solver tolerance.

**Global caveat sentence (abstract or §Limitations).**
> Novelty claims are calibrated to a search of abstracts, titles and patent texts (47 multilingual queries; SPIE and
> IEEE full texts were not accessible); we claim the application and the certificate protocol, not the underlying
> relaxation or moment-theory results.

---

## 4. Residual risks (ranked)
1. An SDP/convex-relaxation bound for the wavefront QCQP inside an IBM patent body or an SPIE paper (2002–2015) not
   indexed by snippet search. Mitigation: full-text search of Google Patents for assignee IBM + "semidefinite" and of
   SPIE DL for "semidefinite" within Optical Microlithography volumes, before submission.
2. The Kohira/Kodama/Takahashi group (Aizu / Tokyo Tech / Kioxia) has both PSD-relaxation and QCQP-mask papers; an
   SDP relaxation of their 0-1 QCQP may exist in a Japanese-language IEICE technical report (one IEICE 2018 report on
   "quadratic programming pixel-based OPC" surfaced; content unknown).
3. A theses-only result (e.g., HKU, Berkeley, Stanford, UCSC) proving convexity of the multi-exposure image set.
4. C2/C3 are mathematically elementary; a reviewer may deny novelty regardless of literature absence. Wording above
   concedes the mathematics and claims only the application.
