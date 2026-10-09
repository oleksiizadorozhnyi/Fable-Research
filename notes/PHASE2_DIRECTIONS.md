# Phase II — Five research directions, evaluation, and selection
(All "state of knowledge" claims below are from WebSearch metadata/abstract snippets; full texts not accessible — see RESEARCH_STATE.md.)

## Direction A — Convex (lifted/SDP) certificates of printability for partially coherent mask synthesis
Question: For a target pattern (bright set B, dark set D), illumination, NA, λ, what is the best achievable dark/bright
intensity ratio over ALL masks (binary / PSM / multi-exposure), and can it be CERTIFIED (not just found by heuristics)?
Motivation: ILT/SMO are nonconvex; practitioners cannot tell whether a pattern is fundamentally unprintable or whether
their optimizer is stuck (reviews PMC12287475, PMC12843865; NVIDIA arXiv 2602.19027 "poor local minima").
State of knowledge: Hopkins image is quadratic in mask -> QCQP (IBM/Rosenbluth: "non-convex quadratic constraints";
step-walk relaxation patents US7587702/7685559, multi-exposure wavefront patent US8495528 use eigen-subspace heuristics).
SDP relaxation of QCQP is textbook (Luo/Ma/So/Ye/Zhang 2010; multicast beamforming Sidiropoulos/Davidson/Luo 2006).
Ma/Arce/Garcia-Frias (Opt. Express 2017, 2018; 2021; 2022; 2025) give information-theoretic "limits" on OPC fidelity
(mutual-information based, not duality certificates).
Gap: no lifted/SDP relaxation, no dual certificates, no multi-exposure convexity statement, no realizability cone, no
measured optimality gap of ILT vs a certified bound, found in the literature reviewed.
Proposed insight: (1) lift X=mm^T; ratio problem becomes an SDP whose dual yields an elementary, independently
checkable certificate (weights mu on dark, lambda on bright with sum mu_d Q_d - sum lambda_b Q_b PSD => no mask beats
ratio sum lambda_b). (2) The relaxation EQUALS the unlimited-multiple-exposure (incoherent superposition) problem, so the
gap between relaxation and best single mask quantifies the value of multi-exposure. (3) In 1D the in-band spectra of
nonnegative (binary-type) masks form exactly the Caratheodory–Toeplitz PSD cone -> exact convex realizability.
(4) Reproduces k1 limits; gives pattern-specific limits for 2D patterns; measures how far gradient ILT is from bound.
Resources: numpy/cvxpy (CLARABEL/SCS). Obstacles: 2D scale (Fourier-domain reduces F ~ (2(1+sigma)NA L/lambda)^2);
relaxation gap may be large; novelty vs IBM eigen-methods must be stated carefully.
Minimal test: 1D dense lines must give ratio bound = 1 (unprintable) below pitch lambda/((1+sigma)NA) and <1 above.
Probability of meaningful result (subjective): ~0.6. Significance if successful: medium-high (tooling for ILT/SMO,
design-rule derivation, "entitlement" analysis).

## Direction B — Information-theoretic capacity of stochastic (photon-shot-noise) EUV pattern transfer
Question: bits per unit area reliably patternable given bandlimited quadratic imaging + Poisson absorption + threshold resist.
State: Ma/Arce (2017-2025) already model lithography as information channel (deterministic imaging; fidelity limits).
Shot-noise numbers known (6.8 photons/nm^2 at 10 mJ/cm^2). Poisson-channel capacity with bandwidth constraints is hard
(Shamai; Lapidoth). Gap: no shot-noise capacity framing found. Obstacles: result would be framework-level, bounds loose,
no experimental validation possible here. Probability ~0.3; significance medium.

## Direction C — Theorem on TCC (Hopkins kernel) eigenvalue decay / number of significant SOCS kernels
Question: prove asymptotics of eigenvalues of W(x,x') = p(x) p*(x') s^(x-x') (source FT times coherent PSF outer product).
State: SOCS (Pati & Kailath 1994; Cobb 1998); Yamazoe stacked pupil shift (TCC = P P^H); 2007 SPIE: "no general guideline
where to cut off". Landau–Widom / Slepian theory for time-frequency limiting gives the template. Gap: no theorem found.
Obstacles: Airy PSF not compactly supported -> polynomial not super-exponential decay; must handle 2D disk geometry.
Probability ~0.5 (numerical evidence easy; proof moderate). Significance medium-low (model-truncation guidelines).

## Direction D — Mathematical theory of superposition (interpretability): optimal geometry / loss bounds
State: Elhage et al 2022; arXiv 2310.06301; "Towards a mathematical theory of superposition" arXiv 2608.27540 (Aug 2026);
June-2026 preprint with upper/lower L2 bounds for power activations; Sharkey et al 2025 open problems.
Gap unclear because recent full texts unreadable here -> high risk of duplicating 2026 results. Probability ~0.25.

## Direction E — Theory of BTSP-based credit assignment (neuroscience)
State: Cone, Clopath, Ponte Costa bioRxiv 2025 (gBTSP); "comprehensive theoretical framework still lacking".
Gap: what target signals/plateau timing make BTSP a valid credit-assignment rule. Obstacles: needs biological data for
validation; theory alone would be a model, not a mechanism. Probability ~0.3; significance high if validated (not possible here).

## Scoring (0-10 each; weights 25/25/20/15/10/5). Subjective, evidence-informed.
| Dir | Originality | Significance | Tractability | Feasibility | Verification | Interdisc | Weighted |
| A | 7 | 7 | 8 | 9 | 9 | 8 | 7.75 |
| B | 5 | 6 | 4 | 6 | 3 | 7 | 5.05 |
| C | 6 | 5 | 6 | 8 | 8 | 4 | 6.15 |
| D | 4 | 7 | 5 | 6 | 4 | 5 | 5.25 |
| E | 6 | 8 | 3 | 4 | 2 | 8 | 5.30 |
Weighted = .25*O + .25*S + .20*T + .15*F + .10*V + .05*I.
Selection: A. Backup: C (shares the same imaging machinery).
