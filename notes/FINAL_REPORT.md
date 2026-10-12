# FINAL SCIENTIFIC DELIVERABLE — v3 (adversarial hardening and clean reproduction, 2026-10-12)
Classification: **VALIDATED RESEARCH CANDIDATE** — publication-ready as a methods paper pending the items in §11.
All numbers below are generated from the result files of one clean execution of `code/run_all_par.sh` (started
12:30 UTC, tests → 21 experiment scripts → independent verification → tables → figures) by `code/make_numbers.py`;
the manuscript `paper/main.tex` contains no hand-typed number.

## 0 Completion-gate status
| Gate | Status | Evidence |
|---|---|---|
| A Precision | PASS | notes/THEORY.md (Theorems 1–9 + session-2/3 corrections); manuscript §3–4 with notation table |
| B Correctness | PASS (model-internal) | proofs reviewed by an independent mathematics review (review_math.md; 8 MAJOR items fixed); every bound re-derived in the clean run; certificate-consistency check (best mask ≥ bound) on every record: 0 violations; 136+ certificates re-verified by 50-digit Cholesky; two solvers; 13 regression tests incl. exact 2D sampling and coherent-mode compression |
| C Originality | PASS with caveats | independent prior-art review (review_priorart.md, 47 queries): SDP/dual certificates, multi-exposure identification, Toeplitz cone and B&B separation not found; nearest neighbours cited; full texts of five key works unread |
| D Significance | PASS (medium) | pattern-specific EL entitlements for DUV/EUV; certified benchmarking of ILT; certified multi-exposure and PSM gains; structure survives vector model and aberrations |
| E Reproducibility | PASS | `code/run_all_par.sh` regenerates every result; `make_numbers.py` regenerates every number/table; logs with exit codes |
| F Adversarial review | PASS with documented limitations | six independent reviews (math, numbers, physics, prior art, editor, code); ADVERSARIAL.md items 1–22; bug log in RESEARCH_STATE.md |

## PART I — Executive research summary
Problem: ILT/SMO are nonconvex; a local optimum says nothing about the achievable optimum. Contribution: for the
Abbe/Hopkins model (scalar or vector, aberrated, multi-focus, immersion) of a periodic cell, the maximum full-range
exposure latitude over all masks at a 10 %-CD edge tolerance is bounded by a small SDP in the diffraction-order domain
whose dual is a one-matrix certificate (Theorem 1), verified independently with a positive margin and a 50-digit
Cholesky. The relaxation is exactly the unlimited multi-exposure problem (Theorem 3); complex phase cannot improve it
under symmetric illumination and symmetric focus sets (Lemma 6); nonnegative 1D masks have an exact Toeplitz
realizability cone (Theorem 4); vector imaging/aberrations/pixelation preserve the structure (Lemmas 7–8); and a
branch-and-bound (Theorem 9; interval and SDP-bounded variants) certifies single-mask optima for small F.
Results (final clean reproduction; manuscript tables): exact k1 cutoffs (DUV, EUV); on 108 non-degenerate random 1D
targets local optimization is certified globally optimal in 73 % (86 % within 1 %); 1D library gaps ≤ 0.4 % except
coherent illumination (10–19 %) and through-focus (1–6 %); all 2D contact arrays within 2.2 %, the 320 nm T-junction
(annular) within 0.25 %, line-end within 0.02 %, quasar T-junction within 1.2 %; certified multi-exposure gains
6.4–16.4 % (F = 5, converged B&B) and 8.9 % (F = 7, 30 M nodes); pixelated nonnegative masks forfeit 3–10 % vs PSM;
pattern-error ILT lands 20–32 % above the certified optimum; SDP bound 1.4–135× tighter than the Rayleigh-quotient
bound; vector TE bounds within 6 % of scalar, TM up to 2× (classical contrast loss); aberrations ≤ 0.1λ move bounds
≤ 2 %. 2D entitlements use the full offset contour including corners (conservative); a 0.25-CD corner-exclusion zone
raises the certified EL by 1.6–4.6× (Table "corner").

## PART II — Five research directions — notes/PHASE2_DIRECTIONS.md (unchanged).
## PART III — Literature and prior art — notes/NOVELTY.md (session-3 calibration), notes/review_priorart.md; manuscript §2.
## PART IV — Formal contribution — notes/THEORY.md; manuscript §4 (Theorems 1–9).
## PART V — Validation — results/RESULTS_TABLES.md, results/figs/fig1–6, results/logs/run_all.log, results/exp15_verify_all.json.
## PART VI — Adversarial analysis — notes/ADVERSARIAL.md items 1–22; notes/review_*.md (six reports).
Session-3 corrections of substance: (1) 1D "conventional" source was a slit → 2D disk/annulus/dipole/quasar quadratures
with convergence study; (2) defocus in air at NA 1.35 → water immersion n = 1.44; (3) EL half-range → full-range;
(4) τ 8 nm → 10 % CD; (5) CRITICAL: 2D edge-offset sampling used a pixel-quantised distance transform (empty rings at the
new τ, biased by one raster step otherwise) → exact geometric offset sampling on rectangle unions with asserted
distances; every 2D number regenerated (contact entitlements dropped from EL ≈ 0.9 to 0.02–0.22); (6) dead orders
dropped; (7) verification script now rebuilds exactly the reported instances; (8) all text numbers macro-generated.
## PART VII — Novelty assessment — notes/NOVELTY.md.
## PART VIII — Manuscript — paper/main.tex → paper/main.pdf (v3, 18 pp.); notes/PAPER_v2.md is the superseded v2 Markdown.
## PART IX — Reproducibility — README.md; code/ (run_all_par.sh, make_numbers.py, 13 tests); results/.
## PART X — Next steps
1. Certify single-mask optima for F ≈ 9–61 (stronger SDP-bounded B&B with better dual recovery); 2. non-periodic
clips and scalable SDP (chordal / Burer–Monteiro); 3. thick-mask (3D) models via boundary-layer linearizations and
rigorous EM validation; 4. NILS/MEEF/EPE and dose-aware constraints; 5. benchmark against open ILT engines on ICCAD-13
clips; 6. external replication; 7. full-text check of Rosenbluth 2002/2007, Liu–Zakhor 1992, Nonaka 2025,
Nouralishahi 2008.

## 11 Publication-readiness assessment
Ready for submission as a methods paper (JM3 / Optics Express / SIAM J. Imaging Sci.) with explicit caveats:
(a) thin-mask periodic scalar/vector model, aerial-image entitlements only; (b) single-mask optima certified only for
F ≤ 7; the quasar T-junction (1.2 %) and the large-F random-2D gaps (≤ 11 %) are unresolved; (c) 2D metric enforces
corner fidelity (conservative; sensitivity reported); (d) literature verified from abstracts/snippets only;
(e) "certified" = exact arithmetic + rounding safeguard, not interval-verified. The LaTeX source compiles cleanly
(pdflatex + bibtex; 0 errors, 0 overfull boxes). Not ready for claims beyond the model (no silicon validation).
