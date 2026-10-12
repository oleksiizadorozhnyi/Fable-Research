# FINAL SCIENTIFIC DELIVERABLE — v2 (hardening pass, 2026-10-12)
Classification: **VALIDATED RESEARCH CANDIDATE** — publication-ready pending the items in §11 (see assessment).

## 0 Completion-gate status
| Gate | Status | Evidence |
|---|---|---|
| A Precision | PASS | notes/THEORY.md §1–2, §6; manuscript §3 |
| B Correctness | PASS (model-internal) | Theorems 1–9 with proofs (elementary / classical); every reported bound re-derived in a single clean run and its certificate re-verified by a full-space eigenvalue test with 1e-11 tolerance and a 50-digit Cholesky; two solvers; 9 regression tests; rigorous B&B certificates for F ≤ 7 |
| C Originality | PASS with caveats | notes/NOVELTY.md (session-2 additions): QCQP framing is IBM prior art; SDP certificates, multi-exposure identification, B&B separation, Toeplitz cone, vector/aberration invariance not found; full texts unread |
| D Significance | PASS (medium) | pattern-specific EL limits; certified benchmarking of ILT; certified multi-exposure and PSM gains; survives vector model |
| E Reproducibility | PASS | code/run_all.sh regenerates everything; README; JSON/logs; figures |
| F Adversarial review | PASS with documented limitations | notes/ADVERSARIAL.md items 1–15; bug log in RESEARCH_STATE.md (4 bugs found and fixed in session 2) |

## PART I — Executive research summary
[[filled after clean run: see manuscript abstract]]

## PART II — Five research directions
notes/PHASE2_DIRECTIONS.md (unchanged).

## PART III — Literature and prior art
notes/NOVELTY.md; manuscript §2 and references. Verification limitation: search snippets only (egress blocks publishers).

## PART IV — Formal contribution
notes/THEORY.md: Theorem 1 (certificate), Lemma 2 (real structure), Theorem 3 (multi-exposure equivalence; Pataki),
Theorem 4 (Toeplitz cone, full proof), Proposition 5 (cutoff), Lemma 6 (complex = real), Lemma 7 (vector/aberrations),
Lemma 8 (pixelation), Theorem 9 (rigorous B&B single-mask bound).

## PART V — Validation
results/RESULTS_TABLES.md (Tables 1–13), results/figs/fig1–6, results/logs/run_all.log.

## PART VI — Adversarial analysis
notes/ADVERSARIAL.md. Session-2 corrections: stale result files (2), multi-focus metric bug, rigorization null-space
flaw, SCS unreliability on degenerate instances, unconverged B&B at F ≥ 9 reported as such.

## PART VII — Novelty assessment — notes/NOVELTY.md.
## PART VIII — Manuscript — notes/PAPER_v2.md.
## PART IX — Reproducibility — README.md; code/; results/.
## PART X — Next steps
1. Interval-arithmetic B&B and/or SDP-bounded B&B to certify single-mask optima for F ≈ 9–61 (would resolve the 320 nm
   2D gaps); 2. non-periodic clips (windowing) and scalable SDP (Burer–Monteiro/chordal); 3. thick-mask (3D) models via
   boundary-layer linearizations (image remains quadratic) and rigorous EM validation; 4. NILS/MEEF/EPE constraints;
   5. benchmark against open ILT engines on ICCAD-13 clips; 6. external replication and submission.

## 11 Publication-readiness assessment
[[filled after clean run]]
