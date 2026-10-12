# Certified bounds for partially coherent mask synthesis (lithography × convex optimization)

Autonomous research package produced 2026-10-09. Everything here is reproducible with Python ≥3.11, numpy, scipy,
cvxpy (CLARABEL/SCS), matplotlib, mpmath.

```
pip install numpy scipy cvxpy matplotlib mpmath
cd code
python3 -W ignore test_certlitho.py            # unit/regression tests (~1 min)
./run_all_1d.sh                                # all 1D experiments (~45 min on 4 cores)
python3 -W ignore exp4_2d.py                   # 2D layouts (~30 min)
python3 -W ignore exp3b_ilt.py exp10_baselines_complex.py exp11_lasserre_pw.py
python3 assemble.py                            # -> results/RESULTS_TABLES.md
python3 -W ignore figs.py fig2d.py             # -> results/figs/
```

Layout
- `code/certlitho.py` — model (Abbe/Hopkins scalar, periodic, diffraction-order domain), SDP relaxations (real/complex/
  pixel/DNN/SMO), dual certificates + rigorization, rounding, local solvers, optimal-face exploration, Lasserre-2.
- `code/pipeline.py` — end-to-end analysis; `code/exp*.py` — experiments; `code/figs.py`, `code/fig2d.py`.
- `notes/` — RESEARCH_STATE.md (persistent state + bug log), PHASE2_DIRECTIONS.md (five directions & scoring),
  THEORY.md (definitions, theorems, proofs), NOVELTY.md, ADVERSARIAL.md, PAPER_DRAFT.md (manuscript), FINAL_REPORT.md.
- `results/` — JSON outputs, RESULTS_TABLES.md, figs/, logs/.

## Manuscript
`paper/main.tex` + `paper/refs.bib` (LaTeX, article class, natbib). Compile with `cd paper && pdflatex main && bibtex main && pdflatex main && pdflatex main`
(figures are read from `../results/figs/`). A Markdown version is kept in `notes/PAPER_v2.md`.
