# Numerical audit of paper/main.tex (and notes/PAPER_v2.md, notes/FINAL_REPORT.md) against results/*.json

Audited 2026-10-12 against results/*.json, results/RESULTS_TABLES.md (assembled 06:56), code/exp*.py, results/logs/, results/figs/.
Severity: CRITICAL = provenance/validity of a headline number; MAJOR = wrong or overstated number/claim; MINOR = rounding, range edge, wording; OK = verified.

## CRITICAL

1. **Clean run did not complete: three scripts were OOM-killed and their JSONs are from the pre-clean session.**
   - Evidence: results/logs/run_all_stdout.log has three `Killed` lines (PIDs 5657, 7382, 6641); results/logs/exp14b_bnb_big.log is 0 bytes; exp7_lasserre.log stops after 2 of 3 candidates (last write 03:53, process "ended" 04:25); exp9_outlier_pw.log stops after trial 46 (no process-window rows; last write 05:37). JSON mtimes: exp7_lasserre.json 2026-10-11 23:30, exp9_outlier_pw.json 2026-10-11 23:33, exp14b_bnb_big.json 2026-10-12 00:17 — all before the run started (03:13).
   - run_all.log's `exit=0` is meaningless: in `echo "=== END $e $(date) exit=$?"` the `$(date)` substitution resets `$?` to 0 before it is read (verified). Fix: `rc=$?; echo "... exit=$rc"`.
   - Affected manuscript numbers: (a) abstract, Table~bnb row "mixed, coherent", §6.4, Discussion, FINAL_REPORT, PAPER_v2 §6.4: `3.0×10^7` nodes, `L=0.4776`, `0.23% gap`, `13.4%` certified gain, "35 min" (exp14b only); (b) §6.4 "each optimal face has rank 2" and "Lasserre-2 single-mask bounds coincide with U to 1e-6" (exp9; the clean-run log does reproduce face rank 2 and the Lasserre values for trials 82/91/110/50 before the kill, but the JSON on disk is the old one and lacks trial 46 / PW rows); (c) §6.4 and Table 2 of RESULTS_TABLES: Lasserre-2 0.4608 / 27% (exp7; clean-run log reproduced both values before the kill; third candidate mixed_conv0.7, F=11, never ran); (d) Fig. 5 left panel (fig5.py reads exp14b_bnb_big.json).
   - Footnote on the title page ("Every number in this manuscript is produced by code/run_all_par.sh from scratch") and §9 ("regenerate every table and figure") are false as stated.
   - Recommended: rerun exp7, exp9, exp14b sequentially (not in the 4-way parallel groups; exp7/exp9 Lasserre-2 and exp14b batch=16384 B&B exceed memory when co-scheduled), fix the exit capture in run_all_par.sh, re-assemble, regenerate fig5; or change the footnote to "...except exp7/exp9/exp14b, whose outputs are from the 2026-10-11 session (logs: exp7.log, exp9.log, exp14b.log)".

2. **Rigorization procedure described in §5 is not the one that produced the tables.**
   - code/certlitho.py was modified at 06:46:47, after `ALL DONE 06:46:32`, adding `rigorize_certificate` v3 (strict margin λmin ≥ +1e-10‖S‖ on the complement of the common null space, null_tol 1e-12) — committed c4d1d4e 06:55. Only exp15 was rerun with v3 (JSON 06:53, after assembly at 06:46; tables re-assembled 06:56).
   - Every other `rigorous` value (exp1, exp2, exp3b, exp4, exp4b, exp5, exp8, exp10, exp11, exp13, exp14, exp16, exp17) was produced by the full-space version at certlitho.py:745 (`λmin ≥ −1e-11‖S‖`, no projection) — i.e. a tolerance that permits a slightly *negative* eigenvalue. FINAL_REPORT §0 ("full-space eigenvalue test with 1e-11 tolerance") describes that version; main.tex §5 describes v3.
   - Consequence: the 320 nm cells (Table 2D rows 6–8), exp3 PSM/DNN bounds, exp16 vector/aberration/pixelation bounds, exp17 (66 bounds), SMO, dense90 ×3 and mixed_conv0.7_PW80 have never passed the strict-margin test or the 50-digit Cholesky. Differences where both exist are ≤1e-8 relative (e.g. isodense_coh 0.36251708581 v2 vs 0.36251708548 v3), so values will not change at 4 d.p., but the text must say which procedure produced which numbers, or run_all_par.sh must be rerun with v3.
   - Recommended text for §5: "Tables were produced with a full-space test λmin(S_D − cS_B) ≥ −10^{-11}‖·‖; the 136 instances of §6.8 were additionally re-rigorized with the strict margin +10^{-10} on N^⊥ and factorized at 50 digits." Or rerun.

## MAJOR

3. **§6.8 vector-model percentages wrong** (also PAPER_v2 abstract and §6.10).
   - Text: "TE bounds are 0–4% below scalar, unpolarized 10–40% above, TM 70–120% above".
   - exp16 vector_1d: TE −2.4/−0.0/−3.9/+0.0%; TM +42.1/+100.7/+63.1/+67.9%; unpol +18.1/+43.8/+25.3/+28.0%. vector_2d: TE +2.3% (contacts120_annular), +10.0% (contacts120_quasar) — *above* scalar; unpol +10.3/+18.1%.
   - Replace: "In 1D, TE bounds are 0–4% below scalar, unpolarized 18–44% above and TM 42–101% above (classical TM contrast loss); for the 2D contacts, unpolarized is 10–18% and TE 2–10% above scalar; tightness is unchanged." Fig. 5 right panel is consistent with the data (it is the text that is wrong).

4. **Rayleigh-quotient factor inconsistent and wrong at the low end.** Abstract: "3–150× tighter"; §6.8: "1–148× weaker"; PAPER_v2 §6.9 "3–150×"; FINAL_REPORT "1–148×". exp10: ratios 5.43, 3.38, 5.24, 148.3, 14.1, 148.3, 2.67, 1.46, 2.67. Table~rayleigh itself lists dense90 σ=0.7 (0.6428 vs 0.4392 = 1.46×). Replace both with "1.5–150×" (or "between 1.5× and 150×").

5. **§6.8 rigorization scale c claim false.** Text: "the rigorization scale c exceeds 0.99999 except on one near-degenerate random instance (bound 1e-3, c=0.835)". exp15: 17 of 136 instances have c < 0.99999; five have c < 0.9999: random_trial16 c=0.8354 (bound 1.25e-3), trial75 0.9897 (3.4e-3), trial87 0.99924 (1.3e-3), trial103 0.99925 (0.142), trial105 0.99964 (5.4e-3). Replace: "c ≥ 0.9999 on 131 of 136 instances; the exceptions are five random instances with c = 0.835, 0.990, 0.9992, 0.9993 and 0.9996 (four with bounds ≤ 6×10^{-3})". PAPER_v2 §6.9 same sentence.

6. **§6.8 solver-agreement sentence stale.** Text: "CLARABEL and SCS agree to median 7×10^{-9}; on five near-degenerate random instances (bounds ≤10^{-2}) SCS is off by up to 0.14 and once negative". exp15: median |Δt| = 4.5e-9; six instances with |Δt| > 1e-5: trial16 (0.0015 vs 0.0002), 75 (0.0034 vs 0.0020), 87 (0.00134 vs 0.00128), 103 (0.1419 vs 0.0004), 105 (0.0054 vs −0.0035), 119 (0.0063 vs 0.0001). trial103's bound is 0.142, not ≤1e-2; it is the "0.14" case. Replace: "median 4.5×10^{-9}; on six random instances SCS deviates by more than 10^{-5}: five near-degenerate ones (bounds ≤6×10^{-3}, once negative) and one with bound 0.142 where SCS returns 4×10^{-4}". PAPER_v2 says "six" and "4.5e-9" (correct); main.tex is the stale one.

7. **§6.7 (2D) null-space sentence stale.** Text: "The weighted dark matrix has a 4-dimensional common null space (orders invisible to every constraint); the certificate is verified on its complement." PAPER_v2 §6.8: "4–16-dimensional null space ... bright matrix vanishes on it to 1e-31". Those numbers come from the 2026-10-09 run (exp4b.log: null_dim 16 quasar, 4 line-end). Current data: exp4b_recert.log and exp4_2d.json give cert_null_dim = 0 for Tjunction320_quasar and lineend320_annular; exp15 has null_dim = 0 on all 136 instances; Tjunction320_annular was never recertified (cert_null_dim absent). Delete the sentence or replace with "Under the v3 criterion no common null space is detected (null_dim = 0) and the certificate is verified on the full space".

8. **§6.5 / abstract / Discussion: "binary masks forfeit a certified 10%/14%/3%" compares two lower bounds.** DNN (binary lower bound) 0.4266/0.4599/0.3777 vs PSM *bound* 0.3879/0.4031/0.3671 — the PSM bound is the multi-exposure optimum (Thm. me), not necessarily attainable by one PSM. Best single PSM found (exp3 `best_psm_pipeline`): 0.3981/0.4036/0.3671, so the certified single-binary vs single-PSM penalty is ≥7.2%/14.0%/2.9%. Replace with "relative to the PSM relaxation (multi-exposure optimum) binary masks forfeit a certified 10/14/3%; relative to the best single PSM found, ≥7/14/3%". Abstract "3–14%" → "3–14%" is still defensible only in the multi-exposure sense; say "at least 3–14% relative to the phase-shift relaxation". PAPER_v2 abstract "≥10% on an iso-dense target" → "≥7%".

9. **§6.3 "every gap above 5% has minimal optimal-face rank 2 or 3".** exp2's `rank` is `rank_eps` of the solver's X (numerical rank), not a minimal face rank; exp9's face exploration gives min_face_rank = 4 for trial16 (gap 14.8%, F=13) in both the stale JSON and the clean-run log. Replace: "the solver's optimal X has numerical rank 2–3 on every instance with gap > 5%; explicit face exploration finds minimal face rank 2 on the four F=5 instances and 4 on trial 16".

10. **Limitations: "the branch-and-bound converges only for F ≤ 7"** and FINAL_REPORT "rigorous B&B certificates for F ≤ 7". exp14/exp14b: converged = True only for F=5; F=7 stopped at 3.0×10^7 nodes unconverged (gap 0.23% > tol 0.2%), F=9/11 bounds vacuous (0.065, ≈0). Replace: "converges only for F = 5; at F = 7 it returns a valid but unconverged bound (0.23% gap after 3×10^7 nodes); at F ≥ 9 the partial bounds are vacuous". PAPER_v2 §5.4 wording ("converges for F = 5 ... valid partial bounds for F = 7–11") is the correct one.

11. **§6.7 "40 random two-rectangle layouts"**: exp17 draws 1–3 rectangles per cell (13/11/16 layouts with 1/2/3 rectangles; sides ≥ 50 nm, τ = 12 nm). "gaps up to 26%" → max 0.265 = 27%. Replace: "40 random layouts of one to three rectangles (sides ≥ 50 nm) in 200/240 nm cells, τ = 12 nm ... up to 27%". ("mostly σ=0.5 and F=37" is correct: all 15 gaps >5% are σ=0.5 (F=13/21) or F=37.)

12. **§6.8 "All 136 reported bounds were recomputed"** overstates coverage. exp15 covers: 6 single-focus 1D library cases + 2 PW80 (isodense conv0.7/dipole0.8) + 109 random + 14 EUV 1D + 5 contact arrays = 136. Not covered: dense90 ×3, mixed_conv0.7_PW80, the three 320 nm cells, exp3 (6 bounds), exp5 SMO (4), exp16 (~30), exp17 (66). The paper reports > 230 bounds. Replace "All 136 reported bounds" with "136 of the reported bounds (the single-focus 1D library, two three-focus cases, all 109 random targets, the 14 EUV cases and the five contact arrays)".

13. **Fig. 4 contradicts Table 2D.** fig2d.py plots exp4 results (6-start, no polish): the printed suptitle reads "best single mask ρ=0.4789; gap 9.17%" while Table 2D/§6.7 report 0.4413 / 0.61% for the same cell, and the caption calls the plotted mask the "best band-limited PSM". Regenerate from results/exp13_Tjunction320_annular_rbest.npy (strong search) or caption it as the exp4 mask and remove "best".

14. **§6.5 "gradient verified by finite differences"**: no finite-difference check exists anywhere in code/ (grep for finite/check_grad/numerical gradient in certlitho.py, pipeline.py, test_certlitho.py, exp3b_ilt.py returns nothing). Remove the clause or add the test.

15. **ILT range inconsistent.** Abstract: "30–67% above the certified optimum"; §6.5: "30–79%"; FINAL_REPORT "30–79%"; PAPER_v2 abstract "30–67%". exp3b: PSM +42/+67/+30%, binary +49/+79/+41%. Use "30–67% (PSM) and 41–79% (binary)" everywhere.

## MINOR

16. §6.4 "converges in 0.3–1.6 s for F=5": exp14 times 1.20, 1.72, 0.32, 0.44 s → "0.3–1.7 s".
17. Abstract & FINAL_REPORT "two of three 320 nm cells are certified within 0.6%": gaps 0.61% and 0.03% → "within 0.7%" or "to 0.61% and 0.03%".
18. §5 Complexity "0.1–1.5 s for F≤15": single-focus 0.12–0.88 s, three-focus 1.35–2.93 s (exp1); "100–650 s for F=61" → 269–649 s. PAPER_v2 §5 "100–200 s" and §6.8 "each SDP solves in < 0.5 s" (contacts: 0.5–1.5 s) are stale.
19. §6.8 "source quadrature converges (0.513–0.526 for 11–21 points)": exp8 0.5161/0.5146/0.5261 → "0.515–0.526".
20. §6 Settings line: EUV 2D contacts (exp17) use annular 0.5/0.9 and quasar 0.5/0.9, not 0.6/0.9; EUV 1D τ = 0.12·pitch (exp5) is never stated ("τ = 8 nm unless stated"); EUV 2D τ = 0.1·pitch stated only in RESULTS_TABLES. Fig. 2 right panel uses τ = 0.1·pitch, 50% duty, raw (un-rigorized) SDP value — caption should say so.
21. §6.7 adversarial: "show 12–37% gaps under coherent illumination (rank 3–5)" holds for 3 of the 6 coherent cases (nearcutoff 12.2%, alternating 14.2%, three_pitches 36.5%); sraf/duty-5% are tight, duty-95% degenerate. Say "three of the six coherent cases". Table 12 row extreme_duty_95pct/dipole shows gap 1592% (bound 1.0e-6, best 1.8e-5): "tight under dipole" needs the caveat "bound and best mask both below 2×10^{-5}".
22. §6.4: the F=11 Lasserre case is `mixed_conv0.7_PW80_edgeonly` (6 bright + 6 dark edge samples, best single 0.9568), not the Table 1 three-focus mixed case (full sample set, best single 0.9742). State "edge-offset samples only".
23. Theorem~me: main.tex "r(r+1)/2 ≤ |D×Z|+|B×Z|"; PAPER_v2 "+1". Align (with t as a free variable and slacks, m = |D×Z|+|B×Z| equalities; main.tex is consistent with its own proof).
24. §9 "9 tests" correct (test_certlitho.py has 9); PAPER_v2 §9 says 8. PAPER_v2 §9 and FINAL_REPORT gate E cite `run_all.sh`; main.tex cites `run_all_par.sh`. README does not mention run_all_par.sh and lists exp11_lasserre_pw.py (the runner executes exp11b).
25. PAPER_v2 §6.8 "0.4386 (annular)" → 0.4387 (0.43865000; Table 2D says 0.4387).
26. Fig. 5 left: 8 bars (exp14 ∪ exp14b) including a duplicate "mixed coh" (6M- and 30M-node runs) and the vacuous F=9/F=11 bars; overlapping x-labels; caption cites Table~bnb (5 rows). Right panel legend covers the isodense/dipole TM bar.
27. RESULTS_TABLES "50-digit certificate check: Cholesky success=False, min eig=−8.25e−12" (exp8 mp_check, full-space matrix at the −1e-11 tolerance) sits next to the paper's "50-digit Cholesky ... 135 succeed" (exp15, projected, strict margin). Not cited in main.tex, but should be explained or dropped from the tables.
28. FINAL_REPORT §0 "1e-11 tolerance" vs main.tex §5 "1e-10 margin": see item 2.
29. Title-page date/footnote: PDF (paper/main.pdf) is modified but uncommitted; README says package "produced 2026-10-09".

## OK (verified against JSON; rounding only)

- §6.1: cutoff 193/(1.7·1.35) = 84.10 nm; k1(20 nm) = 0.244; EUV 22 nm conv0.9 EL 0.362; dipole0.8 bound = 1.0000 at 20 and 22 nm; 14/14 EUV cases gap 0.00%.
- Table 1: all 12 rows (F, foci, bound, EL, best, gap) match exp1; caption τ = 8 nm, dense-90 6 nm matches exp1 (margin=6). Gaps 30.78/13.68/2.50/0.22/0/0/11.22/12.09/16.06/0/0/0 %.
- §6.3: n = 109 (120 trials, 11 skipped), 69.7% / 82.6% / 91.7%, median 1.8e-6, max 17.3%; per-source 85% (conv0.9, n=34), 75% (dipole0.7, n=24), 57% (coh n=28, conv0.5 n=23); 1–4 segments ≥ 40 nm, cells 300/400/500/600, τ = 10 nm, spacing 5 nm — all match exp2.py.
- Table bnb rows 1–4 (F=5): SDP 0.3597/0.6049/0.3464/0.3141, L 0.4061/0.6775/0.4033/0.3343, U 0.4066/0.6782/0.4037/0.3347, nodes 40,383/76,027/17,413/28,517, converged, gains 12.92/12.00/16.44/6.45% → 12.9/12.0/16.4/6.4%; tol 1e-3, budget 6e6, 1+1e-9 inflation and 1e-12 subtraction match bnb_cert.py. Row 5 values match exp14b (subject to item 1).
- Lasserre-2 vs U: 2.3e-8, 3.8e-9, 1.0e-8, 1.4e-7 (≤1e-6 holds); CLARABEL 0.4787045 / SCS 0.4787037 (exp12); 0.4608, 0.9552, statuses optimal_inaccurate, gains 27.1%/9.9% (exp7/exp11).
- §6.5: DNN 0.4266/0.4599/0.3777, PSM 0.3879/0.4031/0.3671, τ = 10 nm (exp3 margin=10); pixel vs Fourier 4.0e-6/4.4e-6/8.1e-7 (≤1e-5); ILT 0.5526/0.6715/0.4768 and 0.5789/0.7218/0.5167; refined 0.4257/0.5302/0.3671; sigmoid resist, threshold 0.3, open-frame units, 1500 Adam steps, 3×2 = 6 settings — match certlitho.ilt_gradient / exp3b.py.
- §6.6: SMO 0.3625 for conv0.7 (7 pts) and 13-point σ=0.9 grid; fixed-source 0.4313–0.5091.
- Table 2D: all 8 rows (F, bound, EL, best, gap, SDP time 1.4/1.5/0.7/0.7/0.5/410/269/649 s) match exp4 + exp13; τ = 12/10/8/14 nm and 12 starts match exp4.py/exp13.py; 1210/1189 samples ≈ 1,200; two-exposure refinement 0.4263 > 0.4135 (did not beat single mask).
- §6.7: random 2D 45.0% / 62.5%, F 13–37; EUV 2D EL 0.131–0.412, gaps ≤ 3e-7, ±40 nm; isolated 20 nm space (extreme_duty_95pct) bounds 1.5e-5 (coh), 1.0e-6 (dipole).
- §6.8: spacing sweep 0.46907 ×5; margin sweep EL 0.094 → 0.737; complex = real on 9/9 cases (max |Δ| 1.1e-5); exp15 n = 136, Cholesky 135 True + 1 degenerate (EUV_p20_dipole0.8), 0 failures; aberrations max +6.9% (Z9, 0.1λ); pixelation 0% for N ≥ 15, +15.3% (50 nm), +46.3% (67 nm); n_resist 1.7, Z5/Z7/Z9 at 0.02/0.05/0.1λ match exp16.py.
- Table rayleigh: 6 rows match exp10.
- §9: run 03:13 → 06:46 = 3 h 33 min on 4 groups ("about 3.5 hours") — but see item 1.
- Figures: fig1–fig6 written 06:46:12–06:46:32, after assembly (06:46:04) and after every JSON they read (exp5 03:18, exp2 04:39, exp4 04:39, exp14 06:42, exp16 04:26, exp17 05:43); fig5 reads the stale exp14b (00:17) (item 1); exp15 (06:53) feeds no figure; RESULTS_TABLES.md (06:56) post-dates the exp15 rerun.
