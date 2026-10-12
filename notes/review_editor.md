# Editorial change log: `paper/main.tex` -> `paper/main_edit.tex`

Scope: professional polish for journal submission (JM3 / Optics Express / SIIMS style). No scientific claim, number, theorem content, `\cite` key or `\label` was changed. `main.tex` is untouched; the edited manuscript is `paper/main_edit.tex` (compiles to `paper/main_edit.pdf`, 15 pages, 0 errors, 0 undefined references/citations, 0 overfull boxes, 0 BibTeX warnings).

Every item below goes beyond typography and should be approved by the authors. Items marked **[APPROVE]** add wording or facts that were not in the original and deserve explicit sign-off.

## 1. Front matter

1. **Author line** replaced by `Author names withheld (autonomous research system with human oversight)`; the `\thanks` footnote (repository and `run_all_par.sh` statement) is kept verbatim. Title and date line unchanged.
2. **Abstract rewritten as a structured abstract** (Purpose / Methods / Results / Conclusions), 249 words including the four labels (original: one 390-word paragraph). All key numbers are retained: $k_1$ cutoffs; 70 % / 83 % tightness on 109 random 1D targets; 0.3 % for 2D contacts; two of three 320 nm cells within 0.6 %; certified multi-exposure gains 6.4-16.4 % and 13.4 %; 3-150x tighter than Rayleigh quotient; binary masks forfeit 3-14 %; ILT 30-67 % above optimum. Dropped from the abstract for length (all remain in the body): the explicit mention of the 50-digit Cholesky / eigenvalue test, the complex-transmission lemma, the "survives moderate aberrations" claim is kept but shortened, and the parenthetical "(incoherent superposition)".

## 2. Structure and headings

3. Section titles made parallel (noun phrases, title case): "Related Work", "Model and Problem Statement", "Theory", "Algorithms", "Experiments", "Discussion", "Limitations". The old "Reproducibility" section was folded into a new unnumbered **"Data and Code Availability"** section (same commands and run-time statement, same content).
4. Section 3 split into subsections "Imaging model", "Synthesis problem and relaxation", and a new **"Notation"** subsection (`sec:notation`).
5. Experiment subsection titles made parallel: "Resolution cutoff" (was "Cutoff recovery"), "Random one-dimensional targets", "Phase-shift versus binary masks and pattern-error ILT", "Adversarial one-dimensional targets" (new label `sec:adv`), "Robustness, solver validation and physical sensitivity" (was "... physics sensitivity"). The experimental settings paragraph now has a `\paragraph{Settings.}` head.
6. Contributions list (C1)-(C4) and the Algorithms section converted from run-in `\emph{}` items to `itemize` / `\paragraph` heads. Related-work paragraphs given `\paragraph` heads; the two short paragraphs on information-theoretic limits and lifting were merged under one head ("Information-theoretic limits and lifting"). Content unchanged.
7. New unnumbered section **"Acknowledgements and Competing Interests"** **[APPROVE]**. Text added: "This manuscript was prepared by an autonomous research system under human oversight; the human overseers reviewed the claims, the verification record and the final text. No external funding was received for this work. The authors declare no competing interests." The funding and competing-interest statements are editorial placeholders and must be confirmed or amended by the authors.

## 3. Notation and symbols

8. New **Notation paragraph + Table 1** (`tab:notation`) listing every symbol used ($\Lambda$, $L_i$, $d$, $m$, $M_k$, $M$, $K$, $F$, $\lambda$, NA, $\sigma$, $s_j$, $w_j$, $P_z$, $z$, $W$, $Z$, $I_{M,z}(x)$, $Q_{x,z}$, $B$, $D$, $\tau$, mask classes, $\rho$, $\rho^\star$, $t_R$, EL, $\mu_d$, $\lambda_b$, $S$, $\tilde S$, $\Pi$, $J$, $U$, $X$, $\kappa$, $L$/$U$ of Theorem 9, $k_1$). The paragraph also fixes the meaning of "bound", "best single mask" and "gap" once, so later sections refer to it.
9. **Symbol clash $\lambda$ (wavelength) vs $\lambda_b$ (bright dual weight)**: not renamed, because the figure legends (Fig. 2, Fig. 5) use "lambda (bright)". Instead the Notation paragraph and table state explicitly that subscripted $\lambda_b$ is the dual weight and unsubscripted $\lambda$ the wavelength. *Authors may prefer renaming the dual weight (e.g. to $\nu_b$) in a future revision together with the figures.*
10. **Rigorization scale renamed $c \to \kappa$** (Algorithms "Rigorization" paragraph, robustness text, Table 7, notation table). Reason: $c$ was also the box center in Theorem 9 and the coefficient vector $c_{j,x,z}$ in Eq. (1). Purely notational.
11. **Common null space renamed $N \to \mathcal N$** in the Rigorization paragraph, because $N$ is the pixel count in Lemma 8. Purely notational.
12. In the proof of Theorem 4 a parenthetical was added: "(here $\sigma^2$ is the innovation variance, unrelated to the partial-coherence factor)". Clarification only.
13. Macros `\NA` and `\EL` introduced for $\mathrm{NA}$ and $\mathrm{EL}$ (typesetting only). The ratio definition received a label `eq:rho` (unused by `\ref`, harmless).

## 4. Theorems and proofs

14. Theorem/lemma/proposition statements are mathematically verbatim. Wording changes only: Theorem 3 title "Relaxation = unlimited multiple exposure" -> "Relaxation equals unlimited multiple exposure"; Theorem 4 title "Nonnegative masks in 1D" -> "Nonnegative masks in one dimension"; "$Q_1..Q_m$" typeset as "$Q_1,\dots,Q_m$"; "conversely if" -> "conversely, if"; "Lemma 6 extends for" -> "extends to"; proof of Theorem 9 items capitalized as sentences. Remark on attainment: "is still certified" -> "is nevertheless certified".

## 5. Tables (all now have self-contained captions with units and definitions; booktabs throughout)

15. **Table 2 (1D library, `tab:1d`)**: caption now defines $F$, foci ($|Z|$, with $\pm 80$ nm meaning $\{-80,0,+80\}$ nm), bound, max EL, best single, gap, and the target names; column heads capitalized. Numbers unchanged.
16. **New Table 3 (`tab:random`)**: the inline list of tightness fractions for the 109 random 1D targets (69.7 % / 82.6 % / 91.7 %, median $1.8\times10^{-6}$, max 17.3 %, per-source 85 / 75 / 57 / 57 %) moved into a table; the sentence in the text is kept.
17. **Table 4 (`tab:bnb`)**: caption defines SDP bound, $L$, $U$, nodes, converged, certified gain $=(L-\text{SDP bound})/\text{SDP bound}$, and the "random $n$" naming. "coh." expanded to "coherent". Numbers unchanged. **[APPROVE]** the stated formula for "certified gain" is the editor's reading of the column (it reproduces 12.9 %, 12.0 %, 16.4 %, 6.4 %, 13.4 % from the listed $L$ and SDP values); please confirm.
18. **New Table 5 (`tab:psm`)**: the inline numbers of the PSM-versus-binary / ILT paragraph (PSM bounds 0.3879/0.4031/0.3671; DNN bounds 0.4266/0.4599/0.3777; penalties 10/14/3 %; ILT 0.553/0.672/0.477 and 0.579/0.722/0.517; refined 0.426/0.530/0.3671) moved into a table. The sentence is retained. **[APPROVE]** the "ILT refined" column assigns 0.426 to iso-dense, 0.530 to mixed and 0.3671 to dense-120, following the iso-dense/mixed/dense-120 order used throughout the paragraph; the original text listed them as "(0.3671) ... (0.426, 0.530)".
19. **ILT percentage clarified [APPROVE]**: the original body text said "30-79 % above the certified optimum" while the abstract said "30-67 %". From the listed numbers, the PSM ILT masks are 30-67 % above the PSM bound and the binary ILT masks up to 79 % above the PSM bound. The body now reads "i.e. 30-67 % (PSM) and up to 79 % (binary) above the certified PSM optimum". No number changed; the referent of each range was made explicit.
20. **Table 6 (2D layouts, `tab:2d`)**: caption defines all columns and units (SDP time in seconds), and expands "strong". Numbers unchanged.
21. **New Table 7 (`tab:sens`)**: the long inline list of the robustness section (interior spacing, $\tau$, source quadrature, complex vs real, Rayleigh quotient, vector TE/unpolarized/TM, aberrations, pixelation, 50-digit Cholesky, $\kappa$, CLARABEL vs SCS) moved into a table; the prose summary is retained in condensed form. All numbers identical to the original sentence.
22. **Table 8 (Rayleigh, `tab:rayleigh`)**: caption now states the conditions (DUV, NA 1.35, $\tau$), that entries are lower bounds on $\rho$ (higher = tighter), and the Lemma 6 coincidence; one row per line in the source.

## 6. Figures

23. **Figures re-ordered to match first reference**: Fig. 1 `fig:k1` (cutoff), Fig. 2 `fig:cert` (1D library; previously never referenced in the text, now cited in Sec. 6.2), Fig. 3 `fig:gaphist` (new, see 24), Fig. 4 `fig:bnb` (multi-exposure gains; its right panel is cited again in Sec. 6.9), Fig. 5 `fig:2d`, Fig. 6 `fig:dataset`. In the original, `fig:dataset` was referenced before `fig:bnb` but numbered after it. Float placement parameters were relaxed and a `\FloatBarrier` inserted before the Discussion so that no float drifts behind the references.
24. **New Figure 3 (`fig:gaphist`) [APPROVE]**: `results/figs/fig3_gap_hist.png` existed but was not used in `main.tex`. It is now included in "Random one-dimensional targets" (histogram of $\log_{10}$ relative gap over the 109 instances). The caption's description ("bimodal; most instances tight to solver precision; instances with gaps $10^{-2}$-$10^{-1}$ have optimal-face rank 2 or 3") restates the text's own claims; please confirm it matches the figure's data.
25. **All captions expanded** to be self-contained: wavelength/NA/illumination, what each panel shows, the meaning of dashed lines and colors, and the relevant table/theorem. Figure 1 caption now says the left dashed line marks $k_1 = 0.25$ **[APPROVE]** (read from the plot) and the right dashed line is $k_1 = 1/(2(1+\sigma))$ (as in the original).
26. **Figure 5 (`fig:2d`) caption flags a figure/table inconsistency [APPROVE]**: the panel title embedded in `fig4_2d_Tjunction.png` reports "best single mask $\rho = 0.4789$; gap 9.17 %", whereas Table 6 reports 0.4413 (0.61 %) for the strong search. The caption now states that the mask shown is from the baseline search and that the strong search reaches 0.4413 against the same certificate. The authors may prefer to regenerate the figure from the strong-search mask instead (figures were not regenerated per instructions).
27. Figure 4 caption notes that for the two non-converged large-$F$ cases in the left panel (iso-dense coherent; mixed, $\sigma = 0.7$, $\pm 80$ nm) the branch-and-bound bound $L$ is not informative, as stated in the text. **[APPROVE]** the panel labels in the PNG ("mixed conv0.7 PW80 coh" etc.) overlap and are hard to read; consider regenerating with rotated labels.

## 7. Wording and consistency (substantive instances only)

28. Hyphenation unified: "phase-shift mask", "line-end", "dark-to-bright", "multi-exposure", "single-mask", "branch-and-bound" (adjective) / "branch and bound" (noun), "band-limited", "rank-one", "edge-offset", "one-/two-dimensional" (spelled out instead of "1D/2D" in headings, theorem names and prose; "1D"/"2D" survive only inside the verbatim figure-panel titles and the `\texttt{}` file names). US spelling adopted consistently ("color", "recognizes", "normalized").
29. The one-sentence Discussion paragraph was split into six sentences with no change of content. The Limitations list was turned into full sentences.
30. "3D mask electromagnetics" -> "three-dimensional mask electromagnetics"; "EPE" expanded to "edge-placement error"; "Prop." -> "Proposition"; "1,500 Adam steps" (thousands separator); "4-dimensional" -> "four-dimensional"; "whereas/and" inserted to resolve run-on sentences in Sec. 6.1, 6.6, 6.7.
31. The Related-Work sentence on Rosenbluth et al. and the Introduction sentence quoting three sources were re-punctuated (semicolons) without changing the quotations or citations.

## 8. Inconsistencies noticed but deliberately NOT changed (for the authors)

- Abstract says the bound is "3-150x tighter" than the Rayleigh-quotient bound; Sec. 6.9 and Table 7 say "1-148x weaker". From Table 8 the ratios are 1.46x (dense 90), 3.4x, 5.2x, 5.4x, 14x and 150x. Please harmonize (likely "1.5-150x").
- Table 6 row "T-junction 320, annular" (0.4413, 0.61 %) vs the embedded title of `fig4_2d_Tjunction.png` (0.4789, 9.17 %): see item 26.
- The "Settings" paragraph lists dipole $\sigma \in \{0.7, 0.8\}$ but only "dipole 0.8" appears in the tables; fine if 0.7 is used elsewhere in `results/`.
- `refs.bib` entry `nvidia2026` has author "Anonymous" and `liu1992` carries a "partially confirmed" note; both are reproduced as is.

## 9. Build

- Packages added: `caption` (small captions, bold labels), `array` (ragged-right paragraph columns in Table 7), `placeins` (`\FloatBarrier`). `siunitx` is not used; the custom `\SI`/`\SIrange` macros are kept.
- Command: `cd paper && pdflatex -interaction=nonstopmode main_edit && bibtex main_edit && pdflatex -interaction=nonstopmode main_edit && pdflatex -interaction=nonstopmode main_edit`.
- Result: 15 pages; `main_edit.log` contains 0 errors, 0 "Overfull", 0 undefined references or citations; `main_edit.blg` contains 0 warnings; all 23 bibliography entries are cited.
