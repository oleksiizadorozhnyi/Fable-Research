import json, glob, numpy as np, os
R = '/home/user/Fable-Research/results/'
def load(n):
    p = R + n
    return json.load(open(p)) if os.path.exists(p) else None
md = []
e1 = load('exp1_1d_features.json')
if e1:
    md.append("### Table 1 — 1D feature library (DUV 193 nm, NA 1.35; edge tolerance 8 nm; dense90: 6 nm)\n| case | F | foci | certified bound ρ | max EL (full range) | best single mask ρ | gap | SDP time (s) |\n|---|---|---|---|---|---|---|---|")
    for r in e1: md.append(f"| {r['name']} | {r['F']} | {len(r['focus'])} | {r['rigorous']:.4f} | {r.get('EL_full_bound', 2*(1-r['rigorous'])/(1+r['rigorous'])):.3f} | {r['best_feasible']:.4f} | {('degenerate' if r.get('gap_rel') is None else f"{r['gap_rel']:.2%}")} | {r['sdp_time']:.2f} |")
e7 = load('exp7_lasserre.json'); e11 = load('exp11_lasserre_pw.json')
rows = (e7 or []) + (e11 or [])
if rows:
    md.append("\n### Table 2 — Certified multi-exposure gain (Lasserre-2 single-mask bound vs SDP multi-exposure bound)\n| case | F | SDP (multi-exposure) bound | Lasserre-2 single-mask bound | best single mask found | certified gain | status |\n|---|---|---|---|---|---|---|")
    for r in rows:
        sb = r.get('rigorous_shor', r.get('shor')); lb = r.get('lasserre2_bound', r.get('lasserre2')); bs = r.get('best_single_mask', r.get('best_single'))
        gain = (lb - sb)/sb if lb else float('nan')
        md.append(f"| {r['name']} | {r['F']} | {sb:.4f} | {lb if lb is None else round(lb,4)} | {bs:.4f} | {gain:.2%} | {r.get('lasserre_status', r.get('status'))} |")
e2 = load('exp2_random_tightness.json')
if e2:
    g = np.array([r['gap_rel'] for r in e2]); md.append(f"\n### Random 1D targets (n={len(e2)}): gap ≤1e-4: {np.mean(g<=1e-4):.1%}; ≤1%: {np.mean(g<=1e-2):.1%}; ≤5%: {np.mean(g<=5e-2):.1%}; median {np.median(g):.1e}; max {g.max():.3f}; worst certificate min-eig {min(r['cert_mineig'] for r in e2):.1e}")
e3 = load('exp3_pixel_binary.json')
if e3:
    md.append("\n### Table 3 — Pixel vs Fourier relaxations, PSM vs binary (DNN), ILT baseline (σ=0.7 conventional)\n| case | F | Fourier PSM bound | pixel PSM bound | pixel nonneg (DNN) bound | ILT-PSM ratio | ILT-PSM refined | ILT-binary ratio | best PSM (pipeline) |\n|---|---|---|---|---|---|---|---|---|")
    for r in e3: md.append(f"| {r['name']} | {r['F']} | {r['t_fourier_psm']:.4f} | {r['t_pixel_psm']:.4f} | {r['t_pixel_nonneg_DNN']:.4f} | {r['ilt_psm_ratio']:.3f} | {r['ilt_psm_refined_on_ratio']:.4f} | {r['ilt_binary_ratio']:.3f} | {r['best_psm_pipeline']:.4f} |")
e3b = load('exp3b_ilt.json')
if e3b:
    md.append("\n### Table 3b — Standard gradient ILT (sigmoid resist, pattern-error loss; best of 6 settings) vs certified bound\n| case | F | certified bound ρ | ILT-PSM best ρ | ILT-PSM then ratio-refined | ILT-binary best ρ | ILT-binary then ratio-refined |\n|---|---|---|---|---|---|---|")
    for r in e3b: md.append(f"| {r['name']} | {r['F']} | {r['certified_bound']:.4f} | {r['ilt_psm_best']:.4f} | {r['ilt_psm_best_refined_on_ratio']:.4f} | {r['ilt_binary_best']:.4f} | {r['ilt_binary_best_refined_on_ratio']:.4f} |")
e5 = load('exp5_smo_euv_toeplitz.json')
if e5:
    md.append("\n### Table 4 — Source–mask entitlement bound (iso-dense, DUV)\n| source | #pts | fixed-source bound | SMO relaxation bound |\n|---|---|---|---|")
    for r in e5['smo']: md.append(f"| {r['src']} | {r['nsrc']} | {r['fixed_source_bound']:.4f} | {r['smo_bound']:.4f} |")
    md.append("\n### Table 5 — EUV (13.5 nm, NA 0.33, ±50 nm focus) dense lines: certified max EL (full range) vs pitch\n| pitch (nm) | k1 | source | F | bound ρ | max EL (full range) | best single mask ρ | gap |\n|---|---|---|---|---|---|---|---|")
    for r in e5['euv_pw']: md.append(f"| {r['pitch']} | {r['k1']:.3f} | {r['src']} | {r['F']} | {r['rigorous']:.4f} | {r.get('EL_full_bound', 2*(1-r['rigorous'])/(1+r['rigorous'])):.3f} | {r['best_feasible']:.4f} | {('degenerate' if r.get('gap_rel') is None else f"{r['gap_rel']:.2%}")} |")
    t4 = e5['theorem4']; md.append(f"\nTheorem 4 numerics: necessity violations {t4['necessity_violations']}/{t4['necessity_trials']}; sufficiency failures {t4['sufficiency_failures']}/{t4['sufficiency_tested']}.")
e8 = load('exp8_robustness.json')
if e8:
    md.append("\n### Robustness (iso-dense, σ=0.7)\n- interior spacing sweep (bound): " + ", ".join(f"{r['nx']} nm: {r['bound']:.5f}" for r in e8['sampling']))
    md.append("- margin sweep (bound / EL): " + ", ".join(f"{r['margin']} nm: {r['bound']:.4f}/{r['EL']:.3f}" for r in e8['margin']))
    md.append("- source quadrature & solver agreement: " + ", ".join(f"{r['nsrc']} pts: CLARABEL {r['clarabel']:.5f} SCS {r['scs']:.5f}" for r in e8['source_quadrature_and_solvers']))
    m = e8['mp_check']; md.append(f"- 50-digit certificate check: Cholesky success={m['cholesky_success']}, min eig={m['mineig_mp50']:.2e}, bound={m['bound']:.6f}")
e4 = load('exp4_2d.json')
if e4:
    md.append("\n### Table 6 — 2D layouts (DUV, NA 1.35), edge tolerance 12 nm (contacts120/100: 10/8 nm; T/line-end: 14 nm)\n| layout | cell (nm) | source pts | F | samples | bound ρ | max EL (full range) | best single mask ρ | gap | rank | cert. null dim | SDP time (s) |\n|---|---|---|---|---|---|---|---|---|---|---|---|")
    for r in e4: md.append(f"| {r['name']} | {r['L']} | {r['nsrc']} | {r['F']} | {r['nB']+r['nD']} | {r['rigorous']:.4f} | {r.get('EL_full_bound', 2*(1-r['rigorous'])/(1+r['rigorous'])):.3f} | {r['best_feasible']:.4f} | {('degenerate' if r.get('gap_rel') is None else f"{r['gap_rel']:.2%}")} | {r['rank']} | {r.get('cert_null_dim', 0)} | {r['sdp_time']:.1f} |")
e9 = load('exp9_outlier_pw.json')
if e9:
    md.append("\n### Outlier / process-window re-examination (strong primal search, face exploration, Lasserre-2 where F≤11)\n| case | F | bound | old best | new best | new gap | min face rank | Lasserre-2 | status |\n|---|---|---|---|---|---|---|---|---|")
    for k, r in e9.items(): md.append(f"| {k} | {r['F']} | {r['bound']:.4f} | {r.get('old_best', float('nan')):.4f} | {r['new_best']:.4f} | {r['gap_new']:.2%} | {r['min_face_rank']} | {r['lasserre2']} | {r['lasserre_status']} |")
e10 = load('exp10_baselines_complex.json')
if e10:
    md.append("\n### Table 7 — Baseline Rayleigh-quotient bound vs SDP; complex-transmission masks\n| case | source | F | SDP real-mask bound | Rayleigh-quotient bound | SDP complex-mask bound |\n|---|---|---|---|---|---|")
    for r in e10: md.append(f"| {r['case']} | {r['src']} | {r['F']} | {r['sdp_real']:.4f} | {r['rayleigh_bound']:.4f} | {r['sdp_complex']:.4f} |")

e13 = load('exp13_2d_strong.json')
if e13:
    md.append("\n### Table 8 — 320 nm 2D layouts: strong single-mask search, optimal-face rank, and two-exposure realization\n| layout | F | bound ρ | min face rank | best single mask ρ | single gap | two-exposure ρ | two-exposure gap |\n|---|---|---|---|---|---|---|---|")
    for r in e13: md.append(f"| {r['name']} | {r['F']} | {r['bound']:.4f} | {r.get('face_rank_trace', r.get('min_face_rank'))} | {r['best_single']:.4f} | {r['gap_single']:.2%} | {r['two_exposure_ratio']:.4f} | {r['gap_two']:.2%} |")
e14 = load('exp14_bnb.json')
if e14:
    md.append("\n### Table 9 — Rigorous branch-and-bound certification of the single-mask optimum (real masks) and certified multi-exposure gain\n| case | F | SDP multi-exposure bound | B&B certified single-mask lower bound L | best single mask U | B&B gap (U-L)/U | nodes | converged | time (s) | certified gain (L-bound)/bound |\n|---|---|---|---|---|---|---|---|---|---|")
    for r in e14: md.append(f"| {r['name']} | {r['F']} | {r['sdp_bound']:.4f} | {r['bnb_L']:.4f} | {r['bnb_U']:.4f} | {(r['bnb_U']-r['bnb_L'])/r['bnb_U']:.2%} | {r['bnb_nodes']} | {r['bnb_converged']} | {r['bnb_time']:.1f} | {r['certified_gain']:.2%} |")
e15 = load('exp15_verify_all.json')
if e15:
    import numpy as _np
    ok = sum(1 for r in e15 if r['cholesky50']); deg = [r['name'] for r in e15 if not r['cholesky50']]
    d = [r['solver_diff'] for r in e15 if r['solver_diff'] is not None]
    md.append(f"\n### Independent certificate verification (exp15): {len(e15)} instances recomputed; 50-digit Cholesky of the reduced certificate matrix succeeded on {ok}; degenerate (S≈0, bound=1, uniform certificate) on {deg}; min rigorization scale c = {min(r['c'] for r in e15):.4f}; CLARABEL/SCS agreement: median |Δt| = {_np.median(d):.1e}, instances with |Δt|>1e-5: " + ", ".join(f"{r['name']} (CLARABEL {r['t_clarabel']:.4f}, SCS {r['t_scs']:.4f})" for r in e15 if r['solver_diff'] is not None and r['solver_diff']>1e-5))
e16 = load('exp16_physics.json')
if e16:
    md.append("\n### Table 10 — Vector (high-NA, n_resist=1.7) vs scalar bounds, 1D\n| case | source | scalar bound | scalar gap | TE bound | TE gap | TM bound | TM gap | unpol bound | unpol gap |\n|---|---|---|---|---|---|---|---|---|---|")
    for r in e16['vector_1d']: md.append(f"| {r['case']} | {r['src']} | {r['scalar_bound']:.4f} | {r['scalar_gap']:.2%} | {r['TE_bound']:.4f} | {r['TE_gap']:.2%} | {r['TM_bound']:.4f} | {r['TM_gap']:.2%} | {r['unpol_bound']:.4f} | {r['unpol_gap']:.2%} |")
    md.append("\n### Table 11 — 2D contacts: vector vs scalar\n| case | scalar bound | scalar gap | unpol bound | unpol gap | TE bound | TE gap |\n|---|---|---|---|---|---|---|")
    for r in e16['vector_2d']: md.append(f"| {r['case']} | {r['scalar_bound']:.4f} | {r['scalar_gap']:.2%} | {r['unpol_bound']:.4f} | {r['unpol_gap']:.2%} | {r['TE_bound']:.4f} | {r['TE_gap']:.2%} |")
    md.append("\n### Aberration sensitivity (iso-dense, σ=0.7): " + "; ".join(f"{r['zernike']} {r['waves']}λ: {r['bound']:.4f}" for r in e16['aberrations']))
    md.append("\n### Pixelation (iso-dense, σ=0.7, F=15): " + "; ".join(f"{r['pixel_nm']:.0f} nm pixels (N={r['N']}): {r['bound']:.4f}" for r in e16['pixelation']) + f" (Fourier {e16['pixelation'][0]['fourier']:.4f})")
e17 = load('exp17_dataset.json')
if e17:
    import numpy as _np
    r2 = e17['random2d']
    if r2:
        g = _np.array([r['gap_rel'] for r in r2]); md.append(f"\n### Random 2D layouts (n={len(r2)}; cells 200/240 nm; annular/quasar/conv0.5; τ=12 nm): gap ≤1e-3: {_np.mean(g<=1e-3):.1%}; ≤1%: {_np.mean(g<=1e-2):.1%}; ≤5%: {_np.mean(g<=5e-2):.1%}; median {_np.median(g):.1e}; max {g.max():.3f}; F range {min(r['F'] for r in r2)}–{max(r['F'] for r in r2)}; median SDP time {_np.median([r['sdp_time'] for r in r2]):.1f} s")
    a1 = e17['adversarial1d']
    if a1:
        md.append("\n### Table 12 — Adversarial 1D targets\n| target | source | F | bound ρ | max EL (full range) | best single ρ | gap | rank |\n|---|---|---|---|---|---|---|---|")
        for r in a1: md.append(f"| {r['name']} | {r['src']} | {r['F']} | {r['rigorous']:.4f} | {r.get('EL_full_bound', 2*(1-r['rigorous'])/(1+r['rigorous'])):.3f} | {r['best_feasible']:.4f} | {('degenerate' if r.get('gap_rel') is None else f"{r['gap_rel']:.2%}")} | {r['rank']} |")
    e2d = e17['euv2d']
    if e2d:
        md.append("\n### Table 13 — EUV 2D contact arrays (13.5 nm, NA 0.33, ±40 nm focus, τ = 0.1 pitch)\n| pitch | k1 | source | F | bound ρ | max EL (full range) | best single ρ | gap | rank |\n|---|---|---|---|---|---|---|---|---|")
        for r in e2d: md.append(f"| {r['pitch']} | {r['k1']:.3f} | {r['src']} | {r['F']} | {r['rigorous']:.4f} | {r.get('EL_full_bound', 2*(1-r['rigorous'])/(1+r['rigorous'])):.3f} | {r['best_feasible']:.4f} | {('degenerate' if r.get('gap_rel') is None else f"{r['gap_rel']:.2%}")} | {r['rank']} |")
open(R + 'RESULTS_TABLES.md', 'w').write("\n".join(md)); print("\n".join(md[-12:]))
