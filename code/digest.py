"""Print every number quoted in the manuscript, straight from results/*.json (run after assemble.py)."""
import json, os, numpy as np
R = '/home/user/Fable-Research/results/'
L = lambda n: json.load(open(R+n)) if os.path.exists(R+n) else None
e1 = L('exp1_1d_features.json'); e2 = L('exp2_random_tightness.json'); e3 = L('exp3_pixel_binary.json'); e3b = L('exp3b_ilt.json')
e5 = L('exp5_smo_euv_toeplitz.json'); e7 = L('exp7_lasserre.json'); e8 = L('exp8_robustness.json'); e9 = L('exp9_outlier_pw.json'); e10 = L('exp10_baselines_complex.json')
e11 = L('exp11_lasserre_pw.json'); e12 = L('exp12_lasserre_xcheck.json'); e13 = L('exp13_2d_strong.json'); e14 = L('exp14_bnb.json'); e14b = L('exp14b_bnb_big.json')
e15 = L('exp15_verify_all.json'); e16 = L('exp16_physics.json'); e17 = L('exp17_dataset.json'); e4 = L('exp4_2d.json')
P = print
P("== Table 1 (1D library)"); [P(f"  {r['name']:26s} bound {r['rigorous']:.4f} EL {r['EL_bound']:.3f} best {r['best_feasible']:.4f} gap {r['gap_rel']:.2%} rank {r['rank']}") for r in e1]
g = np.array([r['gap_rel'] for r in e2]); P(f"== Random 1D n={len(e2)}: <=1e-4 {np.mean(g<=1e-4):.1%}; <=1% {np.mean(g<=1e-2):.1%}; <=5% {np.mean(g<=5e-2):.1%}; median {np.median(g):.1e}; max {g.max():.3f}")
for s_ in sorted(set(r['src'] for r in e2)): P(f"   {s_}: n={sum(r['src']==s_ for r in e2)} exact {np.mean([r['gap_rel']<=1e-4 for r in e2 if r['src']==s_]):.0%}")
P("   gaps>5%:", [(r['trial'], round(r['gap_rel'],3), r['rank']) for r in e2 if r['gap_rel']>0.05])
P("== Lasserre (exp7/exp11):"); [P(f"  {r['name']} shor {r.get('rigorous_shor', r.get('shor')):.4f} las2 {r.get('lasserre2_bound', r.get('lasserre2')):.4f} best {r.get('best_single_mask', r.get('best_single')):.4f} {r.get('lasserre_status', r.get('status'))}") for r in (e7 or [])+(e11 or [])]
P("== Lasserre xcheck:", e12)
P("== B&B:"); [P(f"  {r['name']} F{r['F']} sdp {r['sdp_bound']:.4f} L {r['bnb_L']:.4f} U {r['bnb_U']:.4f} conv {r['bnb_converged']} nodes {r['bnb_nodes']} t {r['bnb_time']:.0f}s gain {r['certified_gain']:.2%}") for r in (e14 or [])+(e14b or [])]
P("== exp9 outliers:"); [P(f"  {k} F{r['F']} bound {r['bound']:.4f} best {r['new_best']:.4f} gap {r['gap_new']:.2%} facerank {r['min_face_rank']} las2 {r['lasserre2']}") for k, r in (e9 or {}).items()]
P("== PSM vs binary (exp3):"); [P(f"  {r['name']} fourier {r['t_fourier_psm']:.4f} pixel {r['t_pixel_psm']:.4f} DNN {r['t_pixel_nonneg_DNN']:.4f} -> binary penalty {(r['t_pixel_nonneg_DNN']/r['t_fourier_psm']-1):.1%}") for r in e3]
P("== ILT (exp3b):"); [P(f"  {r['name']} bound {r['certified_bound']:.4f} ILT psm {r['ilt_psm_best']:.4f} ({r['ilt_psm_best']/r['certified_bound']-1:.0%}) refined {r['ilt_psm_best_refined_on_ratio']:.4f} bin {r['ilt_binary_best']:.4f} ({r['ilt_binary_best']/r['certified_bound']-1:.0%}) refined {r['ilt_binary_best_refined_on_ratio']:.4f}") for r in e3b]
P("== SMO:"); [P(f"  {r['src']} fixed {r['fixed_source_bound']:.4f} smo {r['smo_bound']:.4f}") for r in e5['smo']]
P("== EUV 1D:"); [P(f"  p{r['pitch']} {r['src']} k1 {r['k1']:.3f} bound {r['rigorous']:.4f} EL {r['EL_bound']:.3f} best {r['best_feasible']:.4f} gap {r['gap_rel']:.2%}") for r in e5['euv_pw']]
P("== Theorem 4 numerics:", e5['theorem4'])
P("== Robustness: spacing", [(r['nx'], round(r['bound'],5)) for r in e8['sampling']], "margin", [(r['margin'], round(r['bound'],4), round(r['EL'],3)) for r in e8['margin']], "nsrc", [(r['nsrc'], round(r['clarabel'],4), round(r['diff'],9)) for r in e8['source_quadrature_and_solvers']], "mp", e8['mp_check'])
P("== Rayleigh/complex:"); [P(f"  {r['case']} {r['src']} sdp {r['sdp_real']:.4f} rayleigh {r['rayleigh_bound']:.4f} (x{r['sdp_real']/max(r['rayleigh_bound'],1e-9):.0f}) complex {r['sdp_complex']:.4f}") for r in e10]
P("== 2D (exp4):"); [P(f"  {r['name']} F{r['F']} bound {r['rigorous']:.4f} EL {r['EL_bound']:.3f} best {r['best_feasible']:.4f} gap {r['gap_rel']:.2%} rank {r['rank']} t {r['sdp_time']:.1f}s nulldim {r.get('cert_null_dim')}") for r in e4]
P("== 2D strong (exp13):"); [P(f"  {r['name']} bound {r['bound']:.4f} best {r['best_single']:.4f} gap {r['gap_single']:.2%} two-exp {r['two_exposure_ratio']:.4f} ({r['gap_two']:.2%}) facerank {r.get('face_rank_trace')}") for r in (e13 or [])]
if e15:
    d = [r['solver_diff'] for r in e15 if r['solver_diff'] is not None]
    P(f"== Verify: n={len(e15)} chol ok {sum(r['cholesky50'] for r in e15)} fails {[r['name'] for r in e15 if not r['cholesky50']]} min c {min(r['c'] for r in e15):.4f} median solver diff {np.median(d):.1e} big diffs {[(r['name'], round(r['t_clarabel'],4), round(r['t_scs'],4)) for r in e15 if r['solver_diff'] and r['solver_diff']>1e-5]}")
P("== Vector 1D:"); [P(f"  {r['case']} {r['src']} scalar {r['scalar_bound']:.4f}/{r['scalar_gap']:.2%} TE {r['TE_bound']:.4f}/{r['TE_gap']:.2%} TM {r['TM_bound']:.4f}/{r['TM_gap']:.2%} unpol {r['unpol_bound']:.4f}/{r['unpol_gap']:.2%}") for r in e16['vector_1d']]
P("== Vector 2D:"); [P(f"  {r['case']} scalar {r['scalar_bound']:.4f}/{r['scalar_gap']:.2%} unpol {r['unpol_bound']:.4f}/{r['unpol_gap']:.2%} TE {r['TE_bound']:.4f}/{r['TE_gap']:.2%}") for r in e16['vector_2d']]
P("== Aberr:", [(r['zernike'], r['waves'], round(r['bound'],4)) for r in e16['aberrations']]); P("== Pixelation:", [(round(r['pixel_nm']), round(r['bound'],4)) for r in e16['pixelation']])
if e17:
    r2 = e17['random2d']; g = np.array([r['gap_rel'] for r in r2]); P(f"== Random 2D n={len(r2)}: <=1e-3 {np.mean(g<=1e-3):.1%}; <=1% {np.mean(g<=1e-2):.1%}; <=5% {np.mean(g<=5e-2):.1%}; median {np.median(g):.1e}; max {g.max():.3f}; F {min(r['F'] for r in r2)}-{max(r['F'] for r in r2)}; by src: {{s: (sum(r['src']==s for r in r2), round(float(np.mean([r['gap_rel']<=1e-2 for r in r2 if r['src']==s])),2)) for s in set(r['src'] for r in r2)}}")
    P("== Adversarial 1D:"); [P(f"  {r['name']} {r['src']} F{r['F']} bound {r['rigorous']:.4f} EL {r['EL_bound']:.3f} best {r['best_feasible']:.4f} gap {r['gap_rel']:.2%}") for r in e17['adversarial1d']]
    P("== EUV 2D:"); [P(f"  p{r['pitch']} {r['src']} k1 {r['k1']:.3f} bound {r['rigorous']:.4f} EL {r['EL_bound']:.3f} best {r['best_feasible']:.4f} gap {r['gap_rel']:.2%}") for r in e17['euv2d']]
