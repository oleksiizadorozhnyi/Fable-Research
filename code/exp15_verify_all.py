# Independent high-precision verification of every reported certificate: recompute SDP, rigorize, then verify
# S_c = SD - c SB (restricted to complement of null space of SD) is PSD via mpmath Cholesky at 50 digits; also check
# both solvers agree. Covers exp1 (1D features), exp2 (random), exp5 EUV, exp4 2D contacts.
from pipeline import *
import mpmath as mp
mp.mp.dps = 50
lam, NA = 193.0, 1.35; out = []
def verify(name, rm, B, D):
    try: r = sdp_bound_real(rm, B, D, solvers=("CLARABEL",))
    except Exception as e:                                  # near-degenerate instance: CLARABEL fails; record it, verify the SCS certificate instead
        try: r = sdp_bound_real(rm, B, D, solvers=("SCS",)); r['solver'] = 'SCS (CLARABEL failed)'
        except Exception as e2:
            rec = dict(name=name, F=rm.F, solver_failed=True, error=str(e2)[:120]); print(json.dumps(rec), flush=True); out.append(rec); json.dump(out, open('/home/user/Fable-Research/results/exp15_verify_all.json','w'), indent=1); return
    rig, me, c = rigorize_certificate(rm, r['mu'], r['lam'])
    try: r2 = sdp_bound_real(rm, B, D, solvers=("SCS",)); t2 = r2['t']
    except Exception: t2 = None
    S, ndim = certificate_matrix_projected(rm, r['mu'], r['lam'], c)
    Sr = np.real(S) if np.allclose(S.imag, 0) else S
    degenerate = np.linalg.norm(Sr, 2) < 1e-13
    try:
        if degenerate: chol = None
        else: mp.cholesky(mp.matrix(Sr.tolist())); chol = True
    except Exception: chol = False
    rec = dict(name=name, F=rm.F, t_clarabel=r['t'], t_scs=t2, rigorous=rig, c=c, rel_mineig=float(me), cholesky50=chol, degenerate=bool(degenerate), null_dim=ndim,
               resid_null=float(rigorize_certificate.last_resid_null), solver_diff=(abs(r['t']-t2) if t2 is not None else None), rel_gap_t_vs_rigorous=(r['t']-rig)/max(r['t'],1e-300))
    print(json.dumps(rec), flush=True); out.append(rec); json.dump(out, open('/home/user/Fable-Research/results/exp15_verify_all.json','w'), indent=1)
conv = SRC_DUV['conv0.7']; dip = SRC_DUV['dipole0.8']; coh = SRC_DUV['coh']
for nm, src in [("coh",coh),("conv0.7",conv),("dipole0.8",dip)]:
    for cn, L, segs in [("isodense", 600, [(0,50),(100,150),(200,250),(430,480)]), ("mixed", 500, [(0,40),(100,180),(240,400)])]:
        md, B, D, tgt = build_1d(L, lam, NA, src, segs, TAU_DUV, spacing=4.0); verify(f"{cn}_{nm}", RealModel(md), B, D)
for nm, src in [("conv0.7",conv),("dipole0.8",dip)]:
    md, B, D, tgt = build_1d(600, lam, NA, src, [(0,50),(100,150),(200,250),(430,480)], TAU_DUV, spacing=4.0, focus=(0.0,80.0,-80.0), n_img=NIMG_DUV); verify(f"isodense_{nm}_PW80", RealModel(md), B, D)
md, B, D, tgt = build_1d(500, lam, NA, conv, [(0,40),(100,180),(240,400)], TAU_DUV, spacing=4.0, focus=(0.0,80.0,-80.0), n_img=NIMG_DUV); verify("mixed_conv0.7_PW80", RealModel(md), B, D)
for nm, src in [("conv0.7",conv),("dipole0.8",dip),("coh",coh)]:
    md, B, D, tgt = build_1d(360, lam, NA, src, [(0,45),(90,135),(180,225),(270,315)], 4.5, spacing=4.0); verify(f"dense90_{nm}", RealModel(md), B, D)
e2 = json.load(open('/home/user/Fable-Research/results/exp2_random_tightness.json'))
srcs = {"coh": coh, "conv0.5": SRC_DUV['conv0.5'], "conv0.9": SRC_DUV['conv0.9'], "dipole0.7": SRC_DUV['dipole0.7']}
for rec0 in e2:
    md, B, D, tgt = build_1d(rec0['L'], lam, NA, srcs[rec0['src']], rec0['segs'], TAU_DUV, spacing=5.0); verify(f"random_trial{rec0['trial']}", RealModel(md), B, D)
lamE, NAE = 13.5, 0.33
for pitch in [20, 22, 24, 26, 28, 32, 40]:
    for sname, src in {"conv0.9": SRC_EUV['conv0.9'], "dipole0.8": SRC_EUV['dipole0.8']}.items():
        segsE = [(i*pitch, i*pitch + pitch/2) for i in range(4)]
        md, B, D, tgt = build_1d(4*pitch, lamE, NAE, src, segsE, 0.05*pitch, spacing=pitch/8, focus=(0.0, 50.0, -50.0)); verify(f"EUV_p{pitch}_{sname}", RealModel(md), B, D)
ann = SRC_DUV['annular']; qua = SRC_DUV['quasar']; cv = SRC_DUV['conv0.5']
holes = lambda P, hw: rect_target((P, P), [(P/2 - hw, P/2 - hw, 2*hw, 2*hw)])
Tj = rect_target((320, 320), [(130, 60, 60, 160), (0, 10, 320, 60)]); lineend = rect_target((320, 320), [(75, 30, 50, 180), (195, 270, 50, 200)])
for name, L, tf, src, margin in [("contacts160_annular", (160,160), holes(160,35), ann, 7), ("contacts160_conv0.5", (160,160), holes(160,35), cv, 7), ("contacts120_annular", (120,120), holes(120,25), ann, 5), ("contacts120_quasar", (120,120), holes(120,25), qua, 5), ("contacts100_quasar", (100,100), holes(100,22), qua, 4.4), ("Tjunction320_annular", (320,320), Tj, ann, 6), ("Tjunction320_quasar", (320,320), Tj, qua, 6), ("lineend320_annular", (320,320), lineend, ann, 6)]:
    md, B, D, tgt = build_2d(L, lam, NA, src, tf, margin, spacing=16.0, fine=2.0); verify(name, RealModel(md), B, D)
