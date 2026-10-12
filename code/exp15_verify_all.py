# Independent high-precision verification of every reported certificate: recompute SDP, rigorize, then verify
# S_c = SD - c SB (restricted to complement of null space of SD) is PSD via mpmath Cholesky at 50 digits; also check
# both solvers agree. Covers exp1 (1D features), exp2 (random), exp5 EUV, exp4 2D contacts.
from pipeline import *
import mpmath as mp
mp.mp.dps = 50
lam, NA = 193.0, 1.35; out = []
def verify(name, rm, B, D):
    r = sdp_bound_real(rm, B, D, solvers=("CLARABEL",)); rig, me, c = rigorize_certificate(rm, r['mu'], r['lam'])
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
conv = conventional_source_1d(0.7, NA, 7); dip = dipole_source_1d(0.8, NA); coh = ([0.0],[1.0])
for nm, src in [("coh",coh),("conv0.7",conv),("dipole0.8",dip)]:
    for cn, L, segs in [("isodense", 600, [(0,50),(100,150),(200,250),(430,480)]), ("mixed", 500, [(0,40),(100,180),(240,400)])]:
        md, B, D, tgt = build_1d(L, lam, NA, src, segs, 8, spacing=4.0); verify(f"{cn}_{nm}", RealModel(md), B, D)
for nm, src in [("conv0.7",conv),("dipole0.8",dip)]:
    md, B, D, tgt = build_1d(600, lam, NA, src, [(0,50),(100,150),(200,250),(430,480)], 8, spacing=4.0, focus=(0.0,80.0,-80.0)); verify(f"isodense_{nm}_PW80", RealModel(md), B, D)
e2 = json.load(open('/home/user/Fable-Research/results/exp2_random_tightness.json'))
srcs = {"coh": coh, "conv0.5": conventional_source_1d(0.5, NA, 5), "conv0.9": conventional_source_1d(0.9, NA, 7), "dipole0.7": dipole_source_1d(0.7, NA)}
for rec0 in e2:
    md, B, D, tgt = build_1d(rec0['L'], lam, NA, srcs[rec0['src']], rec0['segs'], 10, spacing=5.0); verify(f"random_trial{rec0['trial']}", RealModel(md), B, D)
lamE, NAE = 13.5, 0.33
for pitch in [20, 22, 24, 26, 28, 32, 40]:
    for sname, src in {"conv0.9": conventional_source_1d(0.9, NAE, 7), "dipole0.8": dipole_source_1d(0.8, NAE)}.items():
        segsE = [(i*pitch, i*pitch + pitch/2) for i in range(4)]
        md, B, D, tgt = build_1d(4*pitch, lamE, NAE, src, segsE, 0.12*pitch, spacing=pitch/8, focus=(0.0, 50.0, -50.0)); verify(f"EUV_p{pitch}_{sname}", RealModel(md), B, D)
ann = annular_source_2d(0.6, 0.9, NA, n=9); qua = quasar_source_2d(0.6, 0.9, NA, n=11); cv = disk_source_2d(0.5, NA, n=5)
holes = lambda P, hw: (lambda a, b: (abs((a+P/2)%P-P/2) < hw) and (abs((b+P/2)%P-P/2) < hw))
for name, L, tf, src, margin in [("contacts160_annular", (160,160), holes(160,35), ann, 12), ("contacts160_conv0.5", (160,160), holes(160,35), cv, 12), ("contacts120_annular", (120,120), holes(120,25), ann, 10), ("contacts120_quasar", (120,120), holes(120,25), qua, 10), ("contacts100_quasar", (100,100), holes(100,22), qua, 8)]:
    md, B, D, tgt = build_2d(L, lam, NA, src, tf, margin, spacing=16.0, fine=2.0); verify(name, RealModel(md), B, D)
