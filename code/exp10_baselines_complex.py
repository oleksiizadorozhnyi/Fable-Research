from pipeline import *
lam, NA = 193.0, 1.35; out = []
srcs = {"coh": SRC_DUV['coh'], "conv0.7": SRC_DUV['conv0.7'], "dipole0.8": SRC_DUV['dipole0.8']}
cases = {"isodense": (600, [(0,50),(100,150),(200,250),(430,480)]), "mixed": (500, [(0,40),(100,180),(240,400)]), "dense90": (360, [(0,45),(90,135),(180,225),(270,315)])}
for cname, (L, segs) in cases.items():
    for sname, src in srcs.items():
        md, B, D, tgt = build_1d(L, lam, NA, src, segs, TAU_DUV if cname != 'dense90' else 4.5, spacing=4.0); rm = RealModel(md)
        r = sdp_bound_real(rm, B, D); rig, _, _ = rigorize_certificate(rm, r['mu'], r['lam'])
        # (a) Rayleigh-quotient (generalized eigenvalue) bound: min_r (avg_D I)/(avg_B I) = lambda_min(QD, QB)  [IBM 'unconstrained' simplification]
        from scipy.linalg import eigh
        QD = sum(rm.Q(d) for d in D)/len(D); QB = sum(rm.Q(b) for b in B)/len(B)
        ev = eigh(QD, QB + 1e-12*np.eye(rm.F), eigvals_only=True); rq = float(ev.min())
        # (b) complex-transmission (arbitrary phase) mask relaxation (Hermitian X): should be <= real bound
        rc = sdp_bound(md, B, D, real_mask=False, solver="CLARABEL")
        rec = dict(case=cname, src=sname, F=rm.F, sdp_real=rig, rayleigh_bound=rq, sdp_complex=rc['t'], rank_complex=rc['rank_eps'])
        print(json.dumps(rec), flush=True); out.append(rec)
json.dump(out, open('/home/user/Fable-Research/results/exp10_baselines_complex.json','w'), indent=1)
