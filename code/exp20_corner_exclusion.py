# Sensitivity of the 2D certified bounds to a corner-exclusion zone: offset samples whose foot point lies within c of a
# boundary corner are dropped (c = 0: full contour incl. corner arcs, the convention of the main tables; c = 0.25 CD and
# 0.15 CD: edge-placement checks away from corners, as in OPC practice).
from pipeline import *
lam, NA = 193.0, 1.35; out = []
holes = lambda P, hw: rect_target((P, P), [(P/2 - hw, P/2 - hw, 2*hw, 2*hw)])
Tj = rect_target((320, 320), [(130, 60, 60, 160), (0, 10, 320, 60)]); lineend = rect_target((320, 320), [(75, 30, 50, 180), (195, 270, 50, 200)])
cases = [("contacts160_annular", (160,160), holes(160,35), SRC_DUV['annular'], 7, 70), ("contacts120_quasar", (120,120), holes(120,25), SRC_DUV['quasar'], 5, 50),
         ("contacts100_quasar", (100,100), holes(100,22), SRC_DUV['quasar'], 4.4, 44), ("Tjunction320_annular", (320,320), Tj, SRC_DUV['annular'], 6, 60), ("lineend320_annular", (320,320), lineend, SRC_DUV['annular'], 6, 50)]
for name, L, tf, src, tau, cd in cases:
    for frac in [0.0, 0.15, 0.25]:
        c = frac * cd; md, B, D, tgt = build_2d(L, lam, NA, src, tf, tau, spacing=16.0, fine=2.0, corner_exclusion=c); rm = RealModel(md)
        t0 = time.time(); r = sdp_bound_real(rm, B, D); rig, _, _ = rigorize_certificate(rm, r['mu'], r['lam'])
        v, rb = strong_primal(rm, B, D, [r['X']], nround=200, nstart=8) if md.F <= 20 else (float('nan'), None)
        rec = dict(name=name, tau=tau, corner_exclusion=c, frac_cd=frac, F=md.F, nB=len(B), nD=len(D), rigorous=rig, EL_full_bound=el_full(rig), best_single=v, time=time.time()-t0)
        print(json.dumps(rec), flush=True); out.append(rec)
        json.dump(out, open('/home/user/Fable-Research/results/exp20_corner_exclusion.json', 'w'), indent=1)
