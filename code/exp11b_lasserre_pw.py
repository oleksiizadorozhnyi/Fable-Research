from pipeline import *
lam, NA = 193.0, 1.35; out = []
conv = SRC_DUV['conv0.7']
for name, L, segs, src, focus in [("mixed_conv0.7_PW80_edgeonly", 500, [(0,40),(100,180),(240,400)], conv, (0.0, 80.0, -80.0))]:
    md, B, D, tgt = build_1d(L, lam, NA, src, segs, 8, spacing=1e9, focus=focus); rm = RealModel(md)   # edge-offset samples only (interior inactive)
    r = sdp_bound_real(rm, B, D); rig, _, _ = rigorize_certificate(rm, r['mu'], r['lam'])
    v, rb = strong_primal(rm, B, D, [r['X']], nround=400, nstart=40)
    print(json.dumps(dict(stage='shor', nB=len(B), nD=len(D), shor=rig, best_single=v)), flush=True)
    t0 = time.time(); las = lasserre2_single_mask(rm, B, D, tmax=1.0, solvers=("CLARABEL","SCS")); dt = time.time()-t0
    rec = dict(name=name, F=rm.F, nB=len(B), nD=len(D), shor=rig, best_single=v, lasserre2=las['t'], status=las['status'], time=dt, M2=las['size_M2'])
    print(json.dumps(rec), flush=True); out.append(rec)
json.dump(out, open('/home/user/Fable-Research/results/exp11_lasserre_pw.json','w'), indent=1)
