from pipeline import *
lam, NA = 193.0, 1.35; out = []
cases = [("isodense_coh", 600, [(0,50),(100,150),(200,250),(430,480)], SRC_DUV['coh'], 8, 96),
         ("mixed_coh", 500, [(0,40),(100,180),(240,400)], SRC_DUV['coh'], 8, 96),
         ("trial7_dipole0.7", 600, [(4.417361850603307, 83.85149016558589), (119.944921485065, 172.99686235814647), (232.48702741810467, 409.4973023849853)], SRC_DUV['dipole0.7'], 10, 80)]
for name, L, segs, (s, w), margin, nx in cases:
    md = Model(L, lam, NA, s, w, nx=nx); tf = lambda v: any(a <= (v % L) < b for a, b in segs)
    B, D, tgt = sets_1d(md, tf, margin); rm = RealModel(md)
    r = sdp_bound_real(rm, B, D); rig, _, _ = rigorize_certificate(rm, r['mu'], r['lam'])
    best, ranks = face_min_rank(rm, B, D, r['t'], ntries=8)
    v, rbest = strong_primal(rm, B, D, [r['X'], best[1]], nround=400, nstart=40)
    rec = dict(name=name, F=rm.F, t_sdp=r['t'], rigorous=rig, rank_ipm=r['rank_eps'], face_ranks=ranks, min_face_rank=best[0], face_eigs=[float(e) for e in best[2][:4]],
               best_feasible=v, gap_rel=(v-rig)/rig)
    print(json.dumps(rec), flush=True); out.append(rec)
json.dump(out, open('/home/user/Fable-Research/results/exp6_gap_check.json','w'), indent=1)
