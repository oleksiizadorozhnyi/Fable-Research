from pipeline import *
lam, NA = 193.0, 1.35; out = {}
e2 = json.load(open('/home/user/Fable-Research/results/exp2_random_tightness.json'))
srcs = {"coh": ([0.0],[1.0]), "conv0.5": conventional_source_1d(0.5, NA, 5), "conv0.9": conventional_source_1d(0.9, NA, 7), "dipole0.7": dipole_source_1d(0.7, NA)}
# (a) outlier trial 16 and the larger-gap trials: strong primal + face exploration
for tr in [16, 82, 91, 110, 50, 46]:
    rec0 = [r for r in e2 if r['trial'] == tr][0]; L = rec0['L']; segs = rec0['segs']; s, w = srcs[rec0['src']]
    md, B, D, tgt = build_1d(L, lam, NA, (s, w), segs, 10, spacing=5.0); rm = RealModel(md)
    r = sdp_bound_real(rm, B, D); rig, _, _ = rigorize_certificate(rm, r['mu'], r['lam'])
    best, ranks = face_min_rank(rm, B, D, r['t'], ntries=6)
    v, rb = strong_primal(rm, B, D, [r['X'], best[1]], nround=600, nstart=60)
    las = lasserre2_single_mask(rm, B, D, tmax=1.0) if rm.F <= 11 else dict(t=None, status='skipped(F>11)')
    rec = dict(trial=tr, src=rec0['src'], F=rm.F, nB=len(B), nD=len(D), segs=segs, bound=rig, old_best=rec0['best_feasible'], new_best=v, gap_new=(v-rig)/rig,
               min_face_rank=best[0], face_eigs=[float(e) for e in best[2][:3]], lasserre2=las['t'], lasserre_status=las['status'])
    print(json.dumps(rec), flush=True); out[f'trial{tr}'] = rec
# (b) process-window cases: strong primal, face rank, Lasserre-2 where feasible
conv = conventional_source_1d(0.7, NA, 7); dip = dipole_source_1d(0.8, NA)
for name, L, segs, (s, w) in [("mixed_conv0.7_PW80", 500, [(0,40),(100,180),(240,400)], conv), ("isodense_conv0.7_PW80", 600, [(0,50),(100,150),(200,250),(430,480)], conv), ("isodense_dipole0.8_PW80", 600, [(0,50),(100,150),(200,250),(430,480)], dip)]:
    md, B, D, tgt = build_1d(L, lam, NA, (s, w), segs, 8, spacing=4.0, focus=(0.0, 80.0, -80.0)); rm = RealModel(md)
    r = sdp_bound_real(rm, B, D); rig, _, _ = rigorize_certificate(rm, r['mu'], r['lam'])
    best, ranks = face_min_rank(rm, B, D, r['t'], ntries=6)
    v, rb = strong_primal(rm, B, D, [r['X'], best[1]], nround=600, nstart=60)
    las = lasserre2_single_mask(rm, B, D, tmax=1.0, solvers=("SCS",)) if rm.F <= 11 else dict(t=None, status='skipped(F>11)')
    rec = dict(name=name, F=rm.F, bound=rig, new_best=v, gap_new=(v-rig)/rig, min_face_rank=best[0], face_eigs=[float(e) for e in best[2][:4]], lasserre2=las['t'], lasserre_status=las['status'])
    print(json.dumps(rec), flush=True); out[name] = rec
    np.save(f'/home/user/Fable-Research/results/exp9_{name}_rbest.npy', rb); np.save(f'/home/user/Fable-Research/results/exp9_{name}_Yface.npy', best[1])
json.dump(out, open('/home/user/Fable-Research/results/exp9_outlier_pw.json','w'), indent=1)
