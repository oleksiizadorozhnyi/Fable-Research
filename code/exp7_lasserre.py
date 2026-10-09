from pipeline import *
lam, NA = 193.0, 1.35; out = []
e1 = json.load(open('/home/user/Fable-Research/results/exp1_1d_features.json'))
cands = [r for r in e1 if r['gap_rel'] > 2e-3 and r['F'] <= 11 and len(r['focus']) == 1]
srcs = {"coh": ([0.0],[1.0]), "conv0.7": conventional_source_1d(0.7, NA, 7), "dipole0.8": dipole_source_1d(0.8, NA)}
for rec0 in cands:
    sname = rec0['name'].split('_')[-1]; src = srcs[sname]
    md, B, D, tgt = build_1d(rec0['L'], lam, NA, src, rec0['segs'], 8 if 'dense90' not in rec0['name'] else 6, spacing=4.0); rm = RealModel(md)
    r = sdp_bound_real(rm, B, D); rig, _, _ = rigorize_certificate(rm, r['mu'], r['lam'])
    best, ranks = face_min_rank(rm, B, D, r['t'], ntries=6)
    v, rbest = strong_primal(rm, B, D, [r['X'], best[1]], nround=400, nstart=40)
    t0 = time.time(); las = lasserre2_single_mask(rm, B, D, tmax=1.0); dt = time.time()-t0
    rec = dict(name=rec0['name'], F=rm.F, shor_bound=r['t'], rigorous_shor=rig, min_face_rank=best[0], lasserre2_bound=las['t'], lasserre_status=las['status'], lasserre_time=dt,
               M2_size=las['size_M2'], best_single_mask=v, multi_exposure_gain_certified=(las['t'] - rig)/rig if las['t'] else None)
    print(json.dumps(rec), flush=True); out.append(rec)
json.dump(out, open('/home/user/Fable-Research/results/exp7_lasserre.json','w'), indent=1)
