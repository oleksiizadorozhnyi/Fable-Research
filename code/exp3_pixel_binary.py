from pipeline import *
lam, NA = 193.0, 1.35; out = []
conv = conventional_source_1d(0.7, NA, 7)
cases = {"isodense": (600, [(0,50),(100,150),(200,250),(430,480)]), "mixed": (500, [(0,40),(100,180),(240,400)]), "dense120": (360, [(0,60),(120,180),(240,300)])}
for name, (L, segs) in cases.items():
    s, w = conv; md, B, D, tgt = build_1d(L, lam, NA, conv, segs, margin=10, spacing=5.0); rm = RealModel(md)
    tf = lambda v: any(a <= (v % L) < b for a, b in segs)
    rF = sdp_bound_real(rm, B, D); N = 120
    rP = sdp_bound_pixel(md, N, B, D, nonneg=False, solver="SCS"); rN = sdp_bound_pixel(md, N, B, D, nonneg=True, solver="SCS")
    tgt_pix = tgt.astype(float)   # target on the model sample grid (where the resist image is evaluated)
    rP_ilt, mP, MP, hP = ilt_gradient(md, tgt_pix, N, B, D, steps=800, lr=0.05, pixel_bounds=(-1,1))
    rB_ilt, mB, MB, hB = ilt_gradient(md, tgt_pix, N, B, D, steps=800, lr=0.05, pixel_bounds=(0,1))
    rP_ref, _ = local_refine_real(rm, rm.from_orders(MP), B, D)
    # best feasible PSM via pipeline
    rec0, _, _, _ = analyze(md, B, D, nround=200, nstart=8)
    # binary (nonneg) feasible: local search in pixel domain directly on ratio with m in [0,1] (projected gradient on log-ratio)
    rec = dict(name=name, F=md.F, nB=len(B), nD=len(D), t_fourier_psm=rF['t'], rank_fourier=rF['rank_eps'], t_pixel_psm=rP['t'], rank_pixel=rP['rank_eps'],
               t_pixel_nonneg_DNN=rN['t'], rank_dnn=rN['rank_eps'], ilt_psm_ratio=rP_ilt, ilt_psm_refined_on_ratio=rP_ref, ilt_binary_ratio=rB_ilt,
               best_psm_pipeline=rec0['best_feasible'], gap_psm=rec0['gap_rel'])
    print(json.dumps(rec), flush=True); out.append(rec)
    np.save(f'/home/user/Fable-Research/results/exp3_{name}_ilt_psm.npy', mP); np.save(f'/home/user/Fable-Research/results/exp3_{name}_ilt_bin.npy', mB)
    json.dump(out, open('/home/user/Fable-Research/results/exp3_pixel_binary.json', 'w'), indent=1)
