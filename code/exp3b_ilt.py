from pipeline import *
lam, NA = 193.0, 1.35; out = []
conv = conventional_source_1d(0.7, NA, 7)
cases = {"isodense": (600, [(0,50),(100,150),(200,250),(430,480)]), "mixed": (500, [(0,40),(100,180),(240,400)]), "dense120": (360, [(0,60),(120,180),(240,300)])}
for name, (L, segs) in cases.items():
    md, B, D, tgt = build_1d(L, lam, NA, conv, segs, margin=10, spacing=5.0); rm = RealModel(md)
    tf = lambda v: any(a <= (v % L) < b for a, b in segs)
    rF = sdp_bound_real(rm, B, D); rig, _, _ = rigorize_certificate(rm, rF['mu'], rF['lam'])
    N = 120; xp = np.arange(N)*L/N; tgt_init = np.array([tf(v) for v in xp]).astype(float)
    res = {}
    for kind, pb in [("psm", (-1, 1)), ("binary", (0, 1))]:
        best = (np.inf, None)
        for a, lr in [(10.0, 0.02), (4.0, 0.05), (25.0, 0.01)]:
            for init in ["target", "random"]:
                r_, m_, M_, h_ = ilt_gradient(md, tgt.astype(float), N, B, D, steps=1500, lr=lr, a=a, pixel_bounds=pb, init_target=(tgt_init if init == "target" else None), rng=np.random.default_rng(0))
                res[f"ilt_{kind}_a{a}_lr{lr}_{init}"] = r_
                if r_ < best[0]: best = (r_, M_)
        res[f"ilt_{kind}_best"] = best[0]
        v, _ = local_refine_real(rm, rm.from_orders(best[1]), B, D); res[f"ilt_{kind}_best_refined_on_ratio"] = v
    rec = dict(name=name, F=rm.F, certified_bound=rig, **res); print(json.dumps(rec), flush=True); out.append(rec)
json.dump(out, open('/home/user/Fable-Research/results/exp3b_ilt.json','w'), indent=1)
