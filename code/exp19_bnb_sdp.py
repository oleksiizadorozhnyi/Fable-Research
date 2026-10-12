from pipeline import *
from bnb_sdp import certify_target
lam, NA = 193.0, 1.35; out = []
cases = [("isodense_coh", 600, [(0,50),(100,150),(200,250),(430,480)], SRC_DUV['coh'], TAU_DUV, 4.0, (0.0,), [0.40, 0.44]),
         ("mixed_conv0.7_PW80_edgeonly", 500, [(0,40),(100,180),(240,400)], SRC_DUV['conv0.7'], TAU_DUV, 1e9, (0.0, 80.0, -80.0), [0.90, 0.93]),
         ("mixed_coh", 500, [(0,40),(100,180),(240,400)], SRC_DUV['coh'], TAU_DUV, 4.0, (0.0,), [0.46])]
for name, L, segs, src, margin, spacing, focus, Ts in cases:
    md, B, D, tgt = build_1d(L, lam, NA, src, segs, margin, spacing=spacing, focus=focus, n_img=NIMG_DUV if len(focus) > 1 else 1.0); rm = RealModel(md)
    r = sdp_bound_real(rm, B, D); rig, _, _ = rigorize_certificate(rm, r['mu'], r['lam'])
    for T in Ts:
        res = certify_target(rm, B, D, T, mu=r['mu'], lam=r['lam'], max_nodes=6000)
        rec = dict(name=name, F=rm.F, sdp_bound=rig, certified_gain_if_ok=(T-rig)/rig, **{k: v for k, v in res.items() if k != 'T'}, T=T); print(json.dumps(rec), flush=True); out.append(rec)
        json.dump(out, open('/home/user/Fable-Research/results/exp19_bnb_sdp.json','w'), indent=1)
        if not res['certified']: break
