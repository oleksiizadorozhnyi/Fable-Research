from pipeline import *
from bnb_cert import bnb_single_mask
lam, NA = 193.0, 1.35; out = []
for name, L, segs, src, margin, spacing, focus, budget in [("mixed_coh", 500, [(0,40),(100,180),(240,400)], SRC_DUV['coh'], TAU_DUV, 4.0, (0.0,), 30_000_000)]:
    md, B, D, tgt = build_1d(L, lam, NA, src, segs, margin, spacing=spacing, focus=focus); rm = RealModel(md)
    r = sdp_bound_real(rm, B, D); rig, _, _ = rigorize_certificate(rm, r['mu'], r['lam']); t_up = ensemble_ratio(rm, r['X'], B, D)
    v, rb = strong_primal(rm, B, D, [r['X']], nround=300, nstart=20)
    res = bnb_single_mask(rm, B, D, tol=2e-3, max_nodes=budget, U0=v, batch=16384)
    rec = dict(name=name, F=rm.F, nB=len(B), nD=len(D), sdp_bound=rig, best_single=v, bnb_L=res['L'], bnb_U=res['U'], bnb_nodes=res['nodes'], bnb_converged=res['converged'], bnb_time=res['time'], certified_gain=(res['L']-t_up)/t_up, sdp_ensemble_ratio=t_up)
    print(json.dumps(rec), flush=True); out.append(rec); json.dump(out, open('/home/user/Fable-Research/results/exp14b_bnb_big.json','w'), indent=1)
