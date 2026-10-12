# Rigorous B&B certification of the single-mask optimum on small-F gap instances -> certified multi-exposure gains.
from pipeline import *
from bnb_cert import bnb_single_mask
lam, NA = 193.0, 1.35; out = []
e2 = json.load(open('/home/user/Fable-Research/results/exp2_random_tightness.json'))
cases = []
for tr in [82, 91, 110, 50]:
    rec0 = [r for r in e2 if r['trial'] == tr][0]; cases.append((f"trial{tr}", rec0['L'], rec0['segs'], SRC_DUV['coh'], 10, 5.0, (0.0,)))
cases.append(("mixed_coh", 500, [(0,40),(100,180),(240,400)], SRC_DUV['coh'], TAU_DUV, 4.0, (0.0,)))
cases.append(("isodense_coh", 600, [(0,50),(100,150),(200,250),(430,480)], SRC_DUV['coh'], TAU_DUV, 4.0, (0.0,)))
cases.append(("mixed_conv0.7_PW80_edgeonly", 500, [(0,40),(100,180),(240,400)], SRC_DUV['conv0.7'], TAU_DUV, 1e9, (0.0, 80.0, -80.0)))
for name, L, segs, src, margin, spacing, focus in cases:
    md, B, D, tgt = build_1d(L, lam, NA, src, segs, margin, spacing=spacing, focus=focus); rm = RealModel(md)
    r = sdp_bound_real(rm, B, D); rig, _, _ = rigorize_certificate(rm, r['mu'], r['lam'])
    v, rb = strong_primal(rm, B, D, [r['X']], nround=300, nstart=20)
    res = bnb_single_mask(rm, B, D, tol=1e-3, max_nodes=6_000_000, U0=v)
    rec = dict(name=name, F=rm.F, nB=len(B), nD=len(D), sdp_bound=rig, best_single=v, bnb_L=res['L'], bnb_U=res['U'], bnb_nodes=res['nodes'], bnb_converged=res['converged'], bnb_time=res['time'],
               certified_gain=(res['L']-rig)/rig)
    print(json.dumps(rec), flush=True); out.append(rec); json.dump(out, open('/home/user/Fable-Research/results/exp14_bnb.json','w'), indent=1)
