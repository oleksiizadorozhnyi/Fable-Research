# Sanity/soundness test of bnb_sdp.certify_target: must certify T below the known optimum region and must NOT certify T above a known feasible mask
import numpy as np, warnings, json, time; warnings.filterwarnings('ignore')
from certlitho import *; from pipeline import *; from bnb_sdp import certify_target, node_certificate; from bnb_cert import bnb_single_mask
lam, NA = 193.0, 1.35
e2 = json.load(open('/home/user/Fable-Research/results/exp2_random_tightness.json')); rec0 = [r for r in e2 if r['trial'] == 82][0]
md, B, D, tgt = build_1d(rec0['L'], lam, NA, SRC_DUV['coh'], rec0['segs'], 10, spacing=5.0); rm = RealModel(md)
r = sdp_bound_real(rm, B, D); rig, _, _ = rigorize_certificate(rm, r['mu'], r['lam'])
v, rv = round_from_Y(rm, r['X'], B, D, nsamp=200); v, rv = local_refine_real(rm, rv, B, D); vp, rp = exact_minimax_polish(rm, rv, B, D); U = min(v, vp); rbest = rp if vp < v else rv
res = bnb_single_mask(rm, B, D, tol=1e-3, max_nodes=2000000, U0=U); print(json.dumps(dict(F=rm.F, sdp=rig, best_single=U, interval_bnb_L=res['L'], interval_bnb_U=res['U'])))
for T, expect in [(0.98*res['L'], True), (1.02*U, False), (0.9*rig, True)]:
    t0 = time.time(); out = certify_target(rm, B, D, T, mu=r['mu'], lam=r['lam'], max_nodes=400)
    out["T_target"]=round(T,5); out["expect_certifiable"]=expect; print(json.dumps(out), flush=True)
    if out['certified'] and T > U: print("  !!! UNSOUND: certified T above a feasible mask ratio", U)
# single-node check: root face certificate with T above U must fail; evaluate node bound vs exact min over random points in the face
i = 0; l = -np.ones(rm.F); u = np.ones(rm.F); l[i] = u[i] = 1.0
Dset = [(z, d) for z in rm.focus for d in D]; Bset = [(z, b) for z in rm.focus for b in B]
A = [rm.Q(d, z) - U*1.01*rm.Q(b, zb) for (z, d) in Dset for (zb, b) in Bset]
LB, ok, sp = node_certificate(rm.F, i, l, u, A); print("root node, T=1.01U: LB", LB, "verified", ok, "s*", sp, "(must not be > 0)")
