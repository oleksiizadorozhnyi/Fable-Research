# bnb_sdp soundness on a tiny box around a known feasible mask (must NOT certify T above its ratio; should certify T below for a local min)
import numpy as np, warnings, json; warnings.filterwarnings('ignore')
from certlitho import *; from pipeline import *; from bnb_sdp import node_certificate
lam, NA = 193.0, 1.35
e2 = json.load(open('/home/user/Fable-Research/results/exp2_random_tightness.json')); rec0 = [r for r in e2 if r['trial'] == 82][0]
md, B, D, tgt = build_1d(rec0['L'], lam, NA, SRC_DUV['coh'], rec0['segs'], 10, spacing=5.0); rm = RealModel(md); F = rm.F
r = sdp_bound_real(rm, B, D)
v, rv = round_from_Y(rm, r['X'], B, D, nsamp=200); v, rv = local_refine_real(rm, rv, B, D); vp, rp = exact_minimax_polish(rm, rv, B, D)
rb = rp if vp < v else rv; U = rm.ratio(rb, B, D); i = int(np.argmax(np.abs(rb))); rb = rb/rb[i]
Dset = [(z, d) for z in rm.focus for d in D]; Bset = [(z, b) for z in rm.focus for b in B]
rng = np.random.default_rng(0)
for hw in [1e-3, 1e-2, 5e-2]:
    l = rb - hw; u = rb + hw; l[i] = u[i] = 1.0; l = np.maximum(l, -1); u = np.minimum(u, 1)
    # exact-ish min of rho over the box by sampling
    rmin = min(rm.ratio(rng.uniform(l, u), B, D) for _ in range(3000)); rmin = min(rmin, U)
    for T in [1.01*U, 0.995*rmin]:
        A = [rm.Q(d, z) - T*rm.Q(b, zb) for (z, d) in Dset for (zb, b) in Bset]
        LB, ok, sp = node_certificate(F, i, l, u, A)
        flag = "UNSOUND!" if (ok and LB > 0 and T > rmin) else ""
        print(f"box hw={hw}: min rho in box~{rmin:.6f}; T={T:.6f} ({'>' if T>rmin else '<'} min rho): LB={LB:.3e} verified={ok} s*={sp:.3e} {flag}")
