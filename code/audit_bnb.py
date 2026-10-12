# Audit (5): bnb_cert interval bounds and final L validity on an F=5 instance (incl. multi-focus F=5)
import numpy as np, warnings, json, heapq; warnings.filterwarnings('ignore')
from certlitho import *; from bnb_cert import bnb_single_mask
lam, NA = 193.0, 1.35; rng = np.random.default_rng(2)
e2 = json.load(open('/home/user/Fable-Research/results/exp2_random_tightness.json'))
rec0 = [r for r in e2 if r['trial'] == 82][0]
cases = [("trial82_coh", build_1d(rec0['L'], lam, NA, ([0.0],[1.0]), rec0['segs'], 10, spacing=5.0)),
         ("F5_coh_test", build_1d(300, lam, NA, ([0.0],[1.0]), [(30,80),(150,230)], 10, spacing=6.0)),
         ("F5_coh_PW", build_1d(300, lam, NA, ([0.0],[1.0]), [(30,80),(150,230)], 10, spacing=6.0, focus=(0.0, 80.0, -80.0)))]
for name, (md, B, D, tgt) in cases:
    rm = RealModel(md); F = rm.F
    # (a) check the box bounds themselves: for random boxes on a face, sample random points in box, q must lie in [lo, hi]
    QD = [np.stack([rm.Q(d, z) for d in D]) for z in rm.focus]
    viol = 0
    for _ in range(300):
        i = rng.integers(F); c = rng.uniform(-1, 1, F); c[i] = 1.0; h = rng.uniform(0, 0.5, F); h[i] = 0.0
        h = np.minimum(h, np.minimum(1-c, 1+c)); h[i] = 0.0   # keep box inside face
        for z in range(len(rm.focus)):
            Q = QD[z]; qc = np.einsum('f,dfg,g->d', c, Q, c); Qc = np.einsum('dfg,g->df', Q, c)
            e1 = 2*np.abs(Qc) @ h; e2 = np.einsum('f,dfg,g->d', h, np.abs(Q), h)
            lo = np.maximum(qc - e1, 0); hi = qc + e1 + e2
            for _ in range(20):
                x = c + rng.uniform(-1, 1, F)*h; q = np.einsum('f,dfg,g->d', x, Q, x)
                viol += int(np.any(q < lo - 1e-12) or np.any(q > hi + 1e-12))
    print(f"{name}: F={F} box-bound violations {viol} / {300*20*len(rm.focus)}")
    # (b) B&B with small budgets (non-converged) and converged; L must be <= rho(r) for all random r and <= U
    r = sdp_bound_real(rm, B, D); rig, _, _ = rigorize_certificate(rm, r['mu'], r['lam'])
    v, rv = round_from_Y(rm, r['X'], B, D, nsamp=200); v, rv = local_refine_real(rm, rv, B, D)
    for mx in [10, 500, 20000, 3000000]:
        res = bnb_single_mask(rm, B, D, tol=1e-3, max_nodes=mx, U0=v)
        L = res['L']; worst = np.inf; nv = 0
        for _ in range(20000):
            x = rng.standard_normal(F) if rng.random() < 0.5 else rv + 0.05*rng.standard_normal(F)
            q = rm.ratio(x, B, D); worst = min(worst, q); nv += int(q < L)
        # also exhaustive-ish: local optimization from many starts
        bestloc = min(local_refine_real(rm, rng.standard_normal(F), B, D, iters=200)[0] for _ in range(10))
        print(f"   max_nodes={mx}: L={L:.6f} U={res['U']:.6f} conv={res['converged']} nodes={res['nodes']} active={res['active']} | random worst={worst:.6f} multistart best={bestloc:.6f} violations={nv} | sdp rig={rig:.6f}")
