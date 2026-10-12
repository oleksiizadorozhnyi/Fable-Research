from pipeline import *
lam, NA = 193.0, 1.35; out = {}
conv = SRC_DUV['conv0.7']; s, w = conv
L = 600; segs = [(0,50),(100,150),(200,250),(430,480)]; tf = lambda v: any(a <= (v % L) < b for a, b in segs)
# (a) sampling density: bound must be nondecreasing in nx (more constraints) and converge
res = []
for nx in [16, 8, 4, 2, 1]:
    md, B, D, tgt = build_1d(L, lam, NA, (s, w), segs, TAU_DUV, spacing=nx); rm = RealModel(md)
    r = sdp_bound_real(rm, B, D); rig, _, _ = rigorize_certificate(rm, r['mu'], r['lam'])
    res.append(dict(nx=nx, nB=len(B), nD=len(D), bound=rig, t=r['time'])); print(json.dumps(res[-1]), flush=True)
out['sampling'] = res
# (b) margin (edge tolerance) sweep at nx=192
res = []
for margin in [2, 2.5, 4, 5, 8, 12, 16, 20]:
    md, B, D, tgt = build_1d(L, lam, NA, (s, w), segs, margin, spacing=4.0); rm = RealModel(md)
    r = sdp_bound_real(rm, B, D); rig, _, _ = rigorize_certificate(rm, r['mu'], r['lam'])
    res.append(dict(margin=margin, nB=len(B), nD=len(D), bound=rig, EL=(1-rig)/(1+rig))); print(json.dumps(res[-1]), flush=True)
out['margin'] = res
# (c) solver cross-check (CLARABEL vs SCS) on 5 instances; (d) number of source points (quadrature) convergence
res = []
for nsrc in [7, 11, 15, 21, 31, 41]:
    s2, w2 = disk_source(0.7, NA, n=nsrc); md, B, D, tgt = build_1d(L, lam, NA, (s2, w2), segs, TAU_DUV, spacing=4.0); rm = RealModel(md)
    rc = sdp_bound_real(rm, B, D, solvers=("CLARABEL",)); rs = sdp_bound_real(rm, B, D, solvers=("SCS",))
    res.append(dict(nsrc=len(w2), grid=nsrc, F=rm.F, clarabel=rc['t'], scs=rs['t'], diff=abs(rc['t']-rs['t']))); print(json.dumps(res[-1]), flush=True)
out['source_quadrature_and_solvers'] = res
# (e) independent high-precision certificate check with mpmath Cholesky on one instance
import mpmath as mp
md, B, D, tgt = build_1d(L, lam, NA, (s, w), segs, 8, spacing=4.0); rm = RealModel(md)
r = sdp_bound_real(rm, B, D); rig, mineig, c = rigorize_certificate(rm, r['mu'], r['lam'])
F = rm.F; SD = np.zeros((F,F)); smu = sum(r['mu'].values()); 
for (z, ix), v in r['mu'].items(): SD += v/smu*rm.Q(ix, z)
SB = np.zeros((F,F)); slam = 0
for (z, ix), v in r['lam'].items(): SB += v/smu*rm.Q(ix, z); slam += v/smu
S = SD - c*SB
mp.mp.dps = 50
Smp = mp.matrix(S.tolist())
try:
    Lc = mp.cholesky(Smp); chol_ok = True
except Exception as e:
    chol_ok = False
# also exact-rational check: min eigenvalue by mpmath
ev = mp.eigsy(Smp)[0]; mineig_mp = min(ev)
out['mp_check'] = dict(bound=c*slam, cholesky_success=chol_ok, mineig_mp50=float(mineig_mp), norm=float(np.linalg.norm(S, 2)))
print(json.dumps(out['mp_check']))
json.dump(out, open('/home/user/Fable-Research/results/exp8_robustness.json','w'), indent=1)
