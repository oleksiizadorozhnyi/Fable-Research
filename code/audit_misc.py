# Audit (6),(9),(3): metric consistency, exact sinc convention, rank_eps with dead orders, ModelX defocus w/ n_img
import numpy as np, warnings, json; warnings.filterwarnings('ignore')
from certlitho import *
lam, NA = 193.0, 1.35; rng = np.random.default_rng(3)
conv = conventional_source_1d(0.7, NA, 7)
# (9) exact analytic check of G in sdp_bound_pixel_int
md, B, D, _ = build_1d(600, lam, NA, conv, [(0,50),(100,150),(200,250),(430,480)], 8, spacing=4.0)
N = 30; n = np.arange(N); k = md.k[:,0]; m = rng.standard_normal(N)
G = (np.sinc(k/N)[:, None]) * np.exp(-2j*np.pi*np.outer(k, n)/N)/N
Mexact_edges = (np.sinc(k/N)[:, None]*np.exp(-2j*np.pi*np.outer(k, n+0.5)/N)/N) @ m   # pixel n = [nL/N,(n+1)L/N)
print(f"pixel_int G: err vs exact(edges at nL/N) {np.abs(G@m-Mexact_edges).max():.2e}; err vs exact(centers at nL/N) {np.abs(G@m-Mexact_edges*np.exp(1j*np.pi*k/N)).max():.2e}")
# does the alignment matter for the bound? compare pixel_int bound with/without half-pixel phase for N=15 (coarse)
import cvxpy as cp
def pix_bound(G):
    X = cp.Variable((N, N), symmetric=True); t = cp.Variable(); cons = [X >> 0]
    for ix in D: Qp = np.real(G.conj().T @ md.Q(ix) @ G); cons.append(cp.trace(0.5*(Qp+Qp.T) @ X) <= t)
    for ix in B: Qp = np.real(G.conj().T @ md.Q(ix) @ G); cons.append(cp.trace(0.5*(Qp+Qp.T) @ X) >= 1)
    cp.Problem(cp.Minimize(t), cons).solve(solver="CLARABEL"); return float(t.value)
for N in [15, 12, 9]:
    n = np.arange(N); G0 = (np.sinc(k/N)[:, None]) * np.exp(-2j*np.pi*np.outer(k, n)/N)/N; G1 = G0*np.exp(-1j*np.pi*k/N)[:, None]
    print(f"  N={N}: bound(centers at nL/N, as coded)={pix_bound(G0):.5f}  bound(edges at nL/N)={pix_bound(G1):.5f}")
# (6) metric consistency: local_refine_real objective at beta->inf equals log rm.ratio; polish feasibility
md, B, D, _ = build_1d(500, lam, NA, conv, [(0,40),(100,180),(240,400)], 8, spacing=8.0, focus=(0.0, 80.0, -80.0)); rm = RealModel(md)
r0 = rng.standard_normal(rm.F); v, r1 = local_refine_real(rm, r0, B, D); vp, rp = exact_minimax_polish(rm, r1, B, D)
dmax = max(rm.image(rp, z)[D].max() for z in rm.focus); bmin = min(rm.image(rp, z)[B].min() for z in rm.focus)
print(f"metric: polish ratio {vp:.6f} == dmax/bmin {dmax/bmin:.6f}; per-focus ratios {[rm.image(rp,z)[D].max()/rm.image(rp,z)[B].min() for z in rm.focus]}")
# (3) rank_eps and dead orders on Tjunction320
Tj = lambda a, b: (abs(a-160) < 30 and 60 < b < 220) or (abs(b-40) < 30)
Y = np.load('/home/user/Fable-Research/results/exp4_Tjunction320_annular_Y.npy')
md4, B4, D4, _ = build_2d((320,320), lam, NA, annular_source_2d(0.6,0.9,NA,9), Tj, 14, spacing=16.0, fine=2.0); rm4 = RealModel(md4)
alive = (np.abs(md4.C[0.0][:, 0, :]) > 0).any(0); dead_r = np.abs(rm4.U[~alive, :]).sum(0) > 1e-12   # real coordinates touching dead orders
ev = np.linalg.eigvalsh(Y)[::-1]; Yd = Y[np.ix_(dead_r, dead_r)]
print(f"Tj320: F={rm4.F} dead real coords={int(dead_r.sum())}; saved Y eigs top6={np.round(ev[:6],5).tolist()}; rank_eps(1e-6)={int(np.sum(ev>1e-6*ev[0]))}; ||Y_dead||={np.linalg.norm(Yd):.3e}, trace(Y)={np.trace(Y):.4f}")
Yk = Y[np.ix_(~dead_r, ~dead_r)]; evk = np.linalg.eigvalsh(Yk)[::-1]; print(f"   rank_eps of Y restricted to live coords: {int(np.sum(evk>1e-6*evk[0]))}, eigs {np.round(evk[:6],5).tolist()}")
# (8) ModelX defocus uses gl=|g|lam clipped, ignoring n_img: show phase for vector model n_img=1.7 vs correct 2 pi n/lam z (sqrt(1-(g lam/n)^2)-1)
segs = [(0,45),(90,135),(180,225),(270,315)]
mdv, *_ = build_1d_x(360, lam, NA, conv, segs, 6, spacing=8.0, vector=True, pol='x', n_img=1.7, focus=(0.0, 80.0))
print("ModelX: defocus phase formula uses sqrt(1-min(|g|lam,1)^2) regardless of n_img (see code line 'gl = np.clip(gn*lamv, 0, 1)')")
# TE/TM for s=0 in 1D conventional source: 'rad' defaults to x-hat
print("ModelX TE/TM: for the on-axis source point s=0, rad := x-hat (arbitrary), azi := y-hat")
