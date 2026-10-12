import sys, json, numpy as np, time; sys.path.insert(0, '/home/user/Fable-Research/code')
import warnings; warnings.filterwarnings('ignore')
from certlitho import *

def analyze(md, B, D, nround=300, nstart=10, iters=400, polish=True, save_prefix=None):
    """Full pipeline on a Model: certified bound (rigorous), rank, rounding, local refine, multistart, SLSQP polish."""
    rm = RealModel(md)
    t0 = time.time(); r = sdp_bound_real(rm, B, D); t_sdp = time.time() - t0
    rig, mineig, c = rigorize_certificate(rm, r['mu'], r['lam'])
    rr, rv = round_from_Y(rm, r['X'], B, D, nsamp=nround)
    rl, rvl = local_refine_real(rm, rv, B, D, iters=iters)
    rmv, rvm = multistart_local_real(rm, B, D, nstart=nstart, iters=iters)
    cands = [(rl, rvl), (rmv, rvm)]
    if polish:
        for v, x in list(cands):
            try:
                vp, xp = exact_minimax_polish(rm, x, B, D)
                if np.isfinite(vp): cands.append((vp, xp))
            except Exception as e: pass
    best, rbest = min(cands, key=lambda z: z[0])
    rec = dict(F=rm.F, nB=int(len(B)), nD=int(len(D)), focus=list(md.focus), status=r['status'], solver=r['solver'], t_sdp=r['t'],
               rigorous=rig, cert_mineig=float(mineig), cert_scale=float(c), rank=r['rank_eps'],
               eig_top=[float(e) for e in r['eig'][:5]], rounded=rr, round_refined=rl, multistart=rmv, best_feasible=best,
               gap_rel=(best - rig) / max(rig, 1e-12), EL_bound=(1 - rig) / (1 + rig), EL_best=(1 - best) / (1 + best),
               sdp_time=t_sdp, total_time=time.time() - t0)
    if save_prefix:
        np.save(save_prefix + '_Y.npy', r['X']); np.save(save_prefix + '_rbest.npy', rbest)
        json.dump(dict(mu={str(k): v for k, v in r['mu'].items()}, lam={str(k): v for k, v in r['lam'].items()}), open(save_prefix + '_dual.json', 'w'))
    return rec, rm, r, rbest

def sets_2d(md, tf, margin):
    pts = md.x; tgt = np.array([tf(a, b) for a, b in pts]); Lx, Ly = md.L
    nx = int(round(np.sqrt(md.nx)))
    def nbrs(i):
        a, b = divmod(i, nx); return [((a+1)%nx)*nx+b, ((a-1)%nx)*nx+b, a*nx+(b+1)%nx, a*nx+(b-1)%nx]
    bnd = np.array([pts[i] for i in range(len(pts)) if any(tgt[i] != tgt[j] for j in nbrs(i))]).reshape(-1, 2)
    def dist(p):
        if len(bnd) == 0: return np.inf
        d = bnd - p; d[:, 0] = (d[:, 0] + Lx/2) % Lx - Lx/2; d[:, 1] = (d[:, 1] + Ly/2) % Ly - Ly/2
        return np.min(np.hypot(d[:, 0], d[:, 1]))
    dd = np.array([dist(p) for p in pts])
    return np.where(tgt & (dd >= margin))[0], np.where(~tgt & (dd >= margin))[0], tgt

# ----------------------------------------------------------------------------- illuminator catalogue (session 3; area-weighted 2-D quadratures)
# DUV: lambda 193 nm, NA 1.35, immersion (defocus in water, n_img = 1.44). EUV: lambda 13.5 nm, NA 0.33, vacuum.
NA_DUV, LAM_DUV, NIMG_DUV = 1.35, 193.0, 1.44
NA_EUV, LAM_EUV, NIMG_EUV = 0.33, 13.5, 1.0
def catalogue(NA):
    return {
        'coh':       (np.array([[0.0, 0.0]]), np.ones(1)),                    # on-axis coherent point
        'conv0.5':   disk_source(0.5, NA, n=21),                              # disk sigma 0.5, 21x21 grid -> 349 pts
        'conv0.7':   disk_source(0.7, NA, n=21),
        'conv0.9':   disk_source(0.9, NA, n=21),
        'dipole0.8': dipole_source(0.8, NA, pole_radius=0.15, n=9),           # two poles r=0.15 at +-0.8 (x axis), 69 pts each
        'dipole0.7': dipole_source(0.7, NA, pole_radius=0.15, n=9),
        'annular':   annular_source(0.6, 0.9, NA, n=25),                      # ring 0.6-0.9, 25x25 grid -> ~280 pts
        'quasar':    quasar_source(0.6, 0.9, NA, opening_deg=30.0, n=31),     # four 30-deg sectors of the 0.6-0.9 ring
        'annularE':  annular_source(0.5, 0.9, NA, n=25),                      # EUV ring 0.5-0.9
        'quasarE':   quasar_source(0.5, 0.9, NA, opening_deg=30.0, n=31),
    }
SRC_DUV = catalogue(NA_DUV); SRC_EUV = catalogue(NA_EUV)
TAU_DUV = 5.0    # edge tolerance for 1-D DUV targets (CD 40-50 nm -> +-10% CD per edge)
def analyze_el(rec):
    """Add full-range exposure-latitude fields (window width / centre dose) and a degenerate flag."""
    rec['EL_full_bound'] = el_full(rec['rigorous']); rec['EL_full_best'] = el_full(rec['best_feasible'])
    rec['degenerate'] = bool(rec['rigorous'] < 1e-6)
    rec['cert_violation'] = bool(rec['best_feasible'] < rec['rigorous'] * (1 - 1e-6) - 1e-12)   # a feasible mask below the bound = invalid certificate
    if rec['cert_violation']: print(f"WARNING: certificate violation: best_feasible {rec['best_feasible']} < rigorous {rec['rigorous']}", flush=True)
    if rec['degenerate']: rec['gap_rel'] = None; rec['gap_abs'] = rec['best_feasible'] - rec['rigorous']
    return rec
_analyze_orig = analyze
def analyze(md, B, D, **kw):
    rec, rm, r, rbest = _analyze_orig(md, B, D, **kw); rec['n_img'] = float(getattr(md, 'n_img', 1.0)); rec['dead_orders_dropped'] = int(getattr(md, 'n_dead_dropped', 0)); return analyze_el(rec), rm, r, rbest
