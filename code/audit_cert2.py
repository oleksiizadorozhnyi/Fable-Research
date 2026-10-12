# Audit (3),(4): dual extraction, rigorize_certificate (last def) validity: rho(r) >= bound for random real masks
import numpy as np, warnings, json; warnings.filterwarnings('ignore')
from certlitho import *; from pipeline import *
lam, NA = 193.0, 1.35; rng = np.random.default_rng(1)
conv = SRC_DUV['conv0.7']; dip = SRC_DUV['dipole0.8']; coh = SRC_DUV['coh']; n = NIMG_DUV
holes = lambda P, hw: (lambda a, b: (abs((a+P/2)%P-P/2) < hw) and (abs((b+P/2)%P-P/2) < hw))
Tj = lambda a, b: (abs(a-160) < 30 and 60 < b < 220) or (abs(b-40) < 30)
inst = []
inst.append(("isodense_conv", build_1d(600, lam, NA, conv, [(0,50),(100,150),(200,250),(430,480)], 8, spacing=4.0, n_img=n)))
inst.append(("mixed_coh", build_1d(500, lam, NA, coh, [(0,40),(100,180),(240,400)], 8, spacing=4.0, n_img=n)))
inst.append(("mixed_conv_PW80", build_1d(500, lam, NA, conv, [(0,40),(100,180),(240,400)], 8, spacing=4.0, focus=(0.0,80.0,-80.0), n_img=n)))
inst.append(("isodense_dip_PW80", build_1d(600, lam, NA, dip, [(0,50),(100,150),(200,250),(430,480)], 8, spacing=4.0, focus=(0.0,80.0,-80.0), n_img=n)))
inst.append(("duty95_coh", build_1d(400, lam, NA, coh, [(0,380)], 6, spacing=4.0, n_img=n)))
inst.append(("contacts120_quasar", build_2d((120,120), lam, NA, SRC_DUV['quasar'], holes(120,25), 10, spacing=16.0, fine=2.0, n_img=n)))
inst.append(("contacts160_annular_F2", build_2d((160,160), lam, NA, SRC_DUV['annular'], holes(160,35), 12, spacing=16.0, fine=2.0, focus=(0.0,60.0), n_img=n)))
#inst.append(("Tj320_annular", build_2d((320,320), lam, NA, SRC_DUV['annular'], Tj, 14, spacing=16.0, fine=2.0, n_img=n)))
for name, (md, B, D, tgt) in inst:
    rm = RealModel(md); r = sdp_bound_real(rm, B, D)
    # dual sign / mapping sanity: sum mu ~ 1 and sum lam ~ t (strong duality)
    smu = sum(r['mu'].values()); slam = sum(r['lam'].values())
    rig, me, c = rigorize_certificate(rm, r['mu'], r['lam'])
    nd = rigorize_certificate.last_null_dim; resid = rigorize_certificate.last_resid_null
    # full-space min eigenvalue of S_c (no projection) as an independent check of the projected test
    F = rm.F; SD = sum(v*rm.Q(ix, z) for (z, ix), v in r['mu'].items())/smu; SB = sum(v*rm.Q(ix, z) for (z, ix), v in r['lam'].items())/smu
    S = SD - c*SB; S = 0.5*(S+S.T); full_rel = np.linalg.eigvalsh(S).min()/np.linalg.norm(S, 2)
    # randomized weak-duality check on random real masks (gaussian, sparse, rounded-from-Y, local-refined)
    worst = np.inf; nviol = 0; N = 3000
    w_, V = np.linalg.eigh(r['X']); Sq = V*np.sqrt(np.clip(w_, 0, None))
    for i in range(N):
        kind = i % 3
        if kind == 0: x = rng.standard_normal(F)
        elif kind == 1: x = rng.standard_normal(F)*(rng.random(F) < 0.4); x = x if np.any(x) else rng.standard_normal(F)
        else: x = Sq @ rng.standard_normal(F)
        v = rm.ratio(x, B, D)
        if v < rig - 1e-9: nviol += 1
        worst = min(worst, v)
    # also local optimum masks
    rr, rv = round_from_Y(rm, r['X'], B, D, nsamp=200); rl, rvl = local_refine_real(rm, rv, B, D)
    worst = min(worst, rl)
    print(json.dumps(dict(name=name, F=F, nB=len(B), nD=len(D), t=r['t'], sum_mu=smu, sum_lam=slam, rig=rig, c=c, rel_mineig_proj=float(me), null_dim=nd, resid_null=resid,
                          full_space_rel_mineig=float(full_rel), worst_random_ratio=worst, best_local=rl, n_violations=nviol, status=r['status'], rank_eps=r['rank_eps'])), flush=True)
