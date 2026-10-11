# Strong single-mask search + optimal-face exploration + 2-exposure realization for the 320 nm 2D layouts.
from pipeline import *
lam, NA = 193.0, 1.35; ann = annular_source_2d(0.6, 0.9, NA, n=9); qua = quasar_source_2d(0.6, 0.9, NA, n=11)
Tj = lambda a, b: (abs(a-160) < 30 and 60 < b < 220) or (abs(b-40) < 30)
lineend = lambda a, b: (abs(a-100) < 25 and 30 < b < 210) or (abs(a-220) < 25 and (b < 150 or b > 270))
out = []
for name, tf, src in [("Tjunction320_annular", Tj, ann), ("Tjunction320_quasar", Tj, qua), ("lineend320_annular", lineend, ann)]:
    md, B, D, tgt = build_2d((320,320), lam, NA, src, tf, 14, spacing=16.0, fine=2.0); rm = RealModel(md)
    t0 = time.time(); r = sdp_bound_real(rm, B, D); rig, me, c = rigorize_certificate(rm, r['mu'], r['lam'])
    # face exploration (few tries; each SDP ~100-300 s)
    best, ranks = face_min_rank(rm, B, D, r['t'], ntries=3)
    # strong primal
    v, rb = strong_primal(rm, B, D, [r['X'], best[1]], nround=300, nstart=24)
    # 2-exposure realization from the lowest-rank face point: take top-2 eigenvectors, refine the PAIR on the summed image
    w, V = np.linalg.eigh(best[1]); idx = np.argsort(w)[::-1]
    def ratio_pair(R):  # R: (F,2)
        worst = 0
        for z in rm.focus:
            I = sum(rm.image(R[:, j], z) for j in range(R.shape[1])); worst = max(worst, I[D].max()/I[B].min())
        return worst
    R0 = V[:, idx[:2]]*np.sqrt(np.clip(w[idx[:2]], 0, None)); v2 = ratio_pair(R0)
    from scipy.optimize import minimize
    from scipy.special import logsumexp, softmax
    def fg(p):
        R = p.reshape(rm.F, 2); f = 0; g = np.zeros_like(R)
        for z in rm.focus:
            C = rm.C[z]; flds = [C @ R[:, j] for j in range(2)]; I = sum(np.einsum('s,sx->x', rm.w, np.abs(f_)**2) for f_ in flds) + 1e-300
            lI = np.log(I); pD = softmax(40*lI[D]); pB = softmax(-40*lI[B]); f += logsumexp(40*lI[D])/40 + logsumexp(-40*lI[B])/40
            for j in range(2):
                G = 2*np.einsum('s,sx,sxk->xk', rm.w, flds[j].conj(), C).real/I[:, None]; g[:, j] += pD @ G[D] - pB @ G[B]
        return f, g.ravel()
    res = minimize(fg, R0.ravel(), jac=True, method='L-BFGS-B', options=dict(maxiter=600)); R2 = res.x.reshape(rm.F, 2); v2r = ratio_pair(R2)
    rec = dict(name=name, F=rm.F, bound=rig, face_ranks=ranks, min_face_rank=best[0], face_eigs=[float(e) for e in best[2][:4]], best_single=v, gap_single=(v-rig)/rig,
               two_exposure_ratio=v2r, gap_two=(v2r-rig)/rig, time=time.time()-t0)
    print(json.dumps(rec), flush=True); out.append(rec); np.save(f'/home/user/Fable-Research/results/exp13_{name}_rbest.npy', rb); np.save(f'/home/user/Fable-Research/results/exp13_{name}_R2.npy', R2)
    json.dump(out, open('/home/user/Fable-Research/results/exp13_2d_strong.json','w'), indent=1)
