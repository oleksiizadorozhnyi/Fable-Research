# Strong single-mask search + optimal-face exploration (CLARABEL only) + explicit 2-exposure realization, 320 nm 2D layouts.
from pipeline import *
from scipy.optimize import minimize
from scipy.special import logsumexp, softmax
lam, NA = 193.0, 1.35; ann = annular_source_2d(0.6, 0.9, NA, n=9); qua = quasar_source_2d(0.6, 0.9, NA, n=11)
Tj = lambda a, b: (abs(a-160) < 30 and 60 < b < 220) or (abs(b-40) < 30)
lineend = lambda a, b: (abs(a-100) < 25 and 30 < b < 210) or (abs(a-220) < 25 and (b < 150 or b > 270))
out = []
def log(*a): print(*a, flush=True)
def refine_pair(rm, R0, B, D, iters=600):
    def fg(p):
        R = p.reshape(rm.F, 2); f = 0; g = np.zeros_like(R)
        for z in rm.focus:
            C = rm.C[z]; flds = [C @ R[:, j] for j in range(2)]; I = sum(np.einsum('s,sx->x', rm.w, np.abs(f_)**2) for f_ in flds) + 1e-300
            lI = np.log(I); pD = softmax(40*lI[D]); pB = softmax(-40*lI[B]); f += logsumexp(40*lI[D])/40 + logsumexp(-40*lI[B])/40
            for j in range(2):
                G = 2*np.einsum('s,sx,sxk->xk', rm.w, flds[j].conj(), C).real/I[:, None]; g[:, j] += pD @ G[D] - pB @ G[B]
        return f, g.ravel()
    res = minimize(fg, R0.ravel(), jac=True, method='L-BFGS-B', options=dict(maxiter=iters)); R = res.x.reshape(rm.F, 2)
    I = sum(rm.image(R[:, j]) for j in range(2)); return I[D].max()/I[B].min(), R
for name, tf, src in [("Tjunction320_annular", Tj, ann), ("Tjunction320_quasar", Tj, qua), ("lineend320_annular", lineend, ann)]:
    t0 = time.time(); md, B, D, tgt = build_2d((320,320), lam, NA, src, tf, 14, spacing=16.0, fine=2.0); rm = RealModel(md)
    r = sdp_bound_real(rm, B, D); rig, me, c = rigorize_certificate(rm, r['mu'], r['lam']); log(name, "SDP", round(r['t'],5), r['solver'], round(time.time()-t0))
    rec = dict(name=name, F=rm.F, bound=rig, ipm_rank=r['rank_eps'], ipm_eigs=[float(e) for e in r['eig'][:4]])
    # face exploration: trace minimization only (1 extra SDP), CLARABEL only
    try:
        best, ranks = face_min_rank(rm, B, D, r['t'], ntries=1, solvers=("CLARABEL",)); rec.update(face_rank_trace=best[0], face_eigs=[float(e) for e in best[2][:4]]); Ys = [r['X'], best[1]]; log(name, "face rank", best[0], round(time.time()-t0))
    except Exception as e: rec.update(face_rank_trace=None); Ys = [r['X']]; log(name, "face failed", e)
    # strong primal (vectorized SLSQP polish inside)
    v, rb = strong_primal(rm, B, D, Ys, nround=200, nstart=12); rec.update(best_single=v, gap_single=(v-rig)/rig); log(name, "best single", round(v,5), round(time.time()-t0))
    # explicit 2-exposure realization from top-2 eigenvectors of the lowest-rank Y
    Y = Ys[-1]; w, V = np.linalg.eigh(Y); idx = np.argsort(w)[::-1]; R0 = V[:, idx[:2]]*np.sqrt(np.clip(w[idx[:2]], 0, None))
    v2, R2 = refine_pair(rm, R0, B, D); rec.update(two_exposure_ratio=v2, gap_two=(v2-rig)/rig, time=time.time()-t0); log(name, "two-exposure", round(v2,5), round(time.time()-t0))
    print(json.dumps(rec), flush=True); out.append(rec); np.save(f'/home/user/Fable-Research/results/exp13_{name}_rbest.npy', rb); np.save(f'/home/user/Fable-Research/results/exp13_{name}_R2.npy', R2)
    json.dump(out, open('/home/user/Fable-Research/results/exp13_2d_strong.json','w'), indent=1)
