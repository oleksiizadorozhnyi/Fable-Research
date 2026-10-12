# Much stronger single-mask search for the unresolved quasar T-junction: many restarts, SLSQP polish of top candidates, basin hopping.
from pipeline import *
lam, NA = 193.0, 1.35; qua = SRC_DUV['quasar']
Tj = lambda a, b: (abs(a-160) < 30 and 60 < b < 220) or (abs(b-40) < 30)
md, B, D, tgt = build_2d((320,320), lam, NA, qua, Tj, 6, spacing=16.0, fine=2.0); rm = RealModel(md)
r = sdp_bound_real(rm, B, D); rig, _, _ = rigorize_certificate(rm, r['mu'], r['lam']); print("bound", rig, flush=True)
rng = np.random.default_rng(11); best = (np.inf, None); t0 = time.time()
rb0 = np.load('/home/user/Fable-Research/results/exp13_Tjunction320_quasar_rbest.npy'); v0 = rm.ratio(rb0, B, D); best = (v0, rb0); print("start from exp13 best", v0, flush=True)
# basin hopping around the incumbent + fresh multistarts
for it in range(60):
    if it % 3 == 0: x0 = rng.standard_normal(rm.F)
    else: x0 = best[1] + 0.3*np.linalg.norm(best[1])*rng.standard_normal(rm.F)/np.sqrt(rm.F)
    v, x = local_refine_real(rm, x0, B, D, iters=600)
    if v < best[0]*1.02:
        try: v2, x2 = exact_minimax_polish(rm, x, B, D, maxiter=400); 
        except Exception: v2, x2 = v, x
        if v2 < v: v, x = v2, x2
    if v < best[0]: best = (v, x); print(f"it {it}: new best {v:.5f} gap {(v-rig)/rig:.2%} t={time.time()-t0:.0f}s", flush=True)
json.dump(dict(bound=rig, best_single=best[0], gap=(best[0]-rig)/rig, iters=60, time=time.time()-t0), open('/home/user/Fable-Research/results/exp18_quasar_strong.json','w'), indent=1)
np.save('/home/user/Fable-Research/results/exp18_quasar_rbest.npy', best[1]); print("done", best[0])
