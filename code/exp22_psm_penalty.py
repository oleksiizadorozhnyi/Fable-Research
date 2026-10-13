# Certified bracket for the nonnegative-mask penalty relative to phase-shift masks, and the ILT suboptimality bracket:
#   rho*_+ / rho*_R - 1  in  [ L_+ / U_R - 1 ,  U_+ / L_R - 1 ],   ILT suboptimality in [ rho_ILT / U_R - 1 , rho_ILT / L_R - 1 ]
# with L_R, L_+ the RIGORIZED Fourier-PSM and pixel-integrated DNN certificates of exp3, both re-verified here exactly
# (50-digit rebuild of SD, SB from the model equations; pixel map G exact; DNN: Z = SD - kappa SB - Nmat >= 1e-30 I, which
# certifies rho >= kappa sum(lam)/sum(mu) over all PSD, entrywise-nonnegative X), U_R the best feasible PSM mask and U_+
# the best feasible nonnegative pixel-integrated mask (bounded L-BFGS from eight starts). Same instances and sample sets as exp3.
from pipeline import *
from exact_verify import *
from scipy.optimize import minimize
from scipy.special import logsumexp, softmax
lam, NA = 193.0, 1.35; out = []
e3 = {r['name']: r for r in json.load(open('/home/user/Fable-Research/results/exp3_pixel_binary.json'))}
conv = SRC_DUV['conv0.7']
cases = {"isodense": (600, [(0,50),(100,150),(200,250),(430,480)]), "mixed": (500, [(0,40),(100,180),(240,400)]), "dense120": (360, [(0,60),(120,180),(240,300)])}
def parse(d): return {(float(k.split(',')[0].strip('( ')), int(k.split(',')[1].strip(') '))): v for k, v in d.items()}
def refine_nonneg(md, B, D, m0, G, iters=600, beta=40.0):
    N = m0.size
    def fg(m):
        M = G @ m; f = 0.0; g = np.zeros(N)
        for z in md.focus:
            I, GM = md.image_grad(M, z); I = I + 1e-300; Gm = np.real(GM @ G)
            lI = np.log(I); pD = softmax(beta*lI[D]); pB = softmax(-beta*lI[B])
            f += logsumexp(beta*lI[D])/beta + logsumexp(-beta*lI[B])/beta; dl = Gm / I[:, None]; g += pD @ dl[D] - pB @ dl[B]
        return f, g
    res = minimize(fg, m0, jac=True, method='L-BFGS-B', bounds=[(0, 1)]*N, options=dict(maxiter=iters))
    m = res.x; M = G @ m; ratio = max(md.image(M, z)[D].max() for z in md.focus) / min(md.image(M, z)[B].min() for z in md.focus)
    return ratio, m
for name, (L, segs) in cases.items():
    md, B, D, tgt = build_1d(L, lam, NA, conv, segs, margin=TAU_DUV, spacing=5.0); rm = RealModel(md); r3 = e3[name]; N = r3['N_pixels']
    du = json.load(open(f'/home/user/Fable-Research/results/exp3_{name}_duals.json')); Nmat = np.load(f'/home/user/Fable-Research/results/exp3_{name}_dnn_Nmat.npy')
    # exact verification of the three certificates
    ver = {}
    for key, basis in [('F', 'U'), ('P', 'U'), ('NN', 'G')]:
        mu, lm, kap = parse(du[key]['mu']), parse(du[key]['lam']), du[key]['kappa']
        S, SB = mp_source_matrices(md, mu, lm); V = real_basis(md) if basis == 'U' else pixel_basis(md, N)
        SDm, SBm = real_sym(S, V), real_sym(SB, V); smu = sum(mu.values())
        if key != 'NN':
            kd, Sr, bis, ok = verify_scale(SDm, SBm, kap)
            ver[key] = dict(kappa=kap, kappa_exact=kd, bound_exact=kd * sum(lm.values()) / smu, bisected=bis, ok=ok)
        else:   # exact Jansson-type bound: t >= kappa*slam / (smu + delta N^2 |D| / lamQ), delta = max(0, -lambda_min(Z)) at 50 digits
            slam = sum(lm.values()); Nm = 0.5*(Nmat + Nmat.T)
            Z = SDm * smu - mp.mpf(kap) * SBm * smu - mp.matrix([[mp.mpf(float(Nm[a, b])) for b in range(N)] for a in range(N)])   # un-normalized (SDm, SBm are divided by smu)
            lminZ = min(mp.eigsy(Z, eigvals_only=True)); delta = max(mp.mpf(0), -lminZ)
            SQ, _ = mp_source_matrices(md, {(z, ix): 1.0 for z in md.focus for ix in D}, {}); lamQ = min(mp.re(v) for v in mp.eighe(SQ, eigvals_only=True)); nDk = len(D) * len(md.focus)
            b_exact = mp.mpf(kap) * slam / (smu + delta * N * N * nDk / lamQ) if lamQ > 0 else mp.mpf(0)
            ver[key] = dict(kappa=kap, kappa_exact=kap, bound_exact=float(b_exact), delta_exact=float(delta), lamQ_exact=float(lamQ), bisected=False, ok=bool(lamQ > 0))
        print(name, key, ver[key], flush=True)
    L_R = ver['F']['bound_exact']; L_plus = ver['NN']['bound_exact']; L_pixpsm = ver['P']['bound_exact']
    U_R = min(r3['best_psm_pipeline'], r3['ilt_psm_refined_on_ratio'])
    G = pixel_map(md, N); xs = (np.arange(N) + 0.5) * L / N; tf = lambda v: any(a <= (v % L) < b for a, b in segs); m_target = np.array([float(tf(v)) for v in xs])
    rng = np.random.default_rng(0); best = (np.inf, None)
    for m0 in [m_target, np.clip(m_target + 0.2*rng.standard_normal(N), 0, 1)] + [rng.uniform(0, 1, N) for _ in range(6)]:
        v, m = refine_nonneg(md, B, D, m0, G); best = min(best, (v, m), key=lambda t: t[0])
    U_plus = best[0]
    assert L_R <= U_R + 1e-12 and L_plus <= U_plus + 1e-12, (L_R, U_R, L_plus, U_plus)
    rec = dict(name=name, F=md.F, N_pixels=N, L_R=L_R, U_R=U_R, L_pixel_psm=L_pixpsm, L_plus=L_plus, U_plus=U_plus, penalty_lo=L_plus/U_R - 1, penalty_hi=U_plus/L_R - 1,
               ilt_psm=r3['ilt_psm_ratio'], exact=ver)
    print(json.dumps(rec), flush=True); out.append(rec); np.save(f'/home/user/Fable-Research/results/masks/exp22_{name}_nonneg_mask.npy', best[1])
    json.dump(out, open('/home/user/Fable-Research/results/exp22_psm_penalty.json', 'w'), indent=1)
