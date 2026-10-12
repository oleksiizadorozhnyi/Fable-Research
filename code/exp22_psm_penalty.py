# Certified bracket for the nonnegative-mask (binary-like) penalty relative to phase-shift masks (referee item 5):
#   rho*_+ / rho*_R - 1  >=  L_+ / U_R - 1      (L_+: DNN lower bound on nonnegative pixel masks; U_R: feasible PSM ratio)
#   rho*_+ / rho*_R - 1  <=  U_+ / L_R - 1      (U_+: feasible nonnegative pixel mask ratio; L_R: certified PSM bound)
# and the ILT suboptimality bracket [rho_ILT/U_R - 1, rho_ILT/L_R - 1]. Same instances and sample sets as exp3.
from pipeline import *
from scipy.optimize import minimize
from scipy.special import logsumexp, softmax
lam, NA = 193.0, 1.35; out = []
e3 = {r['name']: r for r in json.load(open('/home/user/Fable-Research/results/exp3_pixel_binary.json'))}
conv = SRC_DUV['conv0.7']
cases = {"isodense": (600, [(0,50),(100,150),(200,250),(430,480)]), "mixed": (500, [(0,40),(100,180),(240,400)]), "dense120": (360, [(0,60),(120,180),(240,300)])}
def refine_nonneg(md, B, D, m0, iters=600, beta=40.0):
    """L-BFGS-B over pixel values in [0,1] (point-sampled pixel model, as in the DNN bound) on the smoothed log-ratio."""
    N = m0.size; n = np.arange(N); G = np.exp(-2j*np.pi*np.outer(md.k[:, 0], n)/N)/N          # orders_from_pixels as a matrix
    def fg(m):
        M = G @ m; f = 0.0; g = np.zeros(N)
        for z in md.focus:
            I, GM = md.image_grad(M, z); I = I + 1e-300; Gm = np.real(GM @ G)                   # dI/dm = Re(2 dI/dM* . G)
            lI = np.log(I); pD = softmax(beta*lI[D]); pB = softmax(-beta*lI[B])
            f += logsumexp(beta*lI[D])/beta + logsumexp(-beta*lI[B])/beta; dl = Gm / I[:, None]; g += pD @ dl[D] - pB @ dl[B]
        return f, g
    res = minimize(fg, m0, jac=True, method='L-BFGS-B', bounds=[(0, 1)]*N, options=dict(maxiter=iters))
    m = res.x; M = G @ m; ratio = max(md.image(M, z)[D].max() for z in md.focus) / min(md.image(M, z)[B].min() for z in md.focus)
    return ratio, m
for name, (L, segs) in cases.items():
    md, B, D, tgt = build_1d(L, lam, NA, conv, segs, margin=TAU_DUV, spacing=5.0); rm = RealModel(md); N = 120
    r3 = e3[name]; L_R = r3['t_fourier_psm']; L_plus = r3['t_pixel_nonneg_DNN']; U_R = min(r3['best_psm_pipeline'], r3['ilt_psm_refined_on_ratio'])
    xs = (np.arange(N) + 0.5) * L / N; tf = lambda v: any(a <= (v % L) < b for a, b in segs); m_target = np.array([float(tf(v)) for v in xs])
    rng = np.random.default_rng(0); best = (np.inf, None)
    for m0 in [m_target, np.clip(m_target + 0.2*rng.standard_normal(N), 0, 1)] + [rng.uniform(0, 1, N) for _ in range(6)]:
        v, m = refine_nonneg(md, B, D, m0); best = min(best, (v, m), key=lambda t: t[0])
    U_plus = min(best[0], r3['ilt_binary_ratio'])
    rec = dict(name=name, F=md.F, L_R=L_R, U_R=U_R, L_plus=L_plus, U_plus=U_plus, penalty_lo=L_plus/U_R - 1, penalty_hi=U_plus/L_R - 1,
               ilt_psm=r3['ilt_psm_ratio'], ilt_psm_subopt_lo=r3['ilt_psm_ratio']/U_R - 1, ilt_psm_subopt_hi=r3['ilt_psm_ratio']/L_R - 1,
               ilt_bin=r3['ilt_binary_ratio'], ilt_bin_subopt_lo=r3['ilt_binary_ratio']/U_plus - 1, ilt_bin_subopt_hi=r3['ilt_binary_ratio']/L_plus - 1)
    print(json.dumps(rec), flush=True); out.append(rec); np.save(f'/home/user/Fable-Research/results/exp22_{name}_nonneg_mask.npy', best[1])
    json.dump(out, open('/home/user/Fable-Research/results/exp22_psm_penalty.json', 'w'), indent=1)
