from pipeline import *
out = {}
# --- (a) Source-mask entitlement bound vs fixed-source bounds, 1D DUV
lam, NA = 193.0, 1.35
segs = [(0,50),(100,150),(200,250),(430,480)]; L = 600
res = []
for sname, (s, w) in {"coh": SRC_DUV['coh'], "conv0.7": SRC_DUV['conv0.7'], "dipole0.8": SRC_DUV['dipole0.8'], "conv0.9": SRC_DUV['conv0.9']}.items():
    md, B, D, tgt = build_1d(L, lam, NA, (s, w), segs, TAU_DUV, spacing=4.0, compress=False); rm = RealModel(md)
    r = sdp_bound_real(rm, B, D); smo = sdp_bound_smo(rm, B, D)
    res.append(dict(src=sname, nsrc=len(w), fixed_source_bound=r['t'], smo_bound=smo['t'], smo_status=smo['status'], smo_src_weights=[float(x) for x in smo['src_weights']]))
    print(json.dumps(res[-1]), flush=True)
out['smo'] = res
# --- (b) EUV: lam=13.5, NA=0.33, dense lines & iso-dense, pitch sweep, conventional sigma 0.9 and dipole
lam, NA = 13.5, 0.33; res = []
for pitch in [20, 22, 24, 26, 28, 32, 40]:
    for sname, (s, w) in {"conv0.9": SRC_EUV['conv0.9'], "dipole0.8": SRC_EUV['dipole0.8']}.items():
        segsE = [(i*pitch, i*pitch + pitch/2) for i in range(4)]
        md, B, D, tgt = build_1d(4*pitch, lam, NA, (s, w), segsE, 0.05*pitch, spacing=pitch/8, focus=(0.0, 50.0, -50.0))   # +-50nm defocus process window
        rec, rm, r, rb = analyze(md, B, D, nround=200, nstart=6)
        rec.update(pitch=pitch, src=sname, k1=pitch/2*NA/lam); res.append(rec); print(json.dumps({k: rec[k] for k in ['pitch','src','k1','F','t_sdp','rigorous','rank','best_feasible','gap_rel','EL_bound']}), flush=True)
out['euv_pw'] = res
# --- (c) Theorem 4 numerical check: necessity (random nonneg masks -> Toeplitz PSD) and sufficiency (Toeplitz PD -> nonneg realization)
rng = np.random.default_rng(0); Kmax = 6; viol = 0; N = 256
for trial in range(2000):
    m = rng.random(N) * (rng.random(N) < rng.random())      # random nonneg (sparse-ish) mask
    Mf = np.fft.fft(m)/N; M = np.array([Mf[k % N] for k in range(0, Kmax+1)])
    T = np.array([[M[abs(a-b)] if a >= b else np.conj(M[abs(a-b)]) for b in range(Kmax+1)] for a in range(Kmax+1)])
    if np.linalg.eigvalsh(T).min() < -1e-12*np.abs(M[0]): viol += 1
# sufficiency: random Hermitian-symmetric M with T(M) PD -> find nonneg function with these orders via LP (feasibility on fine grid)
from scipy.optimize import linprog
suff_fail = 0; tested = 0
for trial in range(200):
    v = rng.standard_normal(Kmax+1) + 1j*rng.standard_normal(Kmax+1); v[0] = abs(v[0])
    # make a PD Toeplitz by taking autocorrelation of a random vector + small DC boost
    h = rng.standard_normal(Kmax+1) + 1j*rng.standard_normal(Kmax+1)
    M = np.array([np.sum(h[:Kmax+1-k].conj()*h[k:]) for k in range(Kmax+1)]); M[0] *= 1.05
    T = np.array([[M[a-b] if a >= b else np.conj(M[b-a]) for b in range(Kmax+1)] for a in range(Kmax+1)])
    if np.linalg.eigvalsh(T).min() <= 0: continue
    tested += 1
    # LP: find f_n >= 0 on N grid with (1/N) sum f_n e^{-2 pi i k n/N} = M_k, k=0..Kmax  (real & imag parts)
    n = np.arange(N); rows = []; rhs = []
    for k in range(Kmax+1):
        e = np.exp(-2j*np.pi*k*n/N)/N; rows.append(e.real); rhs.append(M[k].real)
        if k > 0: rows.append(e.imag); rhs.append(M[k].imag)
    resLP = linprog(np.zeros(N), A_eq=np.array(rows), b_eq=np.array(rhs), bounds=[(0, None)]*N, method='highs')
    if not resLP.success: suff_fail += 1
out['theorem4'] = dict(necessity_violations=viol, necessity_trials=2000, sufficiency_tested=tested, sufficiency_failures=suff_fail)
print(json.dumps(out['theorem4']))
json.dump(out, open('/home/user/Fable-Research/results/exp5_smo_euv_toeplitz.json','w'), indent=1)
