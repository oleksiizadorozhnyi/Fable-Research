# Exact verification of every reported certificate (referee item 1): the certificate matrix is REBUILT in 50-digit
# arithmetic from the model equations and the instance data (sample coordinates, source points and weights, cell,
# wavelength, NA, index, focus set; all exact as the double-precision numbers that define the instance), using the dual
# weights (mu, lambda) and the rigorization scale kappa found in double precision. S~ = Re(U^H (sum mu Q_d - kappa sum
# lam Q_b) U) is then factorized (Cholesky at 50 digits) after subtracting delta = 1e-30 I, which proves lambda_min(S~)
# >= 1e-30 up to the 50-digit rounding error (~1e-45 F^2 ||S||). Hence kappa*sum(lambda) is a valid lower bound on rho
# for the stated discrete instance. The instance data and duals are archived in results/certs/<name>.npz.
from verify_instances import *
import mpmath as mp, os
mp.mp.dps = 50; out = []; os.makedirs('/home/user/Fable-Research/results/certs', exist_ok=True)
def mp_cert(md, mu, lamb, kappa):
    F = md.F; src = md.src; w = md.w_src if hasattr(md, 'w_src') else md.w; lam = mp.mpf(md.lam); NA = mp.mpf(md.NA); n = mp.mpf(md.n_img)
    L = [mp.mpf(float(v)) for v in md.L]; dim = md.dim
    f = [[mp.mpf(int(kk[i])) / L[i] for i in range(dim)] + ([mp.mpf(0)] if dim == 1 else []) for kk in md.k]   # exact rationals k/L
    S = mp.zeros(F, F)   # complex accumulator
    for z in md.focus:
        zz = mp.mpf(float(z))
        # T_z = sum_s w_s P_s P_s^H  (F x F), pupil and defocus from the model equations
        P = []
        for s_, ws in zip(src, w):
            sx, sy = mp.mpf(float(s_[0])) / lam, mp.mpf(float(s_[1])) / lam; col = []
            for fk in f:
                gx, gy = fk[0] + sx, fk[1] + sy; gn2 = gx*gx + gy*gy
                if gn2 * lam * lam <= NA * NA: col.append(mp.expj(2 * mp.pi * n * zz / lam * (mp.sqrt(max(1 - gn2 * lam * lam / (n * n), 0)) - 1)))
                else: col.append(mp.mpc(0))
            P.append((mp.mpf(float(ws)), col))
        T = mp.zeros(F, F)
        for ws, col in P:
            for a in range(F):
                if col[a] == 0: continue
                ca = ws * col[a]
                for b in range(F):
                    if col[b] != 0: T[a, b] += ca * mp.conj(col[b])
        # constraints: Q_x = D_x T D_x^H with D_x = diag(exp(2 pi i f.x))
        for (zc, ix), v in list(mu.items()) + [((zc, ix), -kappa * v) for (zc, ix), v in lamb.items()]:
            if zc != z or v == 0: continue
            x = md.x[ix]; ph = [mp.expj(2 * mp.pi * sum(f[k][i] * mp.mpf(float(x[i])) for i in range(dim))) for k in range(F)]
            vv = mp.mpf(float(v))
            for a in range(F):
                pa = vv * ph[a]
                for b in range(F): S[a, b] += pa * T[a, b] * mp.conj(ph[b])
    # real structure: U unitary basis (exact entries)
    U = mp.zeros(F, F); col = 0; seen = set(); r2 = 1 / mp.sqrt(2)
    for i in range(F):
        j = int(md.flip[i])
        if i in seen: continue
        if j == i: U[i, col] = 1; col += 1; seen.add(i)
        else: U[i, col] = r2; U[j, col] = r2; col += 1; U[i, col] = mp.mpc(0, 1) * r2; U[j, col] = mp.mpc(0, -1) * r2; col += 1; seen.update([i, j])
    St = U.H * S * U; Sr = mp.zeros(F, F)
    for a in range(F):
        for b in range(F): Sr[a, b] = (mp.re(St[a, b]) + mp.re(St[b, a])) / 2
    return Sr
import sys
LIMIT = int(sys.argv[1]) if len(sys.argv) > 1 else 10**9
for ii, (name, md, B, D) in enumerate(verify_instances()):
    if ii >= LIMIT: break
    rm = RealModel(md); t0 = time.time()
    try: r = sdp_bound_real(rm, B, D, solvers=("CLARABEL",))
    except Exception: r = sdp_bound_real(rm, B, D, solvers=("SCS",))
    rig, me, c = rigorize_certificate(rm, r['mu'], r['lam']); kappa = float(c)
    Sr = mp_cert(md, r['mu'], r['lam'], kappa); F = md.F
    norm = max(abs(Sr[a, b]) for a in range(F) for b in range(F))
    if norm < mp.mpf('1e-40'): status = 'zero certificate (PSD trivially)'; lmin = mp.mpf(0); chol = True
    else:
        try: mp.cholesky(Sr - mp.mpf('1e-30') * mp.eye(F)); chol = True
        except Exception: chol = False
        try: lmin = min(mp.eigsy(Sr, eigvals_only=True)) if F <= 60 else None
        except Exception: lmin = None
        status = 'PD with margin 1e-30' if chol else 'FAILED'
    slam = sum(r['lam'].values()); smu = sum(r['mu'].values()); bound_exact = kappa * slam / smu
    rec = dict(name=name, F=F, nB=len(B), nD=len(D), n_src=int(len(md.src)), kappa=kappa, rigorous=rig, bound_from_duals=bound_exact, cert_norm=float(norm), lambda_min_50digit=(float(lmin) if lmin is not None else None), cholesky_margin_1e30=chol, status=status, time=time.time() - t0)
    np.savez(f'/home/user/Fable-Research/results/certs/{name}.npz', x=md.x, B=np.array(B), D=np.array(D), k=md.k, L=md.L, wavelength=md.lam, NA=md.NA, n_img=md.n_img, focus=np.array(md.focus), src=md.src, w=np.asarray(md.w_src if hasattr(md, 'w_src') else md.w), mu_keys=np.array([(float(z), int(ix)) for (z, ix) in r['mu']]), mu=np.array(list(r['mu'].values())), lam_keys=np.array([(float(z), int(ix)) for (z, ix) in r['lam']]), lam=np.array(list(r['lam'].values())), kappa=kappa, bound=rig)
    print(json.dumps(rec), flush=True); out.append(rec); json.dump(out, open('/home/user/Fable-Research/results/exp21_exact_verify.json', 'w'), indent=1)
