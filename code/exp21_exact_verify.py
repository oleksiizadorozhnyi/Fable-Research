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
    """Exact (50-digit) rebuild of S~ = Re U^H (sum_d mu_d Q_d - kappa sum_b lam_b Q_b) U from the model equations.
       Pupil test: the double-precision margin |g|^2 lam^2 - NA^2 decides every clear case (relative margin > 1e-9 is far
       beyond double rounding); borderline points are decided exactly with rational arithmetic (Fraction of the exact
       binary doubles). At best focus the pupil vectors are 0/1 and T = w * (inside^T inside) is an exact integer count;
       with defocus the phases are evaluated per source point at 50 digits."""
    from fractions import Fraction
    F = md.F; src = np.asarray(md.src, float); w = np.asarray(md.w_src if hasattr(md, 'w_src') else md.w, float); lam = float(md.lam); NA = float(md.NA); n = float(md.n_img)
    Lf = [float(v) for v in md.L]; dim = md.dim; k = np.asarray(md.k, float)
    f2 = np.zeros((F, 2)); f2[:, :dim] = k / np.array(Lf[:dim])
    g = f2[None, :, :] + src[:, None, :] / lam; marg = (g[:, :, 0]**2 + g[:, :, 1]**2) * lam * lam - NA * NA   # (S, F)
    inside = marg <= 0; border = np.abs(marg) < 1e-9 * NA * NA
    if border.any():   # exact rational decision for borderline (source, order) pairs
        L_fr = [Fraction(v) for v in Lf]; lam_fr = Fraction(lam); NA_fr = Fraction(NA)
        for si, ki in zip(*np.nonzero(border)):
            fx = Fraction(int(md.k[ki][0])) / L_fr[0] + Fraction(float(src[si, 0])) / lam_fr
            fy = (Fraction(int(md.k[ki][1])) / L_fr[1] if dim == 2 else Fraction(0)) + Fraction(float(src[si, 1])) / lam_fr
            inside[si, ki] = (fx * fx + fy * fy) * lam_fr * lam_fr <= NA_fr * NA_fr
    assert np.allclose(w, w[0]), "equal-area weights expected"; w0 = mp.mpf(float(w[0]))
    S = mp.zeros(F, F)
    for z in md.focus:
        if float(z) == 0.0:
            N = inside.astype(np.int64).T @ inside.astype(np.int64)                      # exact integer counts (< 2^53)
            T = mp.matrix([[w0 * int(N[a, b]) for b in range(F)] for a in range(F)])
        else:
            zz = mp.mpf(float(z)); lam_m = mp.mpf(lam); n_m = mp.mpf(n); L_m = [mp.mpf(v) for v in Lf]
            T = mp.zeros(F, F)
            fk = [[mp.mpf(int(md.k[a][0])) / L_m[0], (mp.mpf(int(md.k[a][1])) / L_m[1] if dim == 2 else mp.mpf(0))] for a in range(F)]
            for si in range(len(src)):
                live = np.nonzero(inside[si])[0]
                if len(live) < 1: continue
                sx = mp.mpf(float(src[si, 0])) / lam_m; sy = mp.mpf(float(src[si, 1])) / lam_m; ph = {}
                for a in live:
                    gn2 = (fk[a][0] + sx)**2 + (fk[a][1] + sy)**2
                    ph[a] = mp.expj(2 * mp.pi * n_m * zz / lam_m * (mp.sqrt(max(1 - gn2 * lam_m * lam_m / (n_m * n_m), 0)) - 1))
                for a in live:
                    pa = w0 * ph[a]
                    for b in live: T[a, b] += pa * mp.conj(ph[b])
        for (zc, ix), v in list(mu.items()) + [((zc, ix), -kappa * v) for (zc, ix), v in lamb.items()]:
            if zc != z or v == 0: continue
            x = md.x[ix]; phx = [mp.expj(2 * mp.pi * (mp.mpf(int(md.k[a][0])) / mp.mpf(Lf[0]) * mp.mpf(float(x[0])) + ((mp.mpf(int(md.k[a][1])) / mp.mpf(Lf[1]) * mp.mpf(float(x[1]))) if dim == 2 else 0))) for a in range(F)]
            vv = mp.mpf(float(v))
            for a in range(F):
                pa = vv * phx[a]
                for b in range(F): S[a, b] += pa * T[a, b] * mp.conj(phx[b])
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
