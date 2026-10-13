# Exact (50-digit) rebuild of certificate matrices from the model equations and the instance data.
# Pupil test: the double-precision margin |g|^2 lam^2 - NA^2 decides every clear case (relative margin > 1e-9 is far beyond
# double rounding); borderline (source, order) pairs are decided exactly with rational arithmetic on the exact binary
# doubles. At best focus the pupil vectors are 0/1 and T = w * (inside^T inside) is an exact integer count; with defocus
# the phases are evaluated per source point at 50 digits. Returns the complex Hermitian SD, SB (F x F) in the order
# basis; callers apply the real-mask basis U (Fourier certificates) or the pixel map G (pixel certificates).
import numpy as np, mpmath as mp
from fractions import Fraction
mp.mp.dps = 50

def mp_source_matrices(md, mu, lamb):
    F = md.F; src = np.asarray(md.src, float); w = np.asarray(md.w_src if hasattr(md, 'w_src') else md.w, float); lam = float(md.lam); NA = float(md.NA); n = float(md.n_img)
    Lf = [float(v) for v in md.L]; dim = md.dim; k = np.asarray(md.k, float)
    f2 = np.zeros((F, 2)); f2[:, :dim] = k / np.array(Lf[:dim])
    g = f2[None, :, :] + src[:, None, :] / lam; marg = (g[:, :, 0]**2 + g[:, :, 1]**2) * lam * lam - NA * NA
    inside = marg <= 0; border = np.abs(marg) < 1e-9 * NA * NA
    if border.any():
        L_fr = [Fraction(v) for v in Lf]; lam_fr = Fraction(lam); NA_fr = Fraction(NA)
        for si, ki in zip(*np.nonzero(border)):
            fx = Fraction(int(md.k[ki][0])) / L_fr[0] + Fraction(float(src[si, 0])) / lam_fr
            fy = (Fraction(int(md.k[ki][1])) / L_fr[1] if dim == 2 else Fraction(0)) + Fraction(float(src[si, 1])) / lam_fr
            inside[si, ki] = (fx * fx + fy * fy) * lam_fr * lam_fr <= NA_fr * NA_fr
    assert np.allclose(w, w[0]), "equal-area weights expected"; w0 = mp.mpf(float(w[0]))
    S = mp.zeros(F, F); SB = mp.zeros(F, F)
    for z in md.focus:
        if float(z) == 0.0:
            Nc = inside.astype(np.int64).T @ inside.astype(np.int64)
            T = mp.matrix([[w0 * int(Nc[a, b]) for b in range(F)] for a in range(F)])
        else:
            zz = mp.mpf(float(z)); lam_m = mp.mpf(lam); n_m = mp.mpf(n); L_m = [mp.mpf(v) for v in Lf]; T = mp.zeros(F, F)
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
        for (zc, ix), v, tgt in [(kk, v, S) for kk, v in mu.items()] + [(kk, v, SB) for kk, v in lamb.items()]:
            if zc != z or v == 0: continue
            x = md.x[ix]; phx = [mp.expj(2 * mp.pi * (mp.mpf(int(md.k[a][0])) / mp.mpf(Lf[0]) * mp.mpf(float(x[0])) + ((mp.mpf(int(md.k[a][1])) / mp.mpf(Lf[1]) * mp.mpf(float(x[1]))) if dim == 2 else 0))) for a in range(F)]
            vv = mp.mpf(float(v))
            for a in range(F):
                pa = vv * phx[a]
                for b in range(F): tgt[a, b] += pa * T[a, b] * mp.conj(phx[b])
    # model convention: Q_ab = sum_s w_s conj(c_a) c_b  (M^H Q M = sum_s w |c^T M|^2), i.e. the conjugate of D T D^H built above
    for a in range(F):
        for b in range(F): S[a, b] = mp.conj(S[a, b]); SB[a, b] = mp.conj(SB[a, b])
    return S, SB

def real_basis(md):
    F = md.F; U = mp.zeros(F, F); col = 0; seen = set(); r2 = 1 / mp.sqrt(2)
    for i in range(F):
        j = int(md.flip[i])
        if i in seen: continue
        if j == i: U[i, col] = 1; col += 1; seen.add(i)
        else: U[i, col] = r2; U[j, col] = r2; col += 1; U[i, col] = mp.mpc(0, 1) * r2; U[j, col] = mp.mpc(0, -1) * r2; col += 1; seen.update([i, j])
    return U

def pixel_basis(md, N):
    """Exact pixel-integrated map G (F x N): G_{k,n} = sinc(k/N) e^{-2 pi i k (n+1/2)/N} / N."""
    F = md.F; G = mp.zeros(F, N)
    for a in range(F):
        kk = int(md.k[a][0]); sinc = mp.mpf(1) if kk == 0 else mp.sin(mp.pi * kk / N) / (mp.pi * kk / N)
        for nn in range(N): G[a, nn] = sinc * mp.expj(-2 * mp.pi * kk * (nn + mp.mpf(1) / 2) / N) / N
    return G

def real_sym(X, V):
    """Re(V^H X V) symmetrized (exact 50-digit)."""
    Xt = V.H * X * V; n = Xt.rows; R = mp.zeros(n, n)
    for a in range(n):
        for b in range(n): R[a, b] = (mp.re(Xt[a, b]) + mp.re(Xt[b, a])) / 2
    return R

def is_pd(Sr, delta):
    try: mp.cholesky(Sr - delta * mp.eye(Sr.rows)); return True
    except Exception: return False

def double_below(x_mp):
    """Largest double <= the 50-digit value x_mp (so that the SAVED scale is never above the verified one)."""
    d = float(x_mp)
    while mp.mpf(d) > x_mp: d = float(np.nextafter(d, -np.inf))
    return d

def verify_scale(SDm, SBm, kappa, delta=mp.mpf('1e-30'), extra=None):
    """Largest double kappa' <= kappa with SDm - kappa' SBm [- extra] >= delta I at 50 digits (bisection when needed).
       Returns (kappa_double, Sr, bisected, ok). The returned kappa_double is re-verified as a double before returning."""
    def mat(c): return SDm - mp.mpf(c) * SBm - (extra if extra is not None else 0)
    if is_pd(mat(kappa), delta): return kappa, mat(kappa), False, True
    lo, hi = mp.mpf(0), mp.mpf(kappa)
    if not is_pd(mat(0.0), delta): return 0.0, mat(0.0), True, False
    for _ in range(80):
        mid = (lo + hi) / 2
        if is_pd(mat(mid), delta): lo = mid
        else: hi = mid
    kd = double_below(lo)
    while not is_pd(mat(kd), delta): kd = float(np.nextafter(kd, -np.inf))     # guard: saved double must itself verify
    return kd, mat(kd), True, True
