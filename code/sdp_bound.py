"""
Lifted (SDP) relaxation of the dark/bright intensity-ratio mask synthesis problem, Fourier-order domain.

Problem (P):  minimize  t   over masks m (periodic, cell L, pixels N)
              s.t. I_m(x) <= t  for x in D (dark),   I_m(x) >= 1 for x in B (bright)
    I_m(x) = sum_s w_s |sum_k P_s(k) M_k e^{2 pi i f_k x}|^2 = M^H Qx M,   M = in-band orders of m.
Relaxation (SDP):  min t  s.t. Tr(Qx X) <= t (x in D), Tr(Qx X) >= 1 (x in B), X >= 0 (PSD, Hermitian FxF)
   + optional realness (Hermitian-symmetric orders -> X_{k,l} = conj(X_{-k,-l}))
   + optional nonnegativity of mask (1D exact cone: Toeplitz(M) PSD, needs M as variable -> use normalized version)
Dual certificate: mu>=0 on D, lam>=0 on B, sum mu = 1, sum_d mu_d Q_d - sum_b lam_b Q_b PSD  =>  t* >= sum lam_b.
"""
import numpy as np, cvxpy as cp

class FourierLitho1D:
    def __init__(self, L, lam, NA, src_pts, src_w, nx=None, defocus=(0.0,)):
        self.L, self.lam, self.NA = L, lam, NA
        self.src = np.asarray(src_pts, float); self.w = np.asarray(src_w, float)/np.sum(src_w)
        smax = np.max(np.abs(self.src))
        K = int(np.floor((NA + smax)*L/lam + 1e-9))   # |f_k| <= (NA+smax)/lam
        self.k = np.arange(-K, K+1); self.F = len(self.k); self.f = self.k/L
        self.nx = nx or (4*self.F+1)
        self.x = np.arange(self.nx)*L/self.nx
        self.defocus = defocus
        # Qx list per (defocus, x)
        self.Q = {}
        for dz in defocus:
            for ix, x in enumerate(self.x):
                Qx = np.zeros((self.F, self.F), complex)
                for s, ws in zip(self.src, self.w):
                    g = self.f + s/lam; inside = np.abs(g)*lam <= NA
                    gl = np.clip(g*lam, -1, 1)
                    ph = 2*np.pi/lam*dz*(np.sqrt(1-gl**2)-1)
                    a = inside*np.exp(1j*ph)*np.exp(2j*np.pi*self.f*x)   # field(x) = a^T M = a.conj()^H M
                    c = a.conj()
                    Qx += ws*np.outer(c, c.conj())       # M^H Qx M = |a^T M|^2
                self.Q[(dz, ix)] = Qx
    def image_from_orders(self, M, dz=0.0):
        return np.array([np.real(M.conj() @ self.Q[(dz, ix)] @ M) for ix in range(self.nx)])
    def orders_from_pixels(self, m):
        N = len(m); n = np.arange(N)
        return np.array([np.mean(m*np.exp(-2j*np.pi*kk*n/N)) for kk in self.k])

def sdp_ratio_bound(lt, bright_idx, dark_idx, real_mask=True, nonneg_mask=False, defocus=None, solver="CLARABEL", verbose=False):
    F = lt.F
    X = cp.Variable((F, F), hermitian=True)
    t = cp.Variable()
    cons = [X >> 0]
    dzs = defocus if defocus is not None else list(lt.defocus)
    for dz in dzs:
        for ix in dark_idx:  cons.append(cp.real(cp.trace(lt.Q[(dz, ix)] @ X)) <= t)
        for ix in bright_idx: cons.append(cp.real(cp.trace(lt.Q[(dz, ix)] @ X)) >= 1)
    if real_mask:
        # M_{-k} = conj(M_k)  =>  X_{k,l} = M_k conj(M_l) = conj(M_{-k}) M_{-l} = conj(X_{-k,-l})
        idx = lambda kk: kk + (F-1)//2
        for a in lt.k:
            for b in lt.k:
                if (a, b) < (-a, -b):
                    cons.append(X[idx(a), idx(b)] == cp.conj(X[idx(-a), idx(-b)]))
    if nonneg_mask:
        # exact 1D cone for in-band orders of nonnegative masks is Toeplitz(M) PSD; in lifted form we use the
        # necessary condition on X: the "second-moment Toeplitz" -- handled in sdp_ratio_bound_nonneg (needs M var).
        raise NotImplementedError
    prob = cp.Problem(cp.Minimize(t), cons)
    prob.solve(solver=solver, verbose=verbose)
    Xv = X.value
    ev = np.linalg.eigvalsh((Xv+Xv.conj().T)/2)
    return dict(t=t.value, X=Xv, eig=ev[::-1], status=prob.status, prob=prob)

def gaussian_rounding(lt, X, bright_idx, dark_idx, nsamp=2000, real_mask=True, rng=None, defocus=(0.0,)):
    """Sample M ~ CN(0,X) (or real-symmetric structure), rescale so min bright = 1, return best ratio."""
    rng = rng or np.random.default_rng(0)
    F = lt.F
    w, V = np.linalg.eigh((X+X.conj().T)/2); w = np.clip(w, 0, None)
    Sq = V*np.sqrt(w)
    best = (np.inf, None)
    for _ in range(nsamp):
        z = (rng.standard_normal(F)+1j*rng.standard_normal(F))/np.sqrt(2)
        M = Sq @ z
        if real_mask:  # enforce Hermitian symmetry by averaging with conj-flip
            M = 0.5*(M + np.conj(M[::-1]))
        worst = 0.0; ok = True
        for dz in defocus:
            I = lt.image_from_orders(M, dz)
            bmin = I[bright_idx].min(); dmax = I[dark_idx].max()
            if bmin <= 1e-12: ok = False; break
            worst = max(worst, dmax/bmin)
        if ok and worst < best[0]: best = (worst, M)
    return best

def dual_certificate_check(lt, mu, lam, bright_idx, dark_idx, dz=0.0):
    """Returns min eigenvalue of sum mu_d Q_d - sum lam_b Q_b (>= -tol certifies t* >= sum lam)."""
    S = sum(mu[i]*lt.Q[(dz, d)] for i, d in enumerate(dark_idx)) - sum(lam[i]*lt.Q[(dz, b)] for i, b in enumerate(bright_idx))
    return np.linalg.eigvalsh((S+S.conj().T)/2).min(), np.sum(lam)

if __name__ == "__main__":
    lam, NA = 193.0, 1.35
    print("Dense lines, coherent illumination: bound should be 1 (unprintable) for pitch < lam/NA=143nm, <1 above.")
    for pitch in [120, 140, 150, 170, 200, 260]:
        lt = FourierLitho1D(L=pitch, lam=lam, NA=NA, src_pts=[0.0], src_w=[1.0], nx=32)
        # target: line centered at 0 with width pitch/2: bright |x|< pitch/4 (interior 60%), dark elsewhere (interior 60%)
        xc = (lt.x + pitch/2) % pitch - pitch/2
        bright = np.where(np.abs(xc) < 0.3*pitch/2)[0]; dark = np.where(np.abs(xc) > pitch/2 - 0.3*pitch/2)[0]
        r = sdp_ratio_bound(lt, bright, dark)
        print(f"pitch={pitch}: F={lt.F} orders, SDP ratio bound t*={r['t']:.4f}, top eigs={np.round(r['eig'][:3],4)}, status={r['status']}")
