"""
certlitho: certified bounds for partially coherent mask synthesis via lifted (SDP) relaxation.
Scalar Hopkins/Abbe thin-mask model, periodic cell (1D or 2D), diffraction-order (Fourier) domain.

Core objects
  Model: cell (Lx[,Ly]), wavelength lam, NA, source points (direction cosines, |s|<=sigma*NA) with weights,
         focus offsets. In-band orders k with |f_k| <= (NA+smax)/lam. For each sample point x and focus z:
         field_s(x) = c_{s,x}^H M  (c depends on s,x,z);  I(x) = sum_s w_s |c^H M|^2 = M^H Q_x M.
  Ratio problem (P):  t* = min_M  max_{x in D} I(x)  s.t. min_{x in B} I(x) >= 1   (over all foci jointly)
  SDP relaxation (R): t_R = min_X  t  s.t. Tr(Q_x X) <= t (D), Tr(Q_x X) >= 1 (B), X PSD [, symmetry constraints]
  Dual certificate: (mu, lam) >= 0, sum mu = 1, S = sum mu_x Q_x - sum lam_x Q_x PSD  ==>  t* >= sum lam.
"""
import numpy as np, cvxpy as cp, itertools, time

# ----------------------------------------------------------------------------- model
class Model:
    def __init__(self, L, lam, NA, src_pts, src_w, nx, focus=(0.0,), dim=1, xs=None):
        self.dim = dim
        self.L = np.atleast_1d(np.asarray(L, float)) if dim == 2 else np.array([float(L)])
        self.lam, self.NA = lam, NA
        self.src = np.atleast_2d(np.asarray(src_pts, float)) if dim == 2 else np.asarray(src_pts, float).reshape(-1, 1)
        self.w = np.asarray(src_w, float) / np.sum(src_w)
        self.focus = tuple(focus)
        smax = np.max(np.linalg.norm(self.src, axis=1))
        R = (NA + smax) / lam
        if dim == 1:
            K = int(np.floor(R * self.L[0] + 1e-9)); ks = [(k,) for k in range(-K, K + 1)]
        else:
            Kx = int(np.floor(R * self.L[0] + 1e-9)); Ky = int(np.floor(R * self.L[1] + 1e-9))
            ks = [(kx, ky) for kx in range(-Kx, Kx + 1) for ky in range(-Ky, Ky + 1)
                  if np.hypot(kx / self.L[0], ky / self.L[1]) <= R + 1e-12]
        self.k = np.array(ks)                       # (F, dim)
        self.F = len(self.k)
        self.f = self.k / self.L                    # (F, dim) frequencies
        self.kindex = {tuple(kk): i for i, kk in enumerate(self.k)}
        # sample grid
        if dim == 1:
            self.x = (np.arange(nx) * self.L[0] / nx).reshape(-1, 1) if xs is None else np.asarray(xs, float).reshape(-1, 1)
        elif xs is not None:
            self.x = np.asarray(xs, float).reshape(-1, 2)
        else:
            gx = np.arange(nx[0]) * self.L[0] / nx[0]; gy = np.arange(nx[1]) * self.L[1] / nx[1]
            self.x = np.array([(a, b) for a in gx for b in gy])
        self.nx = len(self.x)
        # coherent-field vectors C[z][s] : (nx, F) with field = C @ M
        self.C = {}
        for z in self.focus:
            Cz = []
            for s in self.src:
                g = self.f + s / lam                                    # (F, dim)
                gn = np.linalg.norm(g, axis=1)
                inside = (gn * lam <= NA + 1e-12)
                gl = np.clip(gn * lam, 0, 1)
                ph = 2 * np.pi / lam * z * (np.sqrt(1 - gl ** 2) - 1)
                P = inside * np.exp(1j * ph)                            # pupil incl. defocus
                E = np.exp(2j * np.pi * (self.x @ self.f.T))            # (nx, F)
                Cz.append(E * P[None, :])
            self.C[z] = np.stack(Cz)                                    # (S, nx, F)
        # flip permutation (k -> -k) for real-mask constraint
        self.flip = np.array([self.kindex[tuple(-kk)] for kk in self.k])

    def Q(self, ix, z=0.0):
        Cz = self.C[z][:, ix, :]                                        # (S, F)
        return np.einsum('s,sa,sb->ab', self.w, Cz.conj(), Cz)          # M^H Q M = sum_s w_s |C M|^2

    def image(self, M, z=0.0):
        fld = self.C[z] @ M                                             # (S, nx)
        return np.einsum('s,sx->x', self.w, np.abs(fld) ** 2)

    def orders_from_pixels(self, m):
        """m: pixel array (N,) or (Nx,Ny) on the cell; returns in-band orders (point-sampled mask)."""
        if self.dim == 1:
            Mf = np.fft.fft(m) / m.size
            return np.array([Mf[kk[0] % m.size] for kk in self.k])
        Mf = np.fft.fft2(m) / m.size
        return np.array([Mf[kk[0] % m.shape[0], kk[1] % m.shape[1]] for kk in self.k])

    def pixels_from_orders(self, M, N):
        """bandlimited mask function sampled on N pixels (1D) or (Nx,Ny) grid."""
        if self.dim == 1:
            n = np.arange(N); return np.real(sum(M[i] * np.exp(2j * np.pi * self.k[i, 0] * n / N) for i in range(self.F)))
        gx, gy = np.meshgrid(np.arange(N[0]) / N[0], np.arange(N[1]) / N[1], indexing='ij')
        return np.real(sum(M[i] * np.exp(2j * np.pi * (self.k[i, 0] * gx + self.k[i, 1] * gy)) for i in range(self.F)))

    def ratio(self, M, B, D):
        worst = 0.0
        for z in self.focus:
            I = self.image(M, z); worst = max(worst, I[D].max() / max(I[B].min(), 1e-300))
        return worst

# ----------------------------------------------------------------------------- SDP relaxation
def sdp_bound(model, B, D, real_mask=True, solver="CLARABEL", verbose=False, extra_cons=None):
    F = model.F
    X = cp.Variable((F, F), hermitian=True); t = cp.Variable()
    cons = [X >> 0]
    dual_refs = []   # (kind, z, ix, constraint)
    for z in model.focus:
        for ix in D:
            c = cp.real(cp.trace(model.Q(ix, z) @ X)) <= t; cons.append(c); dual_refs.append(('D', z, ix, c))
        for ix in B:
            c = cp.real(cp.trace(model.Q(ix, z) @ X)) >= 1; cons.append(c); dual_refs.append(('B', z, ix, c))
    if real_mask:
        Pm = np.eye(F)[model.flip]
        cons.append(X == cp.conj(Pm @ X @ Pm.T))  # real-mask (Hermitian-symmetric orders) structure
    if extra_cons: cons += extra_cons(X)
    prob = cp.Problem(cp.Minimize(t), cons)
    t0 = time.time(); prob.solve(solver=solver, verbose=verbose); dt = time.time() - t0
    Xv = 0.5 * (X.value + X.value.conj().T)
    ev = np.linalg.eigvalsh(Xv)[::-1]
    mu, lamb = {}, {}
    for kind, z, ix, c in dual_refs:
        v = float(np.ravel(c.dual_value)[0])
        (mu if kind == 'D' else lamb)[(z, ix)] = abs(v)
    return dict(t=float(t.value), X=Xv, eig=ev, status=prob.status, mu=mu, lam=lamb, time=dt,
                rank_eps=int(np.sum(ev > 1e-6 * ev[0])))

def verify_certificate(model, mu, lamb, tol=1e-7):
    """Independent check of the dual certificate: S = sum mu Q - sum lam Q must be PSD with sum mu = 1.
       Returns (certified_lower_bound, min_eig(S), sum_mu). Valid bound = sum(lam) if min_eig >= -tol."""
    F = model.F
    S = np.zeros((F, F), complex)
    smu = 0.0
    for (z, ix), v in mu.items(): S += v * model.Q(ix, z); smu += v
    # normalize mu to sum 1 (scale lam accordingly) -- certificate is homogeneous
    S /= smu
    slam = 0.0
    for (z, ix), v in lamb.items(): S -= (v / smu) * model.Q(ix, z); slam += v / smu
    mineig = np.linalg.eigvalsh(0.5 * (S + S.conj().T)).min()
    return slam, mineig, smu

def certified_bound_rigorous(model, mu, lamb):
    """Shift certificate to make it rigorously PSD: if min eig = -e < 0, then S + e*I >= 0, i.e.
       sum mu Q - (sum lam Q) + e I >= 0.  Since each Q_b has trace Tr(Q_b)=T_b, we instead use the weaker but valid
       rescaling: for any M, M^H S M >= -e |M|^2. With bright constraint M^H Q_b M >= 1 for all b in B and
       |M|^2 <= 1/(lambda_min(Q_b)) not available in general -> we report the solver bound minus e*kappa where
       kappa = |M|^2 upper bound from sum_b Q_b >= (1/|B|) sum_b Q_b ... (left as diagnostic). Returns raw slam, mineig."""
    return verify_certificate(model, mu, lamb)

# ----------------------------------------------------------------------------- rounding / local solvers
def round_from_X(model, X, B, D, nsamp=500, real_mask=True, rng=None):
    rng = rng or np.random.default_rng(0)
    w, V = np.linalg.eigh(X); w = np.clip(w, 0, None); Sq = V * np.sqrt(w)
    best = (np.inf, None)
    cands = [Sq[:, -1]]  # top eigenvector
    for _ in range(nsamp):
        z = (rng.standard_normal(model.F) + 1j * rng.standard_normal(model.F)) / np.sqrt(2)
        cands.append(Sq @ z)
    for M in cands:
        if real_mask: M = 0.5 * (M + np.conj(M[model.flip]))
        r = model.ratio(M, B, D)
        if r < best[0]: best = (r, M)
    return best

def local_refine(model, M0, B, D, real_mask=True, iters=300):
    """Local minimization of max_D I / min_B I via smooth log-sum-exp surrogate + L-BFGS on real parameters."""
    from scipy.optimize import minimize
    F = model.F
    def unpack(p):
        M = p[:F] + 1j * p[F:]
        if real_mask: M = 0.5 * (M + np.conj(M[model.flip]))
        return M
    beta = 40.0
    def obj(p):
        M = unpack(p); tot = 0.0
        for z in model.focus:
            fld = model.C[z] @ M; I = np.einsum('s,sx->x', model.w, np.abs(fld) ** 2)
            ID, IB = I[D], I[B]
            tot += (np.log(np.sum(np.exp(beta * (ID / ID.max() - 1)))) / beta + np.log(ID.max()))
            tot -= (-np.log(np.sum(np.exp(-beta * (IB / IB.min() - 1)))) / beta + np.log(IB.min()))
        return tot
    p0 = np.concatenate([M0.real, M0.imag])
    res = minimize(obj, p0, method='L-BFGS-B', options=dict(maxiter=iters))
    M = unpack(res.x)
    return model.ratio(M, B, D), M

def multistart_local(model, B, D, nstart=20, real_mask=True, rng=None, iters=300):
    rng = rng or np.random.default_rng(1)
    best = (np.inf, None)
    for _ in range(nstart):
        M0 = rng.standard_normal(model.F) + 1j * rng.standard_normal(model.F)
        r, M = local_refine(model, M0, B, D, real_mask, iters)
        if r < best[0]: best = (r, M)
    return best

# ----------------------------------------------------------------------------- pixel-domain ILT baseline (1D/2D)
def ilt_gradient(model, target_pix, N, B, D, steps=400, lr=0.05, a=10.0, thr=None, rng=None, pixel_bounds=(-1, 1), seed_from=None, init_target=None):
    """Classic ILT: parametrize mask m = tanh(theta) (PSM) or sigmoid (binary), resist = sigmoid(a (I - thr)),
       loss = |resist - target|^2; Adam. Returns ratio metric of the final mask and the mask."""
    rng = rng or np.random.default_rng(0)
    lo, hi = pixel_bounds
    if init_target is not None:   # standard ILT practice: initialize mask from the target layout
        m0 = lo + (hi - lo) * np.clip(init_target.astype(float), 0.02, 0.98)
        theta = np.arctanh(np.clip(2 * (m0 - lo) / (hi - lo) - 1, -0.96, 0.96))
    elif seed_from is None:
        theta = 0.5 * rng.standard_normal(N if np.isscalar(N) else tuple(N))
    else:
        theta = np.arctanh(np.clip((2 * (seed_from - lo) / (hi - lo) - 1), -0.999, 0.999))
    mom, vel = np.zeros_like(theta), np.zeros_like(theta)
    # image via orders: I(x) = sum_s w_s |C M|^2 with M = DFT(m) restricted to band -> gradient via adjoint
    Nn = theta.size
    def fwd(m):
        M = model.orders_from_pixels(m)
        fld = model.C[0.0] @ M
        I = np.einsum('s,sx->x', model.w, np.abs(fld) ** 2)
        return M, fld, I
    def dM_to_dm(gM, shape):
        # M_k = (1/Nn) sum_n m_n e^{-2pi i k n/Nn} -> dm_n = (1/Nn) Re sum_k gM_k^* ... (real m):  grad_m = Re(sum_k conj(dL/dM_k)*... )
        # L real; dL/dm_n = 2 Re( sum_k (dL/dM_k^*) * dM_k^*/dm_n ) handled as: g_n = Re( sum_k gM_k * e^{-2pi i k n /Nn} ) / Nn
        if model.dim == 1:
            n = np.arange(Nn); return np.real(np.array([np.sum(gM * np.exp(+2j * np.pi * model.k[:, 0] * nn / Nn)) for nn in n])) / Nn   # dM_k^*/dm_n = e^{+2 pi i k n/N}/N
        gx, gy = np.meshgrid(np.arange(shape[0]) / shape[0], np.arange(shape[1]) / shape[1], indexing='ij')
        out = np.zeros(shape)
        for i in range(model.F):
            out += np.real(gM[i] * np.exp(+2j * np.pi * (model.k[i, 0] * gx + model.k[i, 1] * gy)))
        return out / Nn
    tgt = target_pix.astype(float)
    hist = []
    for it in range(steps):
        m = lo + (hi - lo) * 0.5 * (1 + np.tanh(theta))
        M, fld, I = fwd(m)
        th = thr if thr is not None else 0.3          # intensity in open-frame units (mask=1 -> I=1); standard fixed threshold
        r = 1 / (1 + np.exp(-a * (I - th)))
        err = r - tgt.ravel()
        dL_dI = 2 * err * a * r * (1 - r)
        # I = sum_s w_s |fld_s|^2, fld_s = C_s M  -> dL/dM^* = sum_s w_s C_s^H (dL_dI * fld_s)
        gM = np.einsum('s,sxk,sx->k', model.w, model.C[0.0].conj(), dL_dI[None, :] * fld)
        gm = 2 * dM_to_dm(gM, theta.shape)
        gtheta = gm * (hi - lo) * 0.5 * (1 - np.tanh(theta) ** 2)
        mom = 0.9 * mom + 0.1 * gtheta; vel = 0.999 * vel + 0.001 * gtheta ** 2
        theta -= lr * mom / (np.sqrt(vel) + 1e-8)
        hist.append(np.mean(err ** 2))
    m = lo + (hi - lo) * 0.5 * (1 + np.tanh(theta))
    M = model.orders_from_pixels(m)
    return model.ratio(M, B, D), m, M, hist

# ----------------------------------------------------------------------------- pixel-domain lifted relaxation (for binary/DNN)
def sdp_bound_pixel(model, N, B, D, nonneg=False, solver="CLARABEL", verbose=False):
    """Lifted relaxation in pixel domain: X_pix = m m^T (real symmetric). Q_pix(x) = G_x^H Q_x G_x, with G = DFT band matrix.
       nonneg=True adds X >= 0 entrywise (valid for masks with m >= 0, e.g. binary/attenuated), i.e. DNN relaxation."""
    assert model.dim == 1
    n = np.arange(N)
    G = np.exp(-2j * np.pi * np.outer(model.k[:, 0], n) / N) / N      # M = G m
    X = cp.Variable((N, N), symmetric=True); t = cp.Variable()
    cons = [X >> 0]
    if nonneg: cons.append(X >= 0)
    for z in model.focus:
        for ix in D:
            Qp = np.real(G.conj().T @ model.Q(ix, z) @ G); Qp = 0.5 * (Qp + Qp.T)
            cons.append(cp.trace(Qp @ X) <= t)
        for ix in B:
            Qp = np.real(G.conj().T @ model.Q(ix, z) @ G); Qp = 0.5 * (Qp + Qp.T)
            cons.append(cp.trace(Qp @ X) >= 1)
    prob = cp.Problem(cp.Minimize(t), cons); prob.solve(solver=solver, verbose=verbose)
    ev = np.linalg.eigvalsh(X.value)[::-1]
    return dict(t=float(t.value), X=X.value, eig=ev, status=prob.status, rank_eps=int(np.sum(ev > 1e-6 * ev[0])))

# ----------------------------------------------------------------------------- helpers
def conventional_source_1d(sigma, NA, n=9):
    s = np.linspace(-sigma * NA, sigma * NA, n); return s, np.ones(n)
def dipole_source_1d(sigma_c, NA):
    return np.array([-sigma_c * NA, sigma_c * NA]), np.ones(2)
def disk_source_2d(sigma, NA, n=7):
    g = np.linspace(-sigma * NA, sigma * NA, n); pts = [(a, b) for a in g for b in g if np.hypot(a, b) <= sigma * NA + 1e-12]
    return np.array(pts), np.ones(len(pts))
def annular_source_2d(s_in, s_out, NA, n=9):
    g = np.linspace(-s_out * NA, s_out * NA, n); pts = [(a, b) for a in g for b in g if s_in * NA <= np.hypot(a, b) <= s_out * NA + 1e-12]
    return np.array(pts), np.ones(len(pts))
def quasar_source_2d(s_in, s_out, NA, n=11):
    g = np.linspace(-s_out * NA, s_out * NA, n)
    pts = [(a, b) for a in g for b in g if s_in * NA <= np.hypot(a, b) <= s_out * NA + 1e-12 and abs(abs(a) - abs(b)) < 0.35 * s_out * NA]
    return np.array(pts), np.ones(len(pts))

def sets_1d(model, target_fn, margin):
    """target_fn(x)->bool bright. Returns indices of bright (B) and dark (D) samples whose distance to the nearest
       target edge (periodic) is >= margin. Edges are located at midpoints between consecutive samples with different target."""
    x = model.x[:, 0]; n = len(x); L = model.L[0]; dx = L / n
    tgt = np.array([bool(target_fn(v)) for v in x])
    edges = np.array([x[i] + 0.5 * dx for i in range(n) if tgt[i] != tgt[(i + 1) % n]])
    if len(edges) == 0:
        d = np.full(n, np.inf)
    else:
        diff = (x[:, None] - edges[None, :] + L / 2) % L - L / 2
        d = np.min(np.abs(diff), axis=1)
    B = np.where(tgt & (d >= margin))[0]; D = np.where(~tgt & (d >= margin))[0]
    return B, D, tgt

# ----------------------------------------------------------------------------- rigorous certificate
def rigorize_certificate(model, mu, lamb, safety=1e-10):
    """Given approximate dual (mu, lam), find the largest c<=1 with  S_c = sum mu Q_d - c*sum lam Q_b  PSD
       (min eigenvalue >= safety*||S_c||), by bisection (min eig is concave in c; c=0 gives a PSD matrix).
       Returns rigorous lower bound c*sum(lam)/sum(mu), achieved min eigenvalue, and c."""
    F = model.F
    SD = np.zeros((F, F), complex); smu = 0.0
    for (z, ix), v in mu.items(): SD += v * model.Q(ix, z); smu += v
    SB = np.zeros((F, F), complex); slam = 0.0
    for (z, ix), v in lamb.items(): SB += v * model.Q(ix, z); slam += v
    SD /= smu; SB /= smu; slam /= smu
    def mineig(c):
        S = SD - c * SB; S = 0.5 * (S + S.conj().T)
        return np.linalg.eigvalsh(S).min() - safety * np.linalg.norm(S, 2)
    if mineig(1.0) >= 0: return slam, mineig(1.0), 1.0
    lo, hi = 0.0, 1.0
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if mineig(mid) >= 0: lo = mid
        else: hi = mid
    return lo * slam, mineig(lo), lo

# ----------------------------------------------------------------------------- real-subspace formulation (real masks)
class RealModel:
    """View of a Model restricted to real masks: orders M = U r with r in R^F, U unitary basis of the
       Hermitian-symmetric subspace {M : M_{-k} = conj(M_k)}. Then I(x) = r^T Qr_x r with Qr_x = Re(U^H Q_x U)."""
    def __init__(self, model):
        self.m = model; F = model.F; self.F = F; self.focus = model.focus; self.w = model.w; self.dim = model.dim
        U = np.zeros((F, F), complex); col = 0; seen = set()
        for i, kk in enumerate(model.k):
            j = model.flip[i]
            if i in seen: continue
            if j == i: U[i, col] = 1.0; col += 1; seen.add(i)
            else:
                U[i, col] = 1/np.sqrt(2); U[j, col] = 1/np.sqrt(2); col += 1
                U[i, col] = 1j/np.sqrt(2); U[j, col] = -1j/np.sqrt(2); col += 1
                seen.update([i, j])
        assert col == F
        self.U = U
        self.C = {z: model.C[z] @ U for z in model.focus}
        self._Q = {}
    def Q(self, ix, z=0.0):
        key = (ix, z)
        if key not in self._Q:
            Qr = np.real(self.U.conj().T @ self.m.Q(ix, z) @ self.U); self._Q[key] = 0.5*(Qr+Qr.T)
        return self._Q[key]
    def image(self, r, z=0.0):
        fld = self.C[z] @ r; return np.einsum('s,sx->x', self.w, np.abs(fld)**2)
    def ratio(self, r, B, D):
        worst = 0.0
        for z in self.focus:
            I = self.image(r, z); worst = max(worst, I[D].max()/I[B].min() if I[B].min() > 0 else np.inf)
        return worst
    def to_orders(self, r): return self.U @ r
    def from_orders(self, M): return np.real(self.U.conj().T @ M)

def sdp_bound_real(rm, B, D, solvers=("CLARABEL", "SCS"), verbose=False, extra_cons=None):
    F = rm.F
    Y = cp.Variable((F, F), symmetric=True); t = cp.Variable()
    cons = [Y >> 0]; refs = []
    for z in rm.focus:
        for ix in D:
            c = cp.trace(rm.Q(ix, z) @ Y) <= t; cons.append(c); refs.append(('D', z, ix, c))
        for ix in B:
            c = cp.trace(rm.Q(ix, z) @ Y) >= 1; cons.append(c); refs.append(('B', z, ix, c))
    if extra_cons: cons += extra_cons(Y)
    prob = cp.Problem(cp.Minimize(t), cons)
    t0 = time.time(); status = None; sv = None
    for sv in solvers:
        try:
            kw = dict(eps=1e-9, max_iters=200000) if sv == "SCS" else {}
            prob.solve(solver=sv, verbose=verbose, **kw); status = prob.status
            if status in ("optimal", "optimal_inaccurate"): break
        except Exception as e:
            status = f"fail:{sv}"
    dt = time.time()-t0
    if Y.value is None: raise RuntimeError("all solvers failed")
    Yv = 0.5*(Y.value+Y.value.T); ev = np.linalg.eigvalsh(Yv)[::-1]
    mu, lamb = {}, {}
    for kind, z, ix, c in refs:
        v = float(np.ravel(c.dual_value)[0]) if c.dual_value is not None else 0.0
        (mu if kind == 'D' else lamb)[(z, ix)] = abs(v)
    return dict(t=float(t.value), X=Yv, eig=ev, status=status, mu=mu, lam=lamb, time=dt, rank_eps=int(np.sum(ev > 1e-6*ev[0])), solver=sv)

def round_from_Y(rm, Y, B, D, nsamp=500, rng=None):
    rng = rng or np.random.default_rng(0)
    w, V = np.linalg.eigh(Y); w = np.clip(w, 0, None); Sq = V*np.sqrt(w)
    best = (np.inf, None)
    for r in [Sq[:, -1]] + [Sq @ rng.standard_normal(rm.F) for _ in range(nsamp)]:
        v = rm.ratio(r, B, D)
        if v < best[0]: best = (v, r)
    return best

def local_refine_real(rm, r0, B, D, iters=400, beta=40.0):
    """Minimize softmax_beta(log I_D) - softmin_beta(log I_B) over foci (smooth surrogate of log ratio), analytic gradient."""
    from scipy.optimize import minimize
    from scipy.special import logsumexp, softmax
    def fg(r):
        f = 0.0; g = np.zeros_like(r)
        for z in rm.focus:
            C = rm.C[z]; fld = C @ r
            I = np.einsum('s,sx->x', rm.w, np.abs(fld)**2) + 1e-300
            G = 2*np.einsum('s,sx,sxk->xk', rm.w, fld.conj(), C).real
            lI = np.log(I); dl = G / I[:, None]
            pD = softmax(beta*lI[D]); pB = softmax(-beta*lI[B])
            f += logsumexp(beta*lI[D])/beta + logsumexp(-beta*lI[B])/beta
            g += pD @ dl[D] - pB @ dl[B]
        return f, g
    res = minimize(fg, r0, jac=True, method='L-BFGS-B', options=dict(maxiter=iters))
    res2 = minimize(fg, res.x, jac=True, method='L-BFGS-B', options=dict(maxiter=iters))
    return rm.ratio(res2.x, B, D), res2.x

def multistart_local_real(rm, B, D, nstart=20, rng=None, iters=400):
    rng = rng or np.random.default_rng(1); best = (np.inf, None)
    for _ in range(nstart):
        v, r = local_refine_real(rm, rng.standard_normal(rm.F), B, D, iters)
        if v < best[0]: best = (v, r)
    return best

def exact_minimax_polish(rm, r0, B, D):
    """Polish by solving the nonconvex minimax as a smooth NLP via SLSQP: min t s.t. I_d <= t, I_b >= 1."""
    from scipy.optimize import minimize
    F = rm.F
    r0 = r0/np.sqrt(min(rm.image(r0, z)[B].min() for z in rm.focus))
    x0 = np.concatenate([r0, [rm.ratio(r0, B, D)]])
    cons = []
    for z in rm.focus:
        for d in D:
            Qd = rm.Q(d, z); cons.append(dict(type='ineq', fun=(lambda x, Qd=Qd: x[-1] - x[:F] @ Qd @ x[:F]), jac=(lambda x, Qd=Qd: np.concatenate([-2*Qd @ x[:F], [1.0]]))))
        for b in B:
            Qb = rm.Q(b, z); cons.append(dict(type='ineq', fun=(lambda x, Qb=Qb: x[:F] @ Qb @ x[:F] - 1), jac=(lambda x, Qb=Qb: np.concatenate([2*Qb @ x[:F], [0.0]]))))
    res = minimize(lambda x: x[-1], x0, jac=lambda x: np.concatenate([np.zeros(F), [1.0]]), constraints=cons, method='SLSQP', options=dict(maxiter=500, ftol=1e-12))
    r = res.x[:F]
    return rm.ratio(r, B, D), r

# ----------------------------------------------------------------------------- optimal-face exploration (rank reduction)
def face_min_rank(rm, B, D, tstar, ntries=6, rng=None, slack=1e-7, solvers=("CLARABEL","SCS")):
    """Explore the optimal face {Y feasible, t <= tstar(1+slack)} by minimizing random linear functionals Tr(G Y),
       G PSD random -> extreme points (Pataki) -> low rank. Returns the lowest-rank solution found and all ranks."""
    rng = rng or np.random.default_rng(0); F = rm.F
    best = None; ranks = []
    for i in range(ntries):
        Y = cp.Variable((F, F), symmetric=True)
        cons = [Y >> 0]
        for z in rm.focus:
            for ix in D: cons.append(cp.trace(rm.Q(ix, z) @ Y) <= tstar*(1+slack))
            for ix in B: cons.append(cp.trace(rm.Q(ix, z) @ Y) >= 1)
        if i == 0: G = np.eye(F)
        else:
            A = rng.standard_normal((F, F)); G = A @ A.T / F
        prob = cp.Problem(cp.Minimize(cp.trace(G @ Y)), cons)
        for sv in solvers:
            try:
                prob.solve(solver=sv, **(dict(eps=1e-9, max_iters=100000) if sv == "SCS" else {}))
                if prob.status in ("optimal", "optimal_inaccurate"): break
            except Exception: pass
        if Y.value is None: continue
        Yv = 0.5*(Y.value+Y.value.T); ev = np.linalg.eigvalsh(Yv)[::-1]
        rk = int(np.sum(ev > 1e-5*ev[0])); ranks.append(rk)
        if best is None or rk < best[0]: best = (rk, Yv, ev)
    return best, ranks

def strong_primal(rm, B, D, Ylist, nround=400, nstart=30, rng=None):
    """Stronger primal search: rounding from several Y's + SLSQP polish of top candidates + multistart."""
    rng = rng or np.random.default_rng(3); cands = []
    for Y in Ylist:
        w, V = np.linalg.eigh(Y); w = np.clip(w, 0, None); Sq = V*np.sqrt(w)
        for r in [Sq[:, -1]] + [Sq @ rng.standard_normal(rm.F) for _ in range(nround)]:
            cands.append((rm.ratio(r, B, D), r))
    cands.sort(key=lambda z: z[0]); top = cands[:8]
    out = []
    for v, r in top:
        v1, r1 = local_refine_real(rm, r, B, D)
        try: v2, r2 = exact_minimax_polish(rm, r1, B, D)
        except Exception: v2, r2 = v1, r1
        out.append(min((v1, r1), (v2, r2), key=lambda z: z[0]))
    vm, rmv = multistart_local_real(rm, B, D, nstart=nstart, rng=rng)
    try: vp, rp = exact_minimax_polish(rm, rmv, B, D); out.append(min((vm, rmv), (vp, rp), key=lambda z: z[0]))
    except Exception: out.append((vm, rmv))
    return min(out, key=lambda z: z[0])

# ----------------------------------------------------------------------------- source-mask (entitlement) relaxation
def sdp_bound_smo(rm, B, D, solvers=("CLARABEL","SCS")):
    """Relaxation of joint source+mask optimization: one PSD block Y_s per source point (its own mask ensemble):
       min t s.t. sum_s Tr(Q_x^{(s)} Y_s) <= t (D), >= 1 (B), Y_s PSD. Equivalent to unlimited multi-exposure with
       coherent point-source exposures. Lower-bounds the ratio achievable by ANY source (supported on the point set) + mask."""
    F = rm.F; S = len(rm.w)
    Ys = [cp.Variable((F, F), symmetric=True) for _ in range(S)]; t = cp.Variable()
    cons = [Y >> 0 for Y in Ys]
    def Qs(ix, z, s):
        c = rm.C[z][s, ix, :]
        Q = np.real(np.outer(c.conj(), c)); return 0.5*(Q+Q.T)
    for z in rm.focus:
        for ix in D: cons.append(sum(cp.trace(Qs(ix, z, s) @ Ys[s]) for s in range(S)) <= t)
        for ix in B: cons.append(sum(cp.trace(Qs(ix, z, s) @ Ys[s]) for s in range(S)) >= 1)
    prob = cp.Problem(cp.Minimize(t), cons)
    for sv in solvers:
        try:
            prob.solve(solver=sv, **(dict(eps=1e-9, max_iters=200000) if sv == "SCS" else {}))
            if prob.status in ("optimal", "optimal_inaccurate"): break
        except Exception: pass
    Yv = [0.5*(Y.value+Y.value.T) for Y in Ys]
    weights = np.array([np.trace(Y) for Y in Yv])
    return dict(t=float(t.value), Ys=Yv, status=prob.status, src_weights=weights/weights.sum())

# ----------------------------------------------------------------------------- Lasserre level-2 bound for the SINGLE-mask problem
def lasserre2_single_mask(rm, B, D, solvers=("CLARABEL","SCS"), tmax=None, verbose=False):
    """Level-2 moment relaxation of  min t s.t. t - r^T Q_d r >= 0 (d in D), r^T Q_b r - 1 >= 0 (b in B), (r,t) in R^{F+1}.
       Returns a lower bound on the single-mask optimum rho*(M_R) that is at least as strong as the Shor/SDP bound."""
    import itertools
    F = rm.F; n = F + 1
    def monos(deg):
        out = []
        for d in range(deg+1):
            for c in itertools.combinations_with_replacement(range(n), d):
                a = [0]*n
                for i in c: a[i] += 1
                out.append(tuple(a))
        return out
    m2 = monos(2); m1 = monos(1); m4 = monos(4)
    yidx = {a: i for i, a in enumerate(m4)}
    y = cp.Variable(len(m4))
    def add(a, b): return tuple(x+z for x, z in zip(a, b))
    def Mmat(basis, poly):
        rows = []
        for a in basis:
            row = []
            for b in basis:
                expr = 0
                for g, cg in poly.items(): expr = expr + cg*y[yidx[add(add(a, b), g)]]
                row.append(expr)
            rows.append(row)
        return cp.bmat(rows)
    one = tuple([0]*n)
    cons = [y[yidx[one]] == 1, Mmat(m2, {one: 1.0}) >> 0]
    def quad_poly(Q, const=0.0, tcoef=0.0):
        p = {}
        if const != 0: p[one] = const
        if tcoef != 0:
            a = [0]*n; a[F] = 1; p[tuple(a)] = tcoef
        for i in range(F):
            for j in range(i, F):
                c = Q[i, j] if i == j else 2*Q[i, j]
                if abs(c) > 1e-14:
                    a = [0]*n; a[i] += 1; a[j] += 1; p[tuple(a)] = p.get(tuple(a), 0) + c
        return p
    for z in rm.focus:
        for d in D: cons.append(Mmat(m1, quad_poly(-rm.Q(d, z), tcoef=1.0)) >> 0)
        for b in B: cons.append(Mmat(m1, quad_poly(rm.Q(b, z), const=-1.0)) >> 0)
    if tmax is not None:
        cons.append(Mmat(m1, {one: tmax, tuple([0]*F+[1]): -1.0}) >> 0)
    tvar = y[yidx[tuple([0]*F+[1])]]
    prob = cp.Problem(cp.Minimize(tvar), cons)
    for sv in solvers:
        try:
            prob.solve(solver=sv, verbose=verbose, **(dict(eps=1e-9, max_iters=200000) if sv == "SCS" else {}))
            if prob.status in ("optimal", "optimal_inaccurate"): break
        except Exception as e:
            if verbose: print("solver fail", sv, e)
    return dict(t=float(tvar.value) if tvar.value is not None else None, status=prob.status, size_M2=len(m2), n_y=len(m4))


# ----------------------------------------------------------------------------- exact edge-offset sampling (1D)
def samples_1d(L, segs, margin, spacing):
    """Deterministic sample set for a 1D periodic target given by bright segments [(a,b),...] in [0,L):
       for every edge e put samples at e-margin and e+margin (exactly), and fill the interiors with a grid of the given
       spacing (only points at distance >= margin from every edge). Returns (x, B, D, tgt)."""
    segs = [(float(a), float(b)) for a, b in segs]
    edges = sorted(set([a for a, b in segs] + [b for a, b in segs]))
    pts = set()
    for e in edges:
        pts.add((e - margin) % L); pts.add((e + margin) % L)
    for v in np.arange(0, L, spacing): pts.add(float(v))
    x = np.array(sorted(pts))
    tf = lambda v: any(a <= (v % L) < b for a, b in segs)
    tgt = np.array([tf(v) for v in x])
    E = np.array(edges)
    d = np.min(np.abs(((x[:, None] - E[None, :] + L/2) % L) - L/2), axis=1)
    keep = d >= margin - 1e-9
    B = np.where(tgt & keep)[0]; D = np.where(~tgt & keep)[0]
    return x, B, D, tgt

def build_1d(L, lam, NA, src, segs, margin, spacing=4.0, focus=(0.0,)):
    s, w = src; x, B, D, tgt = samples_1d(L, segs, margin, spacing)
    md = Model(L, lam, NA, s, w, nx=len(x), focus=focus, xs=x)
    return md, B, D, tgt


# ----------------------------------------------------------------------------- exact edge-offset sampling (2D, periodic)
def samples_2d(L, target_fn, margin, spacing, fine=1.0, band=0.6):
    """2D periodic target given by target_fn(x,y)->bool. Signed distance computed on a fine grid (periodic via tiling);
       sample points: (i) all fine-grid points whose distance to the target boundary lies in [margin, margin+band)
       (offset contours, both sides), (ii) a coarse interior grid with distance >= margin. Returns (x, B, D, tgt, sd_grid)."""
    from scipy.ndimage import distance_transform_edt
    Lx, Ly = float(L[0]), float(L[1]); nx, ny = int(round(Lx/fine)), int(round(Ly/fine))
    gx = (np.arange(nx) + 0.5) * Lx / nx; gy = (np.arange(ny) + 0.5) * Ly / ny
    T = np.array([[bool(target_fn(a, b)) for b in gy] for a in gx])
    T3 = np.tile(T, (3, 3))
    d_out = distance_transform_edt(~T3, sampling=(Lx/nx, Ly/ny))[nx:2*nx, ny:2*ny]   # distance to target for outside pts
    d_in = distance_transform_edt(T3, sampling=(Lx/nx, Ly/ny))[nx:2*nx, ny:2*ny]     # distance to background for inside pts
    sd = np.where(T, d_in, -d_out)   # positive inside
    pts = []; tg = []
    for i in range(nx):
        for j in range(ny):
            a = abs(sd[i, j])
            if margin <= a < margin + band: pts.append((gx[i], gy[j])); tg.append(T[i, j])
    for a in np.arange(0, Lx, spacing):
        for b in np.arange(0, Ly, spacing):
            i = min(int(a / Lx * nx), nx-1); j = min(int(b / Ly * ny), ny-1)
            if abs(sd[i, j]) >= margin + band: pts.append((a, b)); tg.append(T[i, j])
    x = np.array(pts); tgt = np.array(tg)
    return x, np.where(tgt)[0], np.where(~tgt)[0], tgt, sd

def build_2d(L, lam, NA, src, target_fn, margin, spacing=10.0, fine=1.0, focus=(0.0,)):
    s, w = src; x, B, D, tgt, sd = samples_2d(L, target_fn, margin, spacing, fine=fine)
    md = Model(L, lam, NA, s, w, nx=None, focus=focus, dim=2, xs=x)
    return md, B, D, tgt
