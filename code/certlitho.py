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
        """Common-dose dark/bright ratio over the whole focus set: max_{(d,z)} I / min_{(b,z)} I."""
        dmax = max(self.image(M, z)[D].max() for z in self.focus); bmin = min(self.image(M, z)[B].min() for z in self.focus)
        return dmax / bmin if bmin > 0 else np.inf

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
            I = model.image(M, z)
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
        fld = model.fields(M, 0.0)
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
        gM = np.einsum('sxk,sx->k', model.C[0.0].conj(), (model.w[:, None] * dL_dI[None, :]) * fld)
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
        dmax = max(self.image(r, z)[D].max() for z in self.focus); bmin = min(self.image(r, z)[B].min() for z in self.focus)
        return dmax / bmin if bmin > 0 else np.inf
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
        lD, dD, lB, dB = [], [], [], []
        for z in rm.focus:
            I, G = rm.image_grad(r, z); I = I + 1e-300
            lI = np.log(I); dl = G / I[:, None]
            lD.append(lI[D]); dD.append(dl[D]); lB.append(lI[B]); dB.append(dl[B])
        lD = np.concatenate(lD); dD = np.concatenate(dD); lB = np.concatenate(lB); dB = np.concatenate(dB)
        pD = softmax(beta*lD); pB = softmax(-beta*lB)
        f = logsumexp(beta*lD)/beta + logsumexp(-beta*lB)/beta     # = log(max dark) - log(min bright), smoothed, common dose
        g = pD @ dD - pB @ dB
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
    if getattr(rm.m, "compressed", False): raise ValueError("sdp_bound_smo needs per-source channels: build the model with compress=False")
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
def boundary_corners(segs, L):
    """Vertices of the union boundary: segment endpoints at which the boundary turns (an endpoint shared by two collinear
       kept segments, e.g. an over-split straight edge, is not a corner). Returns an (m,2) array in [0,Lx)x[0,Ly)."""
    Lx, Ly = L; ends = {}
    for o, c, lo, hi, sgn in segs:
        for t in (lo, hi):
            p = ((c % Lx, t % Ly) if o == 'v' else (t % Lx, c % Ly)); key = (round(p[0], 6), round(p[1], 6))
            ends.setdefault(key, []).append((o, round(c % (Lx if o == 'v' else Ly), 6), sgn))
    out = []
    for key, inc in ends.items():
        collinear = any(sum(1 for q in inc if q == a) >= 2 for a in set(inc))
        if not collinear: out.append(key)
    return np.array(out, float).reshape(-1, 2)

def samples_2d(L, target_fn, margin, spacing, fine=1.0, band=0.6, corner_exclusion=0.0):
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

# ----------------------------------------------------------------------------- rigorous certificate (null-space aware; overrides)
def rigorize_certificate(model, mu, lamb, safety=1e-10, null_tol=1e-9):
    """Largest c<=1 with S_c = SD - c*SB PSD (SD = sum mu Q_d / sum mu, SB = sum lam Q_b / sum mu), verified by an
       eigenvalue test with a relative safety margin on the orthogonal complement of the near-null space N of SD
       (eigenvalues < null_tol*||SD||). On N both SD and SB must vanish (reported residual resid_null = ||P_N SB P_N||/||SB||;
       a valid certificate has resid_null ~ 0 up to solver tolerance). Returns (bound, min_eig_on_complement, c)."""
    F = model.F
    SD = np.zeros((F, F), complex); smu = 0.0
    for (z, ix), v in mu.items(): SD += v * model.Q(ix, z); smu += v
    SB = np.zeros((F, F), complex); slam = 0.0
    for (z, ix), v in lamb.items(): SB += v * model.Q(ix, z); slam += v
    SD /= smu; SB /= smu; slam /= smu
    SD = 0.5*(SD+SD.conj().T); SB = 0.5*(SB+SB.conj().T)
    w, V = np.linalg.eigh(SD); nD = np.linalg.norm(SD, 2)
    keep = w > null_tol*nD; Vk = V[:, keep]; Vn = V[:, ~keep]
    resid_null = (np.linalg.norm(Vn.conj().T @ SB @ Vn, 2)/max(np.linalg.norm(SB, 2), 1e-300)) if Vn.shape[1] else 0.0
    rigorize_certificate.last_resid_null = float(resid_null); rigorize_certificate.last_null_dim = int(Vn.shape[1])
    def mineig(c):
        S = Vk.conj().T @ (SD - c*SB) @ Vk; S = 0.5*(S+S.conj().T)
        return np.linalg.eigvalsh(S).min() - safety*np.linalg.norm(S, 2)
    if mineig(1.0) >= 0: return slam, mineig(1.0), 1.0
    lo, hi = 0.0, 1.0
    for _ in range(60):
        mid = 0.5*(lo+hi)
        if mineig(mid) >= 0: lo = mid
        else: hi = mid
    return lo*slam, mineig(lo), lo

# ----------------------------------------------------------------------------- vectorized SLSQP minimax polish (overrides)
def exact_minimax_polish(rm, r0, B, D, maxiter=300):
    """min t s.t. I_d(r) <= t (D x Z), I_b(r) >= 1 (B x Z): vectorized constraints/Jacobians for SLSQP."""
    from scipy.optimize import minimize
    F = rm.F
    r0 = r0/np.sqrt(min(rm.image(r0, z)[B].min() for z in rm.focus))
    x0 = np.concatenate([r0, [rm.ratio(r0, B, D)]])
    Cs = [rm.C[z] for z in rm.focus]
    def IG(r):
        Is, Gs = [], []
        for z in rm.focus:
            I, G = rm.image_grad(r, z); Is.append(I); Gs.append(G)
        return Is, Gs
    def cfun(x):
        Is, _ = IG(x[:F]); return np.concatenate([np.concatenate([x[-1] - I[D] for I in Is]), np.concatenate([I[B] - 1 for I in Is])])
    def cjac(x):
        _, Gs = IG(x[:F]); nd = sum(len(D) for _ in Gs); nb = sum(len(B) for _ in Gs)
        J = np.zeros((nd+nb, F+1)); J[:nd, :F] = -np.concatenate([G[D] for G in Gs]); J[:nd, F] = 1.0; J[nd:, :F] = np.concatenate([G[B] for G in Gs])
        return J
    res = minimize(lambda x: x[-1], x0, jac=lambda x: np.concatenate([np.zeros(F), [1.0]]), constraints=[dict(type='ineq', fun=cfun, jac=cjac)], method='SLSQP', options=dict(maxiter=maxiter, ftol=1e-12))
    r = res.x[:F]
    return rm.ratio(r, B, D), r

# ----------------------------------------------------------------------------- physics extensions: aberrations, vector imaging, pixel-integrated masks
def zernike_phase(g, lam, NA, aberr):
    """Pupil phase (radians) from fringe-Zernike coefficients in waves: keys 'Z4' defocus, 'Z5','Z6' astig, 'Z7','Z8' coma,
       'Z9' spherical. g: (F,dim) spatial frequencies; rho = |g| lam/NA, phi = atan2."""
    if not aberr: return np.zeros(len(g))
    gx = g[:, 0]; gy = g[:, 1] if g.shape[1] > 1 else np.zeros(len(g))
    rho = np.hypot(gx, gy)*lam/NA; phi = np.arctan2(gy, gx); W = np.zeros(len(g))
    Z = {'Z4': 2*rho**2 - 1, 'Z5': rho**2*np.cos(2*phi), 'Z6': rho**2*np.sin(2*phi), 'Z7': (3*rho**3 - 2*rho)*np.cos(phi),
         'Z8': (3*rho**3 - 2*rho)*np.sin(phi), 'Z9': 6*rho**4 - 6*rho**2 + 1}
    for k, a in aberr.items(): W += a*Z[k]
    return 2*np.pi*W

class ModelX(Model):
    """Model with optional pupil aberrations (fringe Zernikes, waves) and vector (high-NA) imaging.
       vector=True: image-side plane wave with direction cosines (alpha,beta,gamma) = (lam g / n_img, sqrt(1-...)) carries the
       mask-side polarization p rotated by the minimal rotation taking z to the propagation direction (Rodrigues):
         R p = p - (p.k_t) k_t/(1+gamma) ... explicitly for p=x: (1 - a^2/(1+g), -a b/(1+g), -a); p=y: (-a b/(1+g), 1 - b^2/(1+g), -b).
       Intensity = |E_x|^2+|E_y|^2+|E_z|^2 summed over source points and over the polarization ensemble
       (pol in {'x','y','unpol','TE','TM'}; TE/TM are azimuthal/radial w.r.t. each source point). Optional obliquity
       amplitude factor 1/sqrt(gamma) (aplanatic, 4x reduction). Still quadratic in M: channels = (source, pol, component)."""
    def __init__(self, L, lam, NA, src_pts, src_w, nx, focus=(0.0,), dim=1, xs=None, aberr=None, vector=False, pol='unpol', n_img=1.0, obliquity=True):
        self.aberr = aberr or {}; self.vector = vector; self.pol = pol; self.n_img = n_img; self.obliquity = obliquity
        super().__init__(L, lam, NA, src_pts, src_w, nx, focus, dim, xs)
        # rebuild C with aberration / vector channels
        lamv = lam
        Cs, ws = {z: [] for z in self.focus}, []
        for s, wsrc in zip(self.src, self.w):
            g = self.f + s/lamv; gn = np.linalg.norm(g, axis=1); inside = (gn*lamv <= NA + 1e-12)
            gl = np.clip(gn*lamv, 0, 1); 
            if self.vector:
                a = g[:, 0]*lamv/n_img; b = (g[:, 1] if self.dim == 2 else np.zeros(self.F))*lamv/n_img
                gam = np.sqrt(np.clip(1 - a**2 - b**2, 1e-9, 1))
            E = np.exp(2j*np.pi*(self.x @ self.f.T))
            ab = zernike_phase(g, lamv, NA, self.aberr)
            pols = []
            if not self.vector: pols = [None]
            elif pol == 'x': pols = [np.array([1.0, 0.0])]
            elif pol == 'y': pols = [np.array([0.0, 1.0])]
            elif pol == 'unpol': pols = [np.array([1.0, 0.0]), np.array([0.0, 1.0])]
            elif pol in ('TE', 'TM'):
                sv = np.array([s[0], s[1] if len(s) > 1 else 0.0]); ns = np.linalg.norm(sv)
                rad = sv/ns if ns > 1e-12 else np.array([1.0, 0.0]); azi = np.array([-rad[1], rad[0]])
                pols = [azi if pol == 'TE' else rad]
            for z in self.focus:
                ph = 2*np.pi/lamv*z*(np.sqrt(1 - gl**2) - 1) + ab
                P = inside*np.exp(1j*ph)
                if not self.vector:
                    Cs[z].append(E*P[None, :])
                else:
                    ob = (1/np.sqrt(gam)) if obliquity else np.ones(self.F)
                    for p in pols:
                        ex = p[0]*(1 - a**2/(1+gam)) + p[1]*(-a*b/(1+gam)); ey = p[0]*(-a*b/(1+gam)) + p[1]*(1 - b**2/(1+gam)); ez = -p[0]*a - p[1]*b
                        for comp in (ex, ey, ez): Cs[z].append(E*(P*ob*comp)[None, :])
            nch = 1 if not self.vector else 3*len(pols)
            ws += [wsrc/ (len(pols) if self.vector else 1)]*nch
        self.C = {z: np.stack(Cs[z]) for z in self.focus}; self.w = np.array(ws)
    def Q(self, ix, z=0.0):
        Cz = self.C[z][:, ix, :]; return np.einsum('s,sa,sb->ab', self.w, Cz.conj(), Cz)

def build_1d_x(L, lam, NA, src, segs, margin, spacing=4.0, focus=(0.0,), **kw):
    s, w = src; x, B, D, tgt = samples_1d(L, segs, margin, spacing)
    return ModelX(L, lam, NA, s, w, nx=len(x), focus=focus, xs=x, **kw), B, D, tgt
def build_2d_x(L, lam, NA, src, target_fn, margin, spacing=10.0, fine=1.0, focus=(0.0,), **kw):
    s, w = src; x, B, D, tgt, sd = samples_2d(L, target_fn, margin, spacing, fine=fine)
    return ModelX(L, lam, NA, s, w, nx=None, focus=focus, dim=2, xs=x, **kw), B, D, tgt

def sdp_bound_pixel_int(model, N, B, D, nonneg=False, box=None, solver="CLARABEL"):
    """Pixel-domain lifted relaxation with PIXEL-INTEGRATED transform (mask piecewise constant on N pixels of width L/N):
       M_k = sinc(k/N) * (1/N) sum_n m_n e^{-2 pi i k n /N}. nonneg -> DNN. Returns bound for masks of that pixelation."""
    assert model.dim == 1
    n = np.arange(N); k = model.k[:, 0]
    G = (np.sinc(k/N)[:, None]) * np.exp(-2j*np.pi*np.outer(k, n + 0.5)/N)/N   # pixel n covers [nL/N,(n+1)L/N): centre phase e^{-i pi k/N}
    X = cp.Variable((N, N), symmetric=True); t = cp.Variable(); cons = [X >> 0]
    if nonneg: cons.append(X >= 0)
    for z in model.focus:
        for ix in D:
            Qp = np.real(G.conj().T @ model.Q(ix, z) @ G); cons.append(cp.trace(0.5*(Qp+Qp.T) @ X) <= t)
        for ix in B:
            Qp = np.real(G.conj().T @ model.Q(ix, z) @ G); cons.append(cp.trace(0.5*(Qp+Qp.T) @ X) >= 1)
    prob = cp.Problem(cp.Minimize(t), cons)
    for sv in (solver, "SCS"):
        try:
            prob.solve(solver=sv, **(dict(eps=1e-9, max_iters=200000) if sv == "SCS" else {}))
            if prob.status in ("optimal", "optimal_inaccurate"): break
        except Exception: pass
    ev = np.linalg.eigvalsh(X.value)[::-1]
    return dict(t=float(t.value), status=prob.status, rank_eps=int(np.sum(ev > 1e-6*ev[0])))

# ----------------------------------------------------------------------------- rigorous certificate, full-space version (overrides; session 2 fix)
def rigorize_certificate(model, mu, lamb, tol_rel=1e-11, safety=None):
    """Largest c<=1 such that S_c = SD - c*SB is PSD on the FULL space up to a relative eigenvalue tolerance:
       lambda_min(S_c) >= -tol_rel*||S_c||_2  (floating-point images of exact quantities carry ~1e-16 relative error).
       No null-space projection: if SD has a null direction on which SB > 0, bisection drives c -> 0 (true bound 0).
       Returns (bound = c*sum lam / sum mu, lambda_min(S_c)/||S_c||, c). Diagnostics stored on the function object."""
    F = model.F
    SD = np.zeros((F, F), complex); smu = 0.0
    for (z, ix), v in mu.items(): SD += v * model.Q(ix, z); smu += v
    SB = np.zeros((F, F), complex); slam = 0.0
    for (z, ix), v in lamb.items(): SB += v * model.Q(ix, z); slam += v
    SD /= smu; SB /= smu; slam /= smu
    SD = 0.5*(SD+SD.conj().T); SB = 0.5*(SB+SB.conj().T)
    def rel_mineig(c):
        S = SD - c*SB; n = np.linalg.norm(S, 2)
        return (np.linalg.eigvalsh(S).min()/n) if n > 0 else 0.0
    rigorize_certificate.last_null_dim = 0; rigorize_certificate.last_resid_null = 0.0
    if rel_mineig(1.0) >= -tol_rel: return slam, rel_mineig(1.0), 1.0
    lo, hi = 0.0, 1.0
    for _ in range(60):
        mid = 0.5*(lo+hi)
        if rel_mineig(mid) >= -tol_rel: lo = mid
        else: hi = mid
    return lo*slam, rel_mineig(lo), lo

# ----------------------------------------------------------------------------- rigorous certificate v3 (overrides): strict margin on the complement of the COMMON null space
def rigorize_certificate(model, mu, lamb, margin=1e-10, null_tol=1e-12):
    """Largest c<=1 with S_c = SD - c*SB strictly PSD (lambda_min >= margin*||S_c||) on the orthogonal complement of the
       common null space N = {v : SD v = 0 and SB v = 0} (directions that contribute to no weighted intensity; on N every
       quadratic form in the certificate vanishes identically, so the weak-duality chain is unaffected). N is taken as the
       span of eigenvectors of SD with eigenvalue <= null_tol*||SD|| on which SB also vanishes (||P_N SB P_N|| <= null_tol*||SB||);
       if SB does NOT vanish on the near-null space of SD, no projection is applied and bisection lowers c (to 0 if needed).
       A strictly positive margin makes an exact Cholesky factorisation of the projected matrix possible (exp15).
       Returns (bound = c*sum lam/sum mu, relative min eigenvalue achieved, c); diagnostics on the function object."""
    F = model.F
    SD = np.zeros((F, F), complex); smu = 0.0
    for (z, ix), v in mu.items(): SD += v * model.Q(ix, z); smu += v
    SB = np.zeros((F, F), complex); slam = 0.0
    for (z, ix), v in lamb.items(): SB += v * model.Q(ix, z); slam += v
    SD /= smu; SB /= smu; slam /= smu
    SD = 0.5*(SD+SD.conj().T); SB = 0.5*(SB+SB.conj().T)
    w, V = np.linalg.eigh(SD); nD = np.linalg.norm(SD, 2); nB = max(np.linalg.norm(SB, 2), 1e-300)
    nullmask = w <= null_tol*nD; Vn = V[:, nullmask]
    resid = (np.linalg.norm(Vn.conj().T @ SB @ Vn, 2)/nB) if Vn.shape[1] else 0.0
    common = Vn.shape[1] > 0 and resid <= null_tol
    Vk = V[:, ~nullmask] if common else np.eye(F)
    rigorize_certificate.last_null_dim = int(Vn.shape[1]) if common else 0; rigorize_certificate.last_resid_null = float(resid)
    def rel_mineig(c):
        S = Vk.conj().T @ (SD - c*SB) @ Vk; S = 0.5*(S+S.conj().T); n = np.linalg.norm(S, 2)
        return (np.linalg.eigvalsh(S).min()/n) if n > 0 else 0.0
    if rel_mineig(1.0) >= margin: return slam, rel_mineig(1.0), 1.0
    lo, hi = 0.0, 1.0
    for _ in range(60):
        mid = 0.5*(lo+hi)
        if rel_mineig(mid) >= margin: lo = mid
        else: hi = mid
    rigorize_certificate.last_S_proj = (lambda c: Vk.conj().T @ (SD - c*SB) @ Vk)
    return lo*slam, rel_mineig(lo), lo

def certificate_matrix_projected(model, mu, lamb, c):
    """Return the projected certificate matrix S_c restricted to the complement of the common null space (for exact checks)."""
    F = model.F
    SD = np.zeros((F, F), complex); smu = 0.0
    for (z, ix), v in mu.items(): SD += v * model.Q(ix, z); smu += v
    SB = np.zeros((F, F), complex)
    for (z, ix), v in lamb.items(): SB += v * model.Q(ix, z)
    SD /= smu; SB /= smu; SD = 0.5*(SD+SD.conj().T); SB = 0.5*(SB+SB.conj().T)
    w, V = np.linalg.eigh(SD); nD = np.linalg.norm(SD, 2); nB = max(np.linalg.norm(SB, 2), 1e-300)
    nullmask = w <= 1e-12*nD; Vn = V[:, nullmask]
    common = Vn.shape[1] > 0 and (np.linalg.norm(Vn.conj().T @ SB @ Vn, 2)/nB) <= 1e-12
    Vk = V[:, ~nullmask] if common else np.eye(F)
    S = Vk.conj().T @ (SD - c*SB) @ Vk; return 0.5*(S+S.conj().T), int(Vn.shape[1]) if common else 0

# =============================================================================== Session-3 physics corrections (reviewer F2/F3/F8)
# 1) Source points are ALWAYS 2-D direction-cosine vectors (s_x, s_y); for 1-D (y-invariant) masks the y-component enters
#    the pupil cutoff and the defocus phase exactly as for 2-D masks.  2) Defocus phase is evaluated in the image-side
#    medium of index n_img:  phi = 2 pi n_img z / lambda * ( sqrt(1 - (|g| lambda / n_img)^2) - 1 ).
# 3) Area-weighted quadratures for disk / annular / finite-pole dipole / quasar illuminators.

def compress_channels(model, tol=1e-14, E=None):
    """Exact coherent-mode (TCC eigen) compression of the source sum. For each focus z the channel vectors are E_x * P_s
       (E_x Fourier phases, P_s the pupil/defocus vector of source point s), so
       I(x) = sum_s w_s |(E_x*M)^T P_s|^2 = u^T T_z conj(u),  T_z = sum_s w_s P_s P_s^H  (F x F, PSD, rank <= F).
       With T_z = sum_j lam_j v_j v_j^H the image is sum_j lam_j |(E_x*M)^T v_j|^2: at most F channels replace S source
       points, exactly (up to rounding). Uses model.P[z] (S, F) when present (never forms per-source channel arrays),
       else the per-source channels model.C[z]. Keeps model.w_src (needed by the per-source SMO relaxation)."""
    if getattr(model, 'compressed', False): return model
    if E is None: E = np.exp(2j * np.pi * (model.x @ model.f.T))
    model.w_src = model.w; newC = {}; J = 0
    for z in model.focus:
        P = model.P[z] if getattr(model, 'P', None) and z in model.P else model.C[z][:, 0, :] / E[0][None, :]
        T = (P.T * model.w) @ P.conj(); T = 0.5 * (T + T.conj().T)
        lam_, V = np.linalg.eigh(T); keep = lam_ > tol * max(lam_.max(), 1e-300)
        modes = (V[:, keep] * np.sqrt(lam_[keep])).T
        newC[z] = E[None, :, :] * modes[:, None, :]; J = max(J, modes.shape[0])
    for z in model.focus:
        if newC[z].shape[0] < J: newC[z] = np.concatenate([newC[z], np.zeros((J - newC[z].shape[0],) + newC[z].shape[1:], complex)])
    if hasattr(model, 'C') and model.C: model.C_src = model.C
    model.C = newC; model.w = np.ones(J); model.compressed = True; model.n_channels = J
    return model

class Model(Model):
    def __init__(self, L, lam, NA, src_pts, src_w, nx, focus=(0.0,), dim=1, xs=None, n_img=1.0, aberr=None, compress=True):
        self.dim = dim; self.n_img = n_img; self.aberr = aberr or {}
        self.L = np.atleast_1d(np.asarray(L, float)) if dim == 2 else np.array([float(L)])
        self.lam, self.NA = lam, NA
        src = np.asarray(src_pts, float)
        if src.ndim == 1: src = src.reshape(-1, 1)
        if src.shape[1] == 1: src = np.hstack([src, np.zeros((src.shape[0], 1))])     # legacy 1-D slit points
        self.src = src; self.w = np.asarray(src_w, float) / np.sum(src_w); self.focus = tuple(focus)
        if NA > n_img + 1e-12 and any(abs(z) > 0 for z in focus): raise ValueError(f"NA={NA} exceeds the image-side index n_img={n_img}; defocus needs the immersion index (1.44 for water)")
        smax = np.max(np.linalg.norm(self.src, axis=1)); R = (NA + smax) / lam
        if dim == 1:
            K = int(np.floor(R * self.L[0] + 1e-9)); ks = [(k,) for k in range(-K, K + 1)]
        else:
            Kx = int(np.floor(R * self.L[0] + 1e-9)); Ky = int(np.floor(R * self.L[1] + 1e-9))
            ks = [(kx, ky) for kx in range(-Kx, Kx + 1) for ky in range(-Ky, Ky + 1) if np.hypot(kx / self.L[0], ky / self.L[1]) <= R + 1e-12]
        # drop dead orders: k that no source point ever brings inside the pupil (kept if -k is live, so K = -K)
        kk = np.array(ks, float); ff = kk / self.L; ff2 = ff if dim == 2 else np.hstack([ff, np.zeros((len(ff), 1))])
        live = np.zeros(len(ks), bool)
        for s_ in self.src: live |= (np.linalg.norm(ff2 + s_ / lam, axis=1) * lam <= NA + 1e-12)
        live_set = {tuple(k_) for k_, l in zip(ks, live) if l}
        ks = [k_ for k_ in ks if tuple(k_) in live_set or tuple(-np.array(k_)) in live_set]
        self.n_dead_dropped = int(len(kk) - len(ks))
        self.k = np.array(ks); self.F = len(self.k); self.f = self.k / self.L
        self.f2 = self.f if dim == 2 else np.hstack([self.f, np.zeros((self.F, 1))])
        self.kindex = {tuple(kk): i for i, kk in enumerate(self.k)}
        if dim == 1: self.x = (np.arange(nx) * self.L[0] / nx).reshape(-1, 1) if xs is None else np.asarray(xs, float).reshape(-1, 1)
        elif xs is not None: self.x = np.asarray(xs, float).reshape(-1, 2)
        else:
            gx = np.arange(nx[0]) * self.L[0] / nx[0]; gy = np.arange(nx[1]) * self.L[1] / nx[1]; self.x = np.array([(a, b) for a in gx for b in gy])
        self.nx = len(self.x)
        E = np.exp(2j * np.pi * (self.x @ self.f.T))
        self.flip = np.array([self.kindex[tuple(-kk)] for kk in self.k])
        self.C = {}; self.P = {}
        for z in self.focus:
            g = self.f2[None, :, :] + self.src[:, None, :] / lam                      # (S, F, 2)
            gn = np.linalg.norm(g, axis=2); inside = (gn * lam <= NA + 1e-12)
            arg = np.clip(1 - (gn * lam / n_img) ** 2, 0.0, None)
            ph = 2 * np.pi * n_img * z / lam * (np.sqrt(arg) - 1)
            if self.aberr: ph = ph + np.stack([zernike_phase(gi, lam, NA, self.aberr) for gi in g])
            self.P[z] = inside * np.exp(1j * ph)                                       # (S, F) pupil vectors per source point
            if not compress: self.C[z] = E[None, :, :] * self.P[z][:, None, :]           # (S, nx, F) per-source channels
        if compress: compress_channels(self, E=E)

def disk_source(sigma, NA, n=21, sigma_in=0.0):
    """Area-weighted grid quadrature of a disk (or annulus sigma_in..sigma) of direction cosines, n x n grid over the disk."""
    g = (np.arange(n) + 0.5) / n * 2 - 1; pts = [(a * sigma * NA, b * sigma * NA) for a in g for b in g if sigma_in <= np.hypot(a, b) * sigma <= sigma + 1e-12]
    return np.array(pts), np.ones(len(pts))
def annular_source(s_in, s_out, NA, n=25): return disk_source(s_out, NA, n=n, sigma_in=s_in)
def dipole_source(sigma_c, NA, pole_radius=0.15, n=9, axis='x'):
    """Two finite circular poles of radius pole_radius (in sigma units) centred at (+-sigma_c, 0) [axis x] or (0, +-sigma_c)."""
    pts = []
    for sgn in (+1, -1):
        c = np.array([sgn * sigma_c, 0.0]) if axis == 'x' else np.array([0.0, sgn * sigma_c])
        g = (np.arange(n) + 0.5) / n * 2 - 1
        pts += [((c[0] + a * pole_radius) * NA, (c[1] + b * pole_radius) * NA) for a in g for b in g if np.hypot(a, b) <= 1 + 1e-12]
    return np.array(pts), np.ones(len(pts))
def quasar_source(s_in, s_out, NA, opening_deg=30.0, n=31):
    """Four annular sectors of opening angle `opening_deg` centred on the diagonals (45, 135, 225, 315 deg)."""
    g = (np.arange(n) + 0.5) / n * 2 - 1; pts = []
    for a in g:
        for b in g:
            r = np.hypot(a, b) * s_out; th = np.degrees(np.arctan2(b, a)) % 90
            if s_in <= r <= s_out + 1e-12 and abs(th - 45) <= opening_deg / 2: pts.append((a * s_out * NA, b * s_out * NA))
    return np.array(pts), np.ones(len(pts))

def el_full(rho): return 2 * (1 - rho) / (1 + rho)   # full-range dose latitude (window width / centre dose) for a threshold resist
def two_beam_closed_form(tau, pitch): return (1 - np.sin(2 * np.pi * tau / pitch)) / (1 + np.sin(2 * np.pi * tau / pitch))

def build_1d(L, lam, NA, src, segs, margin, spacing=4.0, focus=(0.0,), n_img=1.0, aberr=None, compress=True):
    s, w = src; x, B, D, tgt = samples_1d(L, segs, margin, spacing)
    return Model(L, lam, NA, s, w, nx=len(x), focus=focus, xs=x, n_img=n_img, aberr=aberr, compress=compress), B, D, tgt
def build_2d(L, lam, NA, src, target_fn, margin, spacing=10.0, fine=1.0, focus=(0.0,), n_img=1.0, aberr=None, corner_exclusion=0.0, compress=True):
    s, w = src; x, B, D, tgt, sd = samples_2d(L, target_fn, margin, spacing, fine=fine, corner_exclusion=corner_exclusion)
    return Model(L, lam, NA, s, w, nx=None, focus=focus, dim=2, xs=x, n_img=n_img, aberr=aberr, compress=compress), B, D, tgt

class ModelX(Model):
    """Vector (high-NA) extension of the corrected Model: channels = (source point, polarization, field component).
       Direction cosines in the image medium: alpha = g_x lambda/n_img, beta = g_y lambda/n_img, gamma = sqrt(1-alpha^2-beta^2).
       Mask-side polarization p (x, y, unpolarized = average of x and y, TE/TM = azimuthal/radial per source point) is rotated
       by the minimal rotation taking z to the propagation direction (Rodrigues); optional obliquity amplitude 1/sqrt(gamma)."""
    def __init__(self, L, lam, NA, src_pts, src_w, nx, focus=(0.0,), dim=1, xs=None, aberr=None, vector=False, pol='unpol', n_img=1.0, obliquity=True, compress=True):
        self.vector = vector; self.pol = pol; self.obliquity = obliquity
        super().__init__(L, lam, NA, src_pts, src_w, nx, focus, dim, xs, n_img=n_img, aberr=aberr, compress=False)
        if not vector:
            if compress: compress_channels(self)
            return
        E = np.exp(2j * np.pi * (self.x @ self.f.T)); self.P = {}; Pz = {z: [] for z in self.focus}; ws = []
        for s, wsrc in zip(self.src, self.w):
            g = self.f2 + s / lam; gn = np.linalg.norm(g, axis=1); inside = (gn * lam <= NA + 1e-12)
            a = g[:, 0] * lam / n_img; b = g[:, 1] * lam / n_img; gam = np.sqrt(np.clip(1 - a**2 - b**2, 1e-9, 1))
            arg = np.clip(1 - (gn * lam / n_img) ** 2, 0.0, None)
            if pol == 'x': pols = [np.array([1.0, 0.0])]
            elif pol == 'y': pols = [np.array([0.0, 1.0])]
            elif pol == 'unpol': pols = [np.array([1.0, 0.0]), np.array([0.0, 1.0])]
            else:
                ns = np.linalg.norm(s); rad = s / ns if ns > 1e-12 else np.array([1.0, 0.0]); azi = np.array([-rad[1], rad[0]]); pols = [azi if pol == 'TE' else rad]
            ob = (1 / np.sqrt(gam)) if obliquity else np.ones(self.F)
            for z in self.focus:
                ph = 2 * np.pi * n_img * z / lam * (np.sqrt(arg) - 1) + zernike_phase(g, lam, NA, self.aberr); P = inside * np.exp(1j * ph)
                for p in pols:
                    ex = p[0] * (1 - a**2 / (1 + gam)) + p[1] * (-a * b / (1 + gam)); ey = p[0] * (-a * b / (1 + gam)) + p[1] * (1 - b**2 / (1 + gam)); ez = -p[0] * a - p[1] * b
                    for comp in (ex, ey, ez): Pz[z].append(P * ob * comp)
            ws += [wsrc / len(pols)] * (3 * len(pols))
        self.w = np.array(ws); self.P = {z: np.array(Pz[z]) for z in self.focus}
        if compress: compress_channels(self, E=E)
        else: self.C = {z: E[None, :, :] * self.P[z][:, None, :] for z in self.focus}

def build_1d_x(L, lam, NA, src, segs, margin, spacing=4.0, focus=(0.0,), **kw):
    s, w = src; x, B, D, tgt = samples_1d(L, segs, margin, spacing)
    return ModelX(L, lam, NA, s, w, nx=len(x), focus=focus, xs=x, **kw), B, D, tgt
def build_2d_x(L, lam, NA, src, target_fn, margin, spacing=10.0, fine=1.0, focus=(0.0,), **kw):
    s, w = src; x, B, D, tgt, sd = samples_2d(L, target_fn, margin, spacing, fine=fine)
    return ModelX(L, lam, NA, s, w, nx=None, focus=focus, dim=2, xs=x, **kw), B, D, tgt

# ----------------------------------------------------------------------------- exact 2D edge-offset sampling (session 3, code audit item 1)
def rect_target(L, rects):
    """Periodic union of axis-aligned rectangles (x0, y0, w, h) in a cell L=(Lx,Ly). Returns a boolean target function
       with attributes .rects and .L so that samples_2d can place samples at EXACT distance tau from the boundary."""
    Lx, Ly = float(L[0]), float(L[1]); rr = [(float(a) % Lx, float(b) % Ly, float(w), float(h)) for a, b, w, h in rects]
    def tf(a, b): return any(((a - x0) % Lx) < w_ and ((b - y0) % Ly) < h_ for x0, y0, w_, h_ in rr)
    tf.rects = rr; tf.L = (Lx, Ly); return tf

def _union_boundary_segments(tf):
    """Axis-aligned boundary segments of the periodic union of tf.rects, each with its outward unit normal.
       Every rectangle edge is split at all candidate coordinates; a piece is kept iff tf is inside on one side and
       outside on the other (so edges buried inside other rectangles, and collinear coincident edges, are discarded)."""
    Lx, Ly = tf.L; eps = 1e-6
    xs_all = sorted({v % Lx for x0, y0, w, h in tf.rects for v in (x0, x0 + w)}); ys_all = sorted({v % Ly for x0, y0, w, h in tf.rects for v in (y0, y0 + h)})
    def cuts(vals, lo, hi, P):
        c = {lo, hi}
        for v in vals:
            for sh in (-P, 0.0, P, 2 * P):
                if lo + eps < v + sh < hi - eps: c.add(v + sh)
        return sorted(c)
    segs = []   # (orientation 'v'|'h', fixed coord, lo, hi, normal sign)
    for x0, y0, w, h in tf.rects:
        for xf in (x0, x0 + w):                       # vertical edges
            for a, b in zip(cuts(ys_all, y0, y0 + h, Ly)[:-1], cuts(ys_all, y0, y0 + h, Ly)[1:]):
                ym = 0.5 * (a + b)
                for sgn in (+1.0, -1.0):
                    if tf(xf + sgn * eps, ym) == False and tf(xf - sgn * eps, ym) == True: segs.append(('v', xf % Lx, a, b, sgn))
        for yf in (y0, y0 + h):                       # horizontal edges
            for a, b in zip(cuts(xs_all, x0, x0 + w, Lx)[:-1], cuts(xs_all, x0, x0 + w, Lx)[1:]):
                xm = 0.5 * (a + b)
                for sgn in (+1.0, -1.0):
                    if tf(xm, yf + sgn * eps) == False and tf(xm, yf - sgn * eps) == True: segs.append(('h', yf % Ly, a, b, sgn))
    return segs

def exact_boundary_distance(pts, segs, L):
    """Exact Euclidean distance (periodic, minimum image over 3x3 cell images) from points (n,2) to the union of segments."""
    Lx, Ly = L; pts = np.asarray(pts, float); d = np.full(len(pts), np.inf)
    for o, c, lo, hi, _ in segs:
        for i in (-1, 0, 1):
            for j in (-1, 0, 1):
                px = pts[:, 0] + i * Lx; py = pts[:, 1] + j * Ly
                if o == 'v': dx = px - c; dy = np.maximum(0.0, np.maximum(lo - py, py - hi))
                else: dy = py - c; dx = np.maximum(0.0, np.maximum(lo - px, px - hi))
                d = np.minimum(d, np.hypot(dx, dy))
    return d

def boundary_corners(segs, L):
    """Vertices of the union boundary: segment endpoints at which the boundary turns (an endpoint shared by two collinear
       kept segments, e.g. an over-split straight edge, is not a corner). Returns an (m,2) array in [0,Lx)x[0,Ly)."""
    Lx, Ly = L; ends = {}
    for o, c, lo, hi, sgn in segs:
        for t in (lo, hi):
            p = ((c % Lx, t % Ly) if o == 'v' else (t % Lx, c % Ly)); key = (round(p[0], 6), round(p[1], 6))
            ends.setdefault(key, []).append((o, round(c % (Lx if o == 'v' else Ly), 6), sgn))
    out = []
    for key, inc in ends.items():
        collinear = any(sum(1 for q in inc if q == a) >= 2 for a in set(inc))
        if not collinear: out.append(key)
    return np.array(out, float).reshape(-1, 2)

def samples_2d(L, target_fn, margin, spacing, fine=1.0, band=0.6, corner_exclusion=0.0):
    """2D periodic target as a union of rectangles (rect_target). Sample points: (i) the offset contours at EXACT distance
       `margin` on both sides of every boundary segment (points every `fine` nm along the contour, plus the circular arcs
       around convex corners), (ii) a coarse interior grid (step `spacing`) at distance > margin. Every returned point is
       at distance >= margin (asserted) and both rings are non-empty (asserted). `band` is kept for signature compatibility.
       Returns (x, B, D, tgt, sd_grid) where sd_grid is the signed distance on a 1 nm raster (for plots)."""
    if not hasattr(target_fn, 'rects'): raise TypeError("samples_2d needs a rect_target (exact geometry); boolean functions are no longer accepted")
    Lx, Ly = float(L[0]), float(L[1]); tau = float(margin); step = float(fine)
    segs = _union_boundary_segments(target_fn); assert segs, "target has no boundary"
    cand = []
    for o, c, lo, hi, sgn in segs:
        ts = np.append(np.arange(lo, hi, step), hi)
        for t in ts:
            if o == 'v': cand += [(c + sgn * tau, t), (c - sgn * tau, t)]
            else: cand += [(t, c + sgn * tau), (t, c - sgn * tau)]
        for t in (lo, hi):                            # corner arcs (both sides; concave-side points are removed by the distance filter)
            cx, cy = (c, t) if o == 'v' else (t, c)
            for ang in np.arange(0, 2 * np.pi, step / tau): cand.append((cx + tau * np.cos(ang), cy + tau * np.sin(ang)))
    cand = np.array(cand); cand[:, 0] %= Lx; cand[:, 1] %= Ly
    cand = np.unique(np.round(cand, 6), axis=0)
    d = exact_boundary_distance(cand, segs, (Lx, Ly)); keep = d >= tau - 1e-7
    if corner_exclusion > 0:                        # drop offset points whose foot point lies within corner_exclusion of a corner
        V = boundary_corners(segs, (Lx, Ly)); dv = np.full(len(cand), np.inf)
        for v in V:
            for i in (-1, 0, 1):
                for j in (-1, 0, 1): dv = np.minimum(dv, np.hypot(cand[:, 0] + i * Lx - v[0], cand[:, 1] + j * Ly - v[1]))
        keep &= dv >= np.hypot(corner_exclusion, tau) - 1e-9
    ring = cand[keep]; tg = np.array([bool(target_fn(a, b)) for a, b in ring])
    pts = [tuple(p) for p in ring]; tgl = list(tg)
    grid = np.array([(a, b) for a in np.arange(0, Lx, spacing) for b in np.arange(0, Ly, spacing)]) if spacing < max(Lx, Ly) else np.zeros((0, 2))
    if len(grid):
        dg = exact_boundary_distance(grid, segs, (Lx, Ly)); sel = dg > tau + 1e-7
        pts += [tuple(p) for p in grid[sel]]; tgl += [bool(target_fn(a, b)) for a, b in grid[sel]]
    x = np.array(pts); tgt = np.array(tgl)
    dall = exact_boundary_distance(x, segs, (Lx, Ly)); assert dall.min() >= tau - 1e-6, dall.min()
    assert tg.any() and (~tg).any(), "empty bright or dark offset ring"
    gx = np.arange(0, Lx, 1.0); gy = np.arange(0, Ly, 1.0); G = np.array([(a, b) for a in gx for b in gy])
    sdg = exact_boundary_distance(G, segs, (Lx, Ly)); ins = np.array([target_fn(a, b) for a, b in G])
    sd = np.where(ins, sdg, -sdg).reshape(len(gx), len(gy))
    return x, np.where(tgt)[0], np.where(~tgt)[0], tgt, sd

def ensemble_ratio(rm, X, B, D):
    """Dark/bright ratio of the summed image of the lifted solution X (real-structured, F x F real symmetric PSD):
       max_d Tr(Qr_d X) / min_b Tr(Qr_b X). A feasible X realizes this ratio with rank(X) exposures (Theorem 3), so the
       value is an UPPER bound on t_R and the right reference for a certified multi-exposure gain."""
    X = np.real(np.asarray(X)); vals = lambda S: [float(np.sum(rm.Q(ix, z) * X)) for z in rm.focus for ix in S]
    return max(vals(D)) / min(vals(B))

# ----------------------------------------------------------------------------- fast field/gradient evaluation (session 3)
# numpy's stacked (S, nx, F) @ (F,) matmul is ~40x slower than the equivalent 2-D BLAS call, and the three-operand einsum
# 's,sx,sxk->xk' ~15x slower than the two-operand form; every primal routine goes through these helpers.
def _fields(self, M, z=0.0):
    """Channel fields (S, nx) for orders/real coefficients M via one 2-D BLAS matvec."""
    C = self.C[z]; S, nx, F = C.shape
    C2 = self.__dict__.setdefault('_C2', {}).get(z)
    if C2 is None or C2.shape[0] != S * nx: C2 = self._C2[z] = C.reshape(S * nx, F)
    return (C2 @ np.asarray(M, dtype=complex)).reshape(S, nx)
def _image_fast(self, M, z=0.0):
    fld = _fields(self, M, z); return np.einsum('s,sx->x', self.w, np.abs(fld) ** 2)
def _image_grad(self, r, z=0.0):
    """(I, G): intensity (nx,) and d I / d r (nx, F) for REAL coefficients r (RealModel) or complex M (Model: Wirtinger 2 dI/dM^*)."""
    fld = _fields(self, r, z); I = np.einsum('s,sx->x', self.w, np.abs(fld) ** 2)
    G = 2 * np.einsum('sx,sxk->xk', (self.w[:, None] * fld).conj(), self.C[z])
    return I, (G.real if np.isrealobj(r) else G)
for _cls in (Model, RealModel, ModelX):
    _cls.fields = _fields; _cls.image = _image_fast; _cls.image_grad = _image_grad
