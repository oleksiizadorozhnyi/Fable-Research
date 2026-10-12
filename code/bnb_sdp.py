"""
SDP-bounded branch-and-bound certificate for the single-mask problem (real masks).
Goal: certify  rho*(M_R) >= T  for a given target T (e.g. T = (1-eps)*U, U = best single mask found).
rho(r) < T  <=>  q_d(r) < T p_b(r) for all (d,b)  (q_d dark forms, p_b bright forms, PSD).
Faces Phi_i = {r_i = 1, |r_j| <= 1} cover all masks up to scale/sign. On a box l <= r <= u (l_i = u_i = 1) we certify that
    {r in box : <Q_d - T P_b, r r^T> <= 0 for all (d,b) in the active sets}  is empty
via the Shor relaxation Z = [[1, r^T],[r, X]] >= 0 with RLT box cuts, by solving  min s  s.t. every constraint <= s  and
verifying an explicit dual certificate: multipliers y >= 0 (sum y = 1), free alpha (Z_00 = 1), beta (X_ii = 1) with
    Lambda := sum_m y_m G_m + alpha E_00 + beta E_ii  >= 0   (checked with a relative margin)
    => for every feasible Z:  s >= sum_m y_m <G_m,Z> + sum_m y_m c_m  >= ... weak duality gives s* >= sum y_m c_m - alpha - beta =: LB.
If LB > 0 the box contains no r with rho(r) < T.  Nodes that cannot be certified are split along the widest coordinate.
Using only a subset of (d,b) pairs is valid (fewer constraints -> relaxation is weaker, never unsound).
"""
import numpy as np, cvxpy as cp, heapq, time

def _constraints_matrices(F, i_fixed, l, u, A_pairs):
    """Build (G_m, c_m) with g_m(Z) = <G_m, Z> + c_m <= s. Z is (F+1)x(F+1) symmetric, index 0 = '1', 1..F = r."""
    Gs, cs = [], []
    n = F + 1
    for A in A_pairs:                       # <A, X> <= 0
        G = np.zeros((n, n)); G[1:, 1:] = A; Gs.append(G); cs.append(0.0)
    def lin(j, coef):                       # coef * r_j  -> symmetric placement in Z[0,j] and Z[j,0]
        G = np.zeros((n, n)); G[0, j+1] = coef/2; G[j+1, 0] = coef/2; return G
    def quad(j, k, coef):                   # coef * X_jk
        G = np.zeros((n, n)); G[j+1, k+1] += coef/2; G[k+1, j+1] += coef/2; return G
    for j in range(F):
        for k in range(j, F):
            if j == i_fixed and k == i_fixed: continue
            # (u_j - r_j)(u_k - r_k) >= 0  ->  -(u_j u_k - u_j r_k - u_k r_j + X_jk) <= 0
            Gs.append(-(lin(k, -u[j]) + lin(j, -u[k]) + quad(j, k, 1.0))); cs.append(-u[j]*u[k])
            # (r_j - l_j)(r_k - l_k) >= 0  ->  -(X_jk - l_k r_j - l_j r_k + l_j l_k) <= 0
            Gs.append(-(quad(j, k, 1.0) + lin(j, -l[k]) + lin(k, -l[j]))); cs.append(-l[j]*l[k])
            # (u_j - r_j)(r_k - l_k) >= 0  ->  -(u_j r_k - u_j l_k - X_jk + l_k r_j) <= 0
            Gs.append(-(lin(k, u[j]) + quad(j, k, -1.0) + lin(j, l[k]))); cs.append(u[j]*l[k])
            if j != k:
                Gs.append(-(lin(j, u[k]) + quad(j, k, -1.0) + lin(k, l[j]))); cs.append(u[k]*l[j])
    return Gs, cs

def node_certificate(F, i_fixed, l, u, A_pairs, margin=1e-9, solver="CLARABEL"):
    """Solve min s and verify the dual certificate. Returns (LB, verified, s_primal)."""
    n = F + 1
    Gs, cs = _constraints_matrices(F, i_fixed, l, u, A_pairs)
    Z = cp.Variable((n, n), symmetric=True); s = cp.Variable()
    cons = [Z >> 0, Z[0, 0] == 1, Z[i_fixed+1, i_fixed+1] == 1]
    lin_cons = [cp.trace(G @ Z) + c <= s for G, c in zip(Gs, cs)]
    prob = cp.Problem(cp.Minimize(s), cons + lin_cons)
    try: prob.solve(solver=solver)
    except Exception:
        try: prob.solve(solver="SCS", eps=1e-9, max_iters=100000)
        except Exception: return -np.inf, False, None
    if s.value is None: return -np.inf, False, None
    y = np.array([max(float(np.ravel(c.dual_value)[0]), 0.0) if c.dual_value is not None else 0.0 for c in lin_cons])
    if y.sum() <= 0: return -np.inf, False, float(s.value)
    y = y / y.sum()
    # Lambda = sum y_m G_m + alpha E00 + beta Eii must be PSD: choose alpha, beta to minimise the violation? They are the
    # equality multipliers from the solver (duals of Z00=1 and Zii=1, sign conventions vary) -> instead fit them: Lambda(alpha,beta)
    # PSD iff ... we take the solver values and verify; LB = sum y c - alpha - beta.
    a = float(np.ravel(cons[1].dual_value)[0]) if cons[1].dual_value is not None else 0.0
    b = float(np.ravel(cons[2].dual_value)[0]) if cons[2].dual_value is not None else 0.0
    S = sum(yy*G for yy, G in zip(y, Gs)); best = (-np.inf, False)
    for sa in (+1, -1):
        for sb in (+1, -1):
            Lam = S.copy(); Lam[0, 0] += sa*a; Lam[i_fixed+1, i_fixed+1] += sb*b
            w = np.linalg.eigvalsh(0.5*(Lam+Lam.T)); nrm = max(np.abs(w).max(), 1e-300)
            if w.min() >= -margin*nrm:
                # rigorous: Lam + |wmin| I >= 0; account by subtracting |min(w,0)| * (Tr Z bound) ... Tr Z <= 1 + sum u_j^2 <= 1+F
                slack = max(-w.min(), 0.0)*(1 + F)
                LB = float(np.dot(y, cs) - sa*a - sb*b) - slack
                if LB > best[0]: best = (LB, True)
    return best[0], best[1], float(s.value)

def certify_target(rm, B, D, T, mu=None, lam=None, max_nodes=20000, wmin=1e-3, verbose=False, active_only=True):
    """Certify rho*(M_R) >= T by B&B over faces. Uses active dark/bright subsets (points with nonzero global dual weight) if given."""
    F = rm.F; zs = list(rm.focus)
    Dset = [(z, d) for z in zs for d in D]; Bset = [(z, b) for z in zs for b in B]
    if active_only and mu is not None:
        Dset = [k for k in Dset if mu.get(k, 0) > 1e-6*max(mu.values())]; Bset = [k for k in Bset if lam.get(k, 0) > 1e-6*max(lam.values())]
    A_pairs = [rm.Q(d, z) - T*rm.Q(b, zb) for (z, d) in Dset for (zb, b) in Bset]
    t0 = time.time(); nodes = 0; failed = []; heap = []; cnt = 0
    for i in range(F):
        l = -np.ones(F); u = np.ones(F); l[i] = u[i] = 1.0; heapq.heappush(heap, (0.0, cnt, i, l, u)); cnt += 1
    while heap:
        if nodes >= max_nodes: break
        _, _, i, l, u = heapq.heappop(heap); nodes += 1
        LB, ok, sp = node_certificate(F, i, l, u, A_pairs)
        if ok and LB > 0: continue                                  # box certified empty
        h = u - l
        if h.max() < wmin: failed.append((i, l.copy(), u.copy(), LB, sp)); continue
        j = int(np.argmax(h)); m = 0.5*(l[j]+u[j])
        l1, u1 = l.copy(), u.copy(); u1[j] = m; l2, u2 = l.copy(), u.copy(); l2[j] = m
        heapq.heappush(heap, (-(h.max()), cnt, i, l1, u1)); cnt += 1; heapq.heappush(heap, (-(h.max()), cnt, i, l2, u2)); cnt += 1
        if verbose and nodes % 50 == 0: print(f"nodes={nodes} open={len(heap)} failed={len(failed)} t={time.time()-t0:.0f}s", flush=True)
    certified = (not heap) and (not failed)
    return dict(T=T, certified=certified, nodes=nodes, open=len(heap), failed=len(failed), time=time.time()-t0, n_pairs=len(A_pairs), nD=len(Dset), nB=len(Bset))
