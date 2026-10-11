"""
Rigorous branch-and-bound lower bound for the SINGLE-mask ratio problem (real masks, r in R^F):
    rho*(M_R) = min_r max_{(d,z)} r^T Q_dz r / min_{(b,z)} r^T Q_bz r.
Scale and sign invariance -> WLOG max_i |r_i| = 1 attained with r_i = +1. We cover R^F \ {0} (up to scale/sign) by the F
faces {r_i = 1, |r_j| <= 1}. Over a box (center c, half-widths h) and a PSD form q(r) = r^T Q r:
    q(c) - 2|Qc|.h  <=  q(r)  <=  q(c) + 2|Qc|.h + h^T|Q|h      (Q PSD => delta^T Q delta >= 0)
Hence  rho(r) >= max_z [ max_d q_d^lo / min_b q_b^hi ]  =: LB(box).  Global B&B with upper bound U from box centers.
Termination: all boxes pruned (LB >= U(1-tol)) or node budget hit; returns certified L = min over leaves of LB (valid
regardless of budget), U, and statistics. Floating-point safety: error terms inflated by (1+1e-9) and LB reduced by 1e-12.
"""
import numpy as np, heapq, time

def bnb_single_mask(rm, B, D, tol=1e-3, max_nodes=2_000_000, batch=4096, verbose=False, U0=None):
    F = rm.F; zs = list(rm.focus)
    QD = [np.stack([rm.Q(d, z) for d in D]) for z in zs]   # (nD,F,F) per focus
    QB = [np.stack([rm.Q(b, z) for b in B]) for z in zs]
    absQD = [np.abs(q) for q in QD]; absQB = [np.abs(q) for q in QB]
    safe = 1 + 1e-9
    def bounds(C, H):   # C,H: (n,F) centers/half-widths  -> LB (n,), ratio at centers (n,)
        n = C.shape[0]; LB = np.zeros(n); UBc = np.zeros(n)
        for z in range(len(zs)):
            qc_d = np.einsum('nf,dfg,ng->nd', C, QD[z], C); qc_b = np.einsum('nf,bfg,ng->nb', C, QB[z], C)
            Qc_d = np.einsum('dfg,ng->ndf', QD[z], C); Qc_b = np.einsum('bfg,ng->nbf', QB[z], C)
            e1_d = 2*np.einsum('ndf,nf->nd', np.abs(Qc_d), H); e1_b = 2*np.einsum('nbf,nf->nb', np.abs(Qc_b), H)
            e2_b = np.einsum('nf,bfg,ng->nb', H, absQB[z], H)
            d_lo = np.maximum(qc_d - safe*e1_d, 0.0); b_hi = qc_b + safe*(e1_b + e2_b)
            lb = d_lo.max(1) / np.maximum(b_hi.min(1), 1e-300)
            LB = np.maximum(LB, lb); UBc = np.maximum(UBc, qc_d.max(1)/np.maximum(qc_b.min(1), 1e-300))
        return LB - 1e-12, UBc
    # initial faces
    C = np.zeros((F, F)); H = np.ones((F, F))
    for i in range(F): C[i, i] = 1.0; H[i, i] = 0.0
    LB, UBc = bounds(C, H); U = min(UBc.min(), U0 if U0 is not None else np.inf)
    cnt = [0]
    def push(lb, c, h):
        cnt[0] += 1; heapq.heappush(active, (lb, cnt[0], c, h))
    active = []
    for i in range(F): push(float(LB[i]), C[i].copy(), H[i].copy())
    leaves_LB = []   # pruned leaves' LB
    nodes = F; t0 = time.time(); Lfinal = None
    while active:
        if nodes >= max_nodes: break
        # pop a batch of lowest-LB boxes
        popped = [heapq.heappop(active) for _ in range(min(batch, len(active)))]
        keep = [p for p in popped if p[0] < U*(1-tol)]
        leaves_LB += [p[0] for p in popped if p[0] >= U*(1-tol)]
        if not keep: continue
        # split each along widest dimension
        Cn = []; Hn = []
        for lb, _, c, h in keep:
            j = int(np.argmax(h)); hh = h.copy(); hh[j] *= 0.5
            c1 = c.copy(); c1[j] -= hh[j]; c2 = c.copy(); c2[j] += hh[j]
            Cn += [c1, c2]; Hn += [hh, hh.copy()]
        Cn = np.array(Cn); Hn = np.array(Hn); LBn, UBn = bounds(Cn, Hn); nodes += len(Cn)
        U = min(U, UBn.min())
        for i in range(len(Cn)):
            if LBn[i] >= U*(1-tol): leaves_LB.append(float(LBn[i]))
            else: push(float(LBn[i]), Cn[i], Hn[i])
        if verbose and nodes % (batch*20) < 2*batch: print(f"nodes={nodes} active={len(active)} U={U:.6f} minLB_active={active[0][0] if active else None:.6f}", flush=True)
    Lactive = active[0][0] if active else np.inf
    L = min(min(leaves_LB) if leaves_LB else np.inf, Lactive)
    return dict(L=float(L), U=float(U), nodes=nodes, active=len(active), converged=(not active), time=time.time()-t0, rel_gap=(U-L)/U)
