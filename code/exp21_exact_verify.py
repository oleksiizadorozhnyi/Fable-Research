# Exact verification of every reported Fourier-domain certificate (referee item 1 of the first report; "saved scale"
# item of the second report): the certificate matrix is REBUILT in 50-digit arithmetic from the model equations and the
# instance data (exact_verify.py) with the duals (mu, lambda) and the double-precision scale kappa; S~ - 1e-30 I is
# Cholesky-factorized at 50 digits. If the double kappa is not safe for the exact matrix, kappa is lowered by an exact
# bisection and the SAVED value is a double <= the verified value that is itself re-verified before being archived.
from verify_instances import *
from exact_verify import *
import os, sys
out = []; os.makedirs('/home/user/Fable-Research/results/certs', exist_ok=True)
LIMIT = int(sys.argv[1]) if len(sys.argv) > 1 and int(sys.argv[1]) > 0 else 10**9
ONLY = set(sys.argv[2].split(',')) if len(sys.argv) > 2 else None
OUTF = '/home/user/Fable-Research/results/exp21_exact_verify.json'
if ONLY: out = [r for r in json.load(open(OUTF)) if r['name'] not in ONLY]
for ii, (name, md, B, D) in enumerate(verify_instances()):
    if ii >= LIMIT: break
    if ONLY and name not in ONLY: continue
    rm = RealModel(md); t0 = time.time()
    try: r = sdp_bound_real(rm, B, D, solvers=("CLARABEL",))
    except Exception: r = sdp_bound_real(rm, B, D, solvers=("SCS",))
    rig, me, c = rigorize_certificate(rm, r['mu'], r['lam']); kappa = float(c); F = md.F
    S, SB = mp_source_matrices(md, r['mu'], r['lam']); U = real_basis(md); SDm, SBm = real_sym(S, U), real_sym(SB, U)
    norm = max(abs((SDm - mp.mpf(kappa) * SBm)[a, b]) for a in range(F) for b in range(F))
    if norm < mp.mpf('1e-40'): kappa_exact, Sr, bisected, ok = kappa, SDm - mp.mpf(kappa) * SBm, False, True; status = 'zero certificate (PSD trivially)'; lmin = mp.mpf(0)
    else:
        kappa_exact, Sr, bisected, ok = verify_scale(SDm, SBm, kappa)
        try: lmin = min(mp.eigsy(Sr, eigvals_only=True)) if F <= 64 else None
        except Exception: lmin = None
        status = ('PD with margin 1e-30' + (' after exact kappa bisection' if bisected else '')) if ok else 'FAILED'
    slam = sum(r['lam'].values()); smu = sum(r['mu'].values()); bound_exact = kappa_exact * slam / smu
    rec = dict(name=name, F=F, nB=len(B), nD=len(D), n_src=int(len(md.src)), kappa=kappa, kappa_exact=kappa_exact, rigorous=rig, bound_exact=bound_exact, cert_norm=float(norm), lambda_min_50digit=(float(lmin) if lmin is not None else None), cholesky_margin_1e30=ok, status=status, time=time.time() - t0)
    np.savez(f'/home/user/Fable-Research/results/certs/{name}.npz', x=md.x, B=np.array(B), D=np.array(D), k=md.k, L=md.L, wavelength=md.lam, NA=md.NA, n_img=md.n_img, focus=np.array(md.focus), src=md.src, w=np.asarray(md.w_src if hasattr(md, 'w_src') else md.w), mu_keys=np.array([(float(z), int(ix)) for (z, ix) in r['mu']]), mu=np.array(list(r['mu'].values())), lam_keys=np.array([(float(z), int(ix)) for (z, ix) in r['lam']]), lam=np.array(list(r['lam'].values())), kappa=kappa, kappa_exact=kappa_exact, bound=rig, bound_exact=bound_exact)
    print(json.dumps(rec), flush=True); out.append(rec); json.dump(sorted(out, key=lambda q: q['name']) if ONLY else out, open(OUTF, 'w'), indent=1)
