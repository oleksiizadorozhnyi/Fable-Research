# Audit (1),(2),(8),(9): Abbe model conventions, RealModel basis, ModelX==Model, DFT roundtrip, defocus/immersion
import numpy as np, warnings; warnings.filterwarnings('ignore')
from certlitho import *
lam, NA = 193.0, 1.35
rng = np.random.default_rng(0)
# ---- (2) RealModel U unitary, spans Hermitian-symmetric subspace, r^T Qr r == M^H Q M
for dim, args in [(1, dict(L=500, src=conventional_source_1d(0.7, NA, 7))), (2, dict(L=(160,160), src=annular_source_2d(0.6,0.9,NA,9)))]:
    if dim == 1:
        md, B, D, _ = build_1d(500, lam, NA, args['src'], [(0,40),(100,180),(240,400)], 8, spacing=4.0, focus=(0.0, 80.0))
    else:
        holes = lambda P, hw: (lambda a, b: (abs((a+P/2)%P-P/2) < hw) and (abs((b+P/2)%P-P/2) < hw))
        md, B, D, _ = build_2d((160,160), lam, NA, args['src'], holes(160,35), 12, spacing=16.0, fine=2.0, focus=(0.0, 60.0))
    rm = RealModel(md); U = rm.U
    print(f"dim={dim} F={md.F} unitary err={np.abs(U.conj().T@U-np.eye(md.F)).max():.2e}")
    worst = 0; worstI = 0
    for _ in range(50):
        r = rng.standard_normal(md.F); M = U @ r
        assert np.abs(M[md.flip] - M.conj()).max() < 1e-12   # Hermitian symmetric
        # real mask -> orders -> should equal M-structure; also random real mask from pixels
        for z in md.focus:
            for ix in list(B[:3]) + list(D[:3]):
                Q = md.Q(ix, z); worst = max(worst, abs(r @ rm.Q(ix, z) @ r - np.real(M.conj() @ Q @ M)))
            worstI = max(worstI, np.abs(rm.image(r, z) - md.image(M, z)).max())
    print(f"   r^T Qr r vs M^H Q M max err {worst:.2e}; image err {worstI:.2e}")
    # real pixel mask -> orders lie in span(U)? (from_orders/to_orders roundtrip)
    for _ in range(20):
        m = rng.standard_normal(256 if dim == 1 else (64,64))
        M = md.orders_from_pixels(m); r = rm.from_orders(M); assert np.abs(rm.to_orders(r) - M).max() < 1e-12
    print("   real pixel masks lie in span(U): OK")
# ---- (1) DFT roundtrip orders_from_pixels / pixels_from_orders
md, B, D, _ = build_1d(500, lam, NA, conventional_source_1d(0.7, NA, 7), [(0,40),(100,180),(240,400)], 8)
rm = RealModel(md); r = rng.standard_normal(md.F); M = rm.to_orders(r)
m = md.pixels_from_orders(M, 128); M2 = md.orders_from_pixels(m); print(f"1D DFT roundtrip err {np.abs(M2-M).max():.2e}")
# check the Fourier convention of orders_from_pixels against the definition M_k = (1/L) int m e^{-2 pi i k x/L}
m = rng.standard_normal(4096); L = md.L[0]; xg = np.arange(4096)*L/4096
Mk = np.array([np.mean(m*np.exp(-2j*np.pi*k*xg/L)) for k in md.k[:,0]]); print(f"1D definition check {np.abs(Mk-md.orders_from_pixels(m)).max():.2e}")
holes = lambda P, hw: (lambda a, b: (abs((a+P/2)%P-P/2) < hw) and (abs((b+P/2)%P-P/2) < hw))
md2, B2, D2, _ = build_2d((160,160), lam, NA, annular_source_2d(0.6,0.9,NA,9), holes(160,35), 12, spacing=16.0, fine=2.0)
rm2 = RealModel(md2); M = rm2.to_orders(rng.standard_normal(md2.F)); m = md2.pixels_from_orders(M, (32,32)); print(f"2D DFT roundtrip err {np.abs(md2.orders_from_pixels(m)-M).max():.2e}")
# image of a shifted mask == shifted image (Abbe convention consistency)
sh = 37.0; Msh = M * np.exp(-2j*np.pi*(md2.f @ np.array([sh, 0.0])))
# compare I_Msh(x) with I_M(x - sh): build second model with shifted sample points
md3 = Model((160,160), lam, NA, *annular_source_2d(0.6,0.9,NA,9), nx=None, dim=2, xs=md2.x - np.array([sh, 0.0]))
print(f"shift covariance err {np.abs(md2.image(Msh) - md3.image(M)).max():.2e}")
# ---- (1) pupil / order set: how many collected orders are 'dead' (never pass the pupil for any source point)?
for name, mdl in [("1D conv0.7 L500", md), ("2D contacts160 annular", md2)]:
    alive = np.zeros(mdl.F, bool)
    for z in mdl.focus: alive |= (np.abs(mdl.C[z][:, 0, :]) > 0).any(0)
    print(f"{name}: F={mdl.F} dead orders={int((~alive).sum())}")
Tj = lambda a, b: (abs(a-160) < 30 and 60 < b < 220) or (abs(b-40) < 30)
md4, B4, D4, _ = build_2d((320,320), lam, NA, annular_source_2d(0.6,0.9,NA,9), Tj, 14, spacing=16.0, fine=2.0)
alive = (np.abs(md4.C[0.0][:, 0, :]) > 0).any(0); print(f"2D Tjunction320 annular: F={md4.F} dead orders={int((~alive).sum())} -> {md4.k[~alive].tolist()}")
# ---- (1) defocus: immersion NA=1.35 -> |g| lam can exceed 1 -> clipped
s, w = conventional_source_1d(0.7, NA, 7); mdf = build_1d(500, lam, NA, (s,w), [(0,40),(100,180),(240,400)], 8, focus=(0.0, 80.0))[0]
cnt = 0; tot = 0
for sp in mdf.src:
    g = mdf.f + sp/lam; gn = np.linalg.norm(g, axis=1)*lam; inside = gn <= NA+1e-12; cnt += int(((gn > 1.0) & inside).sum()); tot += int(inside.sum())
print(f"defocus: passing (order,source) pairs with |g|lam>1 (clipped to 1): {cnt}/{tot}")
# phases actually used vs immersion formula (n=1.44) for a few |g|lam values
for gl in [0.3, 0.7, 1.0, 1.2, 1.35]:
    used = 2*np.pi/lam*80*(np.sqrt(1-min(gl,1)**2)-1); n=1.44; imm = 2*np.pi*n/lam*80*(np.sqrt(1-(gl/n)**2)-1)
    print(f"  |g|lam={gl}: phase used {used:+.4f} rad, immersion(n=1.44) {imm:+.4f} rad")
# ---- (8) ModelX(vector=False) == Model
segs = [(0,45),(90,135),(180,225),(270,315)]
md0, B, D, _ = build_1d(360, lam, NA, (s,w), segs, 6, spacing=8.0, focus=(0.0, 50.0)); mdx, *_ = build_1d_x(360, lam, NA, (s,w), segs, 6, spacing=8.0, focus=(0.0,50.0))
print("ModelX(vector=False) vs Model: C err", max(np.abs(md0.C[z]-mdx.C[z]).max() for z in md0.focus), "w err", np.abs(md0.w-mdx.w).max())
# ModelX vector: per-channel weights, polarization norm
for pol in ['x','y','unpol','TE','TM']:
    mdv, *_ = build_1d_x(360, lam, NA, (s,w), segs, 6, spacing=8.0, vector=True, pol=pol, n_img=1.7, obliquity=False)
    # sum over components of |E_c|^2 per (source,pol) must equal scalar |P|^2 (rotation is unitary)
    S = len(mdv.src); npol = 1 if pol in ('x','y','TE','TM') else 2; C = mdv.C[0.0]; ch = C.reshape(S, npol, 3, mdv.nx, mdv.F)
    nrm = (np.abs(ch)**2).sum(2); nrm0 = np.abs(md0.C[0.0] if md0.nx == mdv.nx else mdx.C[0.0])**2
    print(f"  pol={pol}: sum w = {mdv.w.sum():.3f} (expected 3.0 = 3 comps), max |sum_c|E_c|^2 - |P|^2| = {np.abs(nrm - nrm0[:,None]).max():.2e}")
# vector with n_img=1 and no obliquity, pol x, 1D: Ey=0? E = (1-a^2/(1+g), 0, -a) -> intensity |Ex|^2+|Ez|^2 = 1 -> equals scalar
mdv, *_ = build_1d_x(360, lam, NA, (s,w), segs, 6, spacing=8.0, vector=True, pol='x', n_img=1.0, obliquity=False)
# compare images of random real mask: NOT expected equal (Ex and Ez interfere differently) -- just report
rmv = RealModel(mdv); rm0 = RealModel(md0); r = rng.standard_normal(md0.F)
print("  vector(n=1, no obliq, x-pol) vs scalar image rel diff:", np.abs(rmv.image(r)-rm0.image(r)).max()/rm0.image(r).max())
# ---- (9) sinc convention of sdp_bound_pixel_int: compare G against exact Fourier coefficients of piecewise-constant mask
N = 30; L = md.L[0]; n = np.arange(N); k = md.k[:,0]; m = rng.standard_normal(N)
G = (np.sinc(k/N)[:, None]) * np.exp(-2j*np.pi*np.outer(k, n)/N)/N
xf = np.arange(0, L, L/N/400); mf = m[(xf//(L/N)).astype(int)]   # pixel n covers [nL/N,(n+1)L/N)
Mexact = np.array([np.mean(mf*np.exp(-2j*np.pi*kk*xf/L)) for kk in k])
print(f"pixel_int: |G m - M_exact(edges at nL/N)| = {np.abs(G@m - Mexact).max():.2e};  |G m - M_exact * e^{{+i pi k/N}}| = {np.abs(G@m - Mexact*np.exp(1j*np.pi*k/N)).max():.2e} (i.e. G assumes pixel CENTERS at nL/N)")
