# Re-audit against the session-3 Model: basis, DFT, defocus formula, ModelX==Model, dead orders, new sources
import numpy as np, warnings, json; warnings.filterwarnings('ignore')
from certlitho import *; from pipeline import *
lam, NA, n = LAM_DUV, NA_DUV, NIMG_DUV; rng = np.random.default_rng(0)
# (2) RealModel with 2-D source points for 1-D masks, multi-focus
md, B, D, _ = build_1d(500, lam, NA, SRC_DUV['conv0.7'], [(0,40),(100,180),(240,400)], 5, spacing=4.0, focus=(0.0, 80.0), n_img=n); rm = RealModel(md); U = rm.U
worst = 0
for _ in range(30):
    r = rng.standard_normal(md.F); M = U @ r; assert np.abs(M[md.flip]-M.conj()).max() < 1e-12
    for z in md.focus:
        for ix in list(B[:3])+list(D[:3]): worst = max(worst, abs(r @ rm.Q(ix, z) @ r - np.real(M.conj() @ md.Q(ix, z) @ M)))
print(f"new Model 1D F={md.F} nsrc={len(md.src)}: unitary err {np.abs(U.conj().T@U-np.eye(md.F)).max():.1e}, form err {worst:.1e}")
# defocus phase: check against immersion formula for each (order, source); and against air formula when n_img=1
for nn in [n, 1.0]:
    mdz, *_ = build_1d(500, lam, NA, SRC_DUV['conv0.7'], [(0,40),(100,180),(240,400)], 5, spacing=4.0, focus=(0.0, 80.0), n_img=nn)
    err = 0; nclip = 0
    for si, s in enumerate(mdz.src):
        g = mdz.f2 + s/lam; gn = np.linalg.norm(g, axis=1); inside = gn*lam <= NA+1e-12
        ph_ref = 2*np.pi*nn*80/lam*(np.sqrt(np.clip(1-(gn*lam/nn)**2, 0, None))-1)
        ratio = mdz.C[80.0][si, 0, :]/np.where(inside, mdz.C[0.0][si, 0, :], 1)
        err = max(err, np.abs(np.angle(ratio[inside]) - np.angle(np.exp(1j*ph_ref[inside]))).max()); nclip += int((inside & (gn*lam > nn)).sum())
    print(f"  n_img={nn}: defocus phase matches 2 pi n z/lam (sqrt(1-(g lam/n)^2)-1): err {err:.1e}; passing pairs with |g|lam>n_img (clipped): {nclip}")
# ModelX(vector=False) == Model (trivial now), ModelX vector: rotation norm and weight sum
segs = [(0,45),(90,135),(180,225),(270,315)]
md0, B, D, _ = build_1d(360, lam, NA, SRC_DUV['dipole0.8'], segs, 5, spacing=8.0, focus=(0.0, 50.0), n_img=n); mdx, *_ = build_1d_x(360, lam, NA, SRC_DUV['dipole0.8'], segs, 5, spacing=8.0, focus=(0.0,50.0), n_img=n)
print("ModelX(vector=False) vs Model C err", max(np.abs(md0.C[z]-mdx.C[z]).max() for z in md0.focus))
for pol in ['x','y','unpol','TE','TM']:
    mdv, *_ = build_1d_x(360, lam, NA, SRC_DUV['dipole0.8'], segs, 5, spacing=8.0, vector=True, pol=pol, n_img=n, obliquity=False, focus=(0.0, 50.0))
    S = len(mdv.src); npol = 2 if pol == 'unpol' else 1
    for z in mdv.focus:
        ch = mdv.C[z].reshape(S, npol, 3, mdv.nx, mdv.F); nrm = (np.abs(ch)**2).sum(2); nrm0 = np.abs(md0.C[z])**2
        e = np.abs(nrm - nrm0[:, None]).max()
    # vector defocus consistent with scalar defocus? compare phase of ex channel ratio z/0
    print(f"  pol={pol}: sum w={mdv.w.sum():.3f}, |sum_c|E_c|^2-|P|^2| max={e:.1e}")
# Lemma 6 under the new model with defocus: real == complex relaxation for Z=-Z
rR = sdp_bound_real(RealModel(md0), B, D); rC = sdp_bound(md0, B, D, real_mask=False); print(f"Lemma6 (dipole, foci 0,50 -> Z!=-Z): real {rR['t']:.6f} complex {rC['t']:.6f}")
md1, B1, D1, _ = build_1d(360, lam, NA, SRC_DUV['dipole0.8'], segs, 5, spacing=8.0, focus=(0.0, 50.0, -50.0), n_img=n)
rR = sdp_bound_real(RealModel(md1), B1, D1); rC = sdp_bound(md1, B1, D1, real_mask=False); print(f"Lemma6 (dipole, foci 0,+-50 -> Z=-Z): real {rR['t']:.6f} complex {rC['t']:.6f}")
# dead orders with the new catalogue sources (2D)
holes = lambda P, hw: (lambda a, b: (abs((a+P/2)%P-P/2) < hw) and (abs((b+P/2)%P-P/2) < hw))
for nm in ['annular', 'quasar', 'conv0.5']:
    md2, B2, D2, _ = build_2d((160,160), lam, NA, SRC_DUV[nm], holes(160,35), 12, spacing=16.0, fine=2.0, n_img=n)
    alive = (np.abs(md2.C[0.0][:, 0, :]) > 0).any(0); print(f"2D contacts160 {nm}: F={md2.F} nsrc={len(md2.src)} dead orders={int((~alive).sum())}")
Tj = lambda a, b: (abs(a-160) < 30 and 60 < b < 220) or (abs(b-40) < 30)
md4 = Model((320,320), lam, NA, *SRC_DUV['annular'], nx=None, dim=2, xs=np.array([[0.0,0.0]]), n_img=n); alive = (np.abs(md4.C[0.0][:, 0, :]) > 0).any(0); print(f"2D 320 annular: F={md4.F} dead={int((~alive).sum())}")
# source catalogue sanity: symmetric? sum weights; smax
for k, (s, w) in SRC_DUV.items():
    sym = all(any(np.allclose(-p, q) for q in s) for p in s); print(f"  src {k}: n={len(s)} smax/NA={np.linalg.norm(s,axis=1).max()/NA:.3f} symmetric={sym}")
