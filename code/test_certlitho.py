"""Regression tests for certlitho. Run: python3 -W ignore test_certlitho.py"""
from pipeline import *
lam, NA = 193.0, 1.35
def test_abbe_cutoff():
    s, w = conventional_source_1d(0.7, NA, 7)
    for pitch, expect_unprintable in [(80, True), (84, True), (88, False)]:
        md, B, D, tgt = build_1d(pitch, lam, NA, (s, w), [(pitch/4, 3*pitch/4)], 0.1*pitch, spacing=pitch/16)
        r = sdp_bound_real(RealModel(md), B, D)
        assert (abs(r['t'] - 1) < 1e-6) == expect_unprintable, (pitch, r['t'])
def test_certificate_and_rounding():
    s, w = conventional_source_1d(0.7, NA, 7)
    md, B, D, tgt = build_1d(600, lam, NA, (s, w), [(0,50),(100,150),(200,250),(430,480)], 8, spacing=4.0); rm = RealModel(md)
    r = sdp_bound_real(rm, B, D); rig, mineig, c = rigorize_certificate(rm, r['mu'], r['lam'])
    assert mineig >= -1e-11 and abs(rig - r['t']) < 1e-5 * r['t'] and c > 0.999
    v, rv = round_from_Y(rm, r['X'], B, D, nsamp=100); assert v >= rig - 1e-9   # weak duality: any mask >= bound
    # interior sampling invariance
    md2, B2, D2, _ = build_1d(600, lam, NA, (s, w), [(0,50),(100,150),(200,250),(430,480)], 8, spacing=1.0)
    r2 = sdp_bound_real(RealModel(md2), B2, D2); assert abs(r2['t'] - r['t']) < 1e-5 * r['t'], (r2['t'], r['t'])
def test_pixel_vs_fourier():
    s, w = conventional_source_1d(0.7, NA, 7)
    md, B, D, tgt = build_1d(360, lam, NA, (s, w), [(0,60),(120,180),(240,300)], 10, spacing=5.0)
    rF = sdp_bound_real(RealModel(md), B, D); rP = sdp_bound_pixel(md, 120, B, D, solver="SCS")
    assert abs(rF['t'] - rP['t']) < 2e-4, (rF['t'], rP['t'])
def test_real_equals_complex_symmetric_source():
    s, w = conventional_source_1d(0.7, NA, 7)
    md, B, D, tgt = build_1d(500, lam, NA, (s, w), [(0,40),(100,180),(240,400)], 8, spacing=4.0)
    rR = sdp_bound_real(RealModel(md), B, D); rC = sdp_bound(md, B, D, real_mask=False, solver="CLARABEL")
    assert abs(rR['t'] - rC['t']) < 1e-5, (rR['t'], rC['t'])
def test_ilt_gradient():
    s, w = conventional_source_1d(0.7, NA, 7); segs=[(0,60),(120,180),(240,300)]
    md, B, D, tgt = build_1d(360, lam, NA, (s, w), segs, 10, spacing=5.0); N = 60; rng = np.random.default_rng(0); m = rng.random(N); a, th = 10.0, 0.3
    def loss(m):
        M = md.orders_from_pixels(m); fld = md.C[0.0] @ M; I = np.einsum('s,sx->x', md.w, np.abs(fld)**2); r = 1/(1+np.exp(-a*(I-th))); return np.sum((r-tgt)**2)
    M = md.orders_from_pixels(m); fld = md.C[0.0] @ M; I = np.einsum('s,sx->x', md.w, np.abs(fld)**2); r = 1/(1+np.exp(-a*(I-th))); dL = 2*(r-tgt)*a*r*(1-r)
    gM = np.einsum('s,sxk,sx->k', md.w, md.C[0.0].conj(), dL[None,:]*fld); n = np.arange(N)
    gm = 2*np.real(np.array([np.sum(gM*np.exp(2j*np.pi*md.k[:,0]*nn/N)) for nn in n]))/N
    num = np.array([(loss(m+1e-6*np.eye(N)[i]) - loss(m-1e-6*np.eye(N)[i]))/2e-6 for i in range(N)])
    assert np.linalg.norm(gm-num)/np.linalg.norm(num) < 1e-6
def test_toeplitz_necessity():
    rng = np.random.default_rng(1); N, K = 256, 6
    for _ in range(200):
        m = rng.random(N)*(rng.random(N) < 0.5); Mf = np.fft.fft(m)/N; M = np.array([Mf[k] for k in range(K+1)])
        T = np.array([[M[a-b] if a >= b else np.conj(M[b-a]) for b in range(K+1)] for a in range(K+1)])
        assert np.linalg.eigvalsh(T).min() >= -1e-12*abs(M[0])
def test_bnb_consistency():
    from bnb_cert import bnb_single_mask
    s, w = ([0.0],[1.0]); md, B, D, tgt = build_1d(300, lam, NA, (s, w), [(30,80),(150,230)], 10, spacing=6.0); rm = RealModel(md)   # F=5
    r = sdp_bound_real(rm, B, D); rig, _, _ = rigorize_certificate(rm, r['mu'], r['lam'])
    v, rv = round_from_Y(rm, r['X'], B, D, nsamp=100); v, rv = local_refine_real(rm, rv, B, D)
    res = bnb_single_mask(rm, B, D, tol=5e-3, max_nodes=1500000, U0=v)
    assert res['converged'] and res['L'] <= res['U'] + 1e-12 and res['L'] >= rig - 1e-9, (res, rig)   # single-mask optimum >= multi-exposure bound
def test_vector_model_consistency():
    s, w = conventional_source_1d(0.7, NA, 7); segs = [(0,45),(90,135),(180,225),(270,315)]
    md0, B, D, _ = build_1d(360, lam, NA, (s, w), segs, 6, spacing=8.0); mdx, _, _, _ = build_1d_x(360, lam, NA, (s, w), segs, 6, spacing=8.0)
    r0 = sdp_bound_real(RealModel(md0), B, D); rx = sdp_bound_real(RealModel(mdx), B, D); assert abs(r0['t'] - rx['t']) < 1e-7
    mdv, _, _, _ = build_1d_x(360, lam, NA, (s, w), segs, 6, spacing=8.0, vector=True, pol='x', n_img=1.7)
    rR = sdp_bound_real(RealModel(mdv), B, D); rC = sdp_bound(mdv, B, D, real_mask=False, solver="CLARABEL"); assert abs(rR['t'] - rC['t']) < 1e-5   # Lemma 6 under vector model
    # TM (x-pol for lines along y) must not beat TE
    mdt, _, _, _ = build_1d_x(360, lam, NA, (s, w), segs, 6, spacing=8.0, vector=True, pol='y', n_img=1.7); rT = sdp_bound_real(RealModel(mdt), B, D)
    assert rR['t'] >= rT['t'] - 1e-9

def test_multifocus_metric():
    s, w = conventional_source_1d(0.7, NA, 7)
    md, B, D, tgt = build_1d(500, lam, NA, (s, w), [(0,40),(100,180),(240,400)], 8, spacing=8.0, focus=(0.0, 80.0, -80.0)); rm = RealModel(md)
    r = sdp_bound_real(rm, B, D); rig, _, _ = rigorize_certificate(rm, r['mu'], r['lam'])
    v, rv = round_from_Y(rm, r['X'], B, D, nsamp=100); v, rv = local_refine_real(rm, rv, B, D)
    dmax = max(rm.image(rv, z)[D].max() for z in rm.focus); bmin = min(rm.image(rv, z)[B].min() for z in rm.focus)
    assert abs(rm.ratio(rv, B, D) - dmax/bmin) < 1e-12 and v >= rig - 1e-9, (v, rig)
if __name__ == "__main__":
    for f in [test_abbe_cutoff, test_certificate_and_rounding, test_pixel_vs_fourier, test_real_equals_complex_symmetric_source, test_ilt_gradient, test_toeplitz_necessity, test_bnb_consistency, test_vector_model_consistency, test_multifocus_metric]:
        f(); print("PASS", f.__name__)
