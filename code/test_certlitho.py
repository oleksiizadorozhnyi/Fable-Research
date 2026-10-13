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
    md, B, D, tgt = build_1d(500, lam, NA, (s, w), [(0,40),(100,180),(240,400)], 8, spacing=8.0, focus=(0.0, 80.0, -80.0), n_img=1.44); rm = RealModel(md)
    r = sdp_bound_real(rm, B, D); rig, _, _ = rigorize_certificate(rm, r['mu'], r['lam'])
    v, rv = round_from_Y(rm, r['X'], B, D, nsamp=100); v, rv = local_refine_real(rm, rv, B, D)
    dmax = max(rm.image(rv, z)[D].max() for z in rm.focus); bmin = min(rm.image(rv, z)[B].min() for z in rm.focus)
    assert abs(rm.ratio(rv, B, D) - dmax/bmin) < 1e-12 and v >= rig - 1e-9, (v, rig)
def test_two_beam_closed_form():
    p, tau = 90.0, 6.0; sig = lam/(2*p*NA); s = np.array([[sig*NA, 0.0], [-sig*NA, 0.0]]); w = np.ones(2)
    md, B, D, _ = build_1d(4*p, lam, NA, (s, w), [(i*p, i*p+p/2) for i in range(4)], tau, spacing=4.0)
    r = sdp_bound_real(RealModel(md), B, D); assert abs(r['t'] - two_beam_closed_form(tau, p)) < 2e-5, (r['t'], two_beam_closed_form(tau, p))
def test_catalogue_and_defocus_medium():
    s, w = disk_source(0.7, NA, n=11); segs = [(0,50),(100,150),(200,250),(430,480)]
    md1, B, D, _ = build_1d(600, lam, NA, (s, w), segs, 5, spacing=4.0, focus=(0.0, 80.0, -80.0), n_img=1.44)
    md2, _, _, _ = build_1d(600, lam, NA, (s, w), segs, 5, spacing=4.0, focus=(0.0, 80.0, -80.0), n_img=1.7)
    r1 = sdp_bound_real(RealModel(md1), B, D); r2 = sdp_bound_real(RealModel(md2), B, D); assert r2['t'] < r1['t']   # weaker defocus at higher index
    try: build_1d(600, lam, NA, (s, w), segs, 5, spacing=4.0, focus=(0.0, 80.0, -80.0), n_img=1.0); raise AssertionError("NA > n_img guard did not fire")
    except ValueError: pass
    assert abs(sum(w) - len(w)) < 1e-12 and s.shape[1] == 2
def test_samples_2d_exact():
    from certlitho import _union_boundary_segments, exact_boundary_distance
    Tj = rect_target((320, 320), [(130, 60, 60, 160), (0, 10, 320, 60)]); segs = _union_boundary_segments(Tj)
    assert len(segs) == 8 and abs(sum(hi - lo for _, _, lo, hi, _ in segs) - 940.0) < 1e-9, (len(segs), sum(hi - lo for _, _, lo, hi, _ in segs))
    for tau in (6.0, 7.0, 4.4):
        x, B, D, tgt, sd = samples_2d((320, 320), Tj, tau, 16.0, fine=2.0); d = exact_boundary_distance(x, segs, (320, 320))
        ring = d < tau + 1e-6; assert ring.sum() > 500 and abs(d[ring] - tau).max() < 1e-6 and d.min() >= tau - 1e-6
        assert tgt[B].all() and not tgt[D].any() and len(B) > 100 and len(D) > 100
    # periodic wrap: line-end cell with a rectangle crossing the cell boundary
    le = rect_target((320, 320), [(75, 30, 50, 180), (195, 270, 50, 200)]); x, B, D, tgt, sd = samples_2d((320, 320), le, 6, 16.0, fine=2.0)
    assert exact_boundary_distance(x, _union_boundary_segments(le), (320, 320)).min() >= 6 - 1e-6
    # bound is a valid lower bound for random real masks on a contact array
    holes = rect_target((160, 160), [(45, 45, 70, 70)]); s, w = annular_source(0.6, 0.9, NA, n=9)
    md, B, D, _ = build_2d((160, 160), lam, NA, (s, w), holes, 7, spacing=16.0, fine=2.0); rm = RealModel(md)
    r = sdp_bound_real(rm, B, D); rig, _, _ = rigorize_certificate(rm, r['mu'], r['lam']); rng = np.random.default_rng(0)
    assert all(rm.ratio(rng.standard_normal(rm.F), B, D) >= rig - 1e-9 for _ in range(200))
def test_coherent_mode_compression():
    rng = np.random.default_rng(3); segs = [(0,50),(100,150),(200,250),(430,480)]; s, w = disk_source(0.7, NA, n=11)
    m0, B, D, _ = build_1d(600, lam, NA, (s, w), segs, 5, spacing=4.0, focus=(0.0, 80.0, -80.0), n_img=1.44, compress=False)
    m1, _, _, _ = build_1d(600, lam, NA, (s, w), segs, 5, spacing=4.0, focus=(0.0, 80.0, -80.0), n_img=1.44)
    assert m1.compressed and m1.n_channels <= m1.F < len(m0.w)
    M = rng.standard_normal(m0.F) + 1j * rng.standard_normal(m0.F)
    for z in m0.focus:
        I0, I1 = m0.image(M, z), m1.image(M, z); assert np.abs(I0 - I1).max() <= 1e-10 * I0.max()
        assert np.abs(m0.Q(3, z) - m1.Q(3, z)).max() <= 1e-10 * np.abs(m0.Q(3, z)).max()
    r0 = sdp_bound_real(RealModel(m0), B, D); r1 = sdp_bound_real(RealModel(m1), B, D); assert abs(r0['t'] - r1['t']) < 1e-6
    x0, B, D, _ = build_1d_x(600, lam, NA, (s, w), segs, 5, spacing=4.0, vector=True, pol='TM', n_img=1.7, compress=False)
    x1, _, _, _ = build_1d_x(600, lam, NA, (s, w), segs, 5, spacing=4.0, vector=True, pol='TM', n_img=1.7)
    I0, I1 = x0.image(M), x1.image(M); assert np.abs(I0 - I1).max() <= 1e-10 * I0.max() and x1.n_channels < len(x0.w)
def test_rigorize_no_false_null_space():
    # referee counterexample: a direction where SD ~ 1e-13 and SB ~ 2e-13 caps the true bound at 0.5; treating it as a
    # null space would return 1. Forms: Q_d = diag(1, 1e-13), Q_b = diag(1, 2e-13); rho(M) >= 0.5 is tight at M = e_2.
    class Fake:
        F = 2; focus = (0.0,)
        def Q(self, ix, z=0.0): return np.diag([1.0, 1e-13]) if ix == 0 else np.diag([1.0, 2e-13])
    b, rel, c = rigorize_certificate(Fake(), {(0.0, 0): 1.0}, {(0.0, 1): 1.0})
    assert b <= 0.5 + 1e-9, b
def test_pixel_certificates_exact():
    from exact_verify import mp_source_matrices, real_basis, pixel_basis, real_sym, verify_scale
    s, w = disk_source(0.7, NA, n=21); L = 240.0; segs = [(0, 40), (120, 160)]
    md, B, D, _ = build_1d(L, lam, NA, (s, w), segs, 4, spacing=8.0); N = 48
    rN = sdp_bound_pixel_int(md, N, B, D, nonneg=True, solver="SCS"); rig, rel, c = rigorize_pixel_certificate(md, N, B, D, rN['mu'], rN['lam'], Nmat=rN['Nmat'])
    rP = sdp_bound_pixel_int(md, N, B, D, nonneg=False); rigP, _, cP = rigorize_pixel_certificate(md, N, B, D, rP['mu'], rP['lam'])
    assert rig >= rigP - 1e-9 and rig <= rN['t'] + 1e-6                                # DNN >= PSM pixel bound; rigorized <= solver value
    G = pixel_map(md, N); rng = np.random.default_rng(1)
    for _ in range(100):                                                               # feasibility: random nonnegative pixel masks never beat the DNN bound
        m = rng.uniform(0, 1, N); M = G @ m; I = md.image(M); assert I[D].max() / I[B].min() >= rig - 1e-9
    def parse_keys(d): return d
    S, SB = mp_source_matrices(md, rN['mu'], rN['lam']); V = pixel_basis(md, N); SDm, SBm = real_sym(S, V), real_sym(SB, V)
    smu = sum(rN['mu'].values()); Nm = 0.5*(rN['Nmat'] + rN['Nmat'].T) / smu
    import mpmath as mp; extra = mp.matrix([[mp.mpf(float(Nm[a, b])) for b in range(N)] for a in range(N)])
    kd, Sr, bis, ok = verify_scale(SDm, SBm, c, extra=extra); assert ok and kd <= c
    S2, SB2 = mp_source_matrices(md, rP['mu'], rP['lam']); U = real_basis(md); SDm2, SBm2 = real_sym(S2, U), real_sym(SB2, U); kd2, _, _, ok2 = verify_scale(SDm2, SBm2, cP); assert ok2 and rigP > 0.1
if __name__ == "__main__":
    for f in [test_abbe_cutoff, test_certificate_and_rounding, test_pixel_vs_fourier, test_real_equals_complex_symmetric_source, test_ilt_gradient, test_toeplitz_necessity, test_bnb_consistency, test_vector_model_consistency, test_multifocus_metric, test_two_beam_closed_form, test_catalogue_and_defocus_medium, test_samples_2d_exact, test_coherent_mode_compression, test_rigorize_no_false_null_space, test_pixel_certificates_exact]:
        f(); print("PASS", f.__name__)
