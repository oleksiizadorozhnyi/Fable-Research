# Physics sensitivity: vector imaging (polarizations), aberrations, pixelation; tightness under the vector model.
from pipeline import *
lam, NA = 193.0, 1.35; out = {}
conv = SRC_DUV['conv0.7']; dip = SRC_DUV['dipole0.8']
# (a) vector vs scalar, 1D
res = []
for cname, L, segs, margin in [("dense90", 360, [(0,45),(90,135),(180,225),(270,315)], 6), ("isodense", 600, [(0,50),(100,150),(200,250),(430,480)], 8)]:
    for sname, src in [("conv0.7", conv), ("dipole0.8", dip)]:
        md0, B, D, _ = build_1d(L, lam, NA, src, segs, margin, spacing=4.0); rec0, *_ = analyze(md0, B, D, nround=200, nstart=8)
        row = dict(case=cname, src=sname, scalar_bound=rec0['rigorous'], scalar_best=rec0['best_feasible'], scalar_gap=rec0['gap_rel'])
        for pol in ['TE', 'TM', 'unpol']:
            mdv, B, D, _ = build_1d_x(L, lam, NA, src, segs, margin, spacing=4.0, vector=True, pol=pol, n_img=1.7)
            recv, *_ = analyze(mdv, B, D, nround=200, nstart=8)
            row[f'{pol}_bound'] = recv['rigorous']; row[f'{pol}_best'] = recv['best_feasible']; row[f'{pol}_gap'] = recv['gap_rel']; row[f'{pol}_rank'] = recv['rank']
        print(json.dumps(row), flush=True); res.append(row)
out['vector_1d'] = res
# (b) aberration sensitivity (isodense, conv0.7)
res = []
for key in ['Z5', 'Z7', 'Z9']:
    for amp in [0.0, 0.02, 0.05, 0.1]:
        md, B, D, _ = build_1d_x(600, lam, NA, conv, [(0,50),(100,150),(200,250),(430,480)], TAU_DUV, spacing=4.0, aberr={key: amp})
        r = sdp_bound_real(RealModel(md), B, D); rig, _, _ = rigorize_certificate(RealModel(md), r['mu'], r['lam'])
        res.append(dict(zernike=key, waves=amp, bound=rig)); print(json.dumps(res[-1]), flush=True)
out['aberrations'] = res
# (c) 2D contacts: vector unpolarized / TE vs scalar
res = []
ann = SRC_DUV['annular']; qua = SRC_DUV['quasar']
holes = lambda P, hw: (lambda a, b: (abs((a+P/2)%P-P/2) < hw) and (abs((b+P/2)%P-P/2) < hw))
for name, L, tf, src, margin in [("contacts120_annular", (120,120), holes(120,25), ann, 5), ("contacts120_quasar", (120,120), holes(120,25), qua, 5)]:
    md0, B, D, _ = build_2d(L, lam, NA, src, tf, margin, spacing=16.0, fine=2.0); rec0, *_ = analyze(md0, B, D, nround=150, nstart=4, polish=False)
    row = dict(case=name, scalar_bound=rec0['rigorous'], scalar_gap=rec0['gap_rel'])
    for pol in ['unpol', 'TE']:
        mdv, B, D, _ = build_2d_x(L, lam, NA, src, tf, margin, spacing=16.0, fine=2.0, vector=True, pol=pol, n_img=1.7)
        recv, *_ = analyze(mdv, B, D, nround=150, nstart=4, polish=False); row[f'{pol}_bound'] = recv['rigorous']; row[f'{pol}_gap'] = recv['gap_rel']; row[f'{pol}_rank'] = recv['rank']
    print(json.dumps(row), flush=True); res.append(row)
out['vector_2d'] = res
# (d) pixelation
res = []
md, B, D, _ = build_1d(600, lam, NA, conv, [(0,50),(100,150),(200,250),(430,480)], TAU_DUV, spacing=4.0); rF = sdp_bound_real(RealModel(md), B, D)
for N in [150, 60, 30, 15, 12, 9]:
    rp = sdp_bound_pixel_int(md, N, B, D, solver="SCS"); res.append(dict(N=N, pixel_nm=600/N, bound=rp['t'], fourier=rF['t'], F=md.F)); print(json.dumps(res[-1]), flush=True)
out['pixelation'] = res
json.dump(out, open('/home/user/Fable-Research/results/exp16_physics.json','w'), indent=1)
