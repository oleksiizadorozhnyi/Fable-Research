# PSM (real) versus nonnegative masks on a pixel-integrated 5 nm grid (N = L/5 pixels), and a pattern-error ILT.
# All bounds are RIGORIZED certificates (Fourier PSM: rigorize_certificate v4; pixel PSM and DNN: rigorize_pixel_certificate)
# and the duals are saved for the exact 50-digit verification in exp22.
from pipeline import *
lam, NA = 193.0, 1.35; out = []
conv = SRC_DUV['conv0.7']
cases = {"isodense": (600, [(0,50),(100,150),(200,250),(430,480)]), "mixed": (500, [(0,40),(100,180),(240,400)]), "dense120": (360, [(0,60),(120,180),(240,300)])}
for name, (L, segs) in cases.items():
    md, B, D, tgt = build_1d(L, lam, NA, conv, segs, margin=TAU_DUV, spacing=5.0); rm = RealModel(md); N = int(round(L / 5.0))
    tf = lambda v: any(a <= (v % L) < b for a, b in segs)
    rF = sdp_bound_real(rm, B, D); rigF, _, cF = rigorize_certificate(rm, rF['mu'], rF['lam'])
    rP = sdp_bound_pixel_int(md, N, B, D, nonneg=False, solver="CLARABEL"); rigP, _, cP = rigorize_pixel_certificate(md, N, B, D, rP['mu'], rP['lam'])
    rN = sdp_bound_pixel_int(md, N, B, D, nonneg=True, solver="CLARABEL"); dn = rigorize_dnn_certificate(md, N, B, D, rN['mu'], rN['lam'], rN['Nmat']); rigN, cN = dn['bound'], dn['kappa']
    tgt_pix = tgt.astype(float)
    rP_ilt, mP, MP, hP = ilt_gradient(md, tgt_pix, N, B, D, steps=800, lr=0.05, pixel_bounds=(-1,1))
    rB_ilt, mB, MB, hB = ilt_gradient(md, tgt_pix, N, B, D, steps=800, lr=0.05, pixel_bounds=(0,1))
    rP_ref, _ = local_refine_real(rm, rm.from_orders(MP), B, D)
    rec0, _, _, _ = analyze(md, B, D, nround=200, nstart=8, save_prefix=f'/home/user/Fable-Research/results/masks/exp3_{name}')
    rec = dict(name=name, F=md.F, nB=len(B), nD=len(D), N_pixels=N, pixel_nm=L/N, t_fourier_psm=rF['t'], rig_fourier_psm=rigF, rank_fourier=rF['rank_eps'], t_pixel_psm=rP['t'], rig_pixel_psm=rigP, rank_pixel=rP['rank_eps'],
               t_pixel_nonneg_DNN=rN['t'], rig_pixel_nonneg_DNN=rigN, dnn_dual_infeasibility=dn['delta'], dnn_lamQ=dn['lamQ'], rank_dnn=rN['rank_eps'], status_dnn=rN['status'], ilt_psm_ratio=rP_ilt, ilt_psm_refined_on_ratio=rP_ref, ilt_binary_ratio=rB_ilt,
               best_psm_pipeline=rec0['best_feasible'], gap_psm=(rec0['best_feasible'] - rigF) / rigF)
    print(json.dumps(rec), flush=True); out.append(rec)
    np.save(f'/home/user/Fable-Research/results/exp3_{name}_ilt_psm.npy', mP); np.save(f'/home/user/Fable-Research/results/exp3_{name}_ilt_bin.npy', mB)
    json.dump(dict(F=dict(mu={str(k): v for k, v in rF['mu'].items()}, lam={str(k): v for k, v in rF['lam'].items()}, kappa=cF),
                   P=dict(mu={str(k): v for k, v in rP['mu'].items()}, lam={str(k): v for k, v in rP['lam'].items()}, kappa=cP),
                   NN=dict(mu={str(k): v for k, v in rN['mu'].items()}, lam={str(k): v for k, v in rN['lam'].items()}, kappa=cN), N=N),
              open(f'/home/user/Fable-Research/results/exp3_{name}_duals.json', 'w'))
    np.save(f'/home/user/Fable-Research/results/exp3_{name}_dnn_Nmat.npy', rN['Nmat'])
    json.dump(out, open('/home/user/Fable-Research/results/exp3_pixel_binary.json', 'w'), indent=1)
