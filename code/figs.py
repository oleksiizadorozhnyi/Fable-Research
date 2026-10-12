from pipeline import *
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
lam, NA = 193.0, 1.35
# Fig 1: certificate illustration on isodense conv0.7: target, best image, dual weights
conv = SRC_DUV['conv0.7']; s, w = conv; L = 600; segs = [(0,50),(100,150),(200,250),(430,480)]
md, B, D, tgt = build_1d(L, lam, NA, conv, segs, TAU_DUV, spacing=2.0); rm = RealModel(md)
r = sdp_bound_real(rm, B, D); rig, _, _ = rigorize_certificate(rm, r['mu'], r['lam'])
rr, rv = round_from_Y(rm, r['X'], B, D, nsamp=300); rl, rvl = local_refine_real(rm, rv, B, D)
I = rm.image(rvl); I = I / I[B].min()
x = md.x[:, 0]
fig, ax = plt.subplots(3, 1, figsize=(7, 7), sharex=True)
ax[0].fill_between(x, 0, tgt.astype(float), step='mid', alpha=0.3, color='k', label='target (bright)'); ax[0].plot(x, I, 'b', label='best mask image / min bright')
ax[0].axhline(1, color='g', ls='--', lw=0.8); ax[0].axhline(rig, color='r', ls='--', lw=0.8, label=f'certified bound ratio {rig:.3f}'); ax[0].legend(fontsize=8); ax[0].set_ylabel('intensity')
mu = np.zeros(md.nx); lamv = np.zeros(md.nx)
for (z, ix), v in r['mu'].items(): mu[ix] += v
for (z, ix), v in r['lam'].items(): lamv[ix] += v
ax[1].stem(x[D], mu[D] / mu.sum(), 'r', markerfmt='ro', basefmt=' ', label='dual weights mu (dark)'); ax[1].stem(x[B], lamv[B] / mu.sum(), 'g', markerfmt='go', basefmt=' ', label='dual weights lambda (bright)'); ax[1].legend(fontsize=8); ax[1].set_ylabel('certificate weights')
m = md.pixels_from_orders(rm.to_orders(rvl), 600); ax[2].plot(np.arange(600) * L / 600, m / np.abs(m).max(), 'k'); ax[2].set_ylabel('bandlimited mask amplitude'); ax[2].set_xlabel('x (nm)')
plt.tight_layout(); plt.savefig('/home/user/Fable-Research/results/figs/fig1_certificate_isodense.png', dpi=150); plt.close()
# Fig 2: EUV EL bound vs pitch (from exp5) + DUV k1 sweep
e5 = json.load(open('/home/user/Fable-Research/results/exp5_smo_euv_toeplitz.json'))
fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))
for sname, mk in [("conv0.9", 'o-'), ("dipole0.8", 's-')]:
    rows = [r for r in e5['euv_pw'] if r['src'] == sname]
    ax[0].plot([r['k1'] for r in rows], [r['EL_full_bound'] for r in rows], mk, label=f'{sname}: certified max EL')
    ax[0].plot([r['k1'] for r in rows], [r['EL_full_best'] for r in rows], mk, mfc='none', ls=':', label=f'{sname}: best single mask')
ax[0].axvline(0.25, color='k', ls='--', lw=0.8); ax[0].set_xlabel('k1 = half-pitch·NA/λ (EUV dense lines, ±50 nm focus, τ = 10% CD)'); ax[0].set_ylabel('full-range exposure latitude'); ax[0].legend(fontsize=7)
sig = 0.7; s7, w7 = disk_source(sig, NA, n=21); ks = []; els = []
for pitch in np.arange(80, 240, 8):
    mdp, Bp, Dp, _ = build_1d(float(pitch), lam, NA, (s7, w7), [(pitch/4, 3*pitch/4)], 0.05*pitch, spacing=pitch/16)
    rp = sdp_bound_real(RealModel(mdp), Bp, Dp); ks.append(pitch / 2 * NA / lam); els.append(el_full(rp['t']))
ax[1].plot(ks, els, 'o-'); ax[1].axvline(0.5 / (1 + sig), color='k', ls='--', lw=0.8, label='k1 = 1/(2(1+σ))'); ax[1].set_xlabel('k1 (DUV, σ=0.7 conventional)'); ax[1].set_ylabel('certified max exposure latitude'); ax[1].legend(fontsize=8)
plt.tight_layout(); plt.savefig('/home/user/Fable-Research/results/figs/fig2_k1_limits.png', dpi=150); plt.close()
# Fig 3: gap histogram exp2
e2 = json.load(open('/home/user/Fable-Research/results/exp2_random_tightness.json')); gaps = np.array([min(r['gap_rel'], 1.0) for r in e2])
plt.figure(figsize=(5, 3.5)); plt.hist(np.log10(np.maximum(gaps, 1e-8)), bins=30); plt.xlabel('log10 relative gap (best single mask vs certified multi-exposure bound)'); plt.ylabel('# random 1D instances'); plt.tight_layout()
plt.savefig('/home/user/Fable-Research/results/figs/fig3_gap_hist.png', dpi=150); plt.close()
print("figs written")
