from pipeline import *
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
lam, NA = 193.0, 1.35
e4 = json.load(open('/home/user/Fable-Research/results/exp4_2d.json'))
ann = SRC_DUV['annular']
Tj = rect_target((320, 320), [(130, 60, 60, 160), (0, 10, 320, 60)])
name = "Tjunction320_annular"
if any(r['name'] == name for r in e4):
    md, B, D, tgt = build_2d((320,320), lam, NA, ann, Tj, 6, spacing=16.0, fine=2.0); rm = RealModel(md)
    import os
    rb = np.load(f'/home/user/Fable-Research/results/exp13_{name}_rbest.npy') if os.path.exists(f'/home/user/Fable-Research/results/exp13_{name}_rbest.npy') else np.load(f'/home/user/Fable-Research/results/exp4_{name}_rbest.npy'); dual = json.load(open(f'/home/user/Fable-Research/results/exp4_{name}_dual.json'))
    M = rm.to_orders(rb); N = (128, 128)
    mask = md.pixels_from_orders(M, N)
    gx, gy = np.meshgrid(np.arange(N[0])/N[0]*320, np.arange(N[1])/N[1]*320, indexing='ij')
    pts = np.stack([gx.ravel(), gy.ravel()], 1)
    mdf = Model((320,320), lam, NA, ann[0], ann[1], nx=None, dim=2, xs=pts)
    I = mdf.image(M).reshape(N); Ib = rm.image(rb); I = I / Ib[B].min()
    fig, ax = plt.subplots(1, 4, figsize=(16, 4))
    T = np.array([[Tj(a, b) for b in np.arange(N[1])/N[1]*320] for a in np.arange(N[0])/N[0]*320])
    ax[0].imshow(T.T, origin='lower', extent=[0,320,0,320], cmap='gray'); ax[0].set_title('target')
    ax[1].imshow(mask.T, origin='lower', extent=[0,320,0,320], cmap='RdBu'); ax[1].set_title('best bandlimited mask (PSM)')
    im = ax[2].imshow(I.T, origin='lower', extent=[0,320,0,320], cmap='magma'); ax[2].contour(gx, gy, I, levels=[1.0], colors='c'); plt.colorbar(im, ax=ax[2]); ax[2].set_title('image / min bright; cyan: I=1')
    mu = np.zeros(md.nx); la = np.zeros(md.nx)
    import re
    for k, v in dual['mu'].items(): mu[int(re.findall(r'\d+', k)[-1])] += v
    for k, v in dual['lam'].items(): la[int(re.findall(r'\d+', k)[-1])] += v
    ax[3].scatter(md.x[D,0], md.x[D,1], c=mu[D]/mu.sum(), s=12, cmap='Reds', label='mu (dark)'); ax[3].scatter(md.x[B,0], md.x[B,1], c=la[B]/mu.sum(), s=12, cmap='Greens', marker='s', label='lambda (bright)')
    ax[3].set_xlim(0,320); ax[3].set_ylim(0,320); ax[3].set_aspect('equal'); ax[3].set_title('certificate weights'); ax[3].legend(fontsize=7)
    rec = [r for r in e4 if r['name'] == name][0]; vb = rm.ratio(rb, B, D)
    fig.suptitle(f"{name}: certified bound ρ={rec['rigorous']:.4f} (max full-range EL {rec['EL_full_bound']:.3f}); best single mask ρ={vb:.4f}; gap {(vb-rec['rigorous'])/rec['rigorous']:.2%}")
    plt.tight_layout(); plt.savefig('/home/user/Fable-Research/results/figs/fig4_2d_Tjunction.png', dpi=130); print("fig4 written")
