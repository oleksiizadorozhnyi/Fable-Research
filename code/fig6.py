import json, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
R = '/home/user/Fable-Research/results/'; e17 = json.load(open(R+'exp17_dataset.json'))
fig, ax = plt.subplots(1, 3, figsize=(14, 3.8))
r2 = e17['random2d']
if r2:
    ax[0].scatter([r['F'] for r in r2], [100*max(r['gap_rel'], 1e-6) for r in r2], c=[r['rank'] for r in r2], cmap='viridis', s=25)
    ax[0].set_yscale('log'); ax[0].set_xlabel('collected orders F'); ax[0].set_ylabel('gap best single mask vs bound (%)'); ax[0].set_title(f'random 2D layouts (n={len(r2)}); color = rank of SDP optimum')
a1 = e17['adversarial1d']
if a1:
    names = sorted(set(r['name'] for r in a1)); x = np.arange(len(names)); srcs = ['coh', 'conv0.7', 'dipole0.8']
    for i, sname in enumerate(srcs):
        vals = [next((r['EL_full_bound'] for r in a1 if r['name']==n and r['src']==sname), np.nan) for n in names]
        ax[1].bar(x + (i-1)*0.27, vals, 0.27, label=sname)
    ax[1].set_xticks(x); ax[1].set_xticklabels([n.replace('_','\n') for n in names], fontsize=6); ax[1].set_ylabel('certified max exposure latitude'); ax[1].legend(fontsize=7); ax[1].set_title('adversarial 1D targets')
e2 = e17['euv2d']
if e2:
    for sname, mk in [('annular','o-'), ('quasar','s-')]:
        rows = sorted([r for r in e2 if r['src']==sname], key=lambda r: r['k1'])
        ax[2].plot([r['k1'] for r in rows], [r['EL_full_bound'] for r in rows], mk, label=f'{sname}: certified max EL'); ax[2].plot([r['k1'] for r in rows], [r['EL_full_best'] for r in rows], mk, mfc='none', ls=':', label=f'{sname}: best single mask')
    ax[2].set_xlabel('k1 (EUV contacts, ±40 nm focus)'); ax[2].set_ylabel('exposure latitude'); ax[2].legend(fontsize=7); ax[2].set_title('EUV 2D contact arrays')
plt.tight_layout(); plt.savefig(R+'figs/fig6_dataset.png', dpi=150); print('fig6 written')
