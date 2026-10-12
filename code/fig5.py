import json, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
R = '/home/user/Fable-Research/results/'
e14 = json.load(open(R+'exp14_bnb.json'))
try: e14 += json.load(open(R+'exp14b_bnb_big.json'))
except Exception: pass
fig, ax = plt.subplots(1, 2, figsize=(12, 4.6), gridspec_kw={'width_ratios': [1.25, 1]})
def short(n):
    n = n.replace('_edgeonly', ' (edge)').replace('_PW80', ' ±80nm').replace('conv0.7', 'σ=0.7').replace('_coh', ' coh').replace('isodense', 'iso-dense').replace('trial', 'random ')
    return n.replace('_', ' ')
names = [r['name'] for r in e14]; x = np.arange(len(e14))
ax[0].bar(x-0.25, [r['sdp_bound'] for r in e14], 0.25, label='SDP bound (any number of exposures)')
ax[0].bar(x, [r['bnb_L'] for r in e14], 0.25, label='B&B certified single-mask lower bound L')
ax[0].bar(x+0.25, [r['bnb_U'] for r in e14], 0.25, label='best single mask found U')
ax[0].set_xticks(x); ax[0].set_xticklabels([short(r['name']) + (f" ({r['bnb_nodes']/1e6:.0f}M nodes)" if r.get('bnb_nodes', 0) >= 1e6 else '') for r in e14], rotation=30, ha='right', fontsize=8); ax[0].set_ylabel('dark/bright ratio ρ'); ax[0].legend(fontsize=7); ax[0].set_title('Certified separation of single- and multi-exposure optima', fontsize=10)
e16 = json.load(open(R+'exp16_physics.json'))
rows = e16['vector_1d']; x = np.arange(len(rows)); labels = [f"{r['case'].replace('isodense','iso-dense').replace('_',' ')}\n{r['src'].replace('conv0.7','σ=0.7').replace('dipole0.8','dipole 0.8')}" for r in rows]
for i, (key, lab) in enumerate([('scalar_bound','scalar'), ('TE_bound','vector TE'), ('unpol_bound','vector unpolarized'), ('TM_bound','vector TM')]):
    ax[1].bar(x + (i-1.5)*0.2, [r[key] for r in rows], 0.2, label=lab)
ax[1].set_xticks(x); ax[1].set_xticklabels(labels, fontsize=8); ax[1].set_ylabel('certified bound ρ'); ax[1].legend(fontsize=7); ax[1].set_title('Vector (NA 1.35, n=1.7) versus scalar bounds', fontsize=10)
plt.tight_layout(); plt.savefig(R+'figs/fig5_bnb_vector.png', dpi=150); print('fig5 written')
