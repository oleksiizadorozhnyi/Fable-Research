# Quantify immersion defocus mis-modelling: recompute PW80 bounds with the immersion formula (n=1.44) by rebuilding C
import numpy as np, warnings, json; warnings.filterwarnings('ignore')
from certlitho import *
lam, NA = 193.0, 1.35; n_imm = 1.44
def rebuild_C(md, n):
    for z in md.focus:
        Cz = []
        for s in md.src:
            g = md.f + s/lam; gn = np.linalg.norm(g, axis=1); inside = (gn*lam <= NA + 1e-12)
            if n is None: gl = np.clip(gn*lam, 0, 1); ph = 2*np.pi/lam*z*(np.sqrt(1-gl**2)-1)
            else: ph = 2*np.pi*n/lam*z*(np.sqrt(1-(gn*lam/n)**2)-1)
            E = np.exp(2j*np.pi*(md.x @ md.f.T)); Cz.append(E*(inside*np.exp(1j*ph))[None, :])
        md.C[z] = np.stack(Cz)
conv = conventional_source_1d(0.7, NA, 7); dip = dipole_source_1d(0.8, NA)
for name, L, segs, src in [("isodense_conv0.7_PW80", 600, [(0,50),(100,150),(200,250),(430,480)], conv), ("mixed_conv0.7_PW80", 500, [(0,40),(100,180),(240,400)], conv), ("isodense_dipole0.8_PW80", 600, [(0,50),(100,150),(200,250),(430,480)], dip)]:
    row = dict(name=name)
    for label, n in [("as_coded(n=1,clip)", None), ("immersion_n1.44", n_imm), ("immersion_n1.44_z80*1.44", 'scaled')]:
        md, B, D, tgt = build_1d(L, lam, NA, src, segs, 8, spacing=4.0, focus=(0.0, 80.0, -80.0))
        if n == 'scaled':
            md, B, D, tgt = build_1d(L, lam, NA, src, segs, 8, spacing=4.0, focus=(0.0, 80.0*n_imm, -80.0*n_imm)); rebuild_C(md, n_imm)
        elif n is not None: rebuild_C(md, n)
        rm = RealModel(md); r = sdp_bound_real(rm, B, D); rig, _, _ = rigorize_certificate(rm, r['mu'], r['lam']); row[label] = round(rig, 4)
    print(json.dumps(row), flush=True)
# B&B 'certified gain' definition: ensemble ratio of the SDP primal X* vs t and rig (trial82)
e2 = json.load(open('/home/user/Fable-Research/results/exp2_random_tightness.json')); rec0 = [r for r in e2 if r['trial'] == 82][0]
md, B, D, tgt = build_1d(rec0['L'], lam, NA, ([0.0],[1.0]), rec0['segs'], 10, spacing=5.0); rm = RealModel(md); r = sdp_bound_real(rm, B, D); rig, _, _ = rigorize_certificate(rm, r['mu'], r['lam'])
X = r['X']; w, V = np.linalg.eigh(X); w = np.clip(w, 0, None); X = (V*w) @ V.T
ens = max(np.trace(rm.Q(d) @ X) for d in D)/min(np.trace(rm.Q(b) @ X) for b in B)
print(json.dumps(dict(trial82_t=r['t'], rig=rig, ensemble_ratio_of_Xstar=ens, bnb_L=0.4061467620929617, gain_vs_rig=(0.4061467620929617-rig)/rig, gain_vs_ensemble=(0.4061467620929617-ens)/ens)))
