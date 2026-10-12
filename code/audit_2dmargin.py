# Quantify effect of samples_2d distance bias: exact-distance sampler for square holes vs samples_2d
import numpy as np, warnings, json; warnings.filterwarnings('ignore')
from certlitho import *
lam, NA = 193.0, 1.35
def exact_dist_holes(pts, P, hw):
    a = (pts[:,0] + P/2) % P - P/2; b = (pts[:,1] + P/2) % P - P/2
    dx = np.abs(a) - hw; dy = np.abs(b) - hw
    outside = np.hypot(np.maximum(dx, 0), np.maximum(dy, 0)); inside = -np.minimum(np.maximum(dx, dy), 0)
    ins = (dx < 0) & (dy < 0); return np.where(ins, inside, outside), ins
def samples_2d_exact(P, hw, margin, spacing, fine=2.0, band=0.6):
    nx = int(round(P/fine)); g = (np.arange(nx)+0.5)*P/nx
    G = np.array([(a, b) for a in g for b in g]); d, ins = exact_dist_holes(G, P, hw)
    ring = (d >= margin) & (d < margin+band)
    pts = list(G[ring]); tg = list(ins[ring])
    for a in np.arange(0, P, spacing):
        for b in np.arange(0, P, spacing):
            dd, ii = exact_dist_holes(np.array([[a, b]]), P, hw)
            if dd[0] >= margin+band: pts.append((a, b)); tg.append(ii[0])
    x = np.array(pts); tgt = np.array(tg); return x, np.where(tgt)[0], np.where(~tgt)[0], tgt
holes = lambda P, hw: (lambda a, b: (abs((a+P/2)%P-P/2) < hw) and (abs((b+P/2)%P-P/2) < hw))
ann = annular_source_2d(0.6, 0.9, NA, n=9); qua = quasar_source_2d(0.6, 0.9, NA, n=11); cv = disk_source_2d(0.5, NA, n=5)
for name, P, hw, src, margin in [("contacts160_annular", 160, 35, ann, 12), ("contacts160_conv0.5", 160, 35, cv, 12), ("contacts120_annular", 120, 25, ann, 10), ("contacts120_quasar", 120, 25, qua, 10), ("contacts100_quasar", 100, 22, qua, 8)]:
    md, B, D, tgt = build_2d((P,P), lam, NA, src, holes(P,hw), margin, spacing=16.0, fine=2.0); rm = RealModel(md)
    r = sdp_bound_real(rm, B, D); rig0, _, _ = rigorize_certificate(rm, r['mu'], r['lam'])
    d, _ = exact_dist_holes(md.x, P, hw); dmin = d[np.r_[B, D]].min()
    row = dict(name=name, tau=margin, bound_as_coded=rig0, min_true_dist_as_coded=float(dmin))
    for tau in [margin - 2, margin - 1, margin]:
        x, B2, D2, tgt2 = samples_2d_exact(P, hw, tau, 16.0, fine=2.0)
        md2 = Model((P,P), lam, NA, *src, nx=None, dim=2, xs=x); rm2 = RealModel(md2)
        r2 = sdp_bound_real(rm2, B2, D2); rig2, _, _ = rigorize_certificate(rm2, r2['mu'], r2['lam'])
        row[f'bound_exact_tau{tau}'] = rig2
    print(json.dumps(row), flush=True)
