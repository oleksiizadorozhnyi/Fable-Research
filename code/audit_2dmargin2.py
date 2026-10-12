# Exact offset-contour sampler for square holes (periodic): samples ON the true offset curves at distance tau (step ~2nm) + interior grid
import numpy as np, warnings, json; warnings.filterwarnings('ignore')
from certlitho import *
lam, NA = 193.0, 1.35
def contour_square(hw, tau, step):
    pts = []
    # inside: square half-width hw - tau
    hi = hw - tau
    if hi > 0:
        t = np.arange(-hi, hi, step)
        for v in t: pts += [(v, hi), (v, -hi), (hi, v), (-hi, v)]
    inner = [(a, b, True) for a, b in pts]
    # outside: rounded square, sides at hw+tau, corners quarter circles radius tau centered at (+-hw, +-hw)
    pts = []; ho = hw + tau; t = np.arange(-hw, hw, step)
    for v in t: pts += [(v, ho), (v, -ho), (ho, v), (-ho, v)]
    nth = max(3, int(np.ceil(0.5*np.pi*tau/step)))
    for cx in (-hw, hw):
        for cy in (-hw, hw):
            for th in np.linspace(0, np.pi/2, nth, endpoint=False):
                sx = np.sign(cx); sy = np.sign(cy); pts.append((cx + sx*tau*np.cos(th), cy + sy*tau*np.sin(th)))
    return inner + [(a, b, False) for a, b in pts]
def exact_dist_holes(pts, P, hw):
    a = (pts[:,0] + P/2) % P - P/2; b = (pts[:,1] + P/2) % P - P/2
    dx = np.abs(a) - hw; dy = np.abs(b) - hw
    outside = np.hypot(np.maximum(dx, 0), np.maximum(dy, 0)); inside = -np.minimum(np.maximum(dx, dy), 0)
    ins = (dx < 0) & (dy < 0); return np.where(ins, inside, outside), ins
def samples_exact(P, hw, tau, spacing=16.0, step=2.0):
    c = contour_square(hw, tau, step); x = np.array([(a % P, b % P) for a, b, _ in c]); tg = [i for _, _, i in c]
    pts = list(x); 
    for a in np.arange(0, P, spacing):
        for b in np.arange(0, P, spacing):
            d, ii = exact_dist_holes(np.array([[a, b]]), P, hw)
            if d[0] >= tau + 0.6: pts.append((a, b)); tg.append(bool(ii[0]))
    x = np.array(pts); tgt = np.array(tg); d, ins = exact_dist_holes(x, P, hw); assert np.all(ins == tgt) and d.min() >= tau - 1e-9, (d.min(), tau)
    return x, np.where(tgt)[0], np.where(~tgt)[0]
holes = lambda P, hw: (lambda a, b: (abs((a+P/2)%P-P/2) < hw) and (abs((b+P/2)%P-P/2) < hw))
ann = annular_source_2d(0.6, 0.9, NA, n=9); qua = quasar_source_2d(0.6, 0.9, NA, n=11); cv = disk_source_2d(0.5, NA, n=5)
for name, P, hw, src, margin in [("contacts160_annular", 160, 35, ann, 12), ("contacts120_quasar", 120, 25, qua, 10), ("contacts100_quasar", 100, 22, qua, 8), ("contacts160_conv0.5", 160, 35, cv, 12)]:
    md, B, D, tgt = build_2d((P,P), lam, NA, src, holes(P,hw), margin, spacing=16.0, fine=2.0); rm = RealModel(md)
    r = sdp_bound_real(rm, B, D); rig0, _, _ = rigorize_certificate(rm, r['mu'], r['lam'])
    row = dict(name=name, tau=margin, bound_as_coded=rig0)
    for tau in [margin-2, margin-1.5, margin-1, margin]:
        x, B2, D2 = samples_exact(P, hw, tau)
        md2 = Model((P,P), lam, NA, *src, nx=None, dim=2, xs=x); rm2 = RealModel(md2)
        r2 = sdp_bound_real(rm2, B2, D2); rig2, _, _ = rigorize_certificate(rm2, r2['mu'], r2['lam'])
        row[f'exact_tau{tau}'] = round(rig2, 5)
    print(json.dumps(row), flush=True)
