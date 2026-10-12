# 2D sampling bias at the NEW tolerances (exp4/exp17 session 3): true distance of samples; exact-contour bound vs as-coded
import numpy as np, warnings, json; warnings.filterwarnings('ignore')
from certlitho import *; from pipeline import *
exec(open('audit_2dmargin2.py').read().split("holes = lambda")[0].split("lam, NA = 193.0, 1.35")[1])   # reuse contour_square, exact_dist_holes, samples_exact
holes = lambda P, hw: (lambda a, b: (abs((a+P/2)%P-P/2) < hw) and (abs((b+P/2)%P-P/2) < hw))
rows = []
for name, P, hw, src, tau, fine, lam, NA, n, focus in [("contacts160_annular", 160, 35, SRC_DUV['annular'], 7, 2.0, 193.0, 1.35, 1.44, (0.0,)), ("contacts120_quasar", 120, 25, SRC_DUV['quasar'], 5, 2.0, 193.0, 1.35, 1.44, (0.0,)), ("contacts100_quasar", 100, 22, SRC_DUV['quasar'], 4.4, 2.0, 193.0, 1.35, 1.44, (0.0,)), ("EUV_contacts40_annular", 40, 10, SRC_EUV['annularE'], 2.0, 1.0, 13.5, 0.33, 1.0, (0.0, 40.0, -40.0))]:
    x, B, D, tgt, sd = samples_2d((P,P), holes(P,hw), tau, 16.0 if P > 60 else P/6, fine=fine)
    d, ins = exact_dist_holes(x, P, hw)
    md = Model((P,P), lam, NA, *src, nx=None, dim=2, xs=x, n_img=n, focus=focus); rm = RealModel(md); r = sdp_bound_real(rm, B, D); rig0, _, _ = rigorize_certificate(rm, r['mu'], r['lam'])
    row = dict(name=name, tau=tau, fine=fine, nB=len(B), nD=len(D), min_true_dist_B=round(float(d[B].min()),3), min_true_dist_D=round(float(d[D].min()),3), frac_closer_than_tau=round(float(np.mean(d < tau-1e-9)),3), bound_as_coded=round(rig0,5))
    for t2 in sorted(set([round(d[np.r_[B,D]].min(),2), tau])):
        xe, B2, D2 = samples_exact(P, hw, t2, spacing=16.0 if P > 60 else P/6, step=min(2.0, fine))
        md2 = Model((P,P), lam, NA, *src, nx=None, dim=2, xs=xe, n_img=n, focus=focus); rm2 = RealModel(md2); r2 = sdp_bound_real(rm2, B2, D2); rig2, _, _ = rigorize_certificate(rm2, r2['mu'], r2['lam'])
        row[f'bound_exact_contour_tau{t2}'] = round(rig2, 5)
    print(json.dumps(row), flush=True)
