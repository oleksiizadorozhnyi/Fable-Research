from pipeline import *
lam, NA = 193.0, 1.35; ann = SRC_DUV['annular']; qua = SRC_DUV['quasar']
Tj = lambda a, b: (abs(a-160) < 30 and 60 < b < 220) or (abs(b-40) < 30)
lineend = lambda a, b: (abs(a-100) < 25 and 30 < b < 210) or (abs(a-220) < 25 and (b < 150 or b > 270))
e4 = json.load(open('/home/user/Fable-Research/results/exp4_2d.json'))
for name, tf, src in [("Tjunction320_quasar", Tj, qua), ("lineend320_annular", lineend, ann)]:
    md, B, D, tgt = build_2d((320,320), lam, NA, src, tf, 6, spacing=16.0, fine=2.0); rm = RealModel(md)
    r = sdp_bound_real(rm, B, D); rig, me, c = rigorize_certificate(rm, r['mu'], r['lam'])
    rec = [x for x in e4 if x['name'] == name][0]
    rec.update(t_sdp=r['t'], rigorous=rig, cert_mineig=float(me), cert_scale=float(c), cert_null_dim=rigorize_certificate.last_null_dim, cert_resid_null=rigorize_certificate.last_resid_null,
               gap_rel=(rec['best_feasible']-rig)/rig, EL_bound=(1-rig)/(1+rig), solver=r['solver'], status=r['status'])
    print(json.dumps(dict(name=name, t=r['t'], rigorous=rig, c=c, null_dim=rigorize_certificate.last_null_dim, resid_null=rigorize_certificate.last_resid_null, best=rec['best_feasible'], gap=rec['gap_rel'], solver=r['solver'])), flush=True)
    json.dump(e4, open('/home/user/Fable-Research/results/exp4_2d.json','w'), indent=1)
conv = SRC_DUV['conv0.7']; md, B, D, tgt = build_1d(600, lam, NA, conv, [(0,50),(100,150),(200,250),(430,480)], 8, spacing=4.0); rm = RealModel(md)
r = sdp_bound_real(rm, B, D); print("1D isodense rigorous:", rigorize_certificate(rm, r['mu'], r['lam']), "null_dim", rigorize_certificate.last_null_dim, flush=True)
