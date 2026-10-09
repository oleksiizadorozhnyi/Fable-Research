from pipeline import *
lam, NA = 193.0, 1.35
md, B, D, tgt = build_1d(500, lam, NA, ([0.0],[1.0]), [(0,40),(100,180),(240,400)], 8, spacing=4.0); rm = RealModel(md)
out = {}
for sv in ["CLARABEL", "SCS"]:
    t0 = time.time(); las = lasserre2_single_mask(rm, B, D, tmax=1.0, solvers=(sv,)); out[sv] = dict(t=las['t'], status=las['status'], time=time.time()-t0)
    print(sv, out[sv], flush=True)
json.dump(out, open('/home/user/Fable-Research/results/exp12_lasserre_xcheck.json','w'), indent=1)
