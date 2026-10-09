from pipeline import *
lam, NA = 193.0, 1.35; out = []
def run(name, L, segs, src, focus=(0.0,), margin=8, nx=None):
    md, B, D, tgt = build_1d(L, lam, NA, src, segs, margin, spacing=4.0, focus=focus)
    rec, rm, r, rb = analyze(md, B, D, save_prefix=f'/home/user/Fable-Research/results/exp1_{name}')
    rec.update(name=name, L=L, segs=segs); print(json.dumps(rec), flush=True); out.append(rec)
    json.dump(out, open('/home/user/Fable-Research/results/exp1_1d_features.json', 'w'), indent=1)
conv = conventional_source_1d(0.7, NA, 7); dip = dipole_source_1d(0.8, NA); coh = ([0.0], [1.0])
segs_isodense = [(0,50),(100,150),(200,250),(430,480)]; segs_mixed = [(0,40),(100,180),(240,400)]
for nm, src in [("coh",coh),("conv0.7",conv),("dipole0.8",dip)]:
    run(f"isodense_{nm}", 600, segs_isodense, src); run(f"mixed_{nm}", 500, segs_mixed, src)
run("isodense_conv0.7_PW80", 600, segs_isodense, conv, focus=(0.0, 80.0, -80.0))
run("mixed_conv0.7_PW80", 500, segs_mixed, conv, focus=(0.0, 80.0, -80.0))
run("isodense_dipole0.8_PW80", 600, segs_isodense, dip, focus=(0.0, 80.0, -80.0))
run("dense90_conv0.7", 360, [(0,45),(90,135),(180,225),(270,315)], conv, margin=6)
run("dense90_dipole0.8", 360, [(0,45),(90,135),(180,225),(270,315)], dip, margin=6)
run("dense90_coh", 360, [(0,45),(90,135),(180,225),(270,315)], coh, margin=6)
