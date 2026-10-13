# The list of instances whose certificates are independently verified (exp15: double-precision re-solve + 50-digit
# Cholesky of the computed matrix; exp21: exact 50-digit rebuild of the certificate from the model equations).
from pipeline import *
def verify_instances():
    lam, NA = 193.0, 1.35
    conv = SRC_DUV['conv0.7']; dip = SRC_DUV['dipole0.8']; coh = SRC_DUV['coh']
    for nm, src in [("coh",coh),("conv0.7",conv),("dipole0.8",dip)]:
        for cn, L, segs in [("isodense", 600, [(0,50),(100,150),(200,250),(430,480)]), ("mixed", 500, [(0,40),(100,180),(240,400)])]:
            md, B, D, tgt = build_1d(L, lam, NA, src, segs, TAU_DUV, spacing=4.0); yield f"{cn}_{nm}", md, B, D
    for nm, src in [("conv0.7",conv),("dipole0.8",dip)]:
        md, B, D, tgt = build_1d(600, lam, NA, src, [(0,50),(100,150),(200,250),(430,480)], TAU_DUV, spacing=4.0, focus=(0.0,80.0,-80.0), n_img=NIMG_DUV); yield f"isodense_{nm}_PW80", md, B, D
    md, B, D, tgt = build_1d(500, lam, NA, conv, [(0,40),(100,180),(240,400)], TAU_DUV, spacing=4.0, focus=(0.0,80.0,-80.0), n_img=NIMG_DUV); yield "mixed_conv0.7_PW80", md, B, D
    for nm, src in [("conv0.7",conv),("dipole0.8",dip),("coh",coh)]:
        md, B, D, tgt = build_1d(360, lam, NA, src, [(0,45),(90,135),(180,225),(270,315)], 4.5, spacing=4.0); yield f"dense90_{nm}", md, B, D
    e2 = json.load(open('/home/user/Fable-Research/results/exp2_random_tightness.json'))
    srcs = {"coh": coh, "conv0.5": SRC_DUV['conv0.5'], "conv0.9": SRC_DUV['conv0.9'], "dipole0.7": SRC_DUV['dipole0.7']}
    for rec0 in e2:
        md, B, D, tgt = build_1d(rec0['L'], lam, NA, srcs[rec0['src']], rec0['segs'], TAU_DUV, spacing=5.0); yield f"random_trial{rec0['trial']}", md, B, D
    lamE, NAE = 13.5, 0.33
    for pitch in [20, 22, 24, 26, 28, 32, 40]:
        for sname, src in {"conv0.9": SRC_EUV['conv0.9'], "dipole0.8": SRC_EUV['dipole0.8']}.items():
            segsE = [(i*pitch, i*pitch + pitch/2) for i in range(4)]
            md, B, D, tgt = build_1d(4*pitch, lamE, NAE, src, segsE, 0.05*pitch, spacing=pitch/8, focus=(0.0, 50.0, -50.0)); yield f"EUV_p{pitch}_{sname}", md, B, D
    ann = SRC_DUV['annular']; qua = SRC_DUV['quasar']; cv = SRC_DUV['conv0.5']
    holes = lambda P, hw: rect_target((P, P), [(P/2 - hw, P/2 - hw, 2*hw, 2*hw)])
    Tj = rect_target((320, 320), [(130, 60, 60, 160), (0, 10, 320, 60)]); lineend = rect_target((320, 320), [(75, 30, 50, 180), (195, 270, 50, 200)])
    for name, L, tf, src, margin in [("contacts160_annular", (160,160), holes(160,35), ann, 7), ("contacts160_conv0.5", (160,160), holes(160,35), cv, 7), ("contacts120_annular", (120,120), holes(120,25), ann, 5), ("contacts120_quasar", (120,120), holes(120,25), qua, 5), ("contacts100_quasar", (100,100), holes(100,22), qua, 4.4), ("Tjunction320_annular", (320,320), Tj, ann, 6), ("Tjunction320_quasar", (320,320), Tj, qua, 6), ("lineend320_annular", (320,320), lineend, ann, 6)]:
        md, B, D, tgt = build_2d(L, lam, NA, src, tf, margin, spacing=16.0, fine=2.0); yield name, md, B, D

    # exp17: random 2D layouts, adversarial 1D targets, EUV 2D contacts (instance data read from the result file)
    e17 = json.load(open('/home/user/Fable-Research/results/exp17_dataset.json'))
    srcs2 = {"annular": SRC_DUV['annular'], "quasar": SRC_DUV['quasar'], "conv0.5": SRC_DUV['conv0.5']}
    for r in e17['random2d']:
        Lc = float(r['L']); tf = rect_target((Lc, Lc), r['rects'])
        md, B, D, tgt = build_2d((Lc, Lc), lam, NA, srcs2[r['src']], tf, 5, spacing=16.0, fine=2.0); yield f"random2d_trial{r['trial']}_{r['src']}", md, B, D
    srcs1 = {"coh": coh, "conv0.7": conv, "dipole0.8": dip}
    for r in e17['adversarial1d']:
        md, B, D, tgt = build_1d(r['L'], lam, NA, srcs1[r['src']], [tuple(sg) for sg in r['segs']], r['margin'], spacing=4.0); yield f"adv_{r['name']}_{r['src']}", md, B, D
    holesE = lambda P, hw: rect_target((P, P), [(P/2 - hw, P/2 - hw, 2*hw, 2*hw)])
    for r in e17['euv2d']:
        P = r['pitch']; src = SRC_EUV['annularE'] if r['src'] == 'annular' else SRC_EUV['quasarE']
        md, B, D, tgt = build_2d((P, P), lamE, NAE, src, holesE(P, 0.25*P), 0.05*P, spacing=P/6, fine=P/40, focus=(0.0, 40.0, -40.0)); yield f"euv2d_p{P}_{r['src']}", md, B, D
    # exp20: corner-exclusion variants
    for r in json.load(open('/home/user/Fable-Research/results/exp20_corner_exclusion.json')):
        if r['frac_cd'] == 0.0: continue
        nm = r['name']; cfg = {"contacts160_annular": ((160,160), holes(160,35), ann), "contacts120_quasar": ((120,120), holes(120,25), qua), "contacts100_quasar": ((100,100), holes(100,22), qua), "Tjunction320_annular": ((320,320), Tj, ann), "lineend320_annular": ((320,320), lineend, ann)}[nm]
        md, B, D, tgt = build_2d(cfg[0], lam, NA, cfg[2], cfg[1], r['tau'], spacing=16.0, fine=2.0, corner_exclusion=r['corner_exclusion']); yield f"corner_{nm}_c{r['frac_cd']}", md, B, D
    # exp3 / exp22: PSM instances (Fourier model; the pixel certificates are verified in exp22)
    for name, (L, segs) in {"isodense": (600, [(0,50),(100,150),(200,250),(430,480)]), "mixed": (500, [(0,40),(100,180),(240,400)]), "dense120": (360, [(0,60),(120,180),(240,300)])}.items():
        md, B, D, tgt = build_1d(L, lam, NA, conv, segs, margin=TAU_DUV, spacing=5.0); yield f"psm_{name}", md, B, D
