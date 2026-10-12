from pipeline import *
lam, NA = 193.0, 1.35; out = []
def run(name, L, tf, src, nx=20, margin=12):
    md, B, D, tgt = build_2d(L, lam, NA, src, tf, margin, spacing=16.0, fine=2.0)
    print(name, 'orders', md.F, 'samples', md.nx, 'B', len(B), 'D', len(D), flush=True)
    rec, rm, r, rb = analyze(md, B, D, nround=200, nstart=6, polish=False, save_prefix=f'/home/user/Fable-Research/results/exp4_{name}')
    rec.update(name=name, L=list(L), nsrc=int(len(src[1]))); print(json.dumps(rec), flush=True); out.append(rec)
    json.dump(out, open('/home/user/Fable-Research/results/exp4_2d.json', 'w'), indent=1)
ann = SRC_DUV['annular']; qua = SRC_DUV['quasar']; conv = SRC_DUV['conv0.5']
holes = lambda P, hw: (lambda a, b: (abs((a+P/2)%P-P/2) < hw) and (abs((b+P/2)%P-P/2) < hw))
run("contacts160_annular", (160,160), holes(160, 35), ann, margin=7)
run("contacts160_conv0.5", (160,160), holes(160, 35), conv, margin=7)
run("contacts120_annular", (120,120), holes(120, 25), ann, nx=16, margin=5)
run("contacts120_quasar", (120,120), holes(120, 25), qua, nx=16, margin=5)
run("contacts100_quasar", (100,100), holes(100, 22), qua, nx=16, margin=4.4)
Tj = lambda a, b: (abs(a-160) < 30 and 60 < b < 220) or (abs(b-40) < 30)
run("Tjunction320_annular", (320,320), Tj, ann, nx=24, margin=6)
run("Tjunction320_quasar", (320,320), Tj, qua, nx=24, margin=6)
lineend = lambda a, b: (abs(a-100) < 25 and 30 < b < 210) or (abs(a-220) < 25 and (b < 150 or b > 270))
run("lineend320_annular", (320,320), lineend, ann, nx=24, margin=6)
