from pipeline import *
lam, NA = 193.0, 1.35; rng = np.random.default_rng(42); out = []
srcs = {"coh": SRC_DUV['coh'], "conv0.5": SRC_DUV['conv0.5'], "conv0.9": SRC_DUV['conv0.9'], "dipole0.7": SRC_DUV['dipole0.7']}
for trial in range(120):
    L = float(rng.choice([300, 400, 500, 600])); nseg = int(rng.integers(1, 5))
    cuts = np.sort(rng.uniform(0, L, 2*nseg)); segs = [(cuts[2*i], cuts[2*i+1]) for i in range(nseg) if cuts[2*i+1]-cuts[2*i] >= 40]
    if not segs: continue
    sname = str(rng.choice(list(srcs))); s, w = srcs[sname]
    md, B, D, tgt = build_1d(L, lam, NA, (s, w), segs, margin=TAU_DUV, spacing=5.0)
    if len(B) < 2 or len(D) < 2: continue
    try: rec, rm, r, rb = analyze(md, B, D, nround=200, nstart=8, save_prefix=f'/home/user/Fable-Research/results/masks/exp2_trial{trial}')
    except Exception as e: print("fail", trial, e); continue
    rec.update(trial=trial, L=L, src=sname, segs=[(float(a), float(b)) for a, b in segs], eig_ratio=float(r['eig'][1]/r['eig'][0]))
    print(json.dumps(rec), flush=True); out.append(rec)
    json.dump(out, open('/home/user/Fable-Research/results/exp2_random_tightness.json', 'w'), indent=1)
