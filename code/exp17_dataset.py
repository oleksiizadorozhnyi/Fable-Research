# Expanded dataset: (a) random 2D rectangle layouts in small cells; (b) adversarial 1D instances (near-cutoff mixed pitches,
# SRAF-like narrow spaces, extreme duty cycles); (c) EUV 2D contacts. Records certified bound, best single mask, gap, rank.
from pipeline import *
lam, NA = 193.0, 1.35; rng = np.random.default_rng(7); out = {'random2d': [], 'adversarial1d': [], 'euv2d': []}
ann = annular_source_2d(0.6, 0.9, NA, n=9); qua = quasar_source_2d(0.6, 0.9, NA, n=11); cv = disk_source_2d(0.5, NA, n=5)
srcs2 = {'annular': ann, 'quasar': qua, 'conv0.5': cv}
def save(): json.dump(out, open('/home/user/Fable-Research/results/exp17_dataset.json','w'), indent=1)
# (a) random 2D: cell 200x200 or 240x240, 1-3 random rectangles (>= 50 nm sides), margin 12
for trial in range(40):
    Lc = float(rng.choice([200, 240])); nrect = int(rng.integers(1, 4)); rects = []
    for _ in range(nrect):
        w_ = rng.uniform(50, 0.6*Lc); h_ = rng.uniform(50, 0.6*Lc); x0 = rng.uniform(0, Lc); y0 = rng.uniform(0, Lc); rects.append((x0, y0, w_, h_))
    tf = lambda a, b, rects=rects, Lc=Lc: any(((a - x0) % Lc) < w_ and ((b - y0) % Lc) < h_ for x0, y0, w_, h_ in rects)
    sname = str(rng.choice(list(srcs2)))
    try:
        md, B, D, tgt = build_2d((Lc, Lc), lam, NA, srcs2[sname], tf, 12, spacing=16.0, fine=2.0)
        if len(B) < 4 or len(D) < 4: continue
        rec, rm, r, rb = analyze(md, B, D, nround=150, nstart=4, polish=(md.F <= 21))
        rec.update(trial=trial, L=Lc, src=sname, rects=[[float(v) for v in rr] for rr in rects]); out['random2d'].append(rec); print(json.dumps({k: rec[k] for k in ['trial','src','F','nB','nD','rigorous','best_feasible','gap_rel','rank','sdp_time']}), flush=True); save()
    except Exception as e: print("fail2d", trial, e, flush=True)
# (b) adversarial 1D
conv = conventional_source_1d(0.7, NA, 7); dip = dipole_source_1d(0.8, NA); coh = ([0.0],[1.0])
adv = [("nearcutoff_mixed_p88_p130", 436, [(0,44),(88,132),(176,220),(306,371)], 6),      # dense pitch 88 (just above 84.1 cutoff) + wider feature
       ("sraf_like_narrow_space", 500, [(0,120),(140,160),(180,300)], 8),                  # 20 nm space-feature-space
       ("extreme_duty_5pct", 400, [(0,20)], 6),                                            # isolated 20 nm line
       ("extreme_duty_95pct", 400, [(0,380)], 6),                                          # isolated 20 nm space
       ("alternating_widths", 600, [(0,30),(60,120),(150,180),(210,330),(360,390)], 6),
       ("three_pitches", 720, [(0,45),(90,135),(180,240),(300,360),(420,500),(580,660)], 8)]
for name, L, segs, margin in adv:
    for sname, src in [("coh", coh), ("conv0.7", conv), ("dipole0.8", dip)]:
        try:
            md, B, D, tgt = build_1d(L, lam, NA, src, segs, margin, spacing=4.0)
            if len(B) < 1 or len(D) < 1: continue
            rec, rm, r, rb = analyze(md, B, D, nround=200, nstart=10); rec.update(name=name, src=sname, L=L, segs=segs, margin=margin)
            out['adversarial1d'].append(rec); print(json.dumps({k: rec[k] for k in ['name','src','F','rigorous','best_feasible','gap_rel','rank','EL_bound']}), flush=True); save()
        except Exception as e: print("fail1d", name, sname, e, flush=True)
# (c) EUV 2D contacts: lam 13.5, NA 0.33, pitch 40/36/32 nm, holes 0.5 pitch, annular 0.5-0.9 and quasar, margin 0.1 pitch
lamE, NAE = 13.5, 0.33
annE = annular_source_2d(0.5, 0.9, NAE, n=9); quaE = quasar_source_2d(0.5, 0.9, NAE, n=11)
holes = lambda P, hw: (lambda a, b: (abs((a+P/2)%P-P/2) < hw) and (abs((b+P/2)%P-P/2) < hw))
for P in [44, 40, 36, 32]:
    for sname, src in [("annular", annE), ("quasar", quaE)]:
        try:
            md, B, D, tgt = build_2d((P, P), lamE, NAE, src, holes(P, 0.25*P), 0.1*P, spacing=P/6, fine=P/40, focus=(0.0, 40.0, -40.0))
            rec, rm, r, rb = analyze(md, B, D, nround=150, nstart=6); rec.update(pitch=P, src=sname, k1=0.5*P*NAE/lamE)
            out['euv2d'].append(rec); print(json.dumps({k: rec[k] for k in ['pitch','src','k1','F','rigorous','best_feasible','gap_rel','rank','EL_bound']}), flush=True); save()
        except Exception as e: print("failEUV", P, sname, e, flush=True)
