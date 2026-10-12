# Audit (7): samples_1d / samples_2d edge-offset placement: true distance of every kept sample to the true target edges
import numpy as np, warnings, json; warnings.filterwarnings('ignore')
from certlitho import *
# ---- 1D
for L, segs, margin, spacing in [(500, [(0,40),(100,180),(240,400)], 8, 4.0), (600, [(0,50),(100,150),(200,250),(430,480)], 8, 1.0), (400, [(0,380)], 6, 4.0), (500, [(0,120),(140,160),(180,300)], 8, 4.0), (436, [(0,44),(88,132),(176,220),(306,371)], 6, 4.0)]:
    x, B, D, tgt = samples_1d(L, segs, margin, spacing)
    E = np.array(sorted(set([a for a,b in segs]+[b for a,b in segs])))
    d = np.min(np.abs(((x[:, None] - E[None, :] + L/2) % L) - L/2), axis=1)
    inside = np.array([any(a <= (v % L) < b for a, b in segs) for v in x])
    ok_label = np.all(inside[B]) and not np.any(inside[D])
    nedge_pts = sum(1 for e in E for s in (-1, 1) if np.min(np.abs(((x - ((e + s*margin) % L)) + L/2) % L - L/2)) < 1e-9)
    print(f"1D L={L} segs={segs}: nB={len(B)} nD={len(D)} min dist B={d[B].min():.6f} min dist D={d[D].min():.6f} labels ok={ok_label} edge-offset points present={nedge_pts}/{2*len(E)} ; points exactly at dist=margin: {int(np.sum(np.abs(d[np.r_[B,D]]-margin)<1e-9))}")
# ---- 2D: exact distance to rectangle edges (periodic) for the hole targets and the random-rectangle targets
def exact_dist_holes(pts, P, hw):
    # square hole |x|<hw, |y|<hw centered at 0 (periodic). signed distance to boundary of the square.
    a = (pts[:,0] + P/2) % P - P/2; b = (pts[:,1] + P/2) % P - P/2
    dx = np.abs(a) - hw; dy = np.abs(b) - hw
    outside = np.hypot(np.maximum(dx, 0), np.maximum(dy, 0)); inside = -np.minimum(np.maximum(dx, dy), 0)
    return np.where((dx < 0) & (dy < 0), inside, outside), (dx < 0) & (dy < 0)
holes = lambda P, hw: (lambda a, b: (abs((a+P/2)%P-P/2) < hw) and (abs((b+P/2)%P-P/2) < hw))
for P, hw, margin, fine in [(160, 35, 12, 2.0), (120, 25, 10, 2.0), (100, 22, 8, 2.0), (40, 10, 4.0, 1.0), (160, 35, 12, 1.0), (160, 35, 12, 0.5)]:
    x, B, D, tgt, sd = samples_2d((P,P), holes(P,hw), margin, 16.0, fine=fine)
    dist, ins = exact_dist_holes(x, P, hw)
    print(f"2D holes P={P} hw={hw} margin={margin} fine={fine}: nB={len(B)} nD={len(D)} | true dist: min B={dist[B].min():.3f} min D={dist[D].min():.3f} (margin {margin}) | #samples closer than margin: {int(np.sum(dist < margin - 1e-9))} | label mismatches={int(np.sum(ins != tgt))} | min |sd| used={np.abs(sd).min():.2f}")
# random rectangle targets as in exp17: exact distance via rectangles' edges
rng = np.random.default_rng(7)
def exact_dist_rects(pts, rects, Lc):
    best = np.full(len(pts), np.inf); ins_any = np.zeros(len(pts), bool)
    for x0, y0, w_, h_ in rects:
        a = (pts[:,0] - x0) % Lc; b = (pts[:,1] - y0) % Lc    # in [0,Lc)
        # distance to this rectangle's boundary, periodic: consider shifts
        dmin = np.full(len(pts), np.inf)
        for sa in (-Lc, 0, Lc):
            for sb in (-Lc, 0, Lc):
                aa = a + sa; bb = b + sb; dx = np.maximum(np.maximum(-aa, aa - w_), 0); dy = np.maximum(np.maximum(-bb, bb - h_), 0)
                inside = (aa > 0) & (aa < w_) & (bb > 0) & (bb < h_)
                dd = np.where(inside, np.minimum(np.minimum(aa, w_-aa), np.minimum(bb, h_-bb)), np.hypot(dx, dy))
                dmin = np.minimum(dmin, dd)
        ins = (a < w_) & (b < h_)
        best = np.minimum(best, dmin); ins_any |= ins
    return best, ins_any   # NOTE: for overlapping rects this underestimates distance (edges inside union) -> only use non-overlapping
tot_close = 0; tot = 0; worst = np.inf
for trial in range(12):
    Lc = float(rng.choice([200, 240])); nrect = 1; rects = []
    for _ in range(nrect):
        w_ = rng.uniform(50, 0.6*Lc); h_ = rng.uniform(50, 0.6*Lc); x0 = rng.uniform(0, Lc); y0 = rng.uniform(0, Lc); rects.append((x0, y0, w_, h_))
    tf = lambda a, b, rects=rects, Lc=Lc: any(((a - x0) % Lc) < w_ and ((b - y0) % Lc) < h_ for x0, y0, w_, h_ in rects)
    x, B, D, tgt, sd = samples_2d((Lc, Lc), tf, 12, 16.0, fine=2.0)
    dist, ins = exact_dist_rects(x, rects, Lc)
    close = int(np.sum(dist < 12 - 1e-9)); tot_close += close; tot += len(x); worst = min(worst, dist.min())
    print(f"rand2d trial {trial}: n={len(x)} closer-than-margin={close} min dist={dist.min():.3f} label mismatches={int(np.sum(ins != tgt))}")
print(f"random 2D: {tot_close}/{tot} samples closer than margin=12; worst distance {worst:.3f}")
