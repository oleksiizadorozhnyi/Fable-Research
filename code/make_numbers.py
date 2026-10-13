"""Generate paper/numbers.tex (macros) and paper/tab_*.tex (table bodies) from results/*.json. Run after assemble.py."""
import json, os, numpy as np
R = '/home/user/Fable-Research/results/'; P = '/home/user/Fable-Research/paper/'
L = lambda n: json.load(open(R+n)) if os.path.exists(R+n) else None
e1 = L('exp1_1d_features.json'); e2 = L('exp2_random_tightness.json'); e3 = L('exp3_pixel_binary.json'); e3b = L('exp3b_ilt.json'); e4 = L('exp4_2d.json')
e5 = L('exp5_smo_euv_toeplitz.json'); e7 = L('exp7_lasserre.json'); e8 = L('exp8_robustness.json'); e9 = L('exp9_outlier_pw.json'); e10 = L('exp10_baselines_complex.json')
e11 = L('exp11_lasserre_pw.json'); e12 = L('exp12_lasserre_xcheck.json'); e13 = L('exp13_2d_strong.json'); e14 = L('exp14_bnb.json'); e14b = L('exp14b_bnb_big.json')
e15 = L('exp15_verify_all.json'); e16 = L('exp16_physics.json'); e17 = L('exp17_dataset.json'); e18 = L('exp18_quasar_strong.json'); e19 = L('exp19_bnb_sdp.json'); e20 = L('exp20_corner_exclusion.json')
M = {}
_DIG={'0':'Zero','1':'One','2':'Two','3':'Three','4':'Four','5':'Five','6':'Six','7':'Seven','8':'Eight','9':'Nine'}
def _san(n): return ''.join(_DIG.get(c,c) for c in n if c.isalnum())
def m(name, val): M[_san(name)] = val
def tabbody(rows):
    # rows end with \\; the last row must NOT end with \\ because \input{...} followed by \bottomrule breaks a trailing row terminator
    return "\\\\\n".join(r[:-2] if r.endswith("\\\\") else r for r in rows)+"\n"
def pct(x, d=1):
    v = round(100*x, d); v = 0.0 if abs(v) < 10**(-d) / 2 else v
    return f"{v:.{d}f}\\%"
def f4(x): return f"{x:.4f}"
def elf(rho): return max(2*(1-rho)/(1+rho), 0.0)   # full-range EL, clamped at 0 (rho >= 1: unprintable)
el_full = elf
def gapstr(r):
    g = r.get('gap_rel'); return 'degenerate' if g is None else pct(g, 2 if g < 0.01 else 1)
rows = []
import traceback
def _run(name, fn):
    try: fn()
    except Exception as e: print('block skipped:', name, '->', repr(e)[:120])
def _blk_Table11Dlibrary():
    global M
    names = {'isodense_coh': 'iso-dense, coherent', 'mixed_coh': 'mixed, coherent', 'isodense_conv0.7': 'iso-dense, $\\sigma=0.7$', 'mixed_conv0.7': 'mixed, $\\sigma=0.7$', 'isodense_dipole0.8': 'iso-dense, dipole 0.8', 'mixed_dipole0.8': 'mixed, dipole 0.8',
             'isodense_conv0.7_PW80': 'iso-dense, $\\sigma=0.7$, $\\pm\\SI{80}{nm}$', 'mixed_conv0.7_PW80': 'mixed, $\\sigma=0.7$, $\\pm\\SI{80}{nm}$', 'isodense_dipole0.8_PW80': 'iso-dense, dipole 0.8, $\\pm\\SI{80}{nm}$', 'dense90_conv0.7': 'dense 90, $\\sigma=0.7$', 'dense90_dipole0.8': 'dense 90, dipole 0.8', 'dense90_coh': 'dense 90, coherent'}
    for r in e1: rows.append(f"{names.get(r['name'], r['name'])} & {r['F']} & {len(r['focus'])} & {f4(r['rigorous'])} & {elf(r['rigorous']):.3f} & {f4(r['best_feasible'])} & {elf(r['best_feasible']):.3f} & {gapstr(r)}\\\\")
    open(P+'tab_1d.tex','w').write(tabbody(rows))
    d1 = {r['name']: r for r in e1}
    m('isoCohGap', pct(d1['isodense_coh']['gap_rel'])); m('mixCohGap', pct(d1['mixed_coh']['gap_rel'])); m('isoConvGap', pct(d1['isodense_conv0.7']['gap_rel'])); m('mixConvGap', pct(d1['mixed_conv0.7']['gap_rel']))
    pw = [d1[k]['gap_rel'] for k in ['isodense_conv0.7_PW80','mixed_conv0.7_PW80','isodense_dipole0.8_PW80']]; m('pwGapLo', pct(min(pw),0)); m('pwGapHi', pct(max(pw),0))
    els = [elf(r['rigorous']) for r in e1 if len(r['focus'])==1]; m('elOneLo', f"{min(els):.2f}"); m('elOneHi', f"{max(els):.2f}")
    elp = [elf(r['rigorous']) for r in e1 if len(r['focus'])==3]; m('elPwLo', f"{min(elp):.2f}"); m('elPwHi', f"{max(elp):.2f}")
    m('dipGapMax', pct(max(d1[k]['gap_rel'] for k in ['isodense_dipole0.8','mixed_dipole0.8']), 2)); dense_bounds = [d1[k]['rigorous'] for k in ['dense90_conv0.7','dense90_dipole0.8','dense90_coh']]; m('denseNinetyGapMax', pct(max(d1[k]['gap_rel'] for k in ['dense90_conv0.7','dense90_dipole0.8','dense90_coh']),2))
    m('denseDipoleBound', f"{d1['dense90_dipole0.8']['rigorous']:.5f}")
_run('Table 1 (1D library)', _blk_Table11Dlibrary)
def _blk_random1D():
    global M
    e9 = L('exp9_outlier_pw.json') or {}
    deg = [r for r in e2 if r['rigorous'] < 1e-3]; nd = [r for r in e2 if r['rigorous'] >= 1e-3]
    m('randDegN', str(len(deg))); m('randDegBoundMax', f"{max((r['rigorous'] for r in deg), default=0):.0e}".replace('e-0','\\times10^{-')+'}' if deg else '0')
    m('randDegBestMax', f"{max((min(r['best_feasible'], e9.get('trial%d' % r['trial'], {}).get('new_best', 9)) for r in deg), default=0):.3f}")
    g = np.array([r['gap_rel'] if r.get('gap_rel') is not None else 0.0 for r in nd]); m('randN', str(len(e2))); m('randNonDegN', str(len(nd))); m('randExact', pct(np.mean(g<=1e-4))); m('randOne', pct(np.mean(g<=1e-2))); m('randFive', pct(np.mean(g<=5e-2)))
    m('randMedian', f"{np.median(g):.1e}".replace('e-0','\\times10^{-').replace('e-','\\times10^{-')+'}'); m('randMax', pct(g.max())); m('randMaxFrac', f"{g.max():.3f}")
    bysrc = {s_: np.mean([ (r['gap_rel'] or 0) <=1e-4 for r in nd if r['src']==s_]) for s_ in set(r['src'] for r in nd)}
    for s_, v in bysrc.items(): m('randExact'+s_.replace('.','').replace('0','z'), pct(v,0))
    big = [r for r in nd if (r.get('gap_rel') or 0) > 0.05]; m('randBigN', str(len(big))); m('randBigRankMax', str(max((r['rank'] for r in big), default=0)))
_run('random 1D', _blk_random1D)
def _blk_BB():
    global M
    rows = []; bb = (e14 or []) + (e14b or [])
    lab = {'trial82': 'random 82 (coherent)', 'trial91': 'random 91 (coherent)', 'trial110': 'random 110 (coherent)', 'trial50': 'random 50 (coherent)', 'mixed_coh': 'mixed (coherent)', 'isodense_coh': 'iso-dense (coherent)', 'mixed_conv0.7_PW80_edgeonly': 'mixed, $\\sigma=0.7$, $\\pm\\SI{80}{nm}$'}
    best = {}
    for r in bb:
        if r['name'] in best and r['bnb_L'] <= best[r['name']]['bnb_L']: continue
        best[r['name']] = r
    conv_gains = []; conv_times = []
    for nm, r in best.items():
        conv = 'yes' if r['bnb_converged'] else f"no ({pct((r['bnb_U']-r['bnb_L'])/r['bnb_U'],2)} gap)"
        valid = r['bnb_L'] > r.get('sdp_ensemble_ratio', r['sdp_bound'])   # gain is certified against the realized multi-exposure ratio (upper bound on t_R)
        rows.append(f"{lab.get(nm, nm)} & {r['F']} & {f4(r['sdp_bound'])} & {f4(max(r['bnb_L'],0))} & {f4(r['bnb_U'])} & {r['bnb_nodes']:,} & {conv} & {pct(r['certified_gain']) if valid else 'none'}\\\\")
        if r['bnb_converged']: conv_gains.append(r['certified_gain']); conv_times.append(r['bnb_time'])
    open(P+'tab_bnb.tex','w').write(tabbody(rows))
    m('bnbGainLo', pct(min(conv_gains)) if conv_gains else 'n/a'); m('bnbGainHi', pct(max(conv_gains)) if conv_gains else 'n/a'); m('bnbTimeLo', f"{min(conv_times):.1f}" if conv_times else 'n/a'); m('bnbTimeHi', f"{max(conv_times):.1f}" if conv_times else 'n/a')
    if 'mixed_coh' in best: r = best['mixed_coh']; m('bnbMixGain', pct(r['certified_gain'])); m('bnbMixL', f4(r['bnb_L'])); m('bnbMixU', f4(r['bnb_U'])); m('bnbMixNodes', f"{r['bnb_nodes']/1e6:.0f}"); m('bnbMixGap', pct((r['bnb_U']-r['bnb_L'])/r['bnb_U'],2)); m('bnbMixTime', f"{r['bnb_time']/60:.0f}")
    if e19:
        ok = [r for r in e19 if r['certified']]; m('sdpbnbCertified', "; ".join(f"{lab.get(r['name'], r['name'])}: $T={r['T']}$ (gain {pct(r['certified_gain_if_ok'])}, {r['nodes']} nodes, {r['time']:.0f}\\,s)" for r in ok) or 'none'); m('sdpbnbFailed', "; ".join(f"{lab.get(r['name'], r['name'])}: $T={r['T']}$" for r in e19 if not r['certified']) or 'none')
    else: m('sdpbnbCertified', 'n/a'); m('sdpbnbFailed', 'n/a')
    las = {r['name']: r for r in (e7 or []) + (e11 or [])}
    for k, key in [('lasIso', 'isodense_coh'), ('lasMix', 'mixed_coh'), ('lasPW', 'mixed_conv0.7_PW80_edgeonly')]:
        if key in las:
            r = las[key]; lb = r.get('lasserre2_bound', r.get('lasserre2')); sb = r.get('rigorous_shor', r.get('shor')); m(k+'Bound', f4(lb)); m(k+'Gain', pct((lb-sb)/sb)); m(k+'Best', f4(r.get('best_single_mask', r.get('best_single'))))
    if e12: m('lasXClarabel', f"{e12['CLARABEL']['t']:.7f}"); m('lasXScs', f"{e12['SCS']['t']:.7f}")
    devs = [abs(r['lasserre2_bound'] - r['best_single_mask']) / r['best_single_mask'] for r in e7 if r.get('lasserre2_bound')]
    if e9:
        for k in ('mixed_conv0.7_PW80',): pass
    m('lasMaxDev', f"{max(devs):.0e}".replace('e-0','\\times10^{-').replace('e-','\\times10^{-')+'}')
_run('B&B', _blk_BB)
def _blk_PSMvsbinaryILT():
    global M
    rows = []; e22 = {r['name']: r for r in (L('exp22_psm_penalty.json') or [])}
    plo, phi, ilo, ihi = [], [], [], []
    for r3, rb in zip(e3, e3b):
        q = e22[r3['name']]; ilt = rb['ilt_psm_best']; lo = ilt/q['U_R'] - 1; hi = ilt/q['L_R'] - 1
        plo.append(q['penalty_lo']); phi.append(q['penalty_hi']); ilo.append(lo); ihi.append(hi)
        rows.append(f"{r3['name'].replace('dense120','dense 120').replace('isodense','iso-dense')} & {f4(q['L_R'])} & {f4(q['U_R'])} & {f4(q['L_plus'])} & {f4(q['U_plus'])} & [{pct(q['penalty_lo'])}, {pct(q['penalty_hi'])}] & {ilt:.4f} & [{pct(lo)}, {pct(hi)}]\\\\")
    open(P+'tab_psm.tex','w').write(tabbody(rows))
    m('binPenLo', pct(min(plo),1)); m('binPenHi', pct(max(plo),1)); m('binPenHiLo', pct(min(phi),1)); m('binPenHiHi', pct(max(phi),1))
    m('iltPsmLo', pct(min(ilo),0)); m('iltPsmHi', pct(max(ihi),0)); iltb = [rb['ilt_binary_best']/rb['certified_bound']-1 for rb in e3b]; m('iltBinHi', pct(max(iltb),0))
    m('pixFourierAgree', f"{max(abs(q['L_R']-q['L_pixel_psm'])/q['L_R'] for q in e22.values()):.0e}".replace('e-0','\\times10^{-').replace('e-','\\times10^{-')+'}')
    assert all(q['L_R'] <= q['U_R'] + 1e-12 and q['L_plus'] <= q['U_plus'] + 1e-12 and all(v['ok'] for v in q['exact'].values()) for q in e22.values()), "bracket or exact-verification inconsistency in exp22"
    m('psmPixelN', ", ".join(str(q['N_pixels']) for q in e22.values()))
    dl = [1 - q['L_plus'] / r3['t_pixel_nonneg_DNN'] for q, r3 in zip(e22.values(), e3)]; m('dnnLossLo', pct(min(dl),1)); m('dnnLossHi', pct(max(dl),1))
_run('PSM vs binary, ILT', _blk_PSMvsbinaryILT)
def _blk_SMO():
    global M
    smo = {r['src']: r for r in e5['smo']}; m('smoCoh', f4(smo['coh']['fixed_source_bound'])); m('smoConvFixed', f4(smo['conv0.7']['fixed_source_bound'])); m('smoConvRelax', f4(smo['conv0.7']['smo_bound'])); m('smoNineFixed', f4(smo.get('conv0.9', smo.get('full_grid0.9'))['fixed_source_bound'])); m('smoNineRelax', f4(smo.get('conv0.9', smo.get('full_grid0.9'))['smo_bound'])); m('smoDipFixed', f4(smo['dipole0.8']['fixed_source_bound'])); m('smoDipRelax', f4(smo['dipole0.8']['smo_bound']))
_run('SMO', _blk_SMO)
def _blk_EUV1D():
    global M
    rows = []
    for r in e5['euv_pw']: rows.append(f"{r['pitch']} & {r['k1']:.3f} & {r['src']} & {r['F']} & {f4(r['rigorous'])} & {elf(r['rigorous']):.3f} & {f4(r['best_feasible'])} & {elf(r['best_feasible']):.3f} & {gapstr(r)}\\\\")
    open(P+'tab_euv.tex','w').write(tabbody(rows))
    eu = {(r['pitch'], r['src']): r for r in e5['euv_pw']}
    m('euvTwentyTwoConvEL', f"{elf(eu[(22,'conv0.9')]['rigorous']):.3f}"); m('euvTwentyTwoConvELBest', f"{elf(eu[(22,'conv0.9')]['best_feasible']):.3f}"); _mx = max(e5['euv_pw'], key=lambda r: r.get('gap_rel') or 0); m('euvLoosePitch', str(_mx['pitch'])); m('euvLooseSrc', 'finite-pole dipole' if 'dipole' in _mx['src'] else 'disk'); m('euvLooseELBound', f"{elf(_mx['rigorous']):.3f}"); m('euvLooseELBest', f"{elf(_mx['best_feasible']):.3f}"); m('euvTightN', str(sum(1 for r in e5['euv_pw'] if (r.get('gap_rel') or 0) < 1e-4))); m('euvTwentyTwoDipBound', f4(eu[(22,'dipole0.8')]['rigorous'])); m('euvN', str(len(e5['euv_pw']))); m('euvMaxGap', pct(max((r.get('gap_rel') or 0) for r in e5['euv_pw']),2))
    m('thmFourNec', str(e5['theorem4']['necessity_violations'])); m('thmFourNecN', str(e5['theorem4']['necessity_trials'])); m('thmFourSuf', str(e5['theorem4']['sufficiency_failures'])); m('thmFourSufN', str(e5['theorem4']['sufficiency_tested']))
_run('EUV 1D', _blk_EUV1D)
def _blk_robustness():
    global M
    sp = [r['bound'] for r in e8['sampling']]; m('spacingMin', f"{min(sp):.5f}"); m('spacingMax', f"{max(sp):.5f}")
    mg = e8['margin']; m('marginList', ", ".join(f"$\\tau={r['margin']}$\\,nm: {elf(r['bound']):.2f}" for r in mg))
    q = e8['source_quadrature_and_solvers']; m('quadList', ", ".join(f"${r.get('grid')}^2$ grid ({r.get('nsrc'):,} pts): {r['clarabel']:.5f}" for r in q)); m('quadLastTwoDiff', pct(abs(q[-1]['clarabel']-q[-2]['clarabel'])/q[-1]['clarabel'],2)); _q = [r for r in q if r.get('grid') == 401]; m('quadRefDiff', pct(abs(_q[0]['clarabel']-q[-1]['clarabel'])/q[-1]['clarabel'],3) if _q else 'n/a'); m('quadNconv', f"{(_q[0]['nsrc'] if _q else 0):,}"); m('quadSolverDiffMax', f"{max(r['diff'] for r in q):.0e}".replace('e-0','\\times10^{-').replace('e-','\\times10^{-')+'}')
_run('robustness', _blk_robustness)
def _blk_Rayleighcomplex():
    global M
    rows = []; ratios = []
    for r in e10:
        ratios.append(r['sdp_real']/max(r['rayleigh_bound'],1e-12)); rows.append(f"{r['case'].replace('dense90','dense 90').replace('isodense','iso-dense')} & {r['src'].replace('coh','coherent').replace('conv0.7','$\\sigma=0.7$').replace('dipole0.8','dipole 0.8')} & {r['F']} & {f4(r['sdp_real'])} & {f4(r['rayleigh_bound'])} & {f4(r['sdp_complex'])}\\\\")
    open(P+'tab_rayleigh.tex','w').write(tabbody(rows)); m('rayLo', f"{min(ratios):.1f}"); m('rayHi', f"{max(ratios):.0f}"); m('rayN', str(len(e10))); m('cplxMaxDiff', f"{max(abs(r['sdp_real']-r['sdp_complex']) for r in e10):.0e}".replace('e-0','\\times10^{-').replace('e-','\\times10^{-')+'}')
_run('Rayleigh / complex', _blk_Rayleighcomplex)
def _blk_2D():
    global M
    rows = []; d4 = {r['name']: r for r in e4}; d13 = {r['name']: r for r in (e13 or [])}
    lab2 = {'contacts160_annular': 'contacts 160, annular', 'contacts160_conv0.5': 'contacts 160, $\\sigma=0.5$', 'contacts120_annular': 'contacts 120, annular', 'contacts120_quasar': 'contacts 120, quasar', 'contacts100_quasar': 'contacts 100, quasar', 'Tjunction320_annular': 'T-junction 320, annular', 'Tjunction320_quasar': 'T-junction 320, quasar', 'lineend320_annular': 'line-ends 320, annular'}
    for nm, r in d4.items():
        bst = r['best_feasible']; tag = ''
        if nm in d13: bst = min(bst, d13[nm]['best_single']); tag = ' (strong)'
        if nm == 'Tjunction320_quasar' and e18: bst = min(bst, e18['best_single']); tag = ' (strong)'
        gap = (bst - r['rigorous'])/r['rigorous']
        rows.append(f"{lab2.get(nm, nm)}{tag} & {r['F']} & {f4(r['rigorous'])} & {elf(r['rigorous']):.3f} & {f4(bst)} & {elf(bst):.3f} & {pct(gap,2)} & {r['sdp_time']:.1f}\\\\")
        m('twoD'+nm.replace('_','').replace('.',''), f"{gap:.4f}")
    open(P+'tab_2d.tex','w').write(tabbody(rows))
    cg = [ (min(d4[n]['best_feasible'], d13.get(n, {}).get('best_single', 9)) - d4[n]['rigorous'])/d4[n]['rigorous'] for n in d4 if n.startswith('contacts')]; m('contactsGapMax', pct(max(cg),2))
    def g320(n): b = min(d4[n]['best_feasible'], d13.get(n, {}).get('best_single', 9), (e18['best_single'] if (n=='Tjunction320_quasar' and e18) else 9)); return (b-d4[n]['rigorous'])/d4[n]['rigorous']
    m('tjF', str(d4['Tjunction320_annular']['F'])); m('tjConstraints', f"{d4['Tjunction320_annular']['nB']+d4['Tjunction320_annular']['nD']:,}"); m('contactsRankMax', str(max(d4[n]['rank'] for n in d4 if n.startswith('contacts')))); m('tjAnnGap', pct(g320('Tjunction320_annular'),2)); m('tjQuaGap', pct(g320('Tjunction320_quasar'),1)); m('leGap', pct(g320('lineend320_annular'),2)); m('tjAnnBound', f4(d4['Tjunction320_annular']['rigorous'])); m('tjQuaBound', f4(d4['Tjunction320_quasar']['rigorous'])); m('leBound', f4(d4['lineend320_annular']['rigorous']))
    m('nullDimMax', str(max((r.get('cert_null_dim') or 0) for r in e4))); t61 = [r['sdp_time'] for r in e4 if r['F'] >= 50]; m('sdpTimeBigLo', f"{min(t61):.0f}"); m('sdpTimeBigHi', f"{max(t61):.0f}"); tsm = [r['sdp_time'] for r in e4 if r['F'] <= 15]; m('sdpTimeSmallHi', f"{max(tsm):.1f}")
_run('2D', _blk_2D)
def _blk_expandeddataset():
    global M
    if e17:
        r2 = e17['random2d']; g2 = np.array([r['gap_rel'] if r.get('gap_rel') is not None else 0 for r in r2]); m('rtwoN', str(len(r2))); m('rtwoTight', pct(np.mean(g2<=1e-3))); m('rtwoFive', pct(np.mean(g2<=5e-2))); m('rtwoMax', pct(g2.max(),0)); m('rtwoFlo', str(min(r['F'] for r in r2))); m('rtwoFhi', str(max(r['F'] for r in r2)))
        a1 = e17['adversarial1d']; cohg = [r['gap_rel'] for r in a1 if r['src']=='coh' and r.get('gap_rel') is not None and r['gap_rel']>0.01]; m('advCohGapLo', pct(min(cohg),0) if cohg else '0'); m('advCohGapHi', pct(max(cohg),0) if cohg else '0'); m('advCohBigN', str(len(cohg)))
        eu2 = e17['euv2d']; m('euvTwoDElLo', f"{min(elf(r['rigorous']) for r in eu2):.2f}"); m('euvTwoDElHi', f"{max(elf(r['rigorous']) for r in eu2):.2f}"); m('euvTwoDGapMax', pct(max((r.get('gap_rel') or 0) for r in eu2),2))
_run('expanded dataset', _blk_expandeddataset)
def _blk_verification():
    global M
    if e15:
        d = [r['solver_diff'] for r in e15 if r['solver_diff'] is not None]; m('verN', str(len(e15))); m('verChol', str(sum(1 for r in e15 if r['cholesky50']))); m('verDeg', str(sum(1 for r in e15 if r.get('degenerate')))); m('verFail', str(sum(1 for r in e15 if r['cholesky50'] is False)))
        nde = [r for r in e15 if r.get('rigorous', 0) >= 1e-3]; m('verKappaMin', f"{min(r['c'] for r in nde):.3f}"); m('verKappaBelow', str(sum(1 for r in nde if r['c'] < 0.9999))); m('verNonDegN', str(len(nde))); dg = [r for r in e15 if r.get('rigorous', 0) < 1e-3]; m('verDegKappa', f"{min((r['c'] for r in dg), default=1):.0e}".replace('e-0','\\times10^{-').replace('e-','\\times10^{-')+'}'); m('verSolverMedian', f"{np.median(d):.1e}".replace('e-0','\\times10^{-').replace('e-','\\times10^{-')+'}'); bigd = [r for r in e15 if r['solver_diff'] is not None and r['solver_diff']>1e-5]; m('verBigDiffN', str(len(bigd))); m('verBigDiffMax', f"{max((r['solver_diff'] for r in bigd), default=0):.2f}"); m('verBigDiffBoundMax', f"{max((r['rigorous'] for r in bigd), default=0):.3f}")
_run('verification', _blk_verification)
def _blk_physics():
    global M
    if e16:
        v1 = e16['vector_1d']; te = [r['TE_bound']/r['scalar_bound']-1 for r in v1]; tm = [r['TM_bound']/r['scalar_bound']-1 for r in v1]; up = [r['unpol_bound']/r['scalar_bound']-1 for r in v1]
        m('vecTeLo', pct(min(te),0)); m('vecTeHi', pct(max(te),0)); m('vecTmLo', pct(min(tm),0)); m('vecTmHi', pct(max(tm),0)); m('vecUnLo', pct(min(up),0)); m('vecUnHi', pct(max(up),0)); m('vecGapMax', pct(max(max(r['TE_gap'], r['TM_gap'], r['unpol_gap']) for r in v1),1))
        v2 = e16['vector_2d']; m('vecTwoDUnLo', pct(min(r['unpol_bound']/r['scalar_bound']-1 for r in v2),0)); m('vecTwoDUnHi', pct(max(r['unpol_bound']/r['scalar_bound']-1 for r in v2),0)); m('vecTwoDTeLo', pct(min(r['TE_bound']/r['scalar_bound']-1 for r in v2),0)); m('vecTwoDTeHi', pct(max(r['TE_bound']/r['scalar_bound']-1 for r in v2),0))
        ab = e16['aberrations']; base = [r['bound'] for r in ab if r['waves']==0][0]; m('abMax', pct(max(r['bound']/base-1 for r in ab),0))
        px = e16['pixelation']; F = px[0]['F']; same = [r for r in px if r['N'] >= F]; coarse = [r for r in px if r['N'] < F]; m('pixSameMaxDev', f"{max(abs(r['bound']-r['fourier'])/r['fourier'] for r in same):.0e}".replace('e-0','\\times10^{-').replace('e-','\\times10^{-')+'}'); m('pixCoarse', "; ".join(f"{r['pixel_nm']:.0f}\\,nm: {pct(r['bound']/r['fourier']-1,0)}" for r in coarse))
    pass
_run('physics', _blk_physics)

def _blk_corner():
    global M
    lab = {'contacts160_annular': 'contacts 160, annular', 'contacts120_quasar': 'contacts 120, quasar', 'contacts100_quasar': 'contacts 100, quasar', 'Tjunction320_annular': 'T-junction 320, annular', 'lineend320_annular': 'line-end 320, annular'}
    by = {}
    for r in e20: by.setdefault(r['name'], {})[r['frac_cd']] = r
    rows = []; ratios = []
    for nm, d in by.items():
        r0 = d[0.0]; rows.append(f"{lab.get(nm, nm)} & {r0['tau']:g} & " + " & ".join(f"{f4(d[f]['rigorous'])} & {el_full(d[f]['rigorous']):.3f}" for f in (0.0, 0.15, 0.25)) + "\\\\")
        ratios.append(el_full(d[0.25]['rigorous']) / max(el_full(r0['rigorous']), 1e-12))
    open(P+'tab_corner.tex','w').write(tabbody(rows))
    m('cornerElGainLo', f"{min(ratios):.1f}"); m('cornerElGainHi', f"{max(ratios):.1f}"); m('cornerN', str(len(by)))
    c160 = by['contacts160_annular']; m('cornerContactsOneSixtyEL', f"{el_full(c160[0.0]['rigorous']):.2f}"); m('cornerContactsOneSixtyELQuarter', f"{el_full(c160[0.25]['rigorous']):.2f}")
_run('corner exclusion', _blk_corner)
def _blk_exact():
    global M
    e21 = L('exp21_exact_verify.json') or []
    pd_ = [r for r in e21 if r['status'].startswith('PD')]; zero = [r for r in e21 if r['status'].startswith('zero')]; fail = [r for r in e21 if r['status'] == 'FAILED']
    lm = [r['lambda_min_50digit'] for r in pd_ if r.get('lambda_min_50digit') is not None]
    fmt = lambda v: f"{v:.0e}".replace('e-0','\\times10^{-').replace('e-','\\times10^{-')+'}'
    m('exactBisected', str(sum(1 for r in pd_ if 'bisection' in r['status']))); m('exactN', str(len(e21))); m('exactPD', str(len(pd_))); m('exactZero', str(len(zero))); m('exactFail', str(len(fail))); m('exactLminLo', fmt(min(lm)) if lm else 'n/a'); m('exactLminHi', fmt(max(lm)) if lm else 'n/a')
_run('exact verification', _blk_exact)
def _count_viol():
    n = 0
    for obj in (e1, e2, e4, e17):
        recs = obj if isinstance(obj, list) else [r for v in (obj or {}).values() for r in (v if isinstance(v, list) else [])]
        n += sum(1 for r in recs if isinstance(r, dict) and (r.get('cert_violation') or (r.get('best_feasible') is not None and r.get('rigorous') is not None and r['best_feasible'] < r['rigorous']*(1-1e-6) - 1e-12)))
    return n
_nv = _count_viol(); m('certViolations', str(_nv)); print("certificate violations (feasible mask below its bound):", _nv)
with open(P+'numbers.tex','w') as f:
    for k, v in M.items(): f.write(f"\\newcommand{{\\{k}}}{{{v}}}\n")
print("macros:", len(M)); print("\n".join(f"{k}={v}" for k, v in list(M.items())[:12]))

