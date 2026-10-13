### Table 1 — 1D feature library (DUV 193 nm, NA 1.35; edge tolerance 8 nm; dense90: 6 nm)
| case | F | foci | certified bound ρ | max EL (full range) | best single mask ρ | gap | SDP time (s) |
|---|---|---|---|---|---|---|---|
| isodense_coh | 9 | 1 | 0.5359 | 0.604 | 0.6383 | 19.12% | 0.26 |
| mixed_coh | 7 | 1 | 0.5875 | 0.520 | 0.6439 | 9.59% | 0.62 |
| isodense_conv0.7 | 15 | 1 | 0.7163 | 0.331 | 0.7208 | 0.63% | 0.50 |
| mixed_conv0.7 | 11 | 1 | 0.6812 | 0.379 | 0.6812 | 0.00% | 0.81 |
| isodense_dipole0.8 | 17 | 1 | 0.5658 | 0.555 | 0.5730 | 1.27% | 1.51 |
| mixed_dipole0.8 | 13 | 1 | 0.6020 | 0.497 | 0.6020 | 0.00% | 0.37 |
| isodense_conv0.7_PW80 | 15 | 3 | 0.8305 | 0.185 | 0.8610 | 3.67% | 4.93 |
| mixed_conv0.7_PW80 | 11 | 3 | 0.8370 | 0.177 | 0.8451 | 0.97% | 2.57 |
| isodense_dipole0.8_PW80 | 17 | 3 | 0.6940 | 0.361 | 0.7566 | 9.01% | 4.53 |
| dense90_conv0.7 | 9 | 1 | 0.7632 | 0.269 | 0.7646 | 0.19% | 0.44 |
| dense90_dipole0.8 | 9 | 1 | 0.5279 | 0.618 | 0.5279 | 0.00% | 0.36 |
| dense90_coh | 5 | 1 | 0.5279 | 0.618 | 0.5279 | 0.00% | 0.12 |

### Table 2 — Certified multi-exposure gain (Lasserre-2 single-mask bound vs SDP multi-exposure bound)
| case | F | SDP (multi-exposure) bound | Lasserre-2 single-mask bound | best single mask found | certified gain | status |
|---|---|---|---|---|---|---|
| isodense_coh | 9 | 0.5359 | 0.6383 | 0.6383 | 19.12% | optimal_inaccurate |
| mixed_coh | 7 | 0.5875 | 0.644 | 0.6439 | 9.60% | optimal_inaccurate |
| mixed_conv0.7_PW80_edgeonly | 11 | 0.8370 | 0.8451 | 0.8451 | 0.97% | optimal_inaccurate |

### Random 1D targets (n=109): gap ≤1e-4: 73.4%; ≤1%: 84.4%; ≤5%: 92.7%; median 5.2e-07; max 44.069; worst certificate min-eig 1.0e-10

### Table 3 — Pixel vs Fourier relaxations, PSM vs binary (DNN), ILT baseline (σ=0.7 conventional)
| case | F | Fourier PSM bound | pixel PSM bound | pixel nonneg (DNN) bound | ILT-PSM ratio | ILT-PSM refined | ILT-binary ratio | best PSM (pipeline) |
|---|---|---|---|---|---|---|---|---|
| isodense | 15 | 0.7163 | 0.7168 | 0.7529 | 0.886 | 0.7351 | 1.053 | 0.7208 |
| mixed | 11 | 0.6812 | 0.6820 | 0.7587 | 0.939 | 0.7828 | 0.965 | 0.6812 |
| dense120 | 9 | 0.6735 | 0.6736 | 0.6942 | 0.849 | 0.6810 | 0.857 | 0.6767 |

### Table 3b — Standard gradient ILT (sigmoid resist, pattern-error loss; best of 6 settings) vs certified bound
| case | F | certified bound ρ | ILT-PSM best ρ | ILT-PSM then ratio-refined | ILT-binary best ρ | ILT-binary then ratio-refined |
|---|---|---|---|---|---|---|
| isodense | 15 | 0.7163 | 0.8613 | 0.7351 | 0.9530 | 0.7513 |
| mixed | 11 | 0.6812 | 0.8975 | 0.7828 | 0.9204 | 0.7828 |
| dense120 | 9 | 0.6735 | 0.8261 | 0.6806 | 0.8520 | 0.6806 |

### Table 4 — Source–mask entitlement bound (iso-dense, DUV)
| source | #pts | fixed-source bound | SMO relaxation bound |
|---|---|---|---|
| coh | 1 | 0.5359 | 0.5359 |
| conv0.7 | 349 | 0.7183 | 0.5359 |
| dipole0.8 | 138 | 0.5760 | 0.5359 |
| conv0.9 | 349 | 0.7202 | 0.5359 |

### Table 5 — EUV (13.5 nm, NA 0.33, ±50 nm focus) dense lines: certified max EL (full range) vs pitch
| pitch (nm) | k1 | source | F | bound ρ | max EL (full range) | best single mask ρ | gap |
|---|---|---|---|---|---|---|---|
| 20 | 0.244 | conv0.9 | 7 | 1.0000 | 0.000 | 1.0000 | 0.00% |
| 20 | 0.244 | dipole0.8 | 7 | 1.0000 | 0.000 | 1.0000 | 0.00% |
| 22 | 0.269 | conv0.9 | 9 | 0.9467 | 0.055 | 0.9467 | 0.00% |
| 22 | 0.269 | dipole0.8 | 9 | 0.7404 | 0.298 | 0.7404 | 0.00% |
| 24 | 0.293 | conv0.9 | 9 | 0.8405 | 0.173 | 0.8405 | 0.00% |
| 24 | 0.293 | dipole0.8 | 9 | 0.5714 | 0.546 | 0.5714 | 0.00% |
| 26 | 0.318 | conv0.9 | 9 | 0.7806 | 0.246 | 0.7806 | 0.00% |
| 26 | 0.318 | dipole0.8 | 9 | 0.5457 | 0.588 | 0.5457 | 0.00% |
| 28 | 0.342 | conv0.9 | 11 | 0.7407 | 0.298 | 0.7407 | 0.00% |
| 28 | 0.342 | dipole0.8 | 11 | 0.5551 | 0.572 | 0.5551 | 0.00% |
| 32 | 0.391 | conv0.9 | 11 | 0.6891 | 0.368 | 0.6891 | 0.00% |
| 32 | 0.391 | dipole0.8 | 13 | 0.5902 | 0.515 | 0.5902 | 0.00% |
| 40 | 0.489 | conv0.9 | 15 | 0.6284 | 0.456 | 0.6284 | 0.00% |
| 40 | 0.489 | dipole0.8 | 15 | 0.5740 | 0.541 | 0.5948 | 3.63% |

Theorem 4 numerics: necessity violations 0/2000; sufficiency failures 0/200.

### Robustness (iso-dense, σ=0.7)
- interior spacing sweep (bound): 16 nm: 0.71618, 8 nm: 0.71625, 4 nm: 0.71625, 2 nm: 0.71626, 1 nm: 0.71626
- margin sweep (bound / EL): 2 nm: 0.8747/0.067, 2.5 nm: 0.8460/0.083, 4 nm: 0.7655/0.133, 5 nm: 0.7163/0.165, 8 nm: 0.5878/0.260, 12 nm: 0.4545/0.375, 16 nm: 0.3571/0.474, 20 nm: 0.2901/0.550
- source quadrature & solver agreement: 349 pts: CLARABEL 0.71830 SCS 0.71830, 1313 pts: CLARABEL 0.71443 SCS 0.71443, 2933 pts: CLARABEL 0.71678 SCS 0.71678, 8021 pts: CLARABEL 0.71638 SCS 0.71638, 31757 pts: CLARABEL 0.71624 SCS 0.71624, 126309 pts: CLARABEL 0.71625 SCS 0.71625, 503973 pts: CLARABEL 0.71619 SCS 0.71619
- 50-digit certificate check: Cholesky success=True, min eig=7.15e-11, bound=0.587782

### Table 6 — 2D layouts (DUV, NA 1.35), edge tolerance 12 nm (contacts120/100: 10/8 nm; T/line-end: 14 nm)
| layout | cell (nm) | source pts | F | samples | bound ρ | max EL (full range) | best single mask ρ | gap | rank | cert. null dim | SDP time (s) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| contacts160_annular | [160, 160] | 70192 | 13 | 352 | 0.8040 | 0.217 | 0.8221 | 2.25% | 1 | 0 | 2.4 |
| contacts160_conv0.5 | [160, 160] | 126309 | 9 | 352 | 0.8383 | 0.176 | 0.8384 | 0.01% | 1 | 0 | 1.8 |
| contacts120_annular | [120, 120] | 70192 | 9 | 243 | 0.8749 | 0.133 | 0.8765 | 0.18% | 1 | 0 | 0.9 |
| contacts120_quasar | [120, 120] | 23404 | 9 | 243 | 0.8385 | 0.176 | 0.8394 | 0.11% | 1 | 0 | 1.5 |
| contacts100_quasar | [100, 100] | 23404 | 5 | 211 | 0.9782 | 0.022 | 0.9782 | 0.00% | 1 | 0 | 0.8 |
| Tjunction320_annular | [320, 320] | 70192 | 61 | 1273 | 0.8050 | 0.216 | 0.8663 | 7.62% | 2 | 0 | 2242.7 |
| Tjunction320_quasar | [320, 320] | 23404 | 49 | 1273 | 0.7817 | 0.245 | 0.8501 | 8.75% | 4 | 0 | 794.7 |
| lineend320_annular | [320, 320] | 70192 | 61 | 1288 | 0.8677 | 0.142 | 0.8950 | 3.14% | 2 | 0 | 1323.8 |

### Outlier / process-window re-examination (strong primal search, face exploration, Lasserre-2 where F≤11)
| case | F | bound | old best | new best | new gap | min face rank | Lasserre-2 | status |
|---|---|---|---|---|---|---|---|---|
| trial16 | 13 | 0.0004 | 0.0183 | 0.0183 | 4406.90% | 2 | None | skipped(F>11) |
| trial82 | 5 | 0.6035 | 0.6683 | 0.6683 | 10.75% | 2 | 0.6683098425259854 | optimal_inaccurate |
| trial91 | 5 | 0.7787 | 0.8300 | 0.8300 | 6.59% | 2 | 0.830043033181968 | optimal_inaccurate |
| trial110 | 5 | 0.6268 | 0.8619 | 0.8619 | 37.51% | 2 | 0.8618960894055052 | optimal_inaccurate |
| trial50 | 5 | 0.5572 | 0.5835 | 0.5835 | 4.71% | 2 | 0.5834942690705798 | optimal_inaccurate |
| trial46 | 13 | 0.2854 | 0.2854 | 0.2854 | 0.00% | 1 | None | skipped(F>11) |
| mixed_conv0.7_PW80 | 11 | 0.8370 | nan | 0.8451 | 0.97% | 3 | None | skipped(three foci: Lasserre-2 exceeds memory; see exp11b for the edge-only through-focus instance) |
| isodense_conv0.7_PW80 | 15 | 0.8305 | nan | 0.8610 | 3.67% | 3 | None | skipped(three foci: Lasserre-2 exceeds memory; see exp11b for the edge-only through-focus instance) |
| isodense_dipole0.8_PW80 | 17 | 0.6940 | nan | 0.7566 | 9.01% | 2 | None | skipped(three foci: Lasserre-2 exceeds memory; see exp11b for the edge-only through-focus instance) |

### Table 7 — Baseline Rayleigh-quotient bound vs SDP; complex-transmission masks
| case | source | F | SDP real-mask bound | Rayleigh-quotient bound | SDP complex-mask bound |
|---|---|---|---|---|---|
| isodense | coh | 9 | 0.5359 | 0.1011 | 0.5359 |
| isodense | conv0.7 | 15 | 0.7163 | 0.2526 | 0.7163 |
| isodense | dipole0.8 | 17 | 0.5658 | 0.1205 | 0.5658 |
| mixed | coh | 7 | 0.5875 | 0.0044 | 0.5876 |
| mixed | conv0.7 | 11 | 0.6812 | 0.0370 | 0.6812 |
| mixed | dipole0.8 | 13 | 0.6020 | 0.0141 | 0.6020 |
| dense90 | coh | 5 | 0.5279 | 0.1925 | 0.5279 |
| dense90 | conv0.7 | 9 | 0.7632 | 0.5253 | 0.7632 |
| dense90 | dipole0.8 | 9 | 0.5279 | 0.1925 | 0.5279 |

### Table 8 — 320 nm 2D layouts: strong single-mask search, optimal-face rank, and two-exposure realization
| layout | F | bound ρ | min face rank | best single mask ρ | single gap | two-exposure ρ | two-exposure gap |
|---|---|---|---|---|---|---|---|
| Tjunction320_annular | 61 | 0.8050 | None | 0.8066 | 0.20% | 0.8652 | 7.48% |
| Tjunction320_quasar | 49 | 0.7817 | None | 0.7911 | 1.20% | 0.8368 | 7.05% |
| lineend320_annular | 61 | 0.8677 | None | 0.8681 | 0.04% | 0.8950 | 3.14% |

### Table 9 — Rigorous branch-and-bound certification of the single-mask optimum (real masks) and certified multi-exposure gain
| case | F | SDP multi-exposure bound | B&B certified single-mask lower bound L | best single mask U | B&B gap (U-L)/U | nodes | converged | time (s) | certified gain (L-bound)/bound |
|---|---|---|---|---|---|---|---|---|---|
| trial82 | 5 | 0.3597 | 0.4061 | 0.4066 | 0.10% | 40383 | True | 1.0 | 12.92% |
| trial91 | 5 | 0.6049 | 0.6775 | 0.6782 | 0.10% | 76027 | True | 3.2 | 12.00% |
| trial110 | 5 | 0.3464 | 0.4033 | 0.4037 | 0.10% | 17413 | True | 0.3 | 16.44% |
| trial50 | 5 | 0.3141 | 0.3343 | 0.3347 | 0.10% | 28517 | True | 0.5 | 6.45% |
| mixed_coh | 7 | 0.5875 | 0.6270 | 0.6439 | 2.63% | 6002549 | False | 762.4 | 6.71% |
| isodense_coh | 9 | 0.5359 | 0.1496 | 0.6383 | 76.56% | 6005751 | False | 372.0 | -72.08% |
| mixed_conv0.7_PW80_edgeonly | 11 | 0.8370 | -0.0000 | 0.8451 | 100.00% | 6007797 | False | 172.0 | -100.00% |

### Independent certificate verification (exp15): 143 instances recomputed; 50-digit Cholesky of the reduced certificate matrix succeeded on 143; degenerate (S≈0, bound=1, uniform certificate) on []; min rigorization scale c = 0.0061; CLARABEL/SCS agreement: median |Δt| = 4.8e-09, instances with |Δt|>1e-5: random_trial16 (CLARABEL 0.0664, SCS 0.0712), random_trial75 (CLARABEL 0.0313, SCS 0.0002), random_trial105 (CLARABEL 0.0599, SCS 0.0599), random_trial111 (CLARABEL 0.1082, SCS 0.1082)

### Table 10 — Vector (high-NA, n_resist=1.7) vs scalar bounds, 1D
| case | source | scalar bound | scalar gap | TE bound | TE gap | TM bound | TM gap | unpol bound | unpol gap |
|---|---|---|---|---|---|---|---|---|---|
| dense90 | conv0.7 | 0.6995 | 0.25% | 0.7429 | 0.15% | 0.7756 | 0.00% | 0.7970 | 0.03% |
| dense90 | dipole0.8 | 0.4217 | 0.00% | 0.4258 | 0.00% | 0.8459 | 0.00% | 0.6089 | 0.00% |
| isodense | conv0.7 | 0.5878 | 0.78% | 0.6230 | 0.63% | 0.6631 | 0.37% | 0.6773 | 0.22% |
| isodense | dipole0.8 | 0.4028 | 1.28% | 0.4013 | 1.70% | 0.7212 | 0.00% | 0.5419 | 0.64% |

### Table 11 — 2D contacts: vector vs scalar
| case | scalar bound | scalar gap | unpol bound | unpol gap | TE bound | TE gap |
|---|---|---|---|---|---|---|
| contacts120_annular | 0.8749 | 0.18% | 0.9212 | 1.03% | 0.8808 | 0.16% |
| contacts120_quasar | 0.8385 | 0.11% | 0.8940 | 0.84% | 0.8458 | 0.10% |

### Aberration sensitivity (iso-dense, σ=0.7): Z5 0.0λ: 0.7163; Z5 0.02λ: 0.7165; Z5 0.05λ: 0.7178; Z5 0.1λ: 0.7225; Z7 0.0λ: 0.7163; Z7 0.02λ: 0.7182; Z7 0.05λ: 0.7227; Z7 0.1λ: 0.7327; Z9 0.0λ: 0.7163; Z9 0.02λ: 0.7168; Z9 0.05λ: 0.7197; Z9 0.1λ: 0.7299

### Pixelation (iso-dense, σ=0.7, F=15): 4 nm pixels (N=150): 0.7163; 10 nm pixels (N=60): 0.7163; 20 nm pixels (N=30): 0.7163; 40 nm pixels (N=15): 0.7163; 50 nm pixels (N=12): 0.7196; 67 nm pixels (N=9): 0.7714 (Fourier 0.7163)

### Random 2D layouts (n=24; cells 200/240 nm; annular/quasar/conv0.5; τ=12 nm): gap ≤1e-3: 45.8%; ≤1%: 50.0%; ≤5%: 79.2%; median 1.0e-02; max 0.111; F range 13–37; median SDP time 5.5 s

### Table 12 — Adversarial 1D targets
| target | source | F | bound ρ | max EL (full range) | best single ρ | gap | rank |
|---|---|---|---|---|---|---|---|
| nearcutoff_mixed_p88_p130 | coh | 7 | 0.6183 | 0.472 | 0.6804 | 10.05% | 3 |
| nearcutoff_mixed_p88_p130 | conv0.7 | 11 | 0.7975 | 0.225 | 0.7975 | 0.00% | 2 |
| nearcutoff_mixed_p88_p130 | dipole0.8 | 11 | 0.7600 | 0.273 | 0.7684 | 1.10% | 2 |
| sraf_like_narrow_space | coh | 7 | 0.7126 | 0.336 | 0.8239 | 15.62% | 3 |
| sraf_like_narrow_space | conv0.7 | 11 | 0.9129 | 0.091 | 0.9129 | 0.00% | 1 |
| sraf_like_narrow_space | dipole0.8 | 13 | 0.7699 | 0.260 | 0.7699 | 0.00% | 1 |
| extreme_duty_5pct | coh | 5 | 0.9248 | 0.078 | 0.9248 | 0.00% | 1 |
| extreme_duty_5pct | conv0.7 | 9 | 0.9365 | 0.066 | 0.9384 | 0.20% | 2 |
| extreme_duty_5pct | dipole0.8 | 11 | 0.8989 | 0.106 | 0.8989 | 0.00% | 1 |
| extreme_duty_95pct | coh | 5 | 0.0000 | 2.000 | 0.0077 | degenerate | 3 |
| extreme_duty_95pct | conv0.7 | 9 | 0.5631 | 0.559 | 0.5915 | 5.04% | 2 |
| extreme_duty_95pct | dipole0.8 | 11 | 0.0000 | 2.000 | 0.0015 | degenerate | 5 |
| alternating_widths | coh | 9 | 0.7588 | 0.274 | 1.0518 | 38.61% | 5 |
| alternating_widths | conv0.7 | 15 | 0.9013 | 0.104 | 0.9013 | 0.00% | 1 |
| alternating_widths | dipole0.8 | 17 | 0.8099 | 0.210 | 0.8116 | 0.21% | 2 |
| three_pitches | coh | 11 | 0.5983 | 0.503 | 0.7379 | 23.32% | 5 |
| three_pitches | conv0.7 | 17 | 0.7510 | 0.284 | 0.7583 | 0.98% | 2 |
| three_pitches | dipole0.8 | 19 | 0.6648 | 0.403 | 0.6651 | 0.06% | 3 |

### Table 13 — EUV 2D contact arrays (13.5 nm, NA 0.33, ±40 nm focus, τ = 0.1 pitch)
| pitch | k1 | source | F | bound ρ | max EL (full range) | best single ρ | gap | rank |
|---|---|---|---|---|---|---|---|---|
| 44 | 0.538 | annular | 13 | 0.7698 | 0.260 | 0.7698 | 0.00% | 1 |
| 44 | 0.538 | quasar | 9 | 0.7456 | 0.292 | 0.7456 | 0.00% | 1 |
| 40 | 0.489 | annular | 9 | 0.7722 | 0.257 | 0.7722 | 0.00% | 1 |
| 40 | 0.489 | quasar | 9 | 0.7477 | 0.289 | 0.7477 | 0.00% | 1 |
| 36 | 0.440 | annular | 9 | 0.8171 | 0.201 | 0.8171 | 0.00% | 1 |
| 36 | 0.440 | quasar | 9 | 0.7535 | 0.281 | 0.7535 | 0.00% | 1 |
| 32 | 0.391 | annular | 9 | 0.9099 | 0.094 | 0.9099 | 0.00% | 1 |
| 32 | 0.391 | quasar | 9 | 0.8444 | 0.169 | 0.8444 | 0.00% | 1 |