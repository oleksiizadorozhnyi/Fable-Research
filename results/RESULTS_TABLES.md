### Table 1 — 1D feature library (DUV 193 nm, NA 1.35; edge tolerance 8 nm; dense90: 6 nm)
| case | F | foci | certified bound ρ | max EL (full range) | best single mask ρ | gap | SDP time (s) |
|---|---|---|---|---|---|---|---|
| isodense_coh | 9 | 1 | 0.5359 | 0.604 | 0.6383 | 19.12% | 0.58 |
| mixed_coh | 7 | 1 | 0.5875 | 0.520 | 0.6439 | 9.59% | 0.68 |
| isodense_conv0.7 | 13 | 1 | 0.7183 | 0.328 | 0.7209 | 0.37% | 0.99 |
| mixed_conv0.7 | 11 | 1 | 0.6814 | 0.379 | 0.6814 | 0.00% | 0.69 |
| isodense_dipole0.8 | 17 | 1 | 0.5760 | 0.538 | 0.5779 | 0.34% | 1.88 |
| mixed_dipole0.8 | 13 | 1 | 0.6018 | 0.497 | 0.6018 | 0.00% | 1.32 |
| isodense_conv0.7_PW80 | 13 | 3 | 0.8312 | 0.184 | 0.8737 | 5.12% | 3.84 |
| mixed_conv0.7_PW80 | 11 | 3 | 0.8362 | 0.178 | 0.8467 | 1.26% | 1.47 |
| isodense_dipole0.8_PW80 | 17 | 3 | 0.6973 | 0.357 | 0.7407 | 6.22% | 4.80 |
| dense90_conv0.7 | 9 | 1 | 0.7658 | 0.265 | 0.7669 | 0.15% | 0.51 |
| dense90_dipole0.8 | 9 | 1 | 0.5279 | 0.618 | 0.5279 | 0.00% | 0.48 |
| dense90_coh | 5 | 1 | 0.5279 | 0.618 | 0.5279 | 0.00% | 1.13 |

### Table 2 — Certified multi-exposure gain (Lasserre-2 single-mask bound vs SDP multi-exposure bound)
| case | F | SDP (multi-exposure) bound | Lasserre-2 single-mask bound | best single mask found | certified gain | status |
|---|---|---|---|---|---|---|
| isodense_coh | 9 | 0.5359 | 0.6383 | 0.6383 | 19.12% | optimal_inaccurate |
| mixed_coh | 7 | 0.5875 | 0.6441 | 0.6439 | 9.63% | optimal_inaccurate |
| mixed_conv0.7_PW80_edgeonly | 11 | 0.8362 | 0.8467 | 0.8467 | 1.26% | optimal_inaccurate |

### Random 1D targets (n=109): gap ≤1e-4: 72.5%; ≤1%: 85.3%; ≤5%: 92.7%; median 4.8e-07; max 2033.710; worst certificate min-eig 1.0e-10

### Table 3 — Pixel vs Fourier relaxations, PSM vs binary (DNN), ILT baseline (σ=0.7 conventional)
| case | F | Fourier PSM bound | pixel PSM bound | pixel nonneg (DNN) bound | ILT-PSM ratio | ILT-PSM refined | ILT-binary ratio | best PSM (pipeline) |
|---|---|---|---|---|---|---|---|---|
| isodense | 13 | 0.7183 | 0.7184 | 0.7557 | 0.886 | 0.7370 | 1.049 | 0.7209 |
| mixed | 11 | 0.6814 | 0.6859 | 0.7507 | 0.944 | 0.7777 | 0.963 | 0.6814 |
| dense120 | 9 | 0.6715 | 0.6717 | 0.6943 | 0.847 | 0.6774 | 0.859 | 0.6763 |

### Table 3b — Standard gradient ILT (sigmoid resist, pattern-error loss; best of 6 settings) vs certified bound
| case | F | certified bound ρ | ILT-PSM best ρ | ILT-PSM then ratio-refined | ILT-binary best ρ | ILT-binary then ratio-refined |
|---|---|---|---|---|---|---|
| isodense | 13 | 0.7183 | 0.8620 | 0.7370 | 0.9529 | 0.7543 |
| mixed | 11 | 0.6814 | 0.8977 | 0.7777 | 0.9205 | 0.7777 |
| dense120 | 9 | 0.6715 | 0.8266 | 0.6836 | 0.8534 | 0.6836 |

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
| 22 | 0.269 | conv0.9 | 7 | 0.9462 | 0.055 | 0.9746 | 3.00% |
| 22 | 0.269 | dipole0.8 | 9 | 0.7262 | 0.317 | 0.7262 | 0.00% |
| 24 | 0.293 | conv0.9 | 9 | 0.8347 | 0.180 | 0.8347 | 0.00% |
| 24 | 0.293 | dipole0.8 | 9 | 0.5779 | 0.535 | 0.5779 | 0.00% |
| 26 | 0.318 | conv0.9 | 9 | 0.7768 | 0.251 | 0.7768 | 0.00% |
| 26 | 0.318 | dipole0.8 | 9 | 0.5470 | 0.586 | 0.5470 | 0.00% |
| 28 | 0.342 | conv0.9 | 11 | 0.7424 | 0.296 | 0.7424 | 0.00% |
| 28 | 0.342 | dipole0.8 | 11 | 0.5562 | 0.570 | 0.5562 | 0.00% |
| 32 | 0.391 | conv0.9 | 11 | 0.6899 | 0.367 | 0.6899 | 0.00% |
| 32 | 0.391 | dipole0.8 | 13 | 0.5910 | 0.514 | 0.5910 | 0.00% |
| 40 | 0.489 | conv0.9 | 15 | 0.6295 | 0.455 | 0.6295 | 0.00% |
| 40 | 0.489 | dipole0.8 | 15 | 0.5799 | 0.532 | 0.5979 | 3.10% |

Theorem 4 numerics: necessity violations 0/2000; sufficiency failures 0/200.

### Robustness (iso-dense, σ=0.7)
- interior spacing sweep (bound): 16 nm: 0.71823, 8 nm: 0.71829, 4 nm: 0.71830, 2 nm: 0.71831, 1 nm: 0.71831
- margin sweep (bound / EL): 2 nm: 0.8757/0.066, 2.5 nm: 0.8472/0.083, 4 nm: 0.7672/0.132, 5 nm: 0.7183/0.164, 8 nm: 0.5904/0.258, 12 nm: 0.4573/0.372, 16 nm: 0.3599/0.471, 20 nm: 0.2927/0.547
- source quadrature & solver agreement: 37 pts: CLARABEL 0.70040 SCS 0.70040, 97 pts: CLARABEL 0.70907 SCS 0.70907, 177 pts: CLARABEL 0.71633 SCS 0.71633, 349 pts: CLARABEL 0.71830 SCS 0.71830, 749 pts: CLARABEL 0.71304 SCS 0.71304, 1313 pts: CLARABEL 0.71443 SCS 0.71443
- 50-digit certificate check: Cholesky success=True, min eig=8.40e-11, bound=0.590411

### Table 6 — 2D layouts (DUV, NA 1.35), edge tolerance 12 nm (contacts120/100: 10/8 nm; T/line-end: 14 nm)
| layout | cell (nm) | source pts | F | samples | bound ρ | max EL (full range) | best single mask ρ | gap | rank | cert. null dim | SDP time (s) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| contacts160_annular | [160, 160] | 268 | 13 | 352 | 0.8058 | 0.215 | 0.8236 | 2.21% | 1 | 0 | 1.7 |
| contacts160_conv0.5 | [160, 160] | 349 | 9 | 352 | 0.8401 | 0.174 | 0.8401 | 0.00% | 1 | 0 | 1.7 |
| contacts120_annular | [120, 120] | 268 | 9 | 243 | 0.8750 | 0.133 | 0.8766 | 0.19% | 1 | 0 | 1.2 |
| contacts120_quasar | [120, 120] | 140 | 9 | 243 | 0.8386 | 0.176 | 0.8395 | 0.11% | 1 | 0 | 1.2 |
| contacts100_quasar | [100, 100] | 140 | 5 | 211 | 0.9776 | 0.023 | 0.9776 | 0.00% | 1 | 0 | 1.0 |
| Tjunction320_annular | [320, 320] | 268 | 57 | 1273 | 0.8035 | 0.218 | 0.8652 | 7.67% | 2 | 0 | 1898.9 |
| Tjunction320_quasar | [320, 320] | 140 | 45 | 1273 | 0.7827 | 0.244 | 0.8515 | 8.79% | 4 | 0 | 716.1 |
| lineend320_annular | [320, 320] | 268 | 57 | 1288 | 0.8665 | 0.143 | 0.8935 | 3.11% | 2 | 0 | 1277.2 |

### Outlier / process-window re-examination (strong primal search, face exploration, Lasserre-2 where F≤11)
| case | F | bound | old best | new best | new gap | min face rank | Lasserre-2 | status |
|---|---|---|---|---|---|---|---|---|
| trial16 | 13 | 0.0000 | 0.0730 | 0.0183 | 50974.90% | 2 | None | skipped(F>11) |
| trial82 | 5 | 0.6035 | 0.6683 | 0.6683 | 10.75% | 2 | 0.6683098422692263 | optimal_inaccurate |
| trial91 | 5 | 0.7787 | 0.8300 | 0.8300 | 6.59% | 2 | 0.8300430332529734 | optimal_inaccurate |
| trial110 | 5 | 0.6268 | 0.8619 | 0.8619 | 37.51% | 2 | 0.8618960894106031 | optimal_inaccurate |
| trial50 | 5 | 0.5572 | 0.5835 | 0.5835 | 4.71% | 2 | 0.5834942686814835 | optimal_inaccurate |
| trial46 | 13 | 0.2859 | 0.2859 | 0.2859 | 0.00% | 1 | None | skipped(F>11) |
| mixed_conv0.7_PW80 | 11 | 0.8362 | nan | 0.8467 | 1.26% | 3 | None | skipped(three foci: Lasserre-2 exceeds memory; see exp11b for the edge-only through-focus instance) |
| isodense_conv0.7_PW80 | 13 | 0.8312 | nan | 0.8668 | 4.28% | 3 | None | skipped(three foci: Lasserre-2 exceeds memory; see exp11b for the edge-only through-focus instance) |
| isodense_dipole0.8_PW80 | 17 | 0.6973 | nan | 0.7407 | 6.22% | 3 | None | skipped(three foci: Lasserre-2 exceeds memory; see exp11b for the edge-only through-focus instance) |

### Table 7 — Baseline Rayleigh-quotient bound vs SDP; complex-transmission masks
| case | source | F | SDP real-mask bound | Rayleigh-quotient bound | SDP complex-mask bound |
|---|---|---|---|---|---|
| isodense | coh | 9 | 0.5359 | 0.1011 | 0.5359 |
| isodense | conv0.7 | 13 | 0.7183 | 0.2559 | 0.7183 |
| isodense | dipole0.8 | 17 | 0.5760 | 0.1243 | 0.5760 |
| mixed | coh | 7 | 0.5875 | 0.0044 | 0.5876 |
| mixed | conv0.7 | 11 | 0.6814 | 0.0368 | 0.6814 |
| mixed | dipole0.8 | 13 | 0.6018 | 0.0149 | 0.6018 |
| dense90 | coh | 5 | 0.5279 | 0.1925 | 0.5279 |
| dense90 | conv0.7 | 9 | 0.7658 | 0.5298 | 0.7658 |
| dense90 | dipole0.8 | 9 | 0.5279 | 0.1925 | 0.5279 |

### Table 8 — 320 nm 2D layouts: strong single-mask search, optimal-face rank, and two-exposure realization
| layout | F | bound ρ | min face rank | best single mask ρ | single gap | two-exposure ρ | two-exposure gap |
|---|---|---|---|---|---|---|---|
| Tjunction320_annular | 57 | 0.8035 | None | 0.8055 | 0.25% | 0.8636 | 7.47% |
| Tjunction320_quasar | 45 | 0.7827 | None | 0.7919 | 1.18% | 0.8355 | 6.75% |
| lineend320_annular | 57 | 0.8665 | None | 0.8667 | 0.02% | 0.8935 | 3.12% |

### Table 9 — Rigorous branch-and-bound certification of the single-mask optimum (real masks) and certified multi-exposure gain
| case | F | SDP multi-exposure bound | B&B certified single-mask lower bound L | best single mask U | B&B gap (U-L)/U | nodes | converged | time (s) | certified gain (L-bound)/bound |
|---|---|---|---|---|---|---|---|---|---|
| trial82 | 5 | 0.3597 | 0.4061 | 0.4066 | 0.10% | 40383 | True | 1.0 | 12.92% |
| trial91 | 5 | 0.6049 | 0.6775 | 0.6782 | 0.10% | 76027 | True | 1.7 | 12.00% |
| trial110 | 5 | 0.3464 | 0.4033 | 0.4037 | 0.10% | 17413 | True | 0.4 | 16.44% |
| trial50 | 5 | 0.3141 | 0.3343 | 0.3347 | 0.10% | 28517 | True | 0.5 | 6.45% |
| mixed_coh | 7 | 0.5875 | 0.6270 | 0.6439 | 2.63% | 6002549 | False | 391.2 | 6.71% |
| isodense_coh | 9 | 0.5359 | 0.1496 | 0.6383 | 76.56% | 6005751 | False | 590.6 | -72.08% |
| mixed_conv0.7_PW80_edgeonly | 11 | 0.8362 | -0.0000 | 0.8467 | 100.00% | 6007797 | False | 176.3 | -100.00% |

### Independent certificate verification (exp15): 143 instances recomputed; 50-digit Cholesky of the reduced certificate matrix succeeded on 142; degenerate (S≈0, bound=1, uniform certificate) on ['EUV_p20_conv0.9']; min rigorization scale c = 0.0005; CLARABEL/SCS agreement: median |Δt| = 4.3e-09, instances with |Δt|>1e-5: random_trial75 (CLARABEL 0.0312, SCS 0.0086), random_trial111 (CLARABEL 0.1082, SCS 0.1082)

### Table 10 — Vector (high-NA, n_resist=1.7) vs scalar bounds, 1D
| case | source | scalar bound | scalar gap | TE bound | TE gap | TM bound | TM gap | unpol bound | unpol gap |
|---|---|---|---|---|---|---|---|---|---|
| dense90 | conv0.7 | 0.7028 | 0.19% | 0.7420 | 0.13% | 0.7798 | 0.00% | 0.7990 | 0.02% |
| dense90 | dipole0.8 | 0.4217 | 0.00% | 0.4262 | 0.00% | 0.8459 | 0.00% | 0.6091 | 0.00% |
| isodense | conv0.7 | 0.5904 | 0.62% | 0.6228 | 0.63% | 0.6662 | 0.26% | 0.6792 | 0.39% |
| isodense | dipole0.8 | 0.4118 | 0.21% | 0.4113 | 0.40% | 0.7216 | 0.00% | 0.5474 | 0.09% |

### Table 11 — 2D contacts: vector vs scalar
| case | scalar bound | scalar gap | unpol bound | unpol gap | TE bound | TE gap |
|---|---|---|---|---|---|---|
| contacts120_annular | 0.8750 | 0.19% | 0.9215 | 1.05% | 0.8810 | 0.16% |
| contacts120_quasar | 0.8386 | 0.11% | 0.8943 | 0.84% | 0.8445 | 0.10% |

### Aberration sensitivity (iso-dense, σ=0.7): Z5 0.0λ: 0.7183; Z5 0.02λ: 0.7185; Z5 0.05λ: 0.7199; Z5 0.1λ: 0.7246; Z7 0.0λ: 0.7183; Z7 0.02λ: 0.7204; Z7 0.05λ: 0.7248; Z7 0.1λ: 0.7345; Z9 0.0λ: 0.7183; Z9 0.02λ: 0.7188; Z9 0.05λ: 0.7213; Z9 0.1λ: 0.7309

### Pixelation (iso-dense, σ=0.7, F=15): 4 nm pixels (N=150): 0.7183; 10 nm pixels (N=60): 0.7183; 20 nm pixels (N=30): 0.7183; 40 nm pixels (N=15): 0.7183; 50 nm pixels (N=12): 0.7224; 67 nm pixels (N=9): 0.7732 (Fourier 0.7183)

### Random 2D layouts (n=24; cells 200/240 nm; annular/quasar/conv0.5; τ=12 nm): gap ≤1e-3: 45.8%; ≤1%: 50.0%; ≤5%: 83.3%; median 1.0e-02; max 0.106; F range 13–37; median SDP time 7.1 s

### Table 12 — Adversarial 1D targets
| target | source | F | bound ρ | max EL (full range) | best single ρ | gap | rank |
|---|---|---|---|---|---|---|---|
| nearcutoff_mixed_p88_p130 | coh | 7 | 0.6183 | 0.472 | 0.6804 | 10.05% | 3 |
| nearcutoff_mixed_p88_p130 | conv0.7 | 11 | 0.7910 | 0.233 | 0.7911 | 0.01% | 2 |
| nearcutoff_mixed_p88_p130 | dipole0.8 | 11 | 0.7429 | 0.295 | 0.7615 | 2.51% | 2 |
| sraf_like_narrow_space | coh | 7 | 0.7126 | 0.336 | 0.8239 | 15.62% | 3 |
| sraf_like_narrow_space | conv0.7 | 11 | 0.9151 | 0.089 | 0.9151 | 0.00% | 1 |
| sraf_like_narrow_space | dipole0.8 | 13 | 0.7702 | 0.260 | 0.7702 | 0.00% | 1 |
| extreme_duty_5pct | coh | 5 | 0.9248 | 0.078 | 0.9248 | 0.00% | 1 |
| extreme_duty_5pct | conv0.7 | 9 | 0.9364 | 0.066 | 0.9385 | 0.23% | 2 |
| extreme_duty_5pct | dipole0.8 | 11 | 0.8981 | 0.107 | 0.8981 | 0.00% | 1 |
| extreme_duty_95pct | coh | 5 | 0.0000 | 2.000 | 0.0078 | degenerate | 3 |
| extreme_duty_95pct | conv0.7 | 9 | 0.5703 | 0.547 | 0.5959 | 4.50% | 2 |
| extreme_duty_95pct | dipole0.8 | 11 | 0.0000 | 2.000 | 0.0079 | degenerate | 5 |
| alternating_widths | coh | 9 | 0.7588 | 0.274 | 1.0518 | 38.61% | 5 |
| alternating_widths | conv0.7 | 13 | 0.9029 | 0.102 | 0.9029 | 0.00% | 2 |
| alternating_widths | dipole0.8 | 17 | 0.8155 | 0.203 | 0.8156 | 0.00% | 1 |
| three_pitches | coh | 11 | 0.5983 | 0.503 | 0.7379 | 23.32% | 5 |
| three_pitches | conv0.7 | 17 | 0.7514 | 0.284 | 0.7549 | 0.47% | 2 |
| three_pitches | dipole0.8 | 19 | 0.6644 | 0.403 | 0.6667 | 0.35% | 9 |

### Table 13 — EUV 2D contact arrays (13.5 nm, NA 0.33, ±40 nm focus, τ = 0.1 pitch)
| pitch | k1 | source | F | bound ρ | max EL (full range) | best single ρ | gap | rank |
|---|---|---|---|---|---|---|---|---|
| 44 | 0.538 | annular | 13 | 0.7602 | 0.272 | 0.7602 | 0.00% | 1 |
| 44 | 0.538 | quasar | 9 | 0.7466 | 0.290 | 0.7466 | 0.00% | 1 |
| 40 | 0.489 | annular | 9 | 0.7728 | 0.256 | 0.7728 | 0.00% | 1 |
| 40 | 0.489 | quasar | 9 | 0.7468 | 0.290 | 0.7468 | 0.00% | 1 |
| 36 | 0.440 | annular | 9 | 0.8216 | 0.196 | 0.8216 | 0.00% | 1 |
| 36 | 0.440 | quasar | 9 | 0.7576 | 0.276 | 0.7576 | 0.00% | 1 |
| 32 | 0.391 | annular | 9 | 0.9104 | 0.094 | 0.9104 | 0.00% | 1 |
| 32 | 0.391 | quasar | 9 | 0.8413 | 0.172 | 0.8413 | 0.00% | 1 |