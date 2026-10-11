### Table 1 — 1D feature library (DUV 193 nm, NA 1.35; edge tolerance 8 nm; dense90: 6 nm)
| case | F | foci | certified bound ρ | max EL | best single mask ρ | gap | SDP time (s) |
|---|---|---|---|---|---|---|---|
| isodense_coh | 9 | 1 | 0.3625 | 0.468 | 0.4741 | 30.78% | 0.16 |
| mixed_coh | 7 | 1 | 0.4211 | 0.407 | 0.4787 | 13.68% | 0.16 |
| isodense_conv0.7 | 15 | 1 | 0.4691 | 0.361 | 0.4808 | 2.50% | 0.20 |
| mixed_conv0.7 | 11 | 1 | 0.4825 | 0.349 | 0.4836 | 0.22% | 0.17 |
| isodense_dipole0.8 | 15 | 1 | 0.4313 | 0.397 | 0.4313 | 0.00% | 0.21 |
| mixed_dipole0.8 | 13 | 1 | 0.4211 | 0.407 | 0.4211 | 0.00% | 0.36 |
| isodense_conv0.7_PW80 | 15 | 3 | 0.8131 | 0.103 | 0.9044 | 11.22% | 0.69 |
| mixed_conv0.7_PW80 | 11 | 3 | 0.8691 | 0.070 | 0.9866 | 13.52% | 0.42 |
| isodense_dipole0.8_PW80 | 15 | 3 | 0.7686 | 0.131 | 0.8972 | 16.73% | 0.57 |
| dense90_conv0.7 | 9 | 1 | 0.6428 | 0.217 | 0.6428 | 0.00% | 0.12 |
| dense90_dipole0.8 | 9 | 1 | 0.4217 | 0.407 | 0.4217 | 0.00% | 0.08 |
| dense90_coh | 5 | 1 | 0.4217 | 0.407 | 0.4217 | 0.00% | 0.14 |

### Table 2 — Certified multi-exposure gain (Lasserre-2 single-mask bound vs SDP multi-exposure bound)
| case | F | SDP (multi-exposure) bound | Lasserre-2 single-mask bound | best single mask found | certified gain | status |
|---|---|---|---|---|---|---|
| isodense_coh | 9 | 0.3625 | 0.4608 | 0.4741 | 27.10% | optimal_inaccurate |
| mixed_coh | 7 | 0.4211 | 0.4787 | 0.4787 | 13.68% | optimal_inaccurate |
| mixed_conv0.7_PW80_edgeonly | 11 | 0.8691 | 0.9552 | 0.9568 | 9.91% | optimal_inaccurate |

### Random 1D targets (n=109): gap ≤1e-4: 69.7%; ≤1%: 82.6%; ≤5%: 91.7%; median 1.8e-06; max 0.173; worst certificate min-eig 1.3e-19

### Table 3 — Pixel vs Fourier relaxations, PSM vs binary (DNN), ILT baseline (σ=0.7 conventional)
| case | F | Fourier PSM bound | pixel PSM bound | pixel nonneg (DNN) bound | ILT-PSM ratio | ILT-PSM refined | ILT-binary ratio | best PSM (pipeline) |
|---|---|---|---|---|---|---|---|---|
| isodense | 15 | 0.3879 | 0.3879 | 0.4266 | 6.173 | 0.4332 | 227.751 | 0.3981 |
| mixed | 11 | 0.4031 | 0.4031 | 0.4599 | 29.300 | 0.5174 | 158.776 | 0.4036 |
| dense120 | 9 | 0.3671 | 0.3671 | 0.3777 | 5.956 | 0.3671 | 104.414 | 0.3671 |

### Table 3b — Standard gradient ILT (sigmoid resist, pattern-error loss; best of 6 settings) vs certified bound
| case | F | certified bound ρ | ILT-PSM best ρ | ILT-PSM then ratio-refined | ILT-binary best ρ | ILT-binary then ratio-refined |
|---|---|---|---|---|---|---|
| isodense | 15 | 0.3879 | 0.5526 | 0.4257 | 0.5789 | 0.4257 |
| mixed | 11 | 0.4031 | 0.6715 | 0.5302 | 0.7218 | 0.5301 |
| dense120 | 9 | 0.3671 | 0.4768 | 0.3671 | 0.5167 | 0.3671 |

### Table 4 — Source–mask entitlement bound (iso-dense, DUV)
| source | #pts | fixed-source bound | SMO relaxation bound |
|---|---|---|---|
| coh | 1 | 0.3625 | 0.3625 |
| conv0.7 | 7 | 0.4691 | 0.3625 |
| dipole0.8 | 2 | 0.4313 | 0.4313 |
| full_grid0.9 | 13 | 0.5091 | 0.3625 |

### Table 5 — EUV (13.5 nm, NA 0.33, ±50 nm focus) dense lines: certified max EL vs pitch
| pitch (nm) | k1 | source | F | bound ρ | max EL | best single mask ρ | gap |
|---|---|---|---|---|---|---|---|
| 20 | 0.244 | conv0.9 | 7 | 1.0000 | 0.000 | 1.0000 | 0.00% |
| 20 | 0.244 | dipole0.8 | 7 | 1.0000 | 0.000 | 1.0000 | 0.00% |
| 22 | 0.269 | conv0.9 | 9 | 0.4685 | 0.362 | 0.4685 | 0.00% |
| 22 | 0.269 | dipole0.8 | 7 | 1.0000 | 0.000 | 1.0000 | 0.00% |
| 24 | 0.293 | conv0.9 | 9 | 0.4733 | 0.357 | 0.4733 | 0.00% |
| 24 | 0.293 | dipole0.8 | 9 | 0.2007 | 0.666 | 0.2007 | 0.00% |
| 26 | 0.318 | conv0.9 | 9 | 0.4123 | 0.416 | 0.4123 | 0.00% |
| 26 | 0.318 | dipole0.8 | 9 | 0.1880 | 0.684 | 0.1880 | 0.00% |
| 28 | 0.342 | conv0.9 | 11 | 0.3951 | 0.434 | 0.3951 | 0.00% |
| 28 | 0.342 | dipole0.8 | 9 | 0.2044 | 0.661 | 0.2044 | 0.00% |
| 32 | 0.391 | conv0.9 | 11 | 0.3930 | 0.436 | 0.3930 | 0.00% |
| 32 | 0.391 | dipole0.8 | 11 | 0.2566 | 0.592 | 0.2566 | 0.00% |
| 40 | 0.489 | conv0.9 | 15 | 0.3282 | 0.506 | 0.3282 | 0.00% |
| 40 | 0.489 | dipole0.8 | 15 | 0.1891 | 0.682 | 0.1891 | 0.00% |

Theorem 4 numerics: necessity violations 0/2000; sufficiency failures 0/200.

### Robustness (iso-dense, σ=0.7)
- interior spacing sweep (bound): 16 nm: 0.46907, 8 nm: 0.46907, 4 nm: 0.46907, 2 nm: 0.46907, 1 nm: 0.46907
- margin sweep (bound / EL): 2 nm: 0.8281/0.094, 4 nm: 0.6855/0.187, 8 nm: 0.4691/0.361, 12 nm: 0.3209/0.514, 16 nm: 0.2197/0.640, 20 nm: 0.1516/0.737
- source quadrature & solver agreement: 3 pts: CLARABEL 0.41345 SCS 0.41345, 5 pts: CLARABEL 0.47298 SCS 0.47298, 7 pts: CLARABEL 0.46907 SCS 0.46907, 11 pts: CLARABEL 0.51612 SCS 0.51612, 15 pts: CLARABEL 0.51463 SCS 0.51463, 21 pts: CLARABEL 0.52607 SCS 0.52607
- 50-digit certificate check: Cholesky success=True, min eig=8.25e-11, bound=0.469065

### Table 6 — 2D layouts (DUV, NA 1.35), edge tolerance 12 nm (contacts120/100: 10/8 nm; T/line-end: 14 nm)
| layout | cell (nm) | source pts | F | samples | bound ρ | max EL | best single mask ρ | gap | rank | cert. null dim | SDP time (s) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| contacts160_annular | [160, 160] | 28 | 13 | 320 | 0.5843 | 0.262 | 0.5861 | 0.30% | 1 | 0 | 0.4 |
| contacts160_conv0.5 | [160, 160] | 13 | 9 | 320 | 0.6174 | 0.237 | 0.6177 | 0.04% | 1 | 0 | 0.4 |
| contacts120_annular | [120, 120] | 28 | 9 | 216 | 0.6957 | 0.179 | 0.6962 | 0.07% | 1 | 0 | 0.4 |
| contacts120_quasar | [120, 120] | 20 | 9 | 216 | 0.5740 | 0.271 | 0.5747 | 0.13% | 1 | 0 | 0.3 |
| contacts100_quasar | [100, 100] | 20 | 5 | 197 | 0.8751 | 0.067 | 0.8754 | 0.03% | 1 | 0 | 0.2 |
| Tjunction320_annular | [320, 320] | 28 | 61 | 1210 | 0.4386 | 0.390 | 0.4789 | 9.17% | 6 | 4 | 142.7 |
| Tjunction320_quasar | [320, 320] | 20 | 61 | 1210 | 0.3893 | 0.440 | 0.4444 | 14.17% | 3 | 16 | 104.6 |
| lineend320_annular | [320, 320] | 28 | 61 | 1189 | 0.5194 | 0.316 | 0.5328 | 2.59% | 6 | 4 | 204.3 |

### Outlier / process-window re-examination (strong primal search, face exploration, Lasserre-2 where F≤11)
| case | F | bound | old best | new best | new gap | min face rank | Lasserre-2 | status |
|---|---|---|---|---|---|---|---|---|
| trial16 | 13 | 0.0013 | 0.0014 | 0.0014 | 15.47% | 4 | None | skipped(F>11) |
| trial82 | 5 | 0.3597 | 0.4066 | 0.4066 | 13.03% | 2 | 0.4065533316537036 | optimal_inaccurate |
| trial91 | 5 | 0.6049 | 0.6782 | 0.6782 | 12.11% | 2 | 0.6781623229402676 | optimal_inaccurate |
| trial110 | 5 | 0.3464 | 0.4037 | 0.4037 | 16.56% | 2 | 0.40371765650199054 | optimal |
| trial50 | 5 | 0.3141 | 0.3347 | 0.3347 | 6.55% | 2 | 0.3346778891864531 | optimal_inaccurate |

### Table 7 — Baseline Rayleigh-quotient bound vs SDP; complex-transmission masks
| case | source | F | SDP real-mask bound | Rayleigh-quotient bound | SDP complex-mask bound |
|---|---|---|---|---|---|
| isodense | coh | 9 | 0.3625 | 0.0668 | 0.3625 |
| isodense | conv0.7 | 15 | 0.4691 | 0.1387 | 0.4691 |
| isodense | dipole0.8 | 15 | 0.4313 | 0.0823 | 0.4313 |
| mixed | coh | 7 | 0.4211 | 0.0028 | 0.4211 |
| mixed | conv0.7 | 11 | 0.4825 | 0.0343 | 0.4825 |
| mixed | dipole0.8 | 13 | 0.4211 | 0.0028 | 0.4211 |
| dense90 | coh | 5 | 0.4217 | 0.1582 | 0.4217 |
| dense90 | conv0.7 | 9 | 0.6428 | 0.4392 | 0.6428 |
| dense90 | dipole0.8 | 9 | 0.4217 | 0.1582 | 0.4217 |

### Table 9 — Rigorous branch-and-bound certification of the single-mask optimum (real masks) and certified multi-exposure gain
| case | F | SDP multi-exposure bound | B&B certified single-mask lower bound L | best single mask U | B&B gap (U-L)/U | nodes | converged | time (s) | certified gain (L-bound)/bound |
|---|---|---|---|---|---|---|---|---|---|
| trial82 | 5 | 0.3597 | 0.4061 | 0.4066 | 0.10% | 40383 | True | 0.7 | 12.92% |
| trial91 | 5 | 0.6049 | 0.6775 | 0.6782 | 0.10% | 76027 | True | 1.6 | 12.00% |
| trial110 | 5 | 0.3464 | 0.4033 | 0.4037 | 0.10% | 17413 | True | 0.3 | 16.44% |
| trial50 | 5 | 0.3141 | 0.3343 | 0.3347 | 0.10% | 28517 | True | 0.5 | 6.45% |
| mixed_coh | 7 | 0.4211 | 0.4667 | 0.4787 | 2.51% | 6002541 | False | 288.6 | 10.83% |

### Independent certificate verification (exp15): 136 instances recomputed; 50-digit Cholesky of the reduced certificate matrix succeeded on 134; degenerate (S≈0, bound=1, uniform certificate) on ['EUV_p20_conv0.9', 'EUV_p20_dipole0.8']; min rigorization scale c = 0.8354; CLARABEL/SCS agreement: median |Δt| = 4.5e-09, instances with |Δt|>1e-5: random_trial16 (CLARABEL 0.0015, SCS 0.0002), random_trial75 (CLARABEL 0.0034, SCS 0.0020), random_trial87 (CLARABEL 0.0013, SCS 0.0013), random_trial103 (CLARABEL 0.1419, SCS 0.0004), random_trial105 (CLARABEL 0.0054, SCS -0.0035), random_trial119 (CLARABEL 0.0063, SCS 0.0001)

### Table 10 — Vector (high-NA, n_resist=1.7) vs scalar bounds, 1D
| case | source | scalar bound | scalar gap | TE bound | TE gap | TM bound | TM gap | unpol bound | unpol gap |
|---|---|---|---|---|---|---|---|---|---|
| dense90 | conv0.7 | 0.6428 | 0.00% | 0.6274 | 0.00% | 0.9137 | 0.00% | 0.7590 | 0.00% |
| dense90 | dipole0.8 | 0.4217 | 0.00% | 0.4217 | 0.00% | 0.8465 | 0.00% | 0.6065 | 0.00% |
| isodense | conv0.7 | 0.4691 | 2.50% | 0.4509 | 1.90% | 0.7653 | 0.00% | 0.5878 | 0.04% |
| isodense | dipole0.8 | 0.4313 | 0.00% | 0.4313 | 0.00% | 0.7241 | 0.00% | 0.5522 | 0.00% |

### Table 11 — 2D contacts: vector vs scalar
| case | scalar bound | scalar gap | unpol bound | unpol gap | TE bound | TE gap |
|---|---|---|---|---|---|---|
| contacts120_annular | 0.6957 | 0.07% | 0.7674 | 0.02% | 0.7120 | 0.19% |
| contacts120_quasar | 0.5740 | 0.13% | 0.6781 | 0.01% | 0.6315 | 0.13% |

### Aberration sensitivity (iso-dense, σ=0.7): Z5 0.0λ: 0.4691; Z5 0.02λ: 0.4697; Z5 0.05λ: 0.4733; Z5 0.1λ: 0.4861; Z7 0.0λ: 0.4691; Z7 0.02λ: 0.4692; Z7 0.05λ: 0.4703; Z7 0.1λ: 0.4796; Z9 0.0λ: 0.4691; Z9 0.02λ: 0.4702; Z9 0.05λ: 0.4765; Z9 0.1λ: 0.5013

### Pixelation (iso-dense, σ=0.7, F=15): 4 nm pixels (N=150): 0.4691; 10 nm pixels (N=60): 0.4691; 20 nm pixels (N=30): 0.4691; 40 nm pixels (N=15): 0.4691; 50 nm pixels (N=12): 0.5407; 67 nm pixels (N=9): 0.6864 (Fourier 0.4691)