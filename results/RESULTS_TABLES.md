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
| mixed_coh | 7 | 0.6235 | 0.6408 | 0.6408 | 2.78% | optimal_inaccurate |
| isodense_coh | 9 | 0.5645 | 0.5738 | 0.5738 | 1.63% | optimal_inaccurate |

### Random 1D targets (n=6): gap ≤1e-4: 100.0%; ≤1%: 100.0%; ≤5%: 100.0%; median 2.3e-07; max 0.000; worst certificate min-eig 4.1e-18

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
| layout | cell (nm) | source pts | F | samples | bound ρ | max EL | best single mask ρ | gap | rank | SDP time (s) |
|---|---|---|---|---|---|---|---|---|---|---|
| contacts160_annular | [160, 160] | 28 | 13 | 260 | 0.2324 | 0.623 | 0.2324 | 0.00% | 1 | 0.3 |
| contacts160_conv0.5 | [160, 160] | 13 | 9 | 260 | 0.2661 | 0.580 | 0.2661 | 0.00% | 1 | 0.2 |
| contacts120_annular | [120, 120] | 28 | 9 | 156 | 0.3555 | 0.475 | 0.3555 | 0.00% | 1 | 0.2 |
| contacts120_quasar | [120, 120] | 20 | 9 | 156 | 0.1932 | 0.676 | 0.1932 | 0.00% | 1 | 0.2 |
| contacts100_quasar | [100, 100] | 20 | 5 | 156 | 0.6193 | 0.235 | 0.6193 | 0.00% | 1 | 0.2 |

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