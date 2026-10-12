#!/bin/bash
# Clean reproduction in 4 parallel groups (4 cores), then verification, tables and figures.  cd code && ./run_all_par.sh
cd "$(dirname "$0")"; mkdir -p ../results/logs; L=../results/logs/run_all.log; : > $L
python3 -W ignore test_certlitho.py > ../results/logs/tests.log 2>&1 || { echo "TESTS FAILED" >> $L; exit 1; }
echo "=== tests passed $(date)" >> $L
run() { for e in "$@"; do echo "=== START $e $(date)" >> $L; python3 -W ignore $e.py > ../results/logs/$e.log 2>&1; echo "=== END $e $(date) exit=$?" >> $L; done; }
run exp1_1d_features exp8_robustness exp5_smo_euv_toeplitz exp7_lasserre exp2_random_tightness exp9_outlier_pw exp10_baselines_complex exp11b_lasserre_pw exp12_lasserre_xcheck exp14_bnb &
run exp3_pixel_binary exp3b_ilt exp16_physics &
run exp4_2d exp4b_recert exp13_2d_strong &
run exp17_dataset exp14b_bnb_big &
wait
run exp15_verify_all
python3 assemble.py > /dev/null && echo "=== assembled $(date)" >> $L
run figs fig2d fig5 fig6
echo "=== ALL DONE $(date)" >> $L
