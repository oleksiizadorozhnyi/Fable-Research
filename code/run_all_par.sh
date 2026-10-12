#!/bin/bash
# Clean reproduction in 4 parallel groups (4 cores), then verification, tables and figures.  cd code && ./run_all_par.sh
cd "$(dirname "$0")"; mkdir -p ../results/logs; L=../results/logs/run_all.log; : > $L
python3 -W ignore test_certlitho.py > ../results/logs/tests.log 2>&1 || { echo "TESTS FAILED" >> $L; exit 1; }
echo "=== tests passed $(date)" >> $L
run() { for e in "$@"; do echo "=== START $e $(date)" >> $L; python3 -W ignore $e.py > ../results/logs/$e.log 2>&1; rc=$?; echo "=== END $e $(date) exit=$rc" >> $L; [ $rc -ne 0 ] && echo "FAILED $e rc=$rc" >> $L; done; }
run exp1_1d_features exp7_lasserre exp8_robustness exp5_smo_euv_toeplitz exp3_pixel_binary exp3b_ilt &
run exp2_random_tightness exp9_outlier_pw exp14_bnb exp14b_bnb_big exp19_bnb_sdp exp11b_lasserre_pw exp12_lasserre_xcheck &
run exp4_2d exp4b_recert exp13_2d_strong exp18_quasar_strong exp10_baselines_complex exp20_corner_exclusion &
run exp17_dataset exp16_physics &
wait
run exp15_verify_all
python3 assemble.py > /dev/null && echo "=== assembled $(date)" >> $L
run figs fig2d fig5 fig6
grep -q FAILED $L && echo "=== DONE WITH FAILURES $(date)" >> $L || echo "=== ALL DONE $(date)" >> $L
