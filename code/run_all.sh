#!/bin/bash
# Clean reproduction of every reported number. Runs sequentially (~3-4 h on 4 cores). Usage: cd code && ./run_all.sh
cd "$(dirname "$0")"
mkdir -p ../results/logs
python3 -W ignore test_certlitho.py > ../results/logs/tests.log 2>&1 || { echo "TESTS FAILED"; exit 1; }
for e in exp1_1d_features exp8_robustness exp5_smo_euv_toeplitz exp3_pixel_binary exp3b_ilt exp7_lasserre exp2_random_tightness exp9_outlier_pw exp10_baselines_complex exp11b_lasserre_pw exp12_lasserre_xcheck exp14_bnb exp14b_bnb_big exp16_physics exp4_2d exp4b_recert exp13_2d_strong exp15_verify_all figs fig2d; do
  echo "=== START $e $(date)" >> ../results/logs/run_all.log
  python3 -W ignore $e.py > ../results/logs/$e.log 2>&1
  echo "=== END $e $(date) exit=$?" >> ../results/logs/run_all.log
done
python3 assemble.py > /dev/null && echo "=== assembled $(date)" >> ../results/logs/run_all.log
