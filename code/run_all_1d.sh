#!/bin/bash
cd /home/user/Fable-Research/code
for e in exp1_1d_features exp8_robustness exp5_smo_euv_toeplitz exp3_pixel_binary exp7_lasserre exp2_random_tightness exp9_outlier_pw figs; do
  echo "=== START $e $(date)" >> ../results/logs/run_all_1d.log
  python3 -W ignore $e.py > ../results/logs/$e.log 2>&1
  echo "=== END $e $(date) exit=$?" >> ../results/logs/run_all_1d.log
done
