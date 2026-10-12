#!/bin/bash
cd "$(dirname "$0")"; L=../results/logs/rerun_killed.log; : > $L
for e in exp7_lasserre exp9_outlier_pw exp14b_bnb_big; do echo "=== START $e $(date)" >> $L; python3 -W ignore $e.py > ../results/logs/$e.log 2>&1; rc=$?; echo "=== END $e $(date) exit=$rc" >> $L; done
echo "=== DONE $(date)" >> $L
