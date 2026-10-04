#!/bin/bash
# launch torder.py for the CE variants used in the 2026-10-02 analysis (bash word splitting)
cd "$(dirname "$0")"
mkdir -p mc_logs
while read tag key csel; do
  [ -z "$tag" ] && continue
  nohup python torder.py ce2_${tag}.json ${key} ${csel} --boot 30 > mc_logs/torder_${tag}_${key}_${csel}.log 2>&1 &
done <<LIST
YBMFO E_PM_best_bw_bayes best
YBMFO E_PM_best_bw best
YBMFO E_PM_best ip1+ap+Y
YBMFO E_PM_low best
YBMFO E_PM_clean_bayes best
YBMFO E_PM_all_bayes best
YBCFO E_PM_best_bw_bayes best
YBCFO E_PM_best_bw best
LIST
wait
