#!/bin/bash
cd /data/work/PM/NirvanaHons_Nback
for s in 1 2 3; do
  EMC_LIB=/data/work/EMC2_dev_oo/.Rtmp/claude-1000/-data-work-EMC2-dev-oo/34a015da-344d-4064-b780-dea2402b7c5a/scratchpad/lib_cand Rscript /data/work/EMC2_dev_oo/plans/alloc_prediction_2026-10-15/cont.R cand $s 400 /data/work/EMC2_dev_oo/.Rtmp/claude-1000/-data-work-EMC2-dev-oo/34a015da-344d-4064-b780-dea2402b7c5a/scratchpad/arms/cand_$s.rds > /data/work/EMC2_dev_oo/.Rtmp/claude-1000/-data-work-EMC2-dev-oo/34a015da-344d-4064-b780-dea2402b7c5a/scratchpad/arms/cand_$s.log 2>&1
done
touch /data/work/EMC2_dev_oo/.Rtmp/claude-1000/-data-work-EMC2-dev-oo/34a015da-344d-4064-b780-dea2402b7c5a/scratchpad/arms/cand.done
