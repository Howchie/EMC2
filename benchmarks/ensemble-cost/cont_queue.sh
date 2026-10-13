#!/usr/bin/env bash
# Continuation arms, one fit at a time, 3 chains x 3 cores pinned to 9 CPUs.
cd /data/work/EMC2_ws5
export EMC_LIB=/data/work/EMC2_ws5_libs/p2 R_LIBS=/data/work/EMC2_ws5_libs/p2
out=/data/work/EMC2_ws5_libs/cont; mkdir -p $out
run() { # fixture_name path iter arm
  taskset -c 23-31 Rscript /data/work/EMC2_ws5_libs/h/cont_mr.R $2 $4 $3 $out/${1}_${4}.rds > $out/${1}_${4}.log 2>&1; echo "$1 $4 exit $?"
}
for arm in 32_10 rule 2_2 8_3; do run forstmann /data/work/EMC2_workstreams/fixtures/forstmann_full_fit.rds 400 $arm; done
for arm in 32_10 rule 2_2 8_3; do run T30 /data/work/EMC2_ws5_libs/ckpt/T30_seed1_emc.rds 400 $arm; done
for arm in 32_10 rule 2_2; do run lba7r /data/work/EMC2_workstreams/fixtures/lba7r_sample_stage.rds 200 $arm; done
