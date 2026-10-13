#!/usr/bin/env bash
# Pool size sweep at three rounds (and the old (32,10)), one fit at a time on
# CPUs 23-31, each fit in its own systemd scope so its cgroup CPU is its own.
cd /data/work/EMC2_ws5
out=/data/work/EMC2_ws5_libs/msweep
run() { # name fixture M R iter seed
  systemd-run --user --scope -q env EMC_LIB=/data/work/EMC2_ws5_libs/cpp2 R_LIBS=/data/work/EMC2_ws5_libs/cpp2 taskset -c 23-31 Rscript /data/work/EMC2_ws5_libs/h/cont_m.R $2 $3 $4 $5 $out/${1}_M$3_R$4_s$6.rds $6 > $out/${1}_M$3_R$4_s$6.log 2>&1
  echo "$1 M$3 R$4 seed $6 exit $? $(date -u)"
}
F=/data/work/EMC2_workstreams/fixtures/forstmann_full_fit.rds
T=/data/work/EMC2_ws5_libs/ckpt/T30_seed1_emc.rds
L=/data/work/EMC2_workstreams/fixtures/lba7r_sample_stage.rds
for seed in 1 2 3; do
  for M in 2 8 19 32 64; do run forstmann $F $M 3 400 $seed; done
  run forstmann $F 32 10 400 $seed
done
for seed in 1 2; do
  for M in 2 8 32 64; do run T30 $T $M 3 400 $seed; done
  run T30 $T 32 10 400 $seed
done
for M in 2 32 106; do run lba7r $L $M 3 100 1; done
echo "queue_m exit $(date -u)"
