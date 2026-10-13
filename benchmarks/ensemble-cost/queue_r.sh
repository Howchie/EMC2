#!/usr/bin/env bash
# Rounds at M = 32 on Forstmann, filling each arm to 4 seeds; one fit at a time.
cd /data/work/EMC2_ws5
out=/data/work/EMC2_ws5_libs/msweep
F=/data/work/EMC2_workstreams/fixtures/forstmann_full_fit.rds
run() { # R seed
  systemd-run --user --scope -q env EMC_LIB=/data/work/EMC2_ws5_libs/cpp2 R_LIBS=/data/work/EMC2_ws5_libs/cpp2 taskset -c 23-31 Rscript /data/work/EMC2_ws5_libs/h/cont_m2.R $F 32 $1 400 $out/forstmann_M32_R$1_s$2.rds $2 > $out/forstmann_M32_R$1_s$2.log 2>&1
  echo "R$1 seed $2 exit $? $(date -u)"
}
for seed in 1 2 3 4; do run 5 $seed; run 40 $seed; done
for seed in 3 4; do run 10 $seed; run 20 $seed; done
run 3 4
echo "queue_r exit $(date -u)"
