#!/usr/bin/env bash
# After queue_e: Forstmann ablation of the two plan changes, 3 seeds, one fit
# at a time, 3 chains x 3 cores on CPUs 23-31.
while kill -0 2422681 2>/dev/null; do sleep 60; done
echo "queue_e done $(date -u)"
mkdir -p /data/work/EMC2_ws5_libs/abl_price /data/work/EMC2_ws5_libs/abl_obj
taskset -c 23-31 env MAKEFLAGS=-j9 R CMD INSTALL --no-test-load -l /data/work/EMC2_ws5_libs/abl_price /data/work/EMC2_ws5_libs/abl_price_src > /data/work/EMC2_ws5_libs/tests/install_abl_price.log 2>&1
echo "install abl_price exit $? $(date -u)"
taskset -c 23-31 env MAKEFLAGS=-j9 R CMD INSTALL --no-test-load -l /data/work/EMC2_ws5_libs/abl_obj /data/work/EMC2_ws5_libs/abl_obj_src > /data/work/EMC2_ws5_libs/tests/install_abl_obj.log 2>&1
echo "install abl_obj exit $? $(date -u)"
cd /data/work/EMC2_ws5
fx=/data/work/EMC2_workstreams/fixtures/forstmann_full_fit.rds
out=/data/work/EMC2_ws5_libs/abl
mkdir -p $out
run() { # name lib arm seed
  EMC_LIB=/data/work/EMC2_ws5_libs/$2 R_LIBS=/data/work/EMC2_ws5_libs/$2 taskset -c 23-31 Rscript /data/work/EMC2_ws5_libs/h/cont_mr_v2.R $fx $3 400 $out/forstmann_${1}_s$4.rds $4 > $out/forstmann_${1}_s$4.log 2>&1
  echo "$1 seed $4 exit $? $(date -u)"
}
for seed in 1 2 3; do
  run fixed32_10 p2 32_10 $seed
  run rule_v1 p2 rule $seed
  run price_only abl_price rule $seed
  run objective_only abl_obj rule $seed
  run rule_v2 final rule $seed
done
echo "queue_f exit $(date -u)"
