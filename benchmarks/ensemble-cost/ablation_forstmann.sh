#!/usr/bin/env bash
# After queue_c: install the branch head and the two single-change variants,
# run the Forstmann ablation of the plan changes, then the test suites. One
# job at a time on CPUs 23-31.
while kill -0 1602987 2>/dev/null; do sleep 60; done
echo "queue_c done $(date -u)"
mkdir -p /data/work/EMC2_ws5_libs/final_src /data/work/EMC2_ws5_libs/final /data/work/EMC2_ws5_libs/abl_price /data/work/EMC2_ws5_libs/abl_obj /data/work/EMC2_ws5_libs/abl /data/work/EMC2_ws5_libs/tests
git -C /data/work/EMC2_ws5 archive 141df019 | tar -x -C /data/work/EMC2_ws5_libs/final_src
for v in final abl_price abl_obj; do
  taskset -c 23-31 env MAKEFLAGS=-j9 R CMD INSTALL --no-test-load -l /data/work/EMC2_ws5_libs/$v /data/work/EMC2_ws5_libs/${v}_src > /data/work/EMC2_ws5_libs/tests/install_$v.log 2>&1
  echo "install $v exit $? $(date -u)"
done
cd /data/work/EMC2_ws5
run() { # name lib arm seed
  EMC_LIB=/data/work/EMC2_ws5_libs/$2 R_LIBS=/data/work/EMC2_ws5_libs/$2 taskset -c 23-31 Rscript /data/work/EMC2_ws5_libs/h/cont_mr_v2.R /data/work/EMC2_workstreams/fixtures/forstmann_full_fit.rds $3 400 /data/work/EMC2_ws5_libs/abl/forstmann_${1}_s$4.rds $4 > /data/work/EMC2_ws5_libs/abl/forstmann_${1}_s$4.log 2>&1
  echo "$1 seed $4 exit $? $(date -u)"
}
for seed in 1 2 3; do
  run fixed32_10 p2 32_10 $seed
  run rule_v1 p2 rule $seed
  run price_only abl_price rule $seed
  run objective_only abl_obj rule $seed
  run rule_v2 final rule $seed
done
echo "ablation done $(date -u)"
export NOT_CRAN=true EMC_LIB=/data/work/EMC2_ws5_libs/final R_LIBS=/data/work/EMC2_ws5_libs/final
echo "== full level: ensemble suites"
EMC2_TEST_LEVEL=full taskset -c 23-31 Rscript /data/work/EMC2_workstreams/harness/runtests.R test-ensemble-cpp.R test-ensemble-plan.R test-ensemble-move.R test-ensemble-pool-ll.R test-ensemble-latent-output.R 2>&1 | grep "EMC2 from\|^test-\|^\[\|Error"
echo "== default level: sampler suites"
taskset -c 23-31 Rscript /data/work/EMC2_workstreams/harness/runtests.R test-shear-move.R test-shear-standard-invariance.R test-collapsed-move.R test-subject-hmc.R test-move-allocation.R 2>&1 | grep "EMC2 from\|^test-\|^\[\|Error"
echo "== default level: sampling / fitting suites"
taskset -c 23-31 Rscript /data/work/EMC2_workstreams/harness/runtests.R test-sampling-rejection.R test-sample-proposal-schedule.R test-sampler-guards.R test-sample-store-inplace.R test-fit.R test-refit-alpha.R 2>&1 | grep "EMC2 from\|^test-\|^\[\|Error"
echo "queue_g exit $(date -u)"
