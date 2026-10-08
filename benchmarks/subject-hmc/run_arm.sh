#!/usr/bin/env bash
# Run one fit_arm.R in its own cgroup scope; append the scope's cpu.stat to
# the log and touch <out>.done when the fit has finished.
#   run_arm.sh <lib> <case> <arm> <seed> <outdir>
lib=$1; case=$2; arm=$3; seed=$4; dir=$5
mkdir -p "$dir"
base="$dir/${case}_${arm}_${seed}"
rm -f "$base.done"
systemd-run --user --scope -q bash -c "cd /data/work/EMC2_dev_oo; EMC_LIB=$lib SAVE_EMC=${SAVE_EMC:-0} PROFILE=${PROFILE:-0} CORES_PER_CHAIN=${CORES_PER_CHAIN:-5} timeout ${TIMEOUT:-7200} Rscript benchmarks/subject-hmc/fit_arm.R $case $arm $seed $base.rds > $base.log 2>&1; echo \"exit \$?\" >> $base.log; cat /sys/fs/cgroup\$(cut -d: -f3 /proc/self/cgroup)/cpu.stat >> $base.log; touch $base.done"
