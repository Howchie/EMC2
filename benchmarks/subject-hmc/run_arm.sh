#!/usr/bin/env bash
# Run one fit_arm.R in its own cgroup scope; append the scope's cpu.stat to
# the log and touch <out>.done when the fit has finished. FEATURE, FEATURE_PAR
# and FIXTURE_DIR pass through to fit_arm.R.
#   run_arm.sh <lib> <case> <arm> <seed> <outdir> [cores]
# cores: a taskset list such as 0-14 (default: unpinned). FIT_FILE=1 autosaves
# the chains to <out>_autosave.RData as the fit runs.
lib=$1; case=$2; arm=$3; seed=$4; dir=$5; cores=${6:-}
mkdir -p "$dir"
base="$dir/${case}_${arm}_${seed}"
rm -f "$base.done"
pin=""
if [ -n "$cores" ]; then pin="taskset -c $cores"; fi
fit_file=""
if [ "${FIT_FILE:-0}" = "1" ]; then fit_file="$base"_autosave.RData; fi
systemd-run --user --scope -q bash -c "cd /data/work/EMC2_dev_oo; EMC_LIB=$lib SAVE_EMC=${SAVE_EMC:-0} PROFILE=${PROFILE:-0} CORES_PER_CHAIN=${CORES_PER_CHAIN:-9} FEATURE=${FEATURE:-} FEATURE_PAR=${FEATURE_PAR:-t0} FIXTURE_DIR=${FIXTURE_DIR:-benchmarks/subject-hmc/fixtures} FIT_FILE=$fit_file timeout ${TIMEOUT:-7200} $pin Rscript benchmarks/subject-hmc/fit_arm.R $case $arm $seed $base.rds > $base.log 2>&1; echo \"exit \$?\" >> $base.log; cat /sys/fs/cgroup\$(cut -d: -f3 /proc/self/cgroup)/cpu.stat >> $base.log; touch $base.done"
