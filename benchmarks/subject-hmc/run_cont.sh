#!/usr/bin/env bash
# Run one continue_arm.R in its own cgroup scope, pinned to a core set, and
# append the scope's cpu.stat to the log; touch <out>.done when finished.
#   run_cont.sh <lib> <ckpt.rds> <label> <seed> <iters> <outdir> [cores]
# cores: a taskset list such as 0-14 (default: unpinned). CORES_PER_CHAIN,
# DISCARD, CALIBRATE, N_PARTICLES and ARM_R pass through to continue_arm.R.
lib=$1; ckpt=$2; label=$3; seed=$4; iters=$5; dir=$6; cores=${7:-}
mkdir -p "$dir"
base="$dir/${label}_${seed}"
rm -f "$base.done"
pin=""
if [ -n "$cores" ]; then pin="taskset -c $cores"; fi
# A watchdog logs the scope's memory every 15 s to <out>.mem and kills the run
# if it exceeds MEM_MAX_GB (default 40), so a leak cannot take the machine down.
systemd-run --user --scope -q bash -c "cd /data/work/EMC2_dev_oo; cg=/sys/fs/cgroup\$(cut -d: -f3 /proc/self/cgroup); ( while sleep 15; do m=\$(cat \$cg/memory.current 2>/dev/null || echo 0); echo \"\$(date +%T) \$m\" >> $base.mem; if [ \$m -gt \$(( ${MEM_MAX_GB:-40} * 1024 * 1024 * 1024 )) ]; then echo \"memory cap hit: \$m\" >> $base.log; pkill -9 -f 'continue_arm.R $ckpt $label $seed'; fi; done ) & wd=\$!; EMC_LIB=$lib CORES_PER_CHAIN=${CORES_PER_CHAIN:-9} DISCARD=${DISCARD:-} CALIBRATE=${CALIBRATE:-} N_PARTICLES=${N_PARTICLES:-} ARM_R=${ARM_R:-} timeout ${TIMEOUT:-36000} $pin Rscript benchmarks/subject-hmc/continue_arm.R $ckpt $label $seed $iters $base.rds > $base.log 2>&1; echo \"exit \$?\" >> $base.log; cat \$cg/cpu.stat >> $base.log; kill \$wd 2>/dev/null; touch $base.done"
