#!/usr/bin/env bash
set -euo pipefail
: "${EMC_LIB:?}" "${FIXTURE:?}" "${OUT_DIR:?}"
mkdir -p "$OUT_DIR"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
CORES_PER_CHAIN=${CORES_PER_CHAIN:-2}
CHAINS=${CHAINS:-3}
ITER=${ITER:-60}
STEP_SIZE=${STEP_SIZE:-20}
RESTART_BLOCKS=${RESTART_BLOCKS:-50}
CPU_SET=${CPU_SET:-0-26}
for seed in ${SEEDS:-101 202 303}; do
  for arm in candidate baseline; do
    persistent=1
    if [[ "$arm" == baseline ]]; then persistent=0; fi
    label="worker_retention_cpc${CORES_PER_CHAIN}_${arm}_${seed}"
    trace_file=""
    if [[ "$persistent" == 1 ]]; then
      trace_file="$OUT_DIR/$label.workers"
    fi
    env EMC_LIB="$EMC_LIB" R_LIBS="$EMC_LIB" FIXTURE="$FIXTURE" \
      PERSISTENT="$persistent" OUT="$OUT_DIR/$label.rds" SEED="$seed" \
      CORES_PER_CHAIN="$CORES_PER_CHAIN" CHAINS="$CHAINS" ITER="$ITER" \
      STEP_SIZE="$STEP_SIZE" RESTART_BLOCKS="$RESTART_BLOCKS" \
      WORKER_TRACE_FILE="$trace_file" \
      MEMORY_FILE="$OUT_DIR/$label.memory.csv" \
      /usr/bin/time -v -o "$OUT_DIR/$label.time" \
      taskset -c "$CPU_SET" Rscript benchmarks/dmc2_overhead/chain_persistence.R \
      > "$OUT_DIR/$label.log" 2>&1
  done
done
