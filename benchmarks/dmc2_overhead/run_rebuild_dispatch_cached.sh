#!/usr/bin/env bash
set -euo pipefail
: "${EMC_LIB:?}" "${OUT_DIR:?}"
FIXTURE_DIR=${FIXTURE_DIR:-benchmarks/dmc2_overhead/results}
CPU_SET=${CPU_SET:-0-26}
mkdir -p "$OUT_DIR"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
for seed in ${SEEDS:-101 202 303}; do
  label="EXG_${seed}"
  env EMC_LIB="$EMC_LIB" \
    FIXTURE="$FIXTURE_DIR/${label}_base_default.rds" \
    OUT="$OUT_DIR/${label}.rds" SEED="$seed" \
    /usr/bin/time -v -o "$OUT_DIR/${label}.time" \
    taskset -c "$CPU_SET" Rscript \
      benchmarks/dmc2_overhead/rebuild_dispatch_cached.R \
      > "$OUT_DIR/${label}.log" 2>&1
done
