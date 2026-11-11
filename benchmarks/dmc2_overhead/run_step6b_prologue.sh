#!/usr/bin/env bash
set -euo pipefail
: "${EMC_LIB:?}" "${OUT_DIR:?}"
CPU_SET=${CPU_SET:-0-26}
NP_GRID=${NP_GRID:-1,8,32,128,256}
BATCHES=${BATCHES:-5}
FULL_REPS=${FULL_REPS:-5}
PROLOGUE_REPS=${PROLOGUE_REPS:-100}
mkdir -p "$OUT_DIR"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
for case_name in EXG LNR; do
  echo "Running likelihood prologue benchmark: $case_name"
  env EMC_LIB="$EMC_LIB" R_LIBS="$EMC_LIB" CASE="$case_name" NP_GRID="$NP_GRID" \
    BATCHES="$BATCHES" FULL_REPS="$FULL_REPS" \
    PROLOGUE_REPS="$PROLOGUE_REPS" OUT="$OUT_DIR/$case_name.rds" \
    /usr/bin/time -v -o "$OUT_DIR/$case_name.time" \
    taskset -c "$CPU_SET" Rscript benchmarks/dmc2_overhead/likelihood_prologue.R \
    > "$OUT_DIR/$case_name.log" 2>&1
done
