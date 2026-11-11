#!/usr/bin/env bash
# Run under bench_lock.sh with the candidate install in EMC_LIB.
set -euo pipefail
: "${EMC_LIB:?}" "${OUT_DIR:?}"
CPU_SET=${CPU_SET:-0-26}
SEED=${SEED:-101}
mkdir -p "$OUT_DIR"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
for case_name in EXG LNR; do
  for arm in off on; do
    if [[ "$arm" == off ]]; then context=0; else context=1; fi
    label="${case_name}_${SEED}_${arm}"
    echo "Running $label (native parameter context=$context)"
    env EMC_LIB="$EMC_LIB" R_LIBS="$EMC_LIB" CASE="$case_name" \
      SEED="$SEED" HOST=0 PROFILE=0 NATIVE_CONTEXT="$context" \
      OUT="$OUT_DIR/$label.rds" \
      /usr/bin/time -v -o "$OUT_DIR/$label.time" \
      taskset -c "$CPU_SET" Rscript benchmarks/dmc2_overhead/paper_profile.R \
      > "$OUT_DIR/$label.log" 2>&1
  done
done
