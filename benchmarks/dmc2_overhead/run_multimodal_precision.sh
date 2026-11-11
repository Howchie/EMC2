#!/usr/bin/env bash
set -euo pipefail
: "${LIB:?}" "${OUT_DIR:?}"
SEEDS=${SEEDS:-"101 202 303"}
DRAWS=${DRAWS:-5000}
CPU_SET=${CPU_SET:-0-26}
BASE_REUSES=${BASE_REUSES:-0}
NEW_REUSES=${NEW_REUSES:-Inf}
mkdir -p "$OUT_DIR"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
for seed_value in $SEEDS; do
  for arm in base new; do
    if [[ "$arm" == base ]]; then cache_reuses=$BASE_REUSES; else cache_reuses=$NEW_REUSES; fi
    label="multimodal_${seed_value}_${arm}"
    echo "Running $label with $DRAWS sample draws per chain"
    env EMC_LIB="$LIB" R_LIBS="$LIB" SEED="$seed_value" DRAWS="$DRAWS" \
      CACHE_MAX_REUSES="$cache_reuses" \
      OUT="$OUT_DIR/$label.rds" \
      /usr/bin/time -v -o "$OUT_DIR/$label.time" \
      taskset -c "$CPU_SET" Rscript benchmarks/dmc2_overhead/multimodal_fit.R \
      > "$OUT_DIR/$label.log" 2>&1
  done
done
