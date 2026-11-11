#!/usr/bin/env bash
set -euo pipefail
: "${BASE_LIB:?}" "${NEW_LIB:?}" "${OUT_DIR:?}"
mkdir -p "$OUT_DIR"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
BASE_CACHE_MAX_REUSES=${BASE_CACHE_MAX_REUSES:-default}
NEW_CACHE_MAX_REUSES=${NEW_CACHE_MAX_REUSES:-default}
SEEDS=${SEEDS:-"101 202 303"}
CASES=${CASES:-"EXG LNR"}
for case_name in $CASES; do
  for seed_value in $SEEDS; do
    for arm in base new; do
      if [[ "$arm" == base ]]; then
        lib=$BASE_LIB
        cache_max_reuses=$BASE_CACHE_MAX_REUSES
      else
        lib=$NEW_LIB
        cache_max_reuses=$NEW_CACHE_MAX_REUSES
      fi
      label="${case_name}_${seed_value}_${arm}_precision"
      echo "Running $label"
      env EMC_LIB="$lib" R_LIBS="$lib" CASE="$case_name" SEED="$seed_value" HOST=1 PROFILE=0 OUT="$OUT_DIR/$label.rds" \
        CACHE_MAX_REUSES="$cache_max_reuses" \
        /usr/bin/time -v -o "$OUT_DIR/$label.time" Rscript benchmarks/dmc2_overhead/paper_precision.R > "$OUT_DIR/$label.log" 2>&1
    done
  done
done
