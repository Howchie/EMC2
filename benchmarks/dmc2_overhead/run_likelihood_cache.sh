#!/usr/bin/env bash
# Paired default-host fits for the Step 2 subject-fit cache.
set -euo pipefail
: "${BASE_LIB:?}" "${NEW_LIB:?}" "${OUT_DIR:?}"
SEEDS=${SEEDS:-"101 202 303"}
CASES=${CASES:-"EXG LNR"}
ARMS=${ARMS:-"base new"}
CACHE_MAX_REUSES=${CACHE_MAX_REUSES:-default}
RESULT_TAG=default
if [[ "$CACHE_MAX_REUSES" != "default" ]]; then RESULT_TAG="reuse${CACHE_MAX_REUSES}"; fi
mkdir -p "$OUT_DIR"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
for case_name in $CASES; do
  for seed_value in $SEEDS; do
    for arm in $ARMS; do
      if [[ "$arm" == base ]]; then lib=$BASE_LIB; else lib=$NEW_LIB; fi
      label="${case_name}_${seed_value}_${arm}_${RESULT_TAG}"
      echo "Running $label"
      env EMC_LIB="$lib" R_LIBS="$lib" CASE="$case_name" SEED="$seed_value" \
        HOST=1 PROFILE=0 CACHE_MAX_REUSES="$CACHE_MAX_REUSES" \
        OUT="$OUT_DIR/$label.rds" \
        /usr/bin/time -v -o "$OUT_DIR/$label.time" \
        Rscript benchmarks/dmc2_overhead/paper_default.R > "$OUT_DIR/$label.log" 2>&1
    done
  done
done
