#!/usr/bin/env bash
set -euo pipefail
: "${BASE_LIB:?}" "${CANDIDATE_LIB:?}" "${OUT_DIR:?}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
for arm in base candidate; do
  if [[ "$arm" == base ]]; then lib=$BASE_LIB; else lib=$CANDIDATE_LIB; fi
  echo "Running configurable-bound stop-success cache $arm"
  env EMC_LIB="$lib" R_LIBS="$lib" OUT="$OUT_DIR/ss_cache_$arm.rds" \
    Rscript benchmarks/dmc2_overhead/stop_success_cache_fixture.R \
    > "$OUT_DIR/ss_cache_$arm.log" 2>&1
  echo "Running full SSEXG stop-success likelihood $arm"
  env EMC_LIB="$lib" R_LIBS="$lib" OUT="$OUT_DIR/ss_likelihood_$arm.rds" \
    Rscript benchmarks/dmc2_overhead/stop_success_likelihood_fixture.R \
    > "$OUT_DIR/ss_likelihood_$arm.log" 2>&1
done
