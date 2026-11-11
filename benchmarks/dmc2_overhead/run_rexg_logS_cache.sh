#!/usr/bin/env bash
set -euo pipefail
: "${BASE_LIB:?}" "${CANDIDATE_LIB:?}" "${OUT_DIR:?}"
mkdir -p "$OUT_DIR"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
for arm in base candidate; do
  if [[ "$arm" == base ]]; then lib=$BASE_LIB; else lib=$CANDIDATE_LIB; fi
  env EMC_LIB="$lib" R_LIBS="$lib" \
    OUT="$OUT_DIR/rexg_logS_${arm}.rds" \
    /usr/bin/time -v -o "$OUT_DIR/rexg_logS_${arm}.time" \
    Rscript benchmarks/dmc2_overhead/rexg_logS_cache_fixture.R \
    > "$OUT_DIR/rexg_logS_${arm}.log" 2>&1
done
