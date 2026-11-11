#!/usr/bin/env bash
set -euo pipefail
: "${BASE_LIB:?}" "${CANDIDATE_LIB:?}" "${OUT_DIR:?}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
for arm in base candidate; do
  if [[ "$arm" == base ]]; then lib=$BASE_LIB; else lib=$CANDIDATE_LIB; fi
  echo "Running four-cell REXG normalization cache $arm"
  env EMC_LIB="$lib" R_LIBS="$lib" OUT="$OUT_DIR/rexg_cells_$arm.rds" \
    Rscript benchmarks/dmc2_overhead/kernel_fixture_cells.R > "$OUT_DIR/rexg_cells_$arm.log" 2>&1
done
