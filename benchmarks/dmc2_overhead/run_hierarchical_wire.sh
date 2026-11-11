#!/usr/bin/env bash
set -euo pipefail
: "${NEW_LIB:?}" "${OUT_DIR:?}" "${FIXTURE:?}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
label=hierarchical_wire
EMC_LIB="$NEW_LIB" R_LIBS="$NEW_LIB" CHAINS=3 CORES_PER_CHAIN=2 OUT="$OUT_DIR/$label.rds" \
  Rscript benchmarks/dmc2_overhead/hierarchical_wire.R > "$OUT_DIR/$label.log" 2>&1
