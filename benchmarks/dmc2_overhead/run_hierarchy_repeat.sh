#!/usr/bin/env bash
set -euo pipefail
: "${BASE_LIB:?}" "${OUT_DIR:?}" "${FIXTURE:?}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
for rep in 2 3; do
  label="hierarchical_2_base_rep${rep}"
  echo "Running $label"
  env EMC_LIB="$BASE_LIB" R_LIBS="$BASE_LIB" CHAINS=2 CORES_PER_CHAIN=2 OUT="$OUT_DIR/$label.rds" \
    /usr/bin/time -v -o "$OUT_DIR/$label.time" Rscript benchmarks/dmc2_overhead/hierarchical_fixture.R > "$OUT_DIR/$label.log" 2>&1
done
