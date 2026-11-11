#!/usr/bin/env bash
set -euo pipefail
: "${BASE_LIB:?}" "${NEW_LIB:?}" "${OUT_DIR:?}" "${FIXTURE:?}"
mkdir -p "$OUT_DIR"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
for arm in base new; do
  if [[ "$arm" == base ]]; then lib=$BASE_LIB; else lib=$NEW_LIB; fi
  label="hierarchical_rebuild_${arm}"
  echo "Running $label"
  env EMC_LIB="$lib" R_LIBS="$lib" CHAINS=2 CORES_PER_CHAIN=2 OUT="$OUT_DIR/$label.rds" \
    /usr/bin/time -v -o "$OUT_DIR/$label.time" Rscript benchmarks/dmc2_overhead/hierarchical_rebuild.R > "$OUT_DIR/$label.log" 2>&1
done
