#!/usr/bin/env bash
set -euo pipefail
: "${BASE_LIB:?}" "${NEW_LIB:?}" "${OUT_DIR:?}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
for arm in base new; do
  if [[ "$arm" == base ]]; then lib=$BASE_LIB; else lib=$NEW_LIB; fi
  label="EXG_${arm}_memory"
  echo "Running $label"
  env EMC_LIB="$lib" R_LIBS="$lib" CASE=EXG SEED=101 HOST=1 PROFILE=0 OUT="$OUT_DIR/$label.rds" \
    python3 benchmarks/dmc2_overhead/monitor_tree.py "$OUT_DIR/$label.csv" -- \
    Rscript benchmarks/dmc2_overhead/paper_fixture.R > "$OUT_DIR/$label.log" 2>&1
done
