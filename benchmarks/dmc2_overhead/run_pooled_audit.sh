#!/usr/bin/env bash
set -euo pipefail
: "${BASE_LIB:?}" "${NEW_LIB:?}" "${OUT_DIR:?}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
for mode in base_repeat base_fixed new_fixed; do
  case "$mode" in
    base_repeat) lib=$BASE_LIB; script=pooled_fixture.R;;
    base_fixed) lib=$BASE_LIB; script=pooled_fixed.R;;
    new_fixed) lib=$NEW_LIB; script=pooled_fixed.R;;
  esac
  label="pooled_audit_$mode"
  echo "Running $label"
  env EMC_LIB="$lib" R_LIBS="$lib" OUT="$OUT_DIR/$label.rds" \
    /usr/bin/time -v -o "$OUT_DIR/$label.time" Rscript "benchmarks/dmc2_overhead/$script" > "$OUT_DIR/$label.log" 2>&1
done
