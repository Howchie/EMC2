#!/usr/bin/env bash
set -euo pipefail
: "${BASE_LIB:?}" "${DMC_LIB:?}" "${OUT_DIR:?}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
for case_name in EXG LNR; do
  data="$OUT_DIR/${case_name}_shared_data.rds"
  env EMC_LIB="$BASE_LIB" R_LIBS="$BASE_LIB" CASE="$case_name" OUT="$data" Rscript benchmarks/dmc2_overhead/export_paper_data.R > "$OUT_DIR/${case_name}_export.log" 2>&1
  for arm in base archive; do
    if [[ "$arm" == base ]]; then lib=$BASE_LIB; compat=0; else lib=$DMC_LIB; compat=1; fi
    for seed_value in 101 202 303; do
      label="${case_name}_${seed_value}_${arm}_shared"
      if [[ -f "$OUT_DIR/$label.rds" ]]; then continue; fi
      echo "Running $label"
      env EMC_LIB="$lib" R_LIBS="$lib" CASE="$case_name" SEED="$seed_value" DATA_IN="$data" ARCHIVE_COMPAT="$compat" HOST=1 PROFILE=0 OUT="$OUT_DIR/$label.rds" \
        /usr/bin/time -v -o "$OUT_DIR/$label.time" Rscript benchmarks/dmc2_overhead/paper_shared_legacy.R > "$OUT_DIR/$label.log" 2>&1
    done
  done
done
