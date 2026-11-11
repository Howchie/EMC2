#!/usr/bin/env bash
set -euo pipefail
: "${BASE_LIB:?}" "${NEW_LIB:?}" "${OUT_DIR:?}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
for case_name in EXG LNR; do
  for rep in 2 3; do
    for arm in new base; do
      if [[ "$arm" == base ]]; then lib=$BASE_LIB; else lib=$NEW_LIB; fi
      label="${case_name}_101_${arm}_rep${rep}"
      echo "Running $label"
      env EMC_LIB="$lib" R_LIBS="$lib" CASE="$case_name" SEED=101 HOST=1 PROFILE=0 OUT="$OUT_DIR/$label.rds" \
        /usr/bin/time -v -o "$OUT_DIR/$label.time" Rscript benchmarks/dmc2_overhead/paper_fixture.R > "$OUT_DIR/$label.log" 2>&1
    done
  done
done
