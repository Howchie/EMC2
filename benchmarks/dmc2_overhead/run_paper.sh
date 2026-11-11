#!/usr/bin/env bash
# Run under bench_lock.sh. Libraries must be frozen installs of the two revisions.
set -euo pipefail
: "${BASE_LIB:?}" "${NEW_LIB:?}" "${OUT_DIR:?}"
mkdir -p "$OUT_DIR"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
for case_name in EXG LNR; do
  for seed_value in 101 202 303; do
    for arm in base new; do
      if [[ "$arm" == base ]]; then lib=$BASE_LIB; else lib=$NEW_LIB; fi
      label="${case_name}_${seed_value}_${arm}"
      echo "Running $label"
      env EMC_LIB="$lib" R_LIBS="$lib" CASE="$case_name" SEED="$seed_value" HOST=1 PROFILE=0 OUT="$OUT_DIR/$label.rds" \
        /usr/bin/time -v -o "$OUT_DIR/$label.time" Rscript benchmarks/dmc2_overhead/paper_fixture.R > "$OUT_DIR/$label.log" 2>&1
    done
  done
done
