#!/usr/bin/env bash
set -euo pipefail
: "${BASE_LIB:?}" "${CANDIDATE_LIB:?}" "${OUT_DIR:?}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
for arm in base new; do
  if [[ "$arm" == base ]]; then lib=$BASE_LIB; else lib=$CANDIDATE_LIB; fi
  echo "Running EXG kernel fit parity $arm"
  env EMC_LIB="$lib" R_LIBS="$lib" CASE=EXG SEED=101 HOST=1 PROFILE=0 \
    OUT="$OUT_DIR/rexg_kernel_fit_${arm}.rds" \
    /usr/bin/time -v -o "$OUT_DIR/rexg_kernel_fit_${arm}.time" \
    Rscript benchmarks/dmc2_overhead/paper_fixture.R \
    > "$OUT_DIR/rexg_kernel_fit_${arm}.log" 2>&1
done
