#!/usr/bin/env bash
set -euo pipefail
: "${BASE_LIB:?}" "${CANDIDATE_LIB:?}" "${OUT_DIR:?}"
mkdir -p "$OUT_DIR"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
for seed in ${SEEDS:-101 202 303}; do
  for arm in baseline candidate; do
    if [[ "$arm" == baseline ]]; then lib=$BASE_LIB; else lib=$CANDIDATE_LIB; fi
    label="rexg_logS_${arm}_${seed}"
    env EMC_LIB="$lib" R_LIBS="$lib" SEED="$seed" \
      OUT="$OUT_DIR/$label.rds" \
      /usr/bin/time -v -o "$OUT_DIR/$label.time" \
      Rscript benchmarks/dmc2_overhead/rexg_logS_fit.R \
      > "$OUT_DIR/$label.log" 2>&1
  done
done
