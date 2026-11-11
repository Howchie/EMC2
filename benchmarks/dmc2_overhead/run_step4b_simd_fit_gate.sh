#!/usr/bin/env bash
set -euo pipefail
: "${BASE_LIB:?}" "${CANDIDATE_LIB:?}" "${OUT_DIR:?}"
CPU_SET=${CPU_SET:-0-26}
SEEDS=${SEEDS:-"101 202 303"}
mkdir -p "$OUT_DIR"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
for seed_value in $SEEDS; do
  for arm in base new; do
    if [[ "$arm" == base ]]; then lib=$BASE_LIB; else lib=$CANDIDATE_LIB; fi
    label="EXG_${seed_value}_${arm}"
    echo "Running $label"
    env EMC_LIB="$lib" R_LIBS="$lib" CASE=EXG SEED="$seed_value" \
      HOST=1 PROFILE=0 NATIVE_CONTEXT=0 PARTICLE_WEIGHTS_CPP=0 \
      OUT="$OUT_DIR/$label.rds" \
      /usr/bin/time -v -o "$OUT_DIR/$label.time" \
      taskset -c "$CPU_SET" Rscript benchmarks/dmc2_overhead/paper_profile.R \
      > "$OUT_DIR/$label.log" 2>&1
  done
done
