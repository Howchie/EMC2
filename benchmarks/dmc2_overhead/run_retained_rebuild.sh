#!/usr/bin/env bash
set -euo pipefail
: "${BASE_LIB:?}" "${CANDIDATE_LIB:?}" "${OUT_DIR:?}" "${FIXTURE:?}"
mkdir -p "$OUT_DIR"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
read -r -a seeds <<< "${SEEDS:-101 202 303}"

for seed in "${seeds[@]}"; do
  if (( seed == 202 )); then arms=(candidate baseline); else arms=(baseline candidate); fi
  for arm in "${arms[@]}"; do
    if [[ "$arm" == baseline ]]; then lib=$BASE_LIB; else lib=$CANDIDATE_LIB; fi
    label="retained_rebuild_${seed}_${arm}"
    echo "Running $label"
    env EMC_LIB="$lib" R_LIBS="$lib" SEED="$seed" \
      OUT="$OUT_DIR/$label.rds" \
      REBUILD_TRACE="$OUT_DIR/$label.trace" \
      /usr/bin/time -v -o "$OUT_DIR/$label.time" \
      Rscript benchmarks/dmc2_overhead/retained_rebuild.R \
      > "$OUT_DIR/$label.log" 2>&1
  done
  Rscript benchmarks/dmc2_overhead/analyse_retained_rebuild.R "$OUT_DIR" "$seed"
done
