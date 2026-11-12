#!/usr/bin/env bash
set -euo pipefail
: "${BASE_LIB:?}" "${INPUT_DIR:?}" "${OUT_DIR:?}"
MAX_SAMPLE_DRAWS=${MAX_SAMPLE_DRAWS:-4000}
FORCE=${FORCE:-0}
CPUSET=${CPUSET:-0-26}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
read -r -a seed_values <<< "${POOLED_SEEDS:-1506 1607}"
mkdir -p "$OUT_DIR"
for seed_value in "${seed_values[@]}"; do
  for repeat_arm in base new; do
    label="pooled_${seed_value}_${repeat_arm}"
    input="$INPUT_DIR/$label.rds"
    output="$OUT_DIR/${label}_stop.rds"
    if [[ -f "$output" && "$FORCE" != 1 ]]; then
      echo "Skipping completed same-build repeat $label"
      continue
    fi
    echo "Checking same-build repeat $label against the shared scale-aware target patterned after DMC2 hardfit"
    taskset -c "$CPUSET" env EMC_LIB="$BASE_LIB" R_LIBS="$BASE_LIB" \
      Rscript benchmarks/dmc2_overhead/continue_pooled_scale_gate.R \
      "$input" "$output" "$MAX_SAMPLE_DRAWS" \
      > "$OUT_DIR/${label}_stop.log" 2>&1
  done
done
