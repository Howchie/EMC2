#!/usr/bin/env bash
set -euo pipefail
: "${BASE_LIB:?}" "${NEW_LIB:?}" "${INPUT_DIR:?}" "${OUT_DIR:?}"
MAX_SAMPLE_DRAWS=${MAX_SAMPLE_DRAWS:-4000}
FORCE=${FORCE:-0}
CPUSET=${CPUSET:-0-26}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
read -r -a seed_values <<< "${POOLED_SEEDS:-1506 1607 1708 1809}"
mkdir -p "$OUT_DIR"
for index in "${!seed_values[@]}"; do
  seed_value=${seed_values[$index]}
  if (( index % 2 == 0 )); then arms=(base new); else arms=(new base); fi
  for arm in "${arms[@]}"; do
    if [[ "$arm" == base ]]; then lib=$BASE_LIB; else lib=$NEW_LIB; fi
    label="pooled_${seed_value}_${arm}"
    input="$INPUT_DIR/$label.rds"
    output="$OUT_DIR/${label}_stop.rds"
    if [[ -f "$output" && "$FORCE" != 1 ]]; then
      echo "Skipping completed $label"
      continue
    fi
    echo "Checking $label to the DMC2 hardfit R-hat/ESS subset target"
    taskset -c "$CPUSET" env EMC_LIB="$lib" R_LIBS="$lib" \
      Rscript benchmarks/dmc2_overhead/continue_pooled_scale_gate.R \
      "$input" "$output" "$MAX_SAMPLE_DRAWS" \
      > "$OUT_DIR/${label}_stop.log" 2>&1
  done
done
