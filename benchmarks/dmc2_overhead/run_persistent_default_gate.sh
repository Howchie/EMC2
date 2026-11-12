#!/usr/bin/env bash
set -euo pipefail
: "${EMC_LIB:?}" "${OUT_DIR:?}"
mkdir -p "$OUT_DIR"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
CPU_SET=${CPU_SET:-0-26}
RESTART_BLOCKS=${RESTART_BLOCKS:-25}

for case_name in ${CASES:-EXG LNR}; do
  seed_index=0
  for seed in ${SEEDS:-101 202 303}; do
    seed_index=$((seed_index + 1))
    if (( seed_index % 2 )); then
      arms=(0 1)
    else
      arms=(1 0)
    fi
    for persistent in "${arms[@]}"; do
      mode=fork_per_block
      if [[ "$persistent" == 1 ]]; then mode=persistent; fi
      label="${case_name}_${seed}_${mode}_default"
      echo "Running $label"
      env EMC_LIB="$EMC_LIB" R_LIBS="$EMC_LIB" CASE="$case_name" \
        SEED="$seed" HOST=1 PROFILE=0 PERSISTENT="$persistent" \
        RESTART_BLOCKS="$RESTART_BLOCKS" OUT="$OUT_DIR/$label.rds" \
        /usr/bin/time -v -o "$OUT_DIR/$label.time" \
        taskset -c "$CPU_SET" Rscript benchmarks/dmc2_overhead/paper_default.R \
        > "$OUT_DIR/$label.log" 2>&1
    done
  done
done
