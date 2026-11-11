#!/usr/bin/env bash
# Run under bench_lock.sh with one isolated build containing both code paths.
set -euo pipefail
: "${EMC_LIB:?}" "${OUT_DIR:?}"
CPU_SET=${CPU_SET:-0-26}
EXG_SEEDS=${EXG_SEEDS:-"101 202 303"}
LNR_SEEDS=${LNR_SEEDS:-"101"}
mkdir -p "$OUT_DIR"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
for case_name in EXG LNR; do
  if [[ "$case_name" == EXG ]]; then seeds=$EXG_SEEDS; else seeds=$LNR_SEEDS; fi
  for seed_value in $seeds; do
    for arm in off on; do
      if [[ "$arm" == off ]]; then cpp_weights=0; else cpp_weights=1; fi
      label="${case_name}_${seed_value}_${arm}"
      echo "Running $label (C++ particle weights=$cpp_weights)"
      env EMC_LIB="$EMC_LIB" R_LIBS="$EMC_LIB" CASE="$case_name" \
        SEED="$seed_value" HOST=0 PROFILE=0 NATIVE_CONTEXT=0 \
        PARTICLE_WEIGHTS_CPP="$cpp_weights" OUT="$OUT_DIR/$label.rds" \
        /usr/bin/time -v -o "$OUT_DIR/$label.time" \
        taskset -c "$CPU_SET" Rscript benchmarks/dmc2_overhead/paper_profile.R \
        > "$OUT_DIR/$label.log" 2>&1
    done
  done
done
