#!/usr/bin/env bash
set -euo pipefail
: "${DMC_LIB:?}" "${DMC_SRC:?}" "${OUT_DIR:?}" "${FORTRAN_LIB:?}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 LIBRARY_PATH="$FORTRAN_LIB${LIBRARY_PATH:+:$LIBRARY_PATH}"
R CMD INSTALL -l "$DMC_LIB" "$DMC_SRC" > "$OUT_DIR/dmc-archive-install.log" 2>&1
for script in kernel_archive.R kernel_lnr_archive.R; do
  echo "Running archive $script"
  env EMC_LIB="$DMC_LIB" R_LIBS="$DMC_LIB" PKG=EMC2 Rscript "benchmarks/dmc2_overhead/$script" > "$OUT_DIR/archive_${script%.R}.log" 2>&1
done
for case_name in EXG LNR; do
  for seed_value in 101 202 303; do
    label="${case_name}_${seed_value}_archive_default"
    echo "Running $label"
    env EMC_LIB="$DMC_LIB" R_LIBS="$DMC_LIB" CASE="$case_name" SEED="$seed_value" HOST=1 PROFILE=0 OUT="$OUT_DIR/$label.rds" \
      /usr/bin/time -v -o "$OUT_DIR/$label.time" Rscript benchmarks/dmc2_overhead/paper_default.R > "$OUT_DIR/$label.log" 2>&1
  done
done
