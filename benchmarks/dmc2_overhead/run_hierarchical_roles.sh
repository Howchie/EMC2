#!/usr/bin/env bash
set -euo pipefail
: "${NEW_LIB:?}" "${OUT_DIR:?}" "${FIXTURE:?}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
mkdir -p "$OUT_DIR/hierarchical_roles"
env EMC_LIB="$NEW_LIB" R_LIBS="$NEW_LIB" CHAINS=3 CORES_PER_CHAIN=2 ROLE_DIR="$OUT_DIR/hierarchical_roles" OUT="$OUT_DIR/hierarchical_roles.rds" \
  python3 benchmarks/dmc2_overhead/monitor_tree.py "$OUT_DIR/hierarchical_roles.csv" -- \
  Rscript benchmarks/dmc2_overhead/hierarchical_roles.R > "$OUT_DIR/hierarchical_roles.log" 2>&1
