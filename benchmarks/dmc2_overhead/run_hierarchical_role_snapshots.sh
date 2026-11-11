#!/usr/bin/env bash
set -euo pipefail
: "${NEW_LIB:?}" "${OUT_DIR:?}" "${FIXTURE:?}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
mkdir -p "$OUT_DIR/hierarchical_role_snapshots"
env EMC_LIB="$NEW_LIB" R_LIBS="$NEW_LIB" CHAINS=3 CORES_PER_CHAIN=2 ROLE_DIR="$OUT_DIR/hierarchical_role_snapshots" OUT="$OUT_DIR/hierarchical_role_snapshots.rds" \
  Rscript benchmarks/dmc2_overhead/hierarchical_roles_snapshot.R > "$OUT_DIR/hierarchical_role_snapshots.log" 2>&1
