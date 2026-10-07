#!/usr/bin/env bash
# Warm-ups, then the four production arms per (case, seed); JOBS concurrent
# fits of 3 chains x CORES_PER_CHAIN cores. An arm job waits for its warm-up.
#   EMC_LIB=<lib> JOBS=5 bash benchmarks/production-adaptation/run_all.sh
set -u
cd "$(dirname "$0")/../.."
OUT=${PADAPT_OUT:-benchmarks/production-adaptation/results}
mkdir -p "$OUT/logs"
export PADAPT_OUT=$OUT
JOBS=${JOBS:-5}
CASES=${CASES:-"lba71sim:1 lba71sim:2 jointsim:1 jointsim:2 wide15:1 wide15:2 wide15:3 freeA15:1 freeA15:2 freeA15:3 corr15:1 corr15:2 corr15:3"}
ARMS=${ARMS:-"doubling frozen early extended"}
job() {
  kind=$1; case=$2; seed=$3; arm=${4:-}
  if [ "$kind" = warmup ]; then
    [ -f "$OUT/warmup_${case}_${seed}.rds" ] && exit 0
    CASE=$case SEED=$seed Rscript benchmarks/production-adaptation/warmup.R \
      > "$OUT/logs/warmup_${case}_${seed}.log" 2>&1
  else
    [ -f "$OUT/prod_${case}_${arm}_${seed}.rds" ] && exit 0
    while [ ! -f "$OUT/warmup_${case}_${seed}.rds" ]; do
      grep -q "Error" "$OUT/logs/warmup_${case}_${seed}.log" 2>/dev/null && exit 1
      sleep 30
    done
    CASE=$case SEED=$seed ARM=$arm Rscript benchmarks/production-adaptation/production.R \
      > "$OUT/logs/prod_${case}_${arm}_${seed}.log" 2>&1
  fi
}
export -f job
export OUT
{
  for cs in $CASES; do echo "warmup ${cs%%:*} ${cs##*:}"; done
  # Arms in reverse case order: the fast cases' warm-ups finish first, so their
  # arms take the free slots instead of waiting on a slow warm-up.
  for cs in $(echo $CASES | tr ' ' '\n' | tac); do for a in $ARMS; do echo "prod ${cs%%:*} ${cs##*:} $a"; done; done
} | xargs -P "$JOBS" -L 1 bash -c 'job "$@"' _
echo "all done"
