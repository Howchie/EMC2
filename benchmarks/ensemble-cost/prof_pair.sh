#!/usr/bin/env bash
# Paired per-part ensemble profiles, base and head alternating, pinned to 9 CPUs.
#   benchmarks/ensemble-cost/prof_pair.sh <base lib> <head lib> <out dir>
base=$1; head=$2; out=$3; mkdir -p $out
here=$(cd "$(dirname "$0")" && pwd)
fx=/data/work/EMC2_workstreams/fixtures
fixtures="forstmann:$fx/forstmann_full_fit.rds:200 T30:$fx/T30_seed1_emc.rds:200 lba7r:$fx/lba7r_sample_stage.rds:100"
[ -f $fx/nback31_seed1_emc.rds ] && fixtures="$fixtures nback31:$fx/nback31_seed1_emc.rds:200"
for rep in 1 2; do
  for lib in base head; do
    [ $rep = 2 ] && [ $lib = base ] && continue
    path=${!lib}
    for f in $fixtures; do
      IFS=: read name file iter <<< "$f"
      EMC_LIB=$path R_LIBS=$path EMC_IN=$file ITER=$iter CORES=3 taskset -c 23-31 Rscript $here/prof_ens2.R $out/pair_${lib}_${name}_r$rep.rds > $out/pair_${lib}_${name}_r$rep.log 2>&1
      echo "$lib $name rep $rep exit $?"
    done
  done
done
