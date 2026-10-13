#!/usr/bin/env bash
# Paired per-part ensemble profiles, base and head alternating, pinned to 9 CPUs.
cd /data/work/EMC2_ws5
out=/data/work/EMC2_ws5_libs/prof; mkdir -p $out
fixtures="forstmann:/data/work/EMC2_workstreams/fixtures/forstmann_full_fit.rds:200 T30:/data/work/EMC2_ws5_libs/ckpt/T30_seed1_emc.rds:200 lba7r:/data/work/EMC2_workstreams/fixtures/lba7r_sample_stage.rds:100"
[ -f /data/work/EMC2_ws5_libs/ckpt/nback31_seed1_emc.rds ] && fixtures="$fixtures nback31:/data/work/EMC2_ws5_libs/ckpt/nback31_seed1_emc.rds:200"
for rep in 1 2; do
  for lib in base head; do
    [ $rep = 2 ] && [ $lib = base ] && continue
    for f in $fixtures; do
      IFS=: read name file iter <<< "$f"
      EMC_LIB=/data/work/EMC2_ws5_libs/$lib R_LIBS=/data/work/EMC2_ws5_libs/$lib EMC_IN=$file ITER=$iter CORES=3 taskset -c 23-31 Rscript /data/work/EMC2_ws5_libs/h/prof_ens2.R $out/pair_${lib}_${name}_r$rep.rds > $out/pair_${lib}_${name}_r$rep.log 2>&1
      echo "$lib $name rep $rep exit $?"
    done
  done
done
