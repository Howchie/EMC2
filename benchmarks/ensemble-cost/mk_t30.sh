#!/usr/bin/env bash
# Builds the T30 sample-stage checkpoint with the baseline build.
cd /data/work/EMC2_ws5
export EMC_LIB=/data/work/EMC2_ws5_libs/base R_LIBS=/data/work/EMC2_ws5_libs/base
export CORES_PER_CHAIN=9 SAVE_EMC=1
taskset -c 0-26 Rscript /data/work/EMC2_ws5_libs/h/fit_arm_ws5.R T30 fullc 1 /data/work/EMC2_ws5_libs/ckpt/T30_seed1.rds
