#!/usr/bin/env bash
# Builds the T30 sample-stage checkpoint the profiling runs start from.
#   benchmarks/ensemble-cost/mk_t30.sh <lib>
lib=$1
fixtures=/data/work/EMC2_workstreams/fixtures
cd "$(dirname "$0")/../.."
export EMC_LIB=$lib R_LIBS=$lib
export CORES_PER_CHAIN=9 SAVE_EMC=1
taskset -c 0-26 Rscript benchmarks/subject-hmc/fit_arm.R T30 fullc 1 $fixtures/T30_seed1.rds
