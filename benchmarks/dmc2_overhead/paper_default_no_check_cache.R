# Diagnostic-cache isolation arm: retain mode reuse, recompute every check.
suppressPackageStartupMessages(library(EMC2, lib.loc = Sys.getenv("EMC_LIB")))
utils::assignInNamespace(".emc_laplace_check_lookup",
  function(chain, subject, key) NULL, ns = "EMC2")
source("benchmarks/dmc2_overhead/paper_default.R", local = globalenv())
