# Create a short initialized checkpoint without changing the input fixture.
suppressPackageStartupMessages(library(EMC2, lib.loc = Sys.getenv("EMC_LIB")))
cat("EMC2 from", find.package("EMC2"), "\n")
e <- readRDS(Sys.getenv("EMC_IN"))
stopifnot(inherits(e, "emc"), !isTRUE(e[[1L]]$init))
set.seed(713)
e <- run_emc(e, "preburn", list(iter = 20), step_size = 20, max_tries = 1,
  cores_for_chains = 3, cores_per_chain = 3, verbose = TRUE)
saveRDS(e, Sys.getenv("EMC_OUT"))
