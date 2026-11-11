# Run a fixed number of sample transitions from a saved EMC2 chain state.
# Environment: EMC_LIB, FIXTURE, PERSISTENT, OUT.
suppressPackageStartupMessages(library(EMC2, lib.loc = Sys.getenv("EMC_LIB")))
input <- Sys.getenv("FIXTURE")
output <- Sys.getenv("OUT")
if (!nzchar(input) || !file.exists(input) || !nzchar(output))
  stop("set FIXTURE and OUT to existing/desired RDS paths")

arm <- identical(Sys.getenv("PERSISTENT", "0"), "1")
options(emc2.persistent_chains = arm, emc2.sampler_profile = FALSE)
emc <- EMC2:::restore_custom_kernel_pointers(readRDS(input)$fit, quiet = TRUE)
input_md5 <- unname(tools::md5sum(input))
set.seed(916)
started <- proc.time()
out <- suppressWarnings(EMC2:::run_emc(
  emc, stage = "sample", stop_criteria = list(iter = 60L),
  step_size = 20L, trim = FALSE, verbose = FALSE,
  cores_per_chain = 1L, cores_for_chains = length(emc)))
elapsed <- proc.time() - started
saveRDS(list(fit = out, caller_seed = .Random.seed, elapsed = elapsed,
             persistent = arm, input_md5 = input_md5,
             session = sessionInfo()), output)
cat("persistent", arm, "elapsed", elapsed[["elapsed"]],
    "input_md5", input_md5, "\n")
