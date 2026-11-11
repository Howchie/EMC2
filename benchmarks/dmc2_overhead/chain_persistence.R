# Run a fixed number of sample transitions from a saved EMC2 chain state.
# Environment: EMC_LIB, FIXTURE, PERSISTENT, OUT.
suppressPackageStartupMessages(library(EMC2, lib.loc = Sys.getenv("EMC_LIB")))
input <- Sys.getenv("FIXTURE")
output <- Sys.getenv("OUT")
if (!nzchar(input) || !file.exists(input) || !nzchar(output))
  stop("set FIXTURE and OUT to existing/desired RDS paths")

arm <- identical(Sys.getenv("PERSISTENT", "0"), "1")
memory_file <- Sys.getenv("MEMORY_FILE", "")
if (nzchar(memory_file)) unlink(memory_file)
options(emc2.persistent_chains = arm, emc2.sampler_profile = FALSE,
        emc2.chain_process_memory_file = if (nzchar(memory_file)) memory_file else NULL,
        emc2.persistent_chain_restart_blocks = as.integer(
          Sys.getenv("RESTART_BLOCKS", "50")))
emc <- EMC2:::restore_custom_kernel_pointers(readRDS(input)$fit, quiet = TRUE)
input_md5 <- unname(tools::md5sum(input))
set.seed(916)
started <- proc.time()
out <- suppressWarnings(EMC2:::run_emc(
  emc, stage = "sample",
  stop_criteria = list(iter = as.integer(Sys.getenv("ITER", "60"))),
  step_size = as.integer(Sys.getenv("STEP_SIZE", "20")),
  trim = FALSE, verbose = FALSE,
  cores_per_chain = 1L, cores_for_chains = length(emc)))
elapsed <- proc.time() - started
saveRDS(list(fit = out, caller_seed = .Random.seed, elapsed = elapsed,
             persistent = arm, input_md5 = input_md5,
             session = sessionInfo()), output)
cat("persistent", arm, "elapsed", elapsed[["elapsed"]],
    "input_md5", input_md5, "\n")
