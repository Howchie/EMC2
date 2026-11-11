# Replay a saved hierarchical chain through two sample blocks.
suppressPackageStartupMessages(library(EMC2, lib.loc = Sys.getenv('EMC_LIB')))
input <- Sys.getenv('FIXTURE')
stopifnot(nzchar(input), file.exists(input))
options(emc2.fit_host = FALSE, emc2.sampler_profile = TRUE, emc2.sampler_profile_expensive = TRUE)
emc <- EMC2:::restore_custom_kernel_pointers(readRDS(input), quiet = TRUE)
set.seed(as.integer(Sys.getenv('SEED', '101')))
start <- proc.time()
out <- suppressWarnings(EMC2:::run_emc(emc, stage = 'sample',
  stop_criteria = list(iter = 20), step_size = 10, verbose = FALSE,
  cores_per_chain = as.integer(Sys.getenv('CORES_PER_CHAIN', '2')),
  cores_for_chains = as.integer(Sys.getenv('CHAINS', '2'))))
timing <- proc.time() - start
saveRDS(list(fit = out, seed = .Random.seed, timing = timing,
             blocks = EMC2:::.emc_block_profile(), iterations = EMC2:::.emc_iteration_profile(),
             fixture_md5 = unname(tools::md5sum(input)), session = sessionInfo()), Sys.getenv('OUT'))
cat('timing', timing, '\n')
