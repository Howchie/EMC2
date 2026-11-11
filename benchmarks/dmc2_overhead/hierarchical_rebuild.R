# Replay a saved hierarchical chain through two sample blocks.
suppressPackageStartupMessages(library(EMC2, lib.loc = Sys.getenv('EMC_LIB')))
input <- Sys.getenv('FIXTURE')
stopifnot(nzchar(input), file.exists(input))
options(emc2.fit_host = FALSE, emc2.sampler_profile = TRUE)
.emc_dispatch_log <- list()
trace(EMC2:::.emc_rebuild_lapply, exit = quote({
  decision <- .emc_rebuild_queue$last_decision
  if (!is.null(decision)) {
    log <- get('.emc_dispatch_log', envir = .GlobalEnv)
    assign('.emc_dispatch_log', c(log, list(decision)), envir = .GlobalEnv)
  }
}), print = FALSE)
emc <- EMC2:::restore_custom_kernel_pointers(readRDS(input), quiet = TRUE)
for (j in seq_along(emc)) emc[[j]]$sample_proposals_next <- NULL
set.seed(as.integer(Sys.getenv('SEED', '101')))
start <- proc.time()
out <- suppressWarnings(EMC2:::run_emc(emc, stage = 'sample',
  stop_criteria = list(iter = 200), step_size = 100, verbose = FALSE,
  cores_per_chain = as.integer(Sys.getenv('CORES_PER_CHAIN', '2')),
  cores_for_chains = as.integer(Sys.getenv('CHAINS', '2'))))
timing <- proc.time() - start
saveRDS(list(fit = out, seed = .Random.seed, timing = timing,
             blocks = EMC2:::.emc_block_profile(), iterations = EMC2:::.emc_iteration_profile(),
             rebuild_dispatch = .emc_dispatch_log,
             fixture_md5 = unname(tools::md5sum(input)), session = sessionInfo()), Sys.getenv('OUT'))
cat('timing', timing, '\n')
