pkgload::load_all(".", quiet = TRUE, compile = FALSE)
artifact_dir <- "plans/collapsed_transport_adaptation_2026-10-09"
fit_dir <- "/tmp/collapsed_transport_followup"
source(file.path(artifact_dir, "summarise_ess.R"))
prepared <- readRDS(file.path(fit_dir, "t0_adaptive_allocation_start.rds"))
base <- prepared$start

advance <- function(seed, n_streams) {
  if (is.null(seed)) return(seed)
  for (k in seq_len(n_streams)) seed <- parallel::nextRNGStream(seed)
  seed
}
for (j in seq_along(base)) {
  base[[j]]$rng$chain <- advance(base[[j]]$rng$chain, length(base))
  base[[j]]$rng$gibbs <- advance(base[[j]]$rng$gibbs, length(base))
  base[[j]]$rng$subjects <- lapply(base[[j]]$rng$subjects, advance,
                                  n_streams = length(base))
  base[[j]]$sample_proposals_next <- Inf
}

starts <- list(no_collapsed = unserialize(serialize(base, NULL)),
               t0_adaptive = unserialize(serialize(base, NULL)))
starts$t0_adaptive[[1L]]$collapsed$plan$rounds <- 0.875
for (j in 2:length(starts$t0_adaptive))
  starts$t0_adaptive[[j]]$collapsed$plan$rounds <- 0

results <- list()
for (arm in c("t0_adaptive", "no_collapsed")) {
  started <- proc.time()[["elapsed"]]
  fit <- EMC2::fit(starts[[arm]], stage = "sample",
    stop_criteria = list(sample = list(iter = 200L)), step_size = 100L,
    cores_for_chains = 3L, cores_per_chain = 4L, particle_factor = 10,
    max_tries = 1L, trim = FALSE, verbose = FALSE, verboseProgress = FALSE)
  elapsed <- proc.time()[["elapsed"]] - started
  key <- paste0("t0_policy_repeat2_", arm)
  saveRDS(fit, file.path(fit_dir, paste0(key, ".fit.rds")))
  results[[arm]] <- .summarise_ess(fit, key, 200L,
    prepared$rebuild_seconds, elapsed, 4L)
  cat(key, "seconds", round(elapsed, 1), "group ESS",
      median(results[[arm]]$group_mu_bulk), "participant ESS",
      median(results[[arm]]$subject_52091_bulk), "collapsed acceptance",
      paste(round(results[[arm]]$collapsed_block_acceptance, 3), collapse = ","),
      "\n")
  flush.console()
  saveRDS(results, file.path(artifact_dir,
    "t0_policy_allocation_replicate2_ess.rds"))
  rm(fit)
  gc()
}
