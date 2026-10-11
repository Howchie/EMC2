pkgload::load_all(".", quiet = TRUE, compile = FALSE)
artifact_dir <- "plans/collapsed_transport_adaptation_2026-10-09"
fit_dir <- "/tmp/collapsed_transport_followup"
source(file.path(artifact_dir, "summarise_ess.R"))
prepared <- readRDS(file.path(fit_dir, "length_adaptive_allocation_start.rds"))
base <- prepared$start
options(emc2.subject_hmc = FALSE)

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
  base[[j]]$production_options$emc2.subject_hmc <- FALSE
  base[[j]]$sample_proposals_next <- Inf
}

starts <- list(no_collapsed = unserialize(serialize(base, NULL)),
               allocated = unserialize(serialize(base, NULL)))
for (j in seq_along(starts$no_collapsed))
  starts$no_collapsed[[j]]$collapsed$plan$rounds <- 0

results <- list()
for (arm in c("allocated", "no_collapsed")) {
  started <- proc.time()[["elapsed"]]
  fit <- EMC2::fit(starts[[arm]], stage = "sample",
    stop_criteria = list(sample = list(iter = 200L)), step_size = 100L,
    cores_for_chains = 3L, cores_per_chain = 4L, particle_factor = 10,
    max_tries = 1L, trim = FALSE, verbose = FALSE, verboseProgress = FALSE)
  elapsed <- proc.time()[["elapsed"]] - started
  key <- paste0("length_policy_repeat2_", arm)
  saveRDS(fit, file.path(fit_dir, paste0(key, ".fit.rds")))
  results[[arm]] <- .summarise_ess(fit, key, 200L,
    prepared$rebuild_seconds, elapsed, 4L)
  cat(key, "seconds", round(elapsed, 1), "group ESS/sec",
      median(results[[arm]]$group_mu_bulk /
               (elapsed + prepared$rebuild_seconds)),
      "participant ESS/sec",
      median(results[[arm]]$subject_52091_bulk /
               (elapsed + prepared$rebuild_seconds)),
      "collapsed acceptance",
      paste(round(results[[arm]]$collapsed_block_acceptance, 3), collapse = ","),
      "\n")
  flush.console()
  saveRDS(results,
    file.path(artifact_dir, "length_policy_repeat2_ess_comparison.rds"))
  rm(fit)
  gc()
}
