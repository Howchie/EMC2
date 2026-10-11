pkgload::load_all(".", quiet = TRUE, compile = FALSE)
artifact_dir <- "plans/collapsed_transport_adaptation_2026-10-09"
fit_dir <- "/tmp/collapsed_transport_followup"
source(file.path(artifact_dir, "summarise_ess.R"))
snapshot <- readRDS(file.path(fit_dir, "setup_saturation_start.rds"))
previous <- snapshot$starts$corrected
corrected <- unserialize(serialize(previous, NULL))
for (j in seq_along(corrected)) {
  sol <- EMC2:::.emc_alloc_solve_units(snapshot$captured[[j]]$units)
  corrected[[j]]$collapsed$plan$rounds <- sol$joint
  if (EMC2:::.emc_shear_active(corrected[[j]]))
    corrected[[j]]$shear$plan$rounds <- sol$x[seq_len(2L * corrected[[j]]$n_pars)]
}
starts <- list(previous_budget = previous, incremental_budget = corrected)
for (arm in names(starts))
  cat(arm, "rates", paste(vapply(starts[[arm]], function(ch)
    ch$collapsed$plan$rounds, 0), collapse = ","), "\n")
# Reuse the completed first arm before the incremental-gain fix.
first <- readRDS(file.path(artifact_dir, "setup_saturation_ess.rds"))$corrected_1
stopifnot(!is.null(first))
first$label <- "previous_budget_1"
results <- list(previous_budget_1 = first)
for (replicate in 1:2) {
  order <- if (replicate == 1L) "incremental_budget" else rev(names(starts))
  for (arm in order) {
    start <- starts[[arm]]
    if (replicate > 1L) for (j in seq_along(start)) {
      advance <- function(seed) {
        for (k in seq_len(length(start))) seed <- parallel::nextRNGStream(seed)
        seed
      }
      start[[j]]$rng$chain <- advance(start[[j]]$rng$chain)
      start[[j]]$rng$gibbs <- advance(start[[j]]$rng$gibbs)
      start[[j]]$rng$subjects <- lapply(start[[j]]$rng$subjects, advance)
    }
    # Freeze the calibrated proposals and allocation for this comparison.
    for (j in seq_along(start)) start[[j]]$sample_proposals_next <- Inf
    started <- proc.time()[["elapsed"]]
    fit <- EMC2::fit(start, stage = "sample",
      stop_criteria = list(sample = list(iter = 200L)), step_size = 100L,
      cores_for_chains = 3L, cores_per_chain = 4L, particle_factor = 10,
      max_tries = 1L, trim = FALSE, verbose = FALSE, verboseProgress = FALSE)
    elapsed <- proc.time()[["elapsed"]] - started
    key <- paste(arm, replicate, sep = "_")
    saveRDS(fit, file.path(fit_dir, paste0(key, ".fit.rds")))
    z <- .summarise_ess(fit, key, 200L, first$rebuild_seconds, elapsed, 4L)
    results[[key]] <- z
    saveRDS(results, file.path(artifact_dir, "incremental_budget_ess.rds"))
    cat(key, "seconds", round(elapsed, 1), "group ESS", median(z$group_mu_bulk),
        "subject ESS", median(z$subject_52091_bulk), "max key Rhat", max(z$key_ess$rhat), "\n")
    flush.console()
  }
}
