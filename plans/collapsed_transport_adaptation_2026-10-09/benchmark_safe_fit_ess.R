pkgload::load_all(".", quiet = TRUE, compile = FALSE)
artifact_dir <- "plans/collapsed_transport_adaptation_2026-10-09"
fit_dir <- "/tmp/collapsed_transport_followup"
source(file.path(artifact_dir, "summarise_ess.R"))

starts <- EMC2:::restore_custom_kernel_pointers(
  EMC2:::restore_duplicates(readRDS(
    file.path(fit_dir, "fresh_retuned_calibrated.rds"))), quiet = TRUE)
stopifnot(length(starts) == 3L)
options(emc2.subject_hmc = FALSE)

calibration <- vector("list", length(starts))
rebuild_started <- proc.time()[["elapsed"]]
for (j in seq_along(starts)) {
  ch <- starts[[j]]
  plan <- ch$collapsed$plan
  plan$calibrate <- TRUE
  plan$calibration_screen <- FALSE
  plan$independent_scales <- FALSE
  plan$scale_updates <- 0L
  plan$t0_updates <- 0L
  plan$t0_hyper_updates <- 0L
  plan$length_updates <- 0L
  plan$rounds <- 0
  ch$collapsed$plan <- plan
  ch <- EMC2:::.emc_collapsed_calibrate(ch, n_cores = 4L)
  starts[[j]] <- ch
  calibration[[j]] <- ch$collapsed$plan$calibration
  cat("chain", j, "scale", ch$collapsed$plan$scale,
      "t0 mean", ch$collapsed$plan$scale * ch$collapsed$plan$t0_ratio,
      "t0 hierarchy", ch$collapsed$plan$scale * ch$collapsed$plan$t0_hyper_ratio,
      "acceptance", ch$collapsed$plan$calibration$acceptance,
      "seconds", ch$collapsed$plan$calibration$elapsed, "\n")
}
rebuild_seconds <- proc.time()[["elapsed"]] - rebuild_started
for (j in seq_along(starts)) {
  starts[[j]]$production_options$emc2.subject_hmc <- FALSE
  starts[[j]]$sample_proposals_next <- Inf
  starts[[j]]$collapsed$plan$rounds <- 0
}
saveRDS(list(start = starts, calibration = calibration,
             rebuild_seconds = rebuild_seconds),
  file.path(fit_dir, "safe_fit_ess_start.rds"))

make_forced <- function(start) {
  start[[2L]]$collapsed$plan$rounds <- 4
  start[[2L]]$collapsed$plan$block <-
    EMC2:::.emc_collapsed_record(NULL, start[[2L]]$n_pars)
  start
}
arms <- list(
  no_collapsed = unserialize(serialize(starts, NULL)),
  safe_fit_split_t0 = make_forced(unserialize(serialize(starts, NULL)))
)
results <- list()
for (arm in c("safe_fit_split_t0", "no_collapsed")) {
  started <- proc.time()[["elapsed"]]
  fit <- EMC2::fit(arms[[arm]], stage = "sample",
    stop_criteria = list(sample = list(iter = 200L)), step_size = 100L,
    cores_for_chains = 3L, cores_per_chain = 4L, particle_factor = 10,
    max_tries = 1L, trim = FALSE, verbose = FALSE, verboseProgress = FALSE)
  elapsed <- proc.time()[["elapsed"]] - started
  saveRDS(fit, file.path(fit_dir, paste0("safe_fit_ess_", arm, ".fit.rds")))
  results[[arm]] <- .summarise_ess(fit, paste0("safe_fit_ess_", arm), 200L,
    rebuild_seconds, elapsed, 4L)
  cat(arm, "sample seconds", round(elapsed, 1),
      "total seconds", round(elapsed + rebuild_seconds, 1),
      "group ESS/sec", median(results[[arm]]$group_mu_bulk /
        (elapsed + rebuild_seconds)),
      "t0 group ESS/sec", median(results[[arm]]$group_mu_bulk[
        startsWith(starts[[1L]]$par_names, "t0")] /
        (elapsed + rebuild_seconds)),
      "participant ESS/sec", median(results[[arm]]$subject_52091_bulk /
        (elapsed + rebuild_seconds)),
      "acceptance", paste(round(results[[arm]]$collapsed_block_acceptance, 3),
        collapse = ","), "\n")
  saveRDS(results, file.path(artifact_dir, "safe_fit_ess_comparison.rds"))
  rm(fit)
  gc()
}
