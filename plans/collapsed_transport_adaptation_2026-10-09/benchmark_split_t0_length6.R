pkgload::load_all(".", quiet = TRUE, compile = FALSE)
artifact_dir <- "plans/collapsed_transport_adaptation_2026-10-09"
fit_dir <- "/tmp/collapsed_transport_followup"
source(file.path(artifact_dir, "summarise_ess.R"))
prepared <- readRDS(file.path(fit_dir, "split_t0_allocation_start.rds"))
starts <- unserialize(serialize(prepared$start, NULL))
options(emc2.subject_hmc = FALSE)
for (j in seq_along(starts)) {
  ch <- starts[[j]]
  ch$production_options$emc2.subject_hmc <- FALSE
  ch$sample_proposals_next <- Inf
  ch$collapsed$plan$rounds <- 0
  starts[[j]] <- ch
}

ch <- starts[[2L]]
plan <- ch$collapsed$plan
plan$t0_ratio <- 0.5
plan$t0_mean_ratio <- 0.5
plan$t0_hyper_ratio <- 0.25
plan$length <- 6
plan$rounds <- 4
ctx <- EMC2:::.emc_collapsed_context(ch)
phi <- EMC2:::.emc_collapsed_to_phi(ch$samples$theta_mu[, ch$samples$idx],
  as.matrix(ch$samples$theta_var[, , ch$samples$idx]), ctx$chart)
C <- EMC2:::.emc_collapsed_preconditioner(plan$C, plan$scale,
  plan$t0_ratio, plan$t0_idx, plan$t0_hyper_ratio, plan$t0_mean_idx)
plan$eps <- EMC2:::.emc_collapsed_tune_eps(phi, C, plan$eps, ctx, plan$length)
plan$block <- EMC2:::.emc_collapsed_record(NULL, ch$n_pars)
ch$collapsed$plan <- plan
starts[[2L]] <- ch

started <- proc.time()[["elapsed"]]
fit <- EMC2::fit(starts, stage = "sample",
  stop_criteria = list(sample = list(iter = 200L)), step_size = 100L,
  cores_for_chains = 3L, cores_per_chain = 4L, particle_factor = 10,
  max_tries = 1L, trim = FALSE, verbose = FALSE, verboseProgress = FALSE)
elapsed <- proc.time()[["elapsed"]] - started
saveRDS(fit, file.path(fit_dir, "split_t0_length6.fit.rds"))
result <- .summarise_ess(fit, "split_t0_length6", 200L,
  prepared$rebuild_seconds, elapsed, 4L)
timing <- readRDS(file.path(artifact_dir, "split_t0_counterfactual_ess.rds"))
timing$split_t0_length6 <- result
saveRDS(timing, file.path(artifact_dir, "split_t0_counterfactual_ess.rds"))
t0 <- startsWith(starts[[1L]]$par_names, "t0")
total <- elapsed + prepared$rebuild_seconds
cat("split_t0_length6 sample seconds", round(elapsed, 1),
  "group ESS/sec", median(result$group_mu_bulk) / total,
  "t0 group ESS/sec", median(result$group_mu_bulk[t0]) / total,
  "participant ESS/sec", median(result$subject_52091_bulk) / total,
  "acceptance", paste(round(result$collapsed_block_acceptance, 3),
    collapse = ","), "\n")
