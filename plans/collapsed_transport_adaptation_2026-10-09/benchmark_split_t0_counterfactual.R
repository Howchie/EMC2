pkgload::load_all(".", quiet = TRUE, compile = FALSE)
artifact_dir <- "plans/collapsed_transport_adaptation_2026-10-09"
fit_dir <- "/tmp/collapsed_transport_followup"
source(file.path(artifact_dir, "summarise_ess.R"))
prepared <- readRDS(file.path(fit_dir, "split_t0_allocation_start.rds"))
base <- prepared$start
t0_mu <- startsWith(base[[1L]]$par_names, "t0")
options(emc2.subject_hmc = FALSE)
for (j in seq_along(base)) {
  base[[j]]$production_options$emc2.subject_hmc <- FALSE
  base[[j]]$sample_proposals_next <- Inf
  base[[j]]$collapsed$plan$rounds <- 0
}

make_allocated <- function(start) {
  ch <- start[[2L]]
  plan <- ch$collapsed$plan
  plan$t0_ratio <- 0.5
  plan$t0_mean_ratio <- 0.5
  plan$t0_hyper_ratio <- 0.25
  plan$rounds <- 4
  ctx <- EMC2:::.emc_collapsed_context(ch)
  phi <- EMC2:::.emc_collapsed_to_phi(
    ch$samples$theta_mu[, ch$samples$idx],
    as.matrix(ch$samples$theta_var[, , ch$samples$idx]), ctx$chart)
  C <- EMC2:::.emc_collapsed_preconditioner(plan$C, plan$scale,
    plan$t0_ratio, plan$t0_idx, plan$t0_hyper_ratio, plan$t0_mean_idx)
  plan$eps <- EMC2:::.emc_collapsed_tune_eps(phi, C, plan$eps, ctx,
                                               plan$length)
  plan$block <- EMC2:::.emc_collapsed_record(NULL, ch$n_pars)
  ch$collapsed$plan <- plan
  start[[2L]] <- ch
  start
}

starts <- list(
  no_collapsed = unserialize(serialize(base, NULL)),
  split_t0 = make_allocated(unserialize(serialize(base, NULL))))
results <- list()
for (arm in c("split_t0", "no_collapsed")) {
  started <- proc.time()[["elapsed"]]
  fit <- EMC2::fit(starts[[arm]], stage = "sample",
    stop_criteria = list(sample = list(iter = 200L)), step_size = 100L,
    cores_for_chains = 3L, cores_per_chain = 4L, particle_factor = 10,
    max_tries = 1L, trim = FALSE, verbose = FALSE, verboseProgress = FALSE)
  elapsed <- proc.time()[["elapsed"]] - started
  saveRDS(fit, file.path(fit_dir,
    paste0("split_t0_counterfactual_", arm, ".fit.rds")))
  results[[arm]] <- .summarise_ess(fit,
    paste0("split_t0_counterfactual_", arm), 200L,
    prepared$rebuild_seconds, elapsed, 4L)
  cat(arm, "sample seconds", round(elapsed, 1),
      "group ESS/sec", median(results[[arm]]$group_mu_bulk /
        (elapsed + prepared$rebuild_seconds)),
      "t0 group ESS/sec", median(results[[arm]]$group_mu_bulk[t0_mu] /
        (elapsed + prepared$rebuild_seconds)),
      "participant ESS/sec", median(results[[arm]]$subject_52091_bulk /
        (elapsed + prepared$rebuild_seconds)),
      "acceptance", paste(round(results[[arm]]$collapsed_block_acceptance, 3),
        collapse = ","), "\n")
  saveRDS(results, file.path(artifact_dir,
    "split_t0_counterfactual_ess.rds"))
  rm(fit)
  gc()
}
