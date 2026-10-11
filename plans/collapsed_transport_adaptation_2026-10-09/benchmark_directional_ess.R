pkgload::load_all(".", quiet = TRUE, compile = FALSE)
artifact_dir <- "plans/collapsed_transport_adaptation_2026-10-09"
fit_dir <- "/tmp/collapsed_transport_followup"
source(file.path(artifact_dir, "summarise_ess.R"))
options(emc2.subject_hmc = FALSE)

starts <- EMC2:::restore_custom_kernel_pointers(
  EMC2:::restore_duplicates(readRDS(
    file.path(fit_dir, "safe_fit_ess_start.rds"))$start), quiet = TRUE)
base_chain2_calibration <- starts[[2L]]$collapsed$plan$calibration$elapsed
selected <- readRDS(file.path(fit_dir, "mean_scaled_calibrated.rds"))
stopifnot(length(starts) == 3L,
  identical(starts[[2L]]$samples$idx, selected$samples$idx),
  isTRUE(all.equal(starts[[2L]]$samples$theta_mu,
                   selected$samples$theta_mu, tolerance = 0)),
  isTRUE(all.equal(starts[[2L]]$collapsed$plan$C,
                   selected$collapsed$plan$C, tolerance = 1e-12)))

# The comparison's baseline uses the saved common scale 0.04419, t0 effective
# scales 0.04419 / 0.01105, and four rounds on chain 2. Replace only its chain-2
# plan with the targeted calibration; preserve the same state and RNG streams.
chain <- starts[[2L]]
plan <- selected$collapsed$plan
plan$rounds <- 4
plan$calibrate <- FALSE
plan$calibration_screen <- FALSE
plan$stage <- "sample"
plan$block <- EMC2:::.emc_collapsed_record(NULL, chain$n_pars)
plan$last <- EMC2:::.emc_collapsed_record(NULL, chain$n_pars)
plan$pool <- EMC2:::.emc_collapsed_record(NULL, chain$n_pars)
chain$collapsed$plan <- plan
chain$production_options$emc2.subject_hmc <- FALSE
starts[[2L]] <- chain
for (j in seq_along(starts)) {
  starts[[j]]$production_options$emc2.subject_hmc <- FALSE
  starts[[j]]$sample_proposals_next <- Inf
  starts[[j]]$collapsed$plan$rounds <- if (j == 2L) 4 else 0
}

started <- proc.time()[["elapsed"]]
fit <- EMC2::fit(starts, stage = "sample",
  stop_criteria = list(sample = list(iter = 200L)), step_size = 100L,
  cores_for_chains = 3L, cores_per_chain = 4L, particle_factor = 10,
  max_tries = 1L, trim = FALSE, verbose = FALSE, verboseProgress = FALSE)
sample_seconds <- proc.time()[["elapsed"]] - started
saveRDS(fit, file.path(fit_dir, "safe_fit_ess_directional_200.fit.rds"))

# Existing fixed-scale arm includes the common three-chain calibration cost.
# The targeted arm replaces chain 2's original calibration with the saved
# directional one, whose observed serial cost was 118.152 sec.
fixed <- readRDS(file.path(fit_dir,
  "safe_fit_ess_safe_fit_split_t0.fit.rds"))
fixed_summary <- .summarise_ess(fixed, "fixed_scale", 200L,
  rebuild_seconds = 45.601, sample_seconds = 399.101,
  cores_per_chain = 4L)
directional_summary <- .summarise_ess(fit, "directional_scale", 200L,
  rebuild_seconds = 45.601 - base_chain2_calibration + 118.152,
  sample_seconds = sample_seconds, cores_per_chain = 4L)
comparison <- list(fixed = fixed_summary,
                   directional = directional_summary,
                   base_chain2_calibration_seconds = base_chain2_calibration,
                   sample_seconds = sample_seconds,
                   directional_calibration_seconds = 118.152)
saveRDS(comparison,
  file.path(artifact_dir, "directional_ess_comparison.rds"))
cat("directional sample seconds", round(sample_seconds, 1),
    "total incl shared + directional calibration",
    round(directional_summary$total_seconds, 1), "\n")
cat("fixed median group ESS/sec",
    median(fixed_summary$group_mu_bulk) / fixed_summary$total_seconds,
    "directional",
    median(directional_summary$group_mu_bulk) /
      directional_summary$total_seconds, "\n")
vg_idx <- match("v.g", starts[[1L]]$par_names)
cat("v.g group ESS fixed", fixed_summary$group_mu_bulk[vg_idx],
    "directional", directional_summary$group_mu_bulk[vg_idx],
    "fixed ESS/sec", fixed_summary$group_mu_bulk[vg_idx] /
      fixed_summary$total_seconds, "directional ESS/sec",
    directional_summary$group_mu_bulk[vg_idx] /
      directional_summary$total_seconds, "\n")
cat("collapsed rounds",
    paste(directional_summary$collapsed_rounds, collapse = ","),
    "acceptance",
    paste(round(directional_summary$collapsed_block_acceptance, 3),
      collapse = ","), "\n")
