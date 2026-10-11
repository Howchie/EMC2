pkgload::load_all(".", quiet = TRUE, compile = FALSE)
artifact_dir <- "plans/collapsed_transport_adaptation_2026-10-09"
fit_dir <- "/tmp/collapsed_transport_followup"
source(file.path(artifact_dir, "summarise_ess.R"))
snapshot <- readRDS(file.path(fit_dir, "setup_saturation_start.rds"))
start <- unserialize(serialize(snapshot$starts$corrected, NULL))

# Hold the overall calibrated scale, integrator step, allocation and retained
# state fixed. Shrink only the three t0 group-mean rows of the mass root.
for (j in seq_along(start)) {
  ch <- start[[j]]
  t0_mean <- which(grepl("^t0", ch$par_names))
  if (ch$collapsed$plan$rounds > 0) {
    row_scale <- rep(1, nrow(ch$collapsed$plan$C))
    row_scale[t0_mean] <- .25
    start[[j]]$collapsed$plan$C <- sweep(ch$collapsed$plan$C, 1L,
                                         row_scale, `*`)
    start[[j]]$collapsed$plan$t0_mean_ratio <- .25
  }
  start[[j]]$sample_proposals_next <- Inf
}

started <- proc.time()[["elapsed"]]
fit <- EMC2::fit(start, stage = "sample",
  stop_criteria = list(sample = list(iter = 200L)), step_size = 100L,
  cores_for_chains = 3L, cores_per_chain = 4L, particle_factor = 10,
  max_tries = 1L, trim = FALSE, verbose = FALSE, verboseProgress = FALSE)
elapsed <- proc.time()[["elapsed"]] - started
saveRDS(fit, file.path(fit_dir, "t0_mean_anisotropic_1.fit.rds"))
z <- .summarise_ess(fit, "t0_mean_anisotropic_1", 200L, NA_real_, elapsed, 4L)
saveRDS(z, file.path(artifact_dir, "t0_mean_anisotropic_ess.rds"))
cat("seconds", elapsed, "group mean ESS", median(z$group_mu_bulk),
    "subject ESS", median(z$subject_52091_bulk), "\n")
