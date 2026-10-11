pkgload::load_all(".", quiet = TRUE, compile = FALSE)
fit_dir <- "/tmp/collapsed_transport_followup"
chain <- readRDS(file.path(fit_dir, "mean_scaled_calibrated.rds"))
holder <- structure(list(chain), class = "emc")
chain <- EMC2:::restore_custom_kernel_pointers(holder, quiet = TRUE)[[1L]]

run_one <- function(screen) {
  x <- unserialize(serialize(chain, NULL))
  x <- EMC2:::restore_custom_kernel_pointers(
    structure(list(x), class = "emc"), quiet = TRUE)[[1L]]
  x$collapsed$plan$calibrate <- TRUE
  x$collapsed$plan$calibration_screen <- screen
  started <- proc.time()[["elapsed"]]
  x <- EMC2:::.emc_collapsed_calibrate(x, n_cores = 1L)
  elapsed <- proc.time()[["elapsed"]] - started
  plan <- x$collapsed$plan
  data.frame(mode = plan$calibration$search_mode,
    pilots = plan$calibration$pilots, elapsed_seconds = elapsed,
    scale = plan$scale, t0_mean_scale = plan$scale * plan$t0_ratio,
    t0_hyper_scale = plan$scale * plan$t0_hyper_ratio,
    mean_scale = if (length(plan$mean_idx))
      plan$scale * plan$mean_ratio else NA_real_,
    path_length = plan$length,
    acceptance = plan$calibration$acceptance,
    log_efficiency = plan$calibration$log_score,
    cost_per_round = plan$calibration$cost)
}

full <- run_one(FALSE)
screen <- run_one(TRUE)
out <- rbind(full, screen)
out$paired_screen_rng <- TRUE
write.csv(out,
  "plans/collapsed_transport_adaptation_2026-10-09/paired_rebuild_calibration.csv",
  row.names = FALSE)
print(out, row.names = FALSE)
