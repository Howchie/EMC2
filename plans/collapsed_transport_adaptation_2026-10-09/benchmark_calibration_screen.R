pkgload::load_all(".", quiet = TRUE, compile = FALSE)
fit_dir <- "/tmp/collapsed_transport_followup"
chain <- readRDS(file.path(fit_dir, "mean_scaled_calibrated.rds"))
holder <- structure(list(chain), class = "emc")
chain <- EMC2:::restore_custom_kernel_pointers(holder, quiet = TRUE)[[1L]]
old_calibration <- chain$collapsed$plan$calibration
chain$collapsed$plan$calibrate <- TRUE
chain$collapsed$plan$calibration_screen <- TRUE
started <- proc.time()[["elapsed"]]
screened <- EMC2:::.emc_collapsed_calibrate(chain, n_cores = 1L)
elapsed <- proc.time()[["elapsed"]] - started
out <- data.frame(
  prior_pilots = old_calibration$pilots,
  prior_seconds = old_calibration$elapsed,
  screen_pilots = screened$collapsed$plan$calibration$pilots,
  screen_seconds = elapsed,
  prior_scale = chain$collapsed$plan$scale,
  screen_scale = screened$collapsed$plan$scale,
  prior_acceptance = old_calibration$acceptance,
  screen_acceptance = screened$collapsed$plan$calibration$acceptance
)
write.csv(out,
  "plans/collapsed_transport_adaptation_2026-10-09/calibration_screen_benchmark.csv",
  row.names = FALSE)
print(out, row.names = FALSE)
