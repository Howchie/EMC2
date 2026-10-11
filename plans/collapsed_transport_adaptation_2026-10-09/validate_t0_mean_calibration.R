pkgload::load_all(".", quiet = TRUE, compile = FALSE)
fit_dir <- "/tmp/collapsed_transport_followup"
snapshot <- readRDS(file.path(fit_dir, "setup_saturation_start.rds"))
ch <- snapshot$starts$corrected[[1L]]
C <- ch$collapsed$plan$C
ch <- EMC2:::.emc_collapsed_set_plan(ch, C, "sample")
started <- proc.time()[["elapsed"]]
ch <- EMC2:::.emc_collapsed_calibrate(ch, n_cores = 1L)
elapsed <- proc.time()[["elapsed"]] - started
out <- list(scale = ch$collapsed$plan$scale,
  t0_mean_ratio = ch$collapsed$plan$t0_mean_ratio,
  acceptance = ch$collapsed$plan$calibration$acceptance,
  scale_search = ch$collapsed$plan$calibration$search,
  t0_mean_search = ch$collapsed$plan$calibration$t0_mean_search,
  seconds = elapsed)
print(out)
saveRDS(out, "plans/collapsed_transport_adaptation_2026-10-09/t0_mean_calibration_probe.rds")
