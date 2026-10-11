pkgload::load_all(".", quiet = TRUE, compile = FALSE)
artifact_dir <- "plans/collapsed_transport_adaptation_2026-10-09"
prepared <- readRDS("/tmp/collapsed_transport_followup/length_adaptive_allocation_start.rds")
emc <- EMC2:::restore_custom_kernel_pointers(
  EMC2:::restore_duplicates(prepared$start), quiet = TRUE)
chain <- emc[[1L]]
old <- chain$collapsed$plan
stopifnot(identical(old$stage, "sample"), old$rounds == 0)

chain <- EMC2:::.emc_collapsed_set_plan(chain, old$C, "sample")
stopifnot(isTRUE(chain$collapsed$plan$calibration_screen))
started <- proc.time()[["elapsed"]]
chain <- EMC2:::.emc_collapsed_calibrate(chain, n_cores = 4L)
elapsed <- proc.time()[["elapsed"]] - started
plan <- chain$collapsed$plan
result <- data.frame(mode = plan$calibration$search_mode,
  scale = plan$scale, t0_ratio = plan$t0_ratio, length = plan$length,
  eps = plan$eps, validation_acceptance = plan$calibration$acceptance,
  pilots = plan$calibration$pilots, elapsed_seconds = elapsed)
write.csv(result, file.path(artifact_dir,
  "length_adaptive_calibration_screen.csv"), row.names = FALSE)
print(result, row.names = FALSE)
