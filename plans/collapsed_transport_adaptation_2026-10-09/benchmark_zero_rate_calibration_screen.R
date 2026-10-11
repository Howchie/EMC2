pkgload::load_all(".", quiet = TRUE, compile = FALSE)
prepared <- readRDS("/tmp/collapsed_transport_followup/t0_adaptive_allocation_start.rds")
emc <- EMC2:::restore_custom_kernel_pointers(
  EMC2:::restore_duplicates(prepared$start), quiet = TRUE)
ch <- emc[[1L]]
options(emc2.failure_policy = "strict")
old <- ch$collapsed$plan
stopifnot(is.list(old$calibration), identical(old$stage, "sample"), old$rounds == 0)
ch <- EMC2:::.emc_collapsed_set_plan(ch, old$C, "sample")
stopifnot(isTRUE(ch$collapsed$plan$calibration_screen))
started <- proc.time()[["elapsed"]]
ch <- EMC2:::.emc_collapsed_calibrate(ch, n_cores = 4L)
elapsed <- proc.time()[["elapsed"]] - started
cat("mode", ch$collapsed$plan$calibration$search_mode,
    "pilots", ch$collapsed$plan$calibration$pilots,
    "scale", ch$collapsed$plan$scale,
    "t0 ratio", ch$collapsed$plan$t0_ratio,
    "eps", ch$collapsed$plan$eps,
    "acceptance", ch$collapsed$plan$calibration$acceptance,
    "elapsed", elapsed, "seconds\n")
