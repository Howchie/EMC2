pkgload::load_all(".", quiet = TRUE, compile = FALSE)

start_file <- Sys.getenv("STAGE_CAL_START",
  "/tmp/collapsed_transport_followup/mean_scaled_calibrated.rds")
out_file <- Sys.getenv("STAGE_CAL_OUT",
  "plans/collapsed_transport_adaptation_2026-10-09/stage_transition_calibration.csv")
n_cores <- as.integer(Sys.getenv("STAGE_CAL_CORES", "4"))

saved <- readRDS(start_file)
chain <- if (inherits(saved, "emc")) saved[[1L]] else saved
chain <- EMC2:::restore_custom_kernel_pointers(
  structure(list(chain), class = "emc"), quiet = TRUE)[[1L]]
stopifnot(is.list(chain$collapsed$plan$calibration),
          is.finite(chain$collapsed$plan$scale),
          chain$collapsed$plan$scale > 0)

run_one <- function(screen) {
  x <- unserialize(serialize(chain, NULL))
  x <- EMC2:::restore_custom_kernel_pointers(
    structure(list(x), class = "emc"), quiet = TRUE)[[1L]]
  x$collapsed$plan$stage <- "adapt"
  x <- EMC2:::.emc_collapsed_set_plan(
    x, x$collapsed$plan$C, stage = "sample")
  auto_screen <- x$collapsed$plan$calibration_screen
  stopifnot(isTRUE(auto_screen))
  if (!screen) x$collapsed$plan$calibration_screen <- FALSE
  started <- proc.time()[["elapsed"]]
  x <- EMC2:::.emc_collapsed_calibrate(x, n_cores = n_cores)
  elapsed <- proc.time()[["elapsed"]] - started
  p <- x$collapsed$plan
  data.frame(
    mode = p$calibration$search_mode,
    pilots = p$calibration$pilots,
    elapsed_seconds = elapsed,
    scale = p$scale,
    t0_mean_scale = p$scale * p$t0_ratio,
    t0_hyper_scale = p$scale * p$t0_hyper_ratio,
    mean_index = if (length(p$mean_idx)) p$mean_idx else NA_integer_,
    mean_scale = if (length(p$mean_idx)) p$scale * p$mean_ratio else NA_real_,
    path_length = p$length,
    acceptance = p$calibration$acceptance,
    log_efficiency = p$calibration$log_score,
    fixed_cost = p$calibration$fixed_cost,
    stringsAsFactors = FALSE)
}

full <- run_one(FALSE)
screen <- run_one(TRUE)
out <- rbind(full, screen)
write.csv(out, out_file, row.names = FALSE)
print(out, row.names = FALSE)
