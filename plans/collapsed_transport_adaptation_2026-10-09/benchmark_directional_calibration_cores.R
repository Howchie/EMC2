pkgload::load_all(".", quiet = TRUE, compile = FALSE)
options(emc2.failure_policy = "strict")
fit_dir <- "/tmp/collapsed_transport_followup"
artifact_dir <- "plans/collapsed_transport_adaptation_2026-10-09"
chain <- readRDS(file.path(fit_dir, "safe_fit_ess_start.rds"))$start[[2L]]
plan <- chain$collapsed$plan
for (field in c("scale_updates", "t0_updates", "t0_hyper_updates",
                "mean_updates", "length_updates")) plan[[field]] <- 0L
plan$calibrate <- TRUE
plan$calibration_screen <- FALSE
chain$collapsed$plan <- plan

started <- proc.time()[["elapsed"]]
out <- tryCatch(
  EMC2:::.emc_collapsed_calibrate(chain, n_cores = 4L),
  error = function(e) structure(list(message = conditionMessage(e)),
    class = "directional_calibration_error"))
wall_seconds <- proc.time()[["elapsed"]] - started
if (inherits(out, "directional_calibration_error")) {
  result <- data.frame(status = "error", wall_seconds = wall_seconds,
    calibration_seconds = NA_real_, pilots = NA_integer_, scale = NA_real_,
    mean_idx = NA_integer_, mean_scale = NA_real_, acceptance = NA_real_,
    message = out$message)
} else {
  p <- out$collapsed$plan
  result <- data.frame(status = "ok", wall_seconds = wall_seconds,
    calibration_seconds = p$calibration$elapsed,
    pilots = p$calibration$pilots, scale = p$scale,
    mean_idx = if (length(p$mean_idx)) p$mean_idx else NA_integer_,
    mean_scale = if (length(p$mean_idx)) p$scale * p$mean_ratio else NA_real_,
    acceptance = p$calibration$acceptance, message = "")
}
write.csv(result,
  file.path(artifact_dir, "directional_calibration_4cores.csv"),
  row.names = FALSE)
saveRDS(result,
  file.path(artifact_dir, "directional_calibration_4cores.rds"))
print(result, row.names = FALSE)
