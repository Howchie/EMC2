pkgload::load_all(".", quiet = TRUE, compile = FALSE)
fit <- EMC2:::restore_custom_kernel_pointers(
  EMC2:::restore_duplicates(readRDS("/tmp/lba_sampling_audit/snapshot.rds")),
  quiet = TRUE)
options(emc2.failure_policy = "strict")
results <- vector("list", length(fit))
calibrated <- fit
for (j in seq_along(fit)) {
  ch <- fit[[j]]
  C <- ch$collapsed$plan$C
  ch <- EMC2:::.emc_collapsed_set_plan(ch, C, "sample")
  stopifnot(!isTRUE(ch$collapsed$plan$calibration_screen),
            ch$collapsed$plan$scale == 1,
            ch$collapsed$plan$t0_ratio == 1)
  started <- proc.time()[["elapsed"]]
  ch <- EMC2:::.emc_collapsed_calibrate(ch, n_cores = 4L)
  elapsed <- proc.time()[["elapsed"]] - started
  p <- ch$collapsed$plan
  results[[j]] <- data.frame(chain = j, scale = p$scale,
    t0_mean_ratio = p$t0_ratio, t0_hyper_ratio = p$t0_hyper_ratio,
    effective_t0_mean_scale = p$scale * p$t0_ratio,
    effective_t0_hyper_scale = p$scale * p$t0_hyper_ratio,
    length = p$length, eps = p$eps,
    validation_acceptance = p$calibration$acceptance,
    calibration_seconds = elapsed, pilots = p$calibration$pilots,
    search_mode = p$calibration$search_mode)
  calibrated[[j]] <- ch
  saveRDS(p$calibration,
    sprintf("plans/collapsed_transport_adaptation_2026-10-09/fresh_calibration_chain%d.rds", j))
  cat("chain", j, "scale", p$scale, "t0 mean ratio", p$t0_ratio,
      "t0 hyper ratio", p$t0_hyper_ratio,
      "length", p$length, "eps", p$eps,
      "acceptance", p$calibration$acceptance,
      "seconds", elapsed, "pilots", p$calibration$pilots, "\n")
}
out <- do.call(rbind, results)
write.csv(out,
  "plans/collapsed_transport_adaptation_2026-10-09/fresh_retuned_calibration.csv",
  row.names = FALSE)
saveRDS(calibrated,
  "/tmp/collapsed_transport_followup/fresh_retuned_calibrated.rds")
print(out, row.names = FALSE)
