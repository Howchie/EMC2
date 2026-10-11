artifact_dir <- "plans/collapsed_transport_adaptation_2026-10-09"
fit_dir <- "/tmp/collapsed_transport_followup"
x <- readRDS(file.path(artifact_dir, "directional_ess_comparison.rds"))
starts <- readRDS(file.path(fit_dir, "safe_fit_ess_start.rds"))$start
old_chain2_calibration <- starts[[2L]]$collapsed$plan$calibration$elapsed

# Replace chain 2's old calibration with the measured directional calibration;
# chain 1, chain 3, and all sampling settings are common to both arms.
directional_rebuild <- 45.601 - old_chain2_calibration +
  x$directional_calibration_seconds
x$base_chain2_calibration_seconds <- old_chain2_calibration
x$directional$rebuild_seconds <- directional_rebuild
x$directional$total_seconds <- directional_rebuild + x$sample_seconds
saveRDS(x, file.path(artifact_dir, "directional_ess_comparison.rds"))

vg <- match("v.g", starts[[1L]]$par_names)
tab <- do.call(rbind, lapply(c("fixed", "directional"), function(arm) {
  z <- x[[arm]]
  data.frame(
    arm = arm,
    sample_seconds = z$sample_seconds,
    rebuild_seconds = z$rebuild_seconds,
    total_seconds = z$total_seconds,
    median_group_ess = median(z$group_mu_bulk),
    median_group_ess_per_second = median(z$group_mu_bulk) / z$total_seconds,
    median_group_ess_per_sample_second =
      median(z$group_mu_bulk) / z$sample_seconds,
    median_52091_ess = median(z$subject_52091_bulk),
    median_52091_ess_per_second =
      median(z$subject_52091_bulk) / z$total_seconds,
    v_g_group_ess = z$group_mu_bulk[vg],
    v_g_group_ess_per_second = z$group_mu_bulk[vg] / z$total_seconds,
    median_key_ess = median(z$key_ess$bulk),
    max_key_rhat = max(z$key_ess$rhat)
  )
}))
tab$arm <- factor(tab$arm, levels = c("fixed", "directional"))
tab <- tab[order(tab$arm), ]
tab$arm <- as.character(tab$arm)
write.csv(tab, file.path(artifact_dir, "directional_ess_summary.csv"),
  row.names = FALSE)
print(tab, row.names = FALSE)
cat("baseline chain-2 calibration", old_chain2_calibration,
    "sec; directional calibration", x$directional_calibration_seconds,
    "sec; proposal-pilot core count was 1\n")
