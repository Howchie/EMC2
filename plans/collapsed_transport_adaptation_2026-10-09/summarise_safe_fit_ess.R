artifact_dir <- "plans/collapsed_transport_adaptation_2026-10-09"
fit_dir <- "/tmp/collapsed_transport_followup"
results <- readRDS(file.path(artifact_dir, "safe_fit_ess_comparison.rds"))
starts <- readRDS(file.path(fit_dir, "safe_fit_ess_start.rds"))$start
par_names <- starts[[1L]]$par_names
t0 <- startsWith(par_names, "t0")

summary <- do.call(rbind, lapply(names(results), function(arm) {
  x <- results[[arm]]
  total <- x$total_seconds
  data.frame(
    arm = arm,
    sample_seconds = x$sample_seconds,
    rebuild_seconds = x$rebuild_seconds,
    total_seconds = total,
    collapsed_rounds = paste(x$collapsed_rounds, collapse = ","),
    chain2_acceptance = x$collapsed_block_acceptance[2L],
    group_mean_ess_sec = median(x$group_mu_bulk) / total,
    t0_group_mean_ess_sec = median(x$group_mu_bulk[t0]) / total,
    other_group_mean_ess_sec = median(x$group_mu_bulk[!t0]) / total,
    participant_52091_ess_sec = median(x$subject_52091_bulk) / total,
    key_ess_median = median(x$key_ess$bulk),
    key_max_rhat = max(x$key_ess$rhat, na.rm = TRUE),
    key_rhat_over_1_05 = sum(x$key_ess$rhat > 1.05, na.rm = TRUE)
  )
}))
summary$arm <- factor(summary$arm, levels = c("no_collapsed", "safe_fit_split_t0"))
summary <- summary[order(summary$arm), ]
summary$arm <- as.character(summary$arm)
write.csv(summary, file.path(artifact_dir, "safe_fit_ess_summary.csv"),
          row.names = FALSE)
print(summary, row.names = FALSE)
