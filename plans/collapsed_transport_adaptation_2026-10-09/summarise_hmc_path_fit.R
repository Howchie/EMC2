args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 2L)
inputs <- lapply(args, readRDS)
rows <- list()
for (x in inputs) {
  for (family in c("group_mu_ess", "group_log_sd_ess",
                   "subject_52091_ess", "key_ess")) {
    d <- x[[family]]
    rows[[length(rows) + 1L]] <- data.frame(
      path = x$label, trajectory_length = x$trajectory_length,
      draws_per_chain = x$continuation, family = family,
      min_bulk_ess = min(d$ess_bulk), median_bulk_ess = median(d$ess_bulk),
      max_rhat = max(d$rhat), elapsed_seconds = x$elapsed_seconds,
      cgroup_cpu_seconds = x$cgroup_cpu_seconds,
      min_bulk_ess_per_cpu_second = min(d$ess_bulk) / x$cgroup_cpu_seconds,
      median_bulk_ess_per_cpu_second = median(d$ess_bulk) / x$cgroup_cpu_seconds,
      stringsAsFactors = FALSE)
  }
}
out <- do.call(rbind, rows)
write.csv(out, "plans/collapsed_transport_adaptation_2026-10-09/hmc_path_fit_comparison.csv",
          row.names = FALSE)
print(out, row.names = FALSE)

vars <- c("mu:A", "mu:B", "mu:B_lR_target", "mu:t0_LoadHigh", "mu:t0_LoadLow",
          "52091:A", "52091:B", "52091:B_lR_target",
          "52091:t0_LoadHigh", "52091:t0_LoadLow")
par_rows <- list()
for (x in inputs) {
  d <- x$key_ess[match(vars, x$key_ess$variable), ]
  par_rows[[length(par_rows) + 1L]] <- data.frame(
    path = x$label, trajectory_length = x$trajectory_length,
    variable = d$variable, bulk_ess = d$ess_bulk, rhat = d$rhat)
}
write.csv(do.call(rbind, par_rows),
  "plans/collapsed_transport_adaptation_2026-10-09/hmc_path_fit_key_parameters.csv",
  row.names = FALSE)
