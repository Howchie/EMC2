library(posterior)
artifact_dir <- "plans/collapsed_transport_adaptation_2026-10-09"
fit_dir <- "/tmp/collapsed_transport_followup"
args <- commandArgs(trailingOnly = TRUE)
repeat_id <- if (length(args)) args[[1L]] else ""
fit_prefix <- if (identical(repeat_id, "repeat2"))
  "t0_policy_repeat2_" else "t0_policy_"
timing_file <- if (identical(repeat_id, "repeat2"))
  "t0_policy_allocation_replicate2_ess.rds" else "t0_policy_allocation_ess.rds"
output_suffix <- if (nzchar(repeat_id)) paste0("_", repeat_id) else ""
fits <- list(
  no_collapsed = readRDS(file.path(fit_dir,
    paste0(fit_prefix, "no_collapsed.fit.rds"))),
  t0_adaptive = readRDS(file.path(fit_dir,
    paste0(fit_prefix, "t0_adaptive.fit.rds")))
)
timing <- readRDS(file.path(artifact_dir, timing_file))
summaries <- list()
per_variable <- list()
for (arm in names(fits)) {
  fit <- fits[[arm]]
  p <- fit[[1L]]$n_pars
  n <- fit[[1L]]$n_subjects
  n_draws <- 200L
  par_names <- fit[[1L]]$par_names
  subjects <- dimnames(fit[[1L]]$samples$alpha)[[2L]]
  variables <- c(paste0("mu:", par_names), paste0("lsd:", par_names),
                 paste(rep(subjects, each = p), rep(par_names, n), sep = ":"))
  draws <- array(NA_real_, c(n_draws, length(fit), length(variables)),
                 dimnames = list(NULL, NULL, variables))
  for (j in seq_along(fit)) {
    s <- fit[[j]]$samples
    idx <- tail(which(s$stage[seq_len(s$idx)] == "sample"), n_draws)
    stopifnot(length(idx) == n_draws)
    draws[, j, seq_len(p)] <- t(s$theta_mu[, idx, drop = FALSE])
    draws[, j, p + seq_len(p)] <- t(vapply(idx, function(i)
      log(diag(s$theta_var[, , i])) / 2, numeric(p)))
    draws[, j, 2L * p + seq_len(p * n)] <-
      t(matrix(s$alpha[, , idx], p * n, n_draws))
  }
  diagnostic <- summarise_draws(as_draws_array(draws), ess_bulk, ess_tail, rhat)
  mu <- seq_len(p)
  lsd <- p + seq_len(p)
  alpha <- 2L * p + seq_len(p * n)
  target <- match("52091", subjects)
  target_alpha <- 2L * p + (target - 1L) * p + seq_len(p)
  t0 <- which(startsWith(par_names, "t0"))
  total_time <- timing[[arm]]$total_seconds
  family <- c(rep("group_mean", p), rep("group_log_sd", p),
              rep("participant", p * n))
  per_variable[[arm]] <- data.frame(arm = arm, variable = diagnostic$variable,
    family = family, ess_bulk = diagnostic$ess_bulk,
    ess_tail = diagnostic$ess_tail, rhat = diagnostic$rhat)
  summaries[[arm]] <- data.frame(arm = arm, sample_seconds = timing[[arm]]$sample_seconds,
    rebuild_seconds = timing[[arm]]$rebuild_seconds, total_seconds = total_time,
    group_mu_ess = median(diagnostic$ess_bulk[mu]),
    group_mu_ess_per_second = median(diagnostic$ess_bulk[mu]) / total_time,
    t0_group_mu_ess = median(diagnostic$ess_bulk[t0]),
    t0_group_mu_ess_per_second = median(diagnostic$ess_bulk[t0]) / total_time,
    other_group_mu_ess = median(diagnostic$ess_bulk[setdiff(mu, t0)]),
    group_log_sd_ess = median(diagnostic$ess_bulk[lsd]),
    participant_ess = median(diagnostic$ess_bulk[alpha]),
    participant_ess_per_second = median(diagnostic$ess_bulk[alpha]) / total_time,
    participant_52091_ess = median(diagnostic$ess_bulk[target_alpha]),
    participant_52091_ess_per_second = median(diagnostic$ess_bulk[target_alpha]) / total_time,
    max_rhat = max(diagnostic$rhat), rhat_over_1.01 = sum(diagnostic$rhat > 1.01),
    rhat_over_1.05 = sum(diagnostic$rhat > 1.05))
  print(summaries[[arm]], row.names = FALSE)
  rm(draws, diagnostic)
  gc()
}
summary <- do.call(rbind, summaries)
write.csv(summary, file.path(artifact_dir,
  paste0("t0_policy_allocation", output_suffix, "_family_summary.csv")),
  row.names = FALSE)
long <- do.call(rbind, per_variable)
wide <- reshape(long, idvar = c("variable", "family"), timevar = "arm",
                direction = "wide")
wide$bulk_ess_change <- wide$ess_bulk.t0_adaptive - wide$ess_bulk.no_collapsed
wide$bulk_ess_change_pct <- 100 * wide$bulk_ess_change / wide$ess_bulk.no_collapsed
wide$time_scaled_change <- wide$ess_bulk.t0_adaptive / summary$total_seconds[2L] -
  wide$ess_bulk.no_collapsed / summary$total_seconds[1L]
write.csv(wide[order(wide$time_scaled_change), ], file.path(artifact_dir,
  paste0("t0_policy_allocation", output_suffix, "_per_variable.csv")),
  row.names = FALSE)
