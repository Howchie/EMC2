library(posterior)

artifact_dir <- "plans/collapsed_transport_adaptation_2026-10-09"
fit_dir <- "/tmp/collapsed_transport_followup"
timing <- readRDS(file.path(artifact_dir, "split_t0_counterfactual_ess.rds"))
prepared <- readRDS(file.path(fit_dir, "split_t0_allocation_start.rds"))
arms <- c("no_collapsed", "split_t0")
summaries <- list()
per_variable <- list()

for (arm in arms) {
  fit <- readRDS(file.path(fit_dir,
    paste0("split_t0_counterfactual_", arm, ".fit.rds")))
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
  diag <- summarise_draws(as_draws_array(draws), ess_bulk, ess_tail, rhat)
  family <- c(rep("group_mean", p), rep("group_log_sd", p),
    rep("participant", p * n))
  total_seconds <- timing[[arm]]$sample_seconds + prepared$rebuild_seconds
  per_variable[[arm]] <- data.frame(arm = arm, variable = diag$variable,
    family = family, ess_bulk = diag$ess_bulk, ess_tail = diag$ess_tail,
    rhat = diag$rhat, total_seconds = total_seconds,
    ess_bulk_per_second = diag$ess_bulk / total_seconds)
  t0 <- which(startsWith(par_names, "t0"))
  group_idx <- seq_len(p)
  lsd_idx <- p + seq_len(p)
  participant_idx <- 2L * p + seq_len(p * n)
  summaries[[arm]] <- data.frame(arm = arm,
    sample_seconds = timing[[arm]]$sample_seconds,
    rebuild_seconds = prepared$rebuild_seconds,
    total_seconds = total_seconds,
    group_mean_ess_sec = median(diag$ess_bulk[group_idx]) / total_seconds,
    t0_group_mean_ess_sec = median(diag$ess_bulk[t0]) / total_seconds,
    other_group_mean_ess_sec = median(diag$ess_bulk[setdiff(group_idx, t0)]) /
      total_seconds,
    group_log_sd_ess_sec = median(diag$ess_bulk[lsd_idx]) / total_seconds,
    participant_ess_sec = median(diag$ess_bulk[participant_idx]) / total_seconds,
    sum_variable_ess_sec = sum(diag$ess_bulk, na.rm = TRUE) / total_seconds,
    max_rhat = max(diag$rhat, na.rm = TRUE),
    rhat_over_1.05 = sum(diag$rhat > 1.05, na.rm = TRUE))
  rm(fit, draws, diag)
  gc()
}

summary <- do.call(rbind, summaries)
long <- do.call(rbind, per_variable)
wide <- reshape(long[, c("variable", "family", "arm", "ess_bulk_per_second")],
  idvar = c("variable", "family"), timevar = "arm", direction = "wide")
wide$change_pct <- 100 * (wide$ess_bulk_per_second.split_t0 /
  wide$ess_bulk_per_second.no_collapsed - 1)
write.csv(summary, file.path(artifact_dir,
  "split_t0_counterfactual_family_summary.csv"), row.names = FALSE)
write.csv(long, file.path(artifact_dir,
  "split_t0_counterfactual_per_variable.csv"), row.names = FALSE)
write.csv(wide[order(wide$change_pct), ], file.path(artifact_dir,
  "split_t0_counterfactual_per_variable_change.csv"), row.names = FALSE)
print(summary, row.names = FALSE)
print(aggregate(change_pct ~ family, wide, function(x)
  c(median = median(x, na.rm = TRUE), improved = sum(x > 0, na.rm = TRUE),
    total = sum(is.finite(x)))), row.names = FALSE)
