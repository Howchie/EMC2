library(posterior)
artifact_dir <- "plans/collapsed_transport_adaptation_2026-10-09"
fit_dir <- "/tmp/collapsed_transport_followup"
timing <- readRDS(file.path(artifact_dir, "safe_fit_ess_comparison.rds"))
arms <- c("no_collapsed", "safe_fit_split_t0")
long <- vector("list", length(arms))
names(long) <- arms

for (arm in arms) {
  fit <- readRDS(file.path(fit_dir,
    paste0("safe_fit_ess_", arm, ".fit.rds")))
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
  diag <- summarise_draws(as_draws_array(draws), ess_bulk, rhat)
  family <- c(rep("group_mean", p), rep("group_log_sd", p),
    rep("participant", p * n))
  parameter <- c(par_names, par_names, rep(par_names, n))
  total_seconds <- timing[[arm]]$total_seconds
  long[[arm]] <- data.frame(arm = arm, variable = diag$variable,
    family = family, parameter = parameter, ess_bulk = diag$ess_bulk,
    rhat = diag$rhat, total_seconds = total_seconds,
    ess_bulk_per_second = diag$ess_bulk / total_seconds)
  rm(fit, draws, diag)
  gc()
}

long <- do.call(rbind, long)
write.csv(long,
  file.path(artifact_dir, "safe_fit_ess_per_variable.csv"), row.names = FALSE)
wide <- reshape(long[, c("variable", "family", "parameter", "arm",
  "ess_bulk", "ess_bulk_per_second", "rhat")],
  idvar = c("variable", "family", "parameter"), timevar = "arm",
  direction = "wide")
wide$ess_change_pct <- 100 *
  (wide$ess_bulk.safe_fit_split_t0 / wide$ess_bulk.no_collapsed - 1)
wide$ess_sec_change_pct <- 100 *
  (wide$ess_bulk_per_second.safe_fit_split_t0 /
     wide$ess_bulk_per_second.no_collapsed - 1)
write.csv(wide[order(wide$ess_sec_change_pct), ],
  file.path(artifact_dir, "safe_fit_ess_per_variable_change.csv"),
  row.names = FALSE)

cat("Median variable ESS/sec change by family:\n")
print(aggregate(ess_sec_change_pct ~ family, wide, function(x)
  c(median = median(x, na.rm = TRUE),
    improved = sum(x > 0, na.rm = TRUE),
    total = sum(is.finite(x)))), row.names = FALSE)
cat("Median ESS/sec change by group parameter:\n")
group <- wide[wide$family == "group_mean", ]
print(aggregate(ess_sec_change_pct ~ parameter, group, median,
                na.rm = TRUE), row.names = FALSE)
cat("R-hat counts by family and arm (>1.05):\n")
print(aggregate(rhat ~ family + arm, long, function(x)
  sum(x > 1.05, na.rm = TRUE)), row.names = FALSE)
