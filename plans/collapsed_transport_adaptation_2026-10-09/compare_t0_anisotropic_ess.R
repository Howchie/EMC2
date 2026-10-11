library(posterior)
artifact_dir <- "plans/collapsed_transport_adaptation_2026-10-09"
fits <- list(
  global = readRDS("/tmp/collapsed_transport_followup/setup_saturation_corrected_1.fit.rds"),
  t0_mean_025 = readRDS("/tmp/collapsed_transport_followup/t0_mean_anisotropic_1.fit.rds"))
times <- c(global = 351.673,
           t0_mean_025 = readRDS(file.path(artifact_dir,
             "t0_mean_anisotropic_ess.rds"))$sample_seconds)
summaries <- list()
per_variable <- list()
for (arm in names(fits)) {
  fit <- fits[[arm]]
  p <- fit[[1L]]$n_pars
  n <- fit[[1L]]$n_subjects
  draws_per_chain <- 200L
  par_names <- fit[[1L]]$par_names
  subjects <- dimnames(fit[[1L]]$samples$alpha)[[2L]]
  variables <- c(paste0("mu:", par_names), paste0("lsd:", par_names),
                 paste(rep(subjects, each = p), rep(par_names, n), sep = ":"))
  draws <- array(NA_real_, c(draws_per_chain, length(fit), length(variables)),
                 dimnames = list(NULL, NULL, variables))
  for (j in seq_along(fit)) {
    s <- fit[[j]]$samples
    idx <- tail(which(s$stage[seq_len(s$idx)] == "sample"), draws_per_chain)
    stopifnot(length(idx) == draws_per_chain)
    draws[, j, seq_len(p)] <- t(s$theta_mu[, idx, drop = FALSE])
    draws[, j, p + seq_len(p)] <- t(vapply(idx, function(i)
      log(diag(s$theta_var[, , i])) / 2, numeric(p)))
    draws[, j, 2L * p + seq_len(p * n)] <-
      t(matrix(s$alpha[, , idx], p * n, draws_per_chain))
  }
  diag <- summarise_draws(as_draws_array(draws), ess_bulk, ess_tail, rhat)
  mu <- seq_len(p)
  lsd <- p + seq_len(p)
  alpha <- 2L * p + seq_len(p * n)
  subject_id <- "52091"
  subject_at <- match(subject_id, subjects)
  subject_param <- 2L * p + (subject_at - 1L) * p + seq_len(p)
  t0 <- which(grepl("^t0", par_names))
  family <- c(rep("group_mean", p), rep("group_log_sd", p),
              rep("participant", p * n))
  per_variable[[arm]] <- data.frame(arm = arm, variable = diag$variable,
    family = family, ess_bulk = diag$ess_bulk, ess_tail = diag$ess_tail,
    rhat = diag$rhat)
  summaries[[arm]] <- data.frame(
    arm = arm, sample_seconds = times[arm],
    group_mu_ess = median(diag$ess_bulk[mu]),
    group_mu_ess_per_second = median(diag$ess_bulk[mu]) / times[arm],
    t0_group_mu_ess = median(diag$ess_bulk[t0]),
    t0_group_mu_ess_per_second = median(diag$ess_bulk[t0]) / times[arm],
    other_group_mu_ess = median(diag$ess_bulk[setdiff(mu, t0)]),
    group_lsd_ess = median(diag$ess_bulk[lsd]),
    all_subject_ess = median(diag$ess_bulk[alpha]),
    all_subject_ess_per_second = median(diag$ess_bulk[alpha]) / times[arm],
    subject_52091_ess = median(diag$ess_bulk[subject_param]),
    subject_52091_ess_per_second = median(diag$ess_bulk[subject_param]) / times[arm],
    max_rhat = max(diag$rhat), rhat_over_1.01 = sum(diag$rhat > 1.01),
    rhat_over_1.05 = sum(diag$rhat > 1.05))
  print(summaries[[arm]], row.names = FALSE)
  rm(draws, diag)
  gc()
}
summary <- do.call(rbind, summaries)
write.csv(summary, file.path(artifact_dir,
  "t0_anisotropic_ess_comparison.csv"), row.names = FALSE)
long <- do.call(rbind, per_variable)
wide <- reshape(long, idvar = c("variable", "family"), timevar = "arm",
                direction = "wide")
wide$bulk_ess_change <- wide$ess_bulk.t0_mean_025 - wide$ess_bulk.global
wide$bulk_ess_change_pct <- 100 * wide$bulk_ess_change / wide$ess_bulk.global
write.csv(wide[order(wide$bulk_ess_change_pct), ],
  file.path(artifact_dir, "t0_anisotropic_per_variable.csv"), row.names = FALSE)
