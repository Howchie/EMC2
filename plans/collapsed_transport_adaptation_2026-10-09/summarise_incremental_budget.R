library(posterior)
artifact_dir <- "plans/collapsed_transport_adaptation_2026-10-09"
fit_dir <- "/tmp/collapsed_transport_followup"
results <- readRDS(file.path(artifact_dir, "incremental_budget_ess.rds"))
summary_file <- file.path(artifact_dir, "incremental_budget_summary.csv")
cached <- if (file.exists(summary_file)) read.csv(summary_file) else NULL
summary <- list()
for (key in names(results)) {
  filename <- if (key == "previous_budget_1") "setup_saturation_corrected_1.fit.rds" else
    paste0(key, ".fit.rds")
  diagnostic_file <- file.path(artifact_dir, paste0(key, "_diagnostics.csv"))
  if (!is.null(cached) && file.exists(diagnostic_file)) {
    row <- match(key, paste(cached$arm, cached$replicate, sep = "_"))
    if (!is.na(row) && isTRUE(all.equal(cached$seconds[row], results[[key]]$sample_seconds,
                                       tolerance = 1e-12)) &&
        file.info(diagnostic_file)$mtime >= file.info(file.path(fit_dir, filename))$mtime) {
      summary[[key]] <- cached[row, , drop = FALSE]
      next
    }
  }
  fit <- readRDS(file.path(fit_dir, filename))
  n_draws <- 200L
  p <- fit[[1L]]$n_pars
  n <- fit[[1L]]$n_subjects
  subjects <- dimnames(fit[[1L]]$samples$alpha)[[2L]]
  pars <- fit[[1L]]$par_names
  variables <- c(paste0("mu:", pars), paste0("lsd:", pars),
                 paste(rep(subjects, each = p), rep(pars, n), sep = ":"))
  draws <- array(NA_real_, c(n_draws, length(fit), length(variables)),
                 dimnames = list(NULL, NULL, variables))
  for (j in seq_along(fit)) {
    s <- fit[[j]]$samples
    idx <- tail(which(s$stage[seq_len(s$idx)] == "sample"), n_draws)
    stopifnot(length(idx) == n_draws)
    draws[, j, seq_len(p)] <- t(s$theta_mu[, idx, drop = FALSE])
    draws[, j, p + seq_len(p)] <- t(vapply(idx, function(i)
      log(diag(s$theta_var[, , i])) / 2, numeric(p)))
    draws[, j, 2L * p + seq_len(p * n)] <- t(matrix(s$alpha[, , idx], p * n, n_draws))
  }
  diagnostic <- summarise_draws(as_draws_array(draws), ess_bulk, ess_tail, rhat)
  write.csv(as.data.frame(diagnostic), diagnostic_file,
            row.names = FALSE)
  group <- seq_len(p)
  lsd <- p + seq_len(p)
  alpha <- 2L * p + seq_len(p * n)
  z <- results[[key]]
  summary[[key]] <- data.frame(arm = sub("_[12]$", "", key),
    replicate = as.integer(sub(".*_", "", key)), seconds = z$sample_seconds,
    group_mu_ess = median(diagnostic$ess_bulk[group]),
    group_lsd_ess = median(diagnostic$ess_bulk[lsd]),
    all_subject_ess = median(diagnostic$ess_bulk[alpha]),
    subject_52091_ess = median(z$subject_52091_bulk),
    max_group_rhat = max(diagnostic$rhat[c(group, lsd)]),
    max_subject_rhat = max(diagnostic$rhat[alpha]),
    group_rhat_over_1.01 = sum(diagnostic$rhat[c(group, lsd)] > 1.01),
    subject_rhat_over_1.01 = sum(diagnostic$rhat[alpha] > 1.01))
  rm(fit, draws)
  gc()
}
summary <- do.call(rbind, summary)
rownames(summary) <- NULL
write.csv(summary, summary_file, row.names = FALSE)
print(summary, row.names = FALSE)
