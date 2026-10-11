pkgload::load_all(".", quiet = TRUE, compile = FALSE)
artifact_dir <- "plans/collapsed_transport_adaptation_2026-10-09"
fit_dir <- "/tmp/t30_validation_batch_starts"
n_draws <- 200L
options(emc2.collapsed_move = TRUE, emc2.subject_hmc = FALSE,
        emc2.failure_policy = "strict")

summarise_fit <- function(fit, label, seconds) {
  nch <- length(fit)
  p <- fit[[1L]]$n_pars
  nsub <- fit[[1L]]$n_subjects
  pn <- fit[[1L]]$par_names
  ids <- dimnames(fit[[1L]]$samples$alpha)[[2L]]
  gd <- array(NA_real_, c(n_draws, nch, 2L * p),
    dimnames = list(NULL, NULL, c(paste0("mu:", pn), paste0("logsd:", pn))))
  ad <- array(NA_real_, c(n_draws, nch, p * nsub),
    dimnames = list(NULL, NULL,
      unlist(lapply(ids, function(id) paste0(id, ":", pn)))))
  for (j in seq_len(nch)) {
    ch <- fit[[j]]
    idx <- tail(which(ch$samples$stage[seq_len(ch$samples$idx)] == "sample"), n_draws)
    stopifnot(length(idx) == n_draws)
    gd[, j, seq_len(p)] <- t(ch$samples$theta_mu[, idx, drop = FALSE])
    lsd <- vapply(idx, function(i) 0.5 * log(diag(ch$samples$theta_var[, , i])),
                  numeric(p))
    gd[, j, p + seq_len(p)] <- t(lsd)
    ad[, j, ] <- t(matrix(ch$samples$alpha[, , idx], nrow = p * nsub))
  }
  group <- posterior::summarise_draws(posterior::as_draws_array(gd),
    posterior::ess_bulk, posterior::rhat)
  subject <- posterior::summarise_draws(posterior::as_draws_array(ad),
    posterior::ess_bulk, posterior::rhat)
  ess_col <- grep("ess_bulk", names(group), value = TRUE)[1L]
  rhat_col <- grep("rhat", names(group), value = TRUE)[1L]
  ae <- grep("ess_bulk", names(subject), value = TRUE)[1L]
  ar <- grep("rhat", names(subject), value = TRUE)[1L]
  data.frame(
    arm = label, sample_seconds = seconds,
    rounds = paste(vapply(fit, function(ch) ch$collapsed$plan$rounds, 0.0),
                   collapse = ","),
    acceptance = paste(vapply(fit, function(ch) {
      b <- ch$collapsed$plan$block
      if (!b$n) NA_real_ else b$acc / b$n
    }, 0.0), collapse = ","),
    group_mu_ess_sec = median(group[[ess_col]][startsWith(group$variable, "mu:")]) / seconds,
    group_logsd_ess_sec = median(group[[ess_col]][startsWith(group$variable, "logsd:")]) / seconds,
    t0_group_mu_ess_sec = median(group[[ess_col]][startsWith(group$variable, "mu:t0")]) / seconds,
    all_subject_ess_sec = median(subject[[ae]]) / seconds,
    subject_52091_ess_sec = median(subject[[ae]][startsWith(subject$variable,
      "52091:")]) / seconds,
    max_group_rhat = max(group[[rhat_col]], na.rm = TRUE),
    max_subject_rhat = max(subject[[ar]], na.rm = TRUE),
    stringsAsFactors = FALSE)
}

results <- list()
order <- c("screen8", "screen16")
for (arm in order) {
  start <- EMC2:::restore_custom_kernel_pointers(
    EMC2:::restore_duplicates(readRDS(file.path(fit_dir, paste0(arm, ".rds")))),
    quiet = TRUE)
  tic <- proc.time()[["elapsed"]]
  fit <- EMC2::fit(start, stage = "sample",
    stop_criteria = list(sample = list(iter = n_draws)), step_size = n_draws,
    cores_for_chains = 3L, cores_per_chain = 4L, particle_factor = 10,
    max_tries = 1L, trim = FALSE, verbose = FALSE, verboseProgress = FALSE)
  seconds <- proc.time()[["elapsed"]] - tic
  results[[arm]] <- summarise_fit(fit, arm, seconds)
  saveRDS(fit, file.path(fit_dir, paste0(arm, "_fit.rds")))
  cat(arm, "seconds", round(seconds, 2), "rounds",
      results[[arm]]$rounds, "acceptance", results[[arm]]$acceptance, "\n")
  rm(fit, start)
  gc()
}
summary <- do.call(rbind, results)
write.csv(summary, file.path(artifact_dir,
  "validation_batch_ess_comparison.csv"), row.names = FALSE)
print(summary, row.names = FALSE)
