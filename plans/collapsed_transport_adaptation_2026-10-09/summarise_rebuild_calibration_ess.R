pkgload::load_all(".", quiet = TRUE, compile = FALSE)
artifact_dir <- "plans/collapsed_transport_adaptation_2026-10-09"
fit_dir <- "/tmp/collapsed_transport_followup"
repeat_id <- as.integer(Sys.getenv("EMC2_CALIBRATION_REPEAT", "1"))
stopifnot(length(repeat_id) == 1L, is.finite(repeat_id), repeat_id >= 1L)
timing_path <- if (repeat_id == 1L) "rebuild_calibration_ess.rds" else
  paste0("rebuild_calibration_ess_repeat", repeat_id, ".rds")
timing <- readRDS(file.path(artifact_dir, timing_path))

summarise_one <- function(label) {
  fit <- readRDS(file.path(fit_dir,
    if (repeat_id == 1L)
      paste0("rebuild_calibration_", label, "_ess.fit.rds") else
        paste0("rebuild_calibration_", label, "_ess_rep", repeat_id, ".fit.rds")))
  chain <- fit[[1L]]
  pnames <- dimnames(chain$samples$theta_mu)[[1L]]
  idx <- tail(which(chain$samples$stage[seq_len(chain$samples$idx)] == "sample"),
              200L)
  stopifnot(length(idx) == 200L)
  as_ess <- function(x, names) {
    a <- array(NA_real_, dim = c(nrow(x), 1L, ncol(x)),
      dimnames = list(NULL, NULL, names))
    for (j in seq_len(ncol(x))) a[, 1L, j] <- x[, j]
    posterior::summarise_draws(posterior::as_draws_array(a),
      posterior::ess_bulk)[[2L]]
  }
  group_draws <- t(chain$samples$theta_mu[, idx, drop = FALSE])
  s <- match("52091", dimnames(chain$samples$alpha)[[2L]])
  stopifnot(!is.na(s))
  subject_draws <- t(matrix(chain$samples$alpha[, s, idx],
    nrow = chain$n_pars, ncol = length(idx)))
  group_ess <- as.numeric(as_ess(group_draws, pnames))
  subject_ess <- as.numeric(as_ess(subject_draws, pnames))
  subject_names <- dimnames(chain$samples$alpha)[[2L]]
  all_subject_draws <- matrix(NA_real_, length(idx),
    chain$n_pars * chain$n_subjects)
  all_subject_names <- character(ncol(all_subject_draws))
  for (subject_idx in seq_len(chain$n_subjects)) {
    cols <- (subject_idx - 1L) * chain$n_pars + seq_len(chain$n_pars)
    all_subject_draws[, cols] <- t(matrix(
      chain$samples$alpha[, subject_idx, idx],
      nrow = chain$n_pars, ncol = length(idx)))
    all_subject_names[cols] <- paste(subject_names[subject_idx], pnames,
                                      sep = ":")
  }
  all_subject_ess <- as.numeric(as_ess(all_subject_draws, all_subject_names))
  logsd_draws <- matrix(NA_real_, length(idx), chain$n_pars,
    dimnames = list(NULL, pnames))
  for (k in seq_along(idx))
    logsd_draws[k, ] <- log(sqrt(diag(chain$samples$theta_var[, , idx[k]])))
  logsd_ess <- as.numeric(as_ess(logsd_draws, paste0("logsd:", pnames)))
  plan <- chain$collapsed$plan
  calibration_seconds <- timing[[label]]$calibration_seconds
  sample_seconds <- timing[[label]]$sample_seconds
  total_seconds <- calibration_seconds + sample_seconds
  acceptance <- if (plan$block$n > 0L) plan$block$acc / plan$block$n else NA_real_
  summary <- data.frame(
    arm = label,
    calibration_seconds = calibration_seconds,
    sample_seconds = sample_seconds,
    total_seconds = total_seconds,
    collapsed_acceptance = acceptance,
    median_group_mu_ess = median(group_ess),
    median_group_mu_ess_per_second = median(group_ess) / total_seconds,
    median_group_logsd_ess = median(logsd_ess),
    median_group_logsd_ess_per_second = median(logsd_ess) / total_seconds,
    median_subject_ess = median(subject_ess),
    median_subject_ess_per_second = median(subject_ess) / total_seconds,
    median_all_subject_ess = median(all_subject_ess),
    median_all_subject_ess_per_second = median(all_subject_ess) / total_seconds,
    mean_all_subject_ess = mean(all_subject_ess),
    fraction_all_subject_ess_below_10 = mean(all_subject_ess < 10),
    v_g_ess = group_ess[match("v.g", pnames)],
    v_g_ess_per_second = group_ess[match("v.g", pnames)] / total_seconds,
    slow_mean = pnames[plan$mean_idx],
    slow_mean_ess = group_ess[plan$mean_idx],
    slow_mean_ess_per_second = group_ess[plan$mean_idx] / total_seconds,
    scale = plan$scale,
    t0_mean_scale = plan$scale * plan$t0_ratio,
    t0_hyper_scale = plan$scale * plan$t0_hyper_ratio,
    mean_scale = if (length(plan$mean_idx))
      plan$scale * plan$mean_ratio else NA_real_,
    path_length = plan$length,
    eps = plan$eps,
    calibration_pilots = plan$calibration$pilots
  )
  per_variable <- rbind(
    data.frame(arm = label, family = "group_mu", variable = pnames,
      ess = group_ess, ess_per_total_second = group_ess / total_seconds),
    data.frame(arm = label, family = "group_logsd",
      variable = paste0("logsd:", pnames), ess = logsd_ess,
      ess_per_total_second = logsd_ess / total_seconds),
    data.frame(arm = label, family = "subject_52091", variable = pnames,
      ess = subject_ess, ess_per_total_second = subject_ess / total_seconds),
    data.frame(arm = label, family = "all_subject", variable = all_subject_names,
      ess = all_subject_ess,
      ess_per_total_second = all_subject_ess / total_seconds)
  )
  list(summary = summary, per_variable = per_variable,
       calibration_seconds = calibration_seconds,
       sample_seconds = sample_seconds,
       calibration = timing[[label]]$calibration)
}

full <- summarise_one("full")
screen <- summarise_one("screen")
summary <- rbind(full$summary, screen$summary)
per_variable <- rbind(full$per_variable, screen$per_variable)
summary_stub <- if (repeat_id == 1L) "rebuild_calibration_ess" else
  paste0("rebuild_calibration_ess_repeat", repeat_id)
write.csv(summary, file.path(artifact_dir,
  paste0(summary_stub, "_summary.csv")), row.names = FALSE)
write.csv(per_variable, file.path(artifact_dir,
  paste0(summary_stub, "_per_variable.csv")), row.names = FALSE)
saveRDS(list(full = full, screen = screen), file.path(artifact_dir,
  paste0(summary_stub, ".rds")))
print(summary, row.names = FALSE)
