pkgload::load_all(".", quiet = TRUE, compile = FALSE)
artifact_dir <- "plans/collapsed_transport_adaptation_2026-10-09"
fit_dir <- "/tmp/collapsed_transport_followup"
repeat_id <- as.integer(Sys.getenv("EMC2_CALIBRATION_REPEAT", "1"))
stopifnot(length(repeat_id) == 1L, is.finite(repeat_id), repeat_id >= 1L)
base <- readRDS(file.path(fit_dir, "mean_scaled_calibrated.rds"))
base <- EMC2:::restore_custom_kernel_pointers(
  structure(list(base), class = "emc"), quiet = TRUE)[[1L]]
if (repeat_id > 1L) {
  advance <- function(seed) if (is.null(seed)) NULL else parallel::nextRNGStream(seed)
  for (i in seq_len(repeat_id - 1L)) {
    base$rng$chain <- advance(base$rng$chain)
    base$rng$gibbs <- advance(base$rng$gibbs)
    base$rng$subjects <- lapply(base$rng$subjects, advance)
  }
}
base_rng <- base$rng$gibbs
options(emc2.subject_hmc = FALSE)

summarise_one <- function(fit, label, calibration_seconds, sample_seconds) {
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
  subject <- match("52091", dimnames(chain$samples$alpha)[[2L]])
  subject_draws <- t(matrix(chain$samples$alpha[, subject, idx],
    nrow = chain$n_pars, ncol = length(idx)))
  group_ess <- as_ess(group_draws, pnames)
  subject_ess <- as_ess(subject_draws, pnames)
  logsd_draws <- matrix(NA_real_, length(idx), chain$n_pars,
    dimnames = list(NULL, pnames))
  for (k in seq_along(idx))
    logsd_draws[k, ] <- log(sqrt(diag(chain$samples$theta_var[, , idx[k]])))
  logsd_ess <- as_ess(logsd_draws, paste0("logsd:", pnames))
  plan <- chain$collapsed$plan
  acceptance <- if (plan$block$n > 0L) plan$block$acc / plan$block$n else NA_real_
  summary <- data.frame(
    arm = label,
    calibration_seconds = calibration_seconds,
    sample_seconds = sample_seconds,
    total_seconds = calibration_seconds + sample_seconds,
    collapsed_acceptance = acceptance,
    median_group_mu_ess = median(group_ess),
    median_group_mu_ess_per_second = median(group_ess) /
      (calibration_seconds + sample_seconds),
    median_group_logsd_ess = median(logsd_ess),
    median_group_logsd_ess_per_second = median(logsd_ess) /
      (calibration_seconds + sample_seconds),
    median_subject_ess = median(subject_ess),
    median_subject_ess_per_second = median(subject_ess) /
      (calibration_seconds + sample_seconds),
    v_g_ess = group_ess[match("v.g", pnames)],
    v_g_ess_per_second = group_ess[match("v.g", pnames)] /
      (calibration_seconds + sample_seconds),
    slow_mean = pnames[plan$mean_idx],
    slow_mean_ess = group_ess[plan$mean_idx],
    slow_mean_ess_per_second = group_ess[plan$mean_idx] /
      (calibration_seconds + sample_seconds),
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
      ess = as.numeric(group_ess), ess_per_total_second =
        as.numeric(group_ess) / summary$total_seconds),
    data.frame(arm = label, family = "group_logsd",
      variable = paste0("logsd:", pnames), ess = as.numeric(logsd_ess),
      ess_per_total_second = as.numeric(logsd_ess) / summary$total_seconds),
    data.frame(arm = label, family = "subject_52091", variable = pnames,
      ess = as.numeric(subject_ess), ess_per_total_second =
        as.numeric(subject_ess) / summary$total_seconds)
  )
  list(summary = summary, per_variable = per_variable)
}

run_one <- function(label, screen) {
  chain <- unserialize(serialize(base, NULL))
  chain <- EMC2:::restore_custom_kernel_pointers(
    structure(list(chain), class = "emc"), quiet = TRUE)[[1L]]
  chain$collapsed$plan$calibrate <- TRUE
  chain$collapsed$plan$calibration_screen <- screen
  started <- proc.time()[["elapsed"]]
  chain <- EMC2:::.emc_collapsed_calibrate(chain, n_cores = 1L)
  calibration_seconds <- proc.time()[["elapsed"]] - started
  chain$rng$gibbs <- base_rng # common production stream; pilots are discarded
  plan <- chain$collapsed$plan
  plan$rounds <- 4
  plan$calibrate <- FALSE
  plan$calibration_screen <- FALSE
  plan$last <- EMC2:::.emc_collapsed_record(NULL, chain$n_pars)
  plan$pool <- EMC2:::.emc_collapsed_record(NULL, chain$n_pars)
  plan$block <- EMC2:::.emc_collapsed_record(NULL, chain$n_pars)
  chain$collapsed$plan <- plan
  chain$sample_proposals_next <- Inf
  chain$production_options$emc2.subject_hmc <- FALSE
  started <- proc.time()[["elapsed"]]
  fit <- EMC2::fit(structure(list(chain), class = "emc"), stage = "sample",
    stop_criteria = list(iter = 200L), step_size = 100L,
    cores_for_chains = 1L, cores_per_chain = 4L, particle_factor = 10,
    max_tries = 1L, trim = FALSE, verbose = FALSE,
    verboseProgress = FALSE)
  sample_seconds <- proc.time()[["elapsed"]] - started
  saveRDS(fit, file.path(fit_dir,
    paste0("rebuild_calibration_", label, "_ess_rep", repeat_id, ".fit.rds")))
  result <- summarise_one(fit, label, calibration_seconds, sample_seconds)
  list(summary = result$summary, per_variable = result$per_variable,
    calibration_seconds = calibration_seconds, sample_seconds = sample_seconds,
    calibration = fit[[1L]]$collapsed$plan$calibration)
}

arm_order <- if (repeat_id %% 2L == 1L) c("full", "screen") else c("screen", "full")
results <- list()
for (label in arm_order)
  results[[label]] <- run_one(label, screen = identical(label, "screen"))
full <- results$full
screen <- results$screen
summary <- rbind(full$summary, screen$summary)
per_variable <- rbind(full$per_variable, screen$per_variable)
write.csv(summary, file.path(artifact_dir,
  paste0("rebuild_calibration_ess_repeat", repeat_id, "_summary.csv")),
  row.names = FALSE)
write.csv(per_variable, file.path(artifact_dir,
  paste0("rebuild_calibration_ess_repeat", repeat_id, "_per_variable.csv")),
  row.names = FALSE)
saveRDS(results, file.path(artifact_dir,
  paste0("rebuild_calibration_ess_repeat", repeat_id, ".rds")))
print(summary, row.names = FALSE)
