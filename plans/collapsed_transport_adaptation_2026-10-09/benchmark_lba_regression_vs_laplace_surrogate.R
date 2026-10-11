pkgload::load_all(".", quiet = TRUE, compile = FALSE)

artifact_dir <- "plans/collapsed_transport_adaptation_2026-10-09"
target <- Sys.getenv("EMC2_SURROGATE_TARGET", "LBA")
if (identical(target, "LBA")) {
  prepared <- readRDS("/tmp/collapsed_transport_followup/safe_fit_ess_start.rds")
  chain <- EMC2:::restore_custom_kernel_pointers(
    EMC2:::restore_duplicates(prepared$start), quiet = TRUE)[[2L]]
} else if (identical(target, "T30")) {
  archive <- Sys.getenv("EMC2_FULLC_RDS",
    ".Rtmp/claude-1000/-data-work-EMC2-dev-oo/3baace33-74ca-4333-ba0b-800d175b2295/scratchpad/runs/final/T30_fullc_1_emc.rds")
  chain <- EMC2:::restore_custom_kernel_pointers(
    EMC2:::restore_duplicates(readRDS(archive)), quiet = TRUE)[[1L]]
} else stop("EMC2_SURROGATE_TARGET must be LBA or T30")
prefix <- tolower(target)
n_valid <- as.integer(Sys.getenv("EMC2_SURROGATE_VALID", "8"))
max_attempts <- as.integer(Sys.getenv("EMC2_SURROGATE_MAX", "500"))

idx <- chain$samples$idx
p <- chain$n_pars
n <- chain$n_subjects
alpha0 <- matrix(chain$samples$alpha[, , idx], p)
ll0 <- as.numeric(chain$samples$subj_ll[, idx])
mu0 <- chain$samples$theta_mu[, idx]
V0 <- as.matrix(chain$samples$theta_var[, , idx])
base_plan <- chain$collapsed$plan
chart <- EMC2:::.emc_collapsed_chart(EMC2:::.emc_shear_v_spec(chain), p)
phi0 <- EMC2:::.emc_collapsed_to_phi(mu0, V0, chart)
base_C <- base_plan$C
t0_idx <- EMC2:::.emc_collapsed_t0_chart_idx(chain$par_names, chart)
t0_mean_idx <- EMC2:::.emc_collapsed_t0_mean_chart_idx(chain$par_names)

# Freeze the t0 chart outputs and give every other row its unshrunk learned
# metric. The exact target still enforces all participant-specific RT bounds.
plan <- base_plan
if (identical(target, "LBA")) plan$C[t0_idx, ] <- 0
plan$scale <- 1
plan$t0_ratio <- 1
plan$t0_mean_ratio <- 1
plan$t0_hyper_ratio <- 1
plan$mean_idx <- integer()
plan$mean_ratio <- 1
plan$length <- if (isTRUE(is.finite(base_plan$length)) &&
    base_plan$length > 0) base_plan$length else EMC2:::.EMC_COLLAPSED_LENGTH
plan$eps <- NA_real_

make_context <- function(base_fits) {
  base_fits <- lapply(base_fits, EMC2:::.emc_collapsed_safe_fit)
  pooled_value <- if (is.null(chain$pooled)) NULL else chain$pooled$value
  fits <- lapply(EMC2:::.emc_collapsed_condition(base_fits,
    pooled_value), EMC2:::.emc_collapsed_safe_fit)
  spec <- EMC2:::.emc_shear_v_spec(chain)
  list(fits = fits, base_fits = base_fits, spec = spec,
    prior = chain$prior, chart = chart,
    cpp = EMC2:::.emc_collapsed_cpp_model(chart, fits, chain$prior, spec))
}

regression_ctx <- EMC2:::.emc_collapsed_context(chain)
laplace_ctx <- make_context(EMC2:::.emc_shear_fits(chain))
contexts <- list(regression = regression_ctx, laplace = laplace_ctx)

check_t0 <- function(rows) {
  lower <- upper <- logical(length(rows))
  for (s in seq_along(rows)) {
    pars <- EMC2:::get_pars_oo(rows[[s]], chain$data[[s]], chain$model)
    t0 <- pars[, "t0"]
    rt <- chain$data[[s]]$rt
    stopifnot(length(t0) == length(rt))
    lower[s] <- any(is.finite(t0) & t0 < 0.05)
    upper[s] <- any(is.finite(t0) & is.finite(rt) & t0 >= rt)
  }
  c(lower = sum(lower), upper = sum(upper))
}

run_arm <- function(label, ctx, seed) {
  cur <- EMC2:::.emc_shear_state(mu0, V0, ctx$fits, chain$prior, ctx$spec)
  stopifnot(!is.null(cur), is.finite(cur$lp))
  x0 <- list(cur = cur, alpha = alpha0, ll = ll0, rounds = 0L,
    wpool = NULL, wpool_part = NULL)
  C_candidate <- EMC2:::.emc_collapsed_preconditioner(plan$C, plan$scale,
    plan$t0_ratio, t0_idx, plan$t0_hyper_ratio, t0_mean_idx,
    plan$mean_idx, plan$mean_ratio)
  plan$eps <- EMC2:::.emc_collapsed_tune_eps(phi0, C_candidate,
    base_plan$eps, ctx, plan$length)
  rows_out <- vector("list", n_valid)
  subject_out <- vector("list", n_valid)
  attempts <- valid <- lower_invalid <- upper_invalid <- 0L
  capture <- NULL
  ns <- asNamespace("EMC2")
  binding <- ".emc_wpool_subject_ll"
  old_ll <- get(binding, envir = ns, inherits = FALSE)
  unlockBinding(binding, ns)
  assign(binding, function(pool, rows, ctx, part = NULL, source) {
    if (identical(source, "collapsed")) {
      bounds <- check_t0(rows)
      capture <<- list(rows = rows, bounds = bounds,
        valid = !any(bounds > 0))
      if (!capture$valid)
        return(list(ll = rep(-Inf, length(rows)), pool = pool))
      result <- old_ll(pool, rows, ctx, part, source)
      capture$ll <<- as.numeric(result$ll)
      return(result)
    }
    old_ll(pool, rows, ctx, part, source)
  }, envir = ns)
  lockBinding(binding, ns)
  on.exit({
    unlockBinding(binding, ns)
    assign(binding, old_ll, envir = ns)
    lockBinding(binding, ns)
  }, add = TRUE)
  started <- proc.time()[["elapsed"]]

  set.seed(seed)
  while (attempts < max_attempts && valid < n_valid) {
    attempts <- attempts + 1L
    capture <- NULL
    x <- EMC2:::.emc_collapsed_round(x0, plan, chain, ctx)
    if (is.null(capture)) next
    bounds <- capture$bounds
    if (any(bounds > 0)) {
      lower_invalid <- lower_invalid + as.integer(bounds["lower"] > 0)
      upper_invalid <- upper_invalid + as.integer(bounds["upper"] > 0)
      next
    }
    rows <- capture$rows
    ll_new <- capture$ll
    alpha_new <- matrix(vapply(seq_len(n), function(s)
      as.numeric(rows[[s]][1L, chain$par_names]), numeric(p)), p)
    delta <- (ll_new - vapply(seq_len(n), function(s)
      EMC2:::.emc_collapsed_log_lhat(alpha_new[, s, drop = FALSE],
        ctx$fits[s]), numeric(1))) -
      (ll0 - EMC2:::.emc_collapsed_log_lhat(alpha0, ctx$fits))
    log_correction <- sum(delta)
    acc <- x$acc
    metric_jump2 <- sum(forwardsolve(base_C, x$phi_jump)^2)
    mu_jump <- x$jump[seq_len(p)]
    lsd_jump <- x$jump[p + seq_len(p)]
    valid <- valid + 1L
    row <- data.frame(variant = label, proposal = valid,
      attempts = attempts, exact_acceptance = acc,
      log_correction = log_correction,
      metric_esjd = acc * metric_jump2,
      lower_bound_violations = bounds["lower"],
      rt_bound_violations = bounds["upper"])
    names_mu <- paste0("mu_esjd_", make.names(chain$par_names))
    names_sd <- paste0("logsd_esjd_", make.names(chain$par_names))
    row[names_mu] <- as.list(acc * mu_jump^2)
    row[names_sd] <- as.list(acc * lsd_jump^2)
    rows_out[[valid]] <- row
    subject_out[[valid]] <- data.frame(variant = label,
      proposal = valid, subject = seq_len(n), correction = delta)
  }
  elapsed <- proc.time()[["elapsed"]] - started
  detail <- if (valid) do.call(rbind, rows_out[seq_len(valid)]) else
    data.frame()
  mean_cols <- paste0("mu_esjd_", make.names(chain$par_names))
  logsd_cols <- paste0("logsd_esjd_", make.names(chain$par_names))
  mean_esjd <- if (valid) sum(colMeans(detail[, mean_cols, drop = FALSE])) else NA_real_
  logsd_esjd <- if (valid) sum(colMeans(detail[, logsd_cols, drop = FALSE])) else NA_real_
  summary <- data.frame(variant = label, n_subjects = n,
    chart_dim = chart$dim, scale = 1, eps = plan$eps,
    proposals = valid, total_attempts = attempts,
    fraction_invalid_by_t0 = (lower_invalid + upper_invalid) / attempts,
    lower_bound_invalid_attempts = lower_invalid,
    rt_bound_invalid_attempts = upper_invalid,
    mean_expected_acceptance = if (valid) mean(detail$exact_acceptance) else NA_real_,
    median_log_correction = if (valid) median(detail$log_correction) else NA_real_,
    mean_metric_esjd = if (valid) mean(detail$metric_esjd) else NA_real_,
    metric_esjd_per_second = if (valid)
      mean(detail$metric_esjd) * valid / elapsed else NA_real_,
    mean_group_mean_esjd = mean_esjd,
    group_mean_esjd_per_second = mean_esjd * valid / elapsed,
    mean_logsd_esjd = logsd_esjd,
    logsd_esjd_per_second = logsd_esjd * valid / elapsed,
    seconds = elapsed, exact_likelihood_proposals = valid)
  list(summary = summary, detail = detail,
    subject = if (valid) do.call(rbind, subject_out[seq_len(valid)]) else
      data.frame(variant = character(), proposal = integer(),
        subject = integer(), correction = numeric()))
}

results <- list(
  run_arm("regression", regression_ctx, 20261031L),
  run_arm("laplace", laplace_ctx, 20261032L))
summary <- do.call(rbind, lapply(results, `[[`, "summary"))
detail <- do.call(rbind, lapply(results, `[[`, "detail"))
subject <- do.call(rbind, lapply(results, `[[`, "subject"))
write.csv(summary, file.path(artifact_dir,
  paste0(prefix, "_regression_vs_laplace_summary.csv")), row.names = FALSE)
write.csv(detail, file.path(artifact_dir,
  paste0(prefix, "_regression_vs_laplace_proposals.csv")), row.names = FALSE)
write.csv(subject, file.path(artifact_dir,
  paste0(prefix, "_regression_vs_laplace_subject_corrections.csv")),
  row.names = FALSE)
print(summary, row.names = FALSE)
