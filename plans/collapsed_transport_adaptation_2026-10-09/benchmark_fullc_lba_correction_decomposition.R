pkgload::load_all(".", quiet = TRUE, compile = FALSE)

archive <- Sys.getenv("EMC2_FULLC_RDS",
  "/data/work/EMC2_dev_oo/.Rtmp/claude-1000/-data-work-EMC2-dev-oo/3baace33-74ca-4333-ba0b-800d175b2295/scratchpad/runs/final/T30_fullc_1_emc.rds")
lba_start <- Sys.getenv("EMC2_LBA_START",
  "/tmp/collapsed_transport_followup/safe_fit_ess_start.rds")
artifact_dir <- "plans/collapsed_transport_adaptation_2026-10-09"
n_valid <- as.integer(Sys.getenv("EMC2_DECOMP_VALID", "24"))
max_attempts <- as.integer(Sys.getenv("EMC2_DECOMP_MAX", "1200"))

make_state <- function(chain) {
  i <- chain$samples$idx
  ctx <- EMC2:::.emc_collapsed_context(chain)
  cur <- EMC2:::.emc_shear_state(chain$samples$theta_mu[, i],
    as.matrix(chain$samples$theta_var[, , i]), ctx$fits,
    chain$prior, ctx$spec)
  stopifnot(!is.null(cur), is.finite(cur$lp))
  x <- list(cur = cur,
    alpha = matrix(chain$samples$alpha[, , i], chain$n_pars),
    ll = chain$samples$subj_ll[, i], rounds = 0L,
    wpool = NULL, wpool_part = NULL)
  list(chain = chain, ctx = ctx, x = x)
}

check_t0 <- function(rows, chain) {
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

run_target <- function(label, state, plan, seed) {
  chain <- state$chain; ctx <- state$ctx; x0 <- state$x
  n <- chain$n_subjects; p <- chain$n_pars
  subject_ids <- dimnames(chain$samples$alpha)[[2L]]
  if (is.null(subject_ids)) subject_ids <- as.character(seq_len(n))
  plan$scale <- 1
  plan$t0_ratio <- 1
  plan$t0_mean_ratio <- 1
  plan$t0_hyper_ratio <- 1
  plan$mean_idx <- integer()
  plan$mean_ratio <- 1

  if (identical(label, "LBA")) {
    # Match the original-scale counterfactual: group-chart t0 rows fixed,
    # while all non-t0 chart rows retain the full metric.
    plan$C[plan$t0_idx, ] <- 0
    phi0 <- EMC2:::.emc_collapsed_to_phi(x0$cur$mu, x0$cur$V, ctx$chart)
    C <- EMC2:::.emc_collapsed_preconditioner(plan$C, 1, 1,
      plan$t0_idx, 1, plan$t0_mean_idx, integer(), 1)
    plan$eps <- EMC2:::.emc_collapsed_tune_eps(phi0, C, plan$eps, ctx,
      if (isTRUE(is.finite(plan$length))) plan$length else
        EMC2:::.EMC_COLLAPSED_LENGTH)
  }

  records <- subject_records <- list(); capture <- NULL
  ns <- asNamespace("EMC2"); binding <- ".emc_wpool_subject_ll"
  old_ll <- get(binding, envir = ns, inherits = FALSE)
  unlockBinding(binding, ns)
  assign(binding, function(pool, rows, ctx, part = NULL, source) {
    if (identical(source, "collapsed")) {
      b <- check_t0(rows, chain)
      capture <<- list(rows = rows, bounds = b, valid = !any(b > 0))
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

  set.seed(seed)
  started <- proc.time()[["elapsed"]]
  attempts <- valid_count <- 0L
  while (attempts < max_attempts && valid_count < n_valid) {
    attempts <- attempts + 1L; capture <- NULL
    x <- EMC2:::.emc_collapsed_round(x0, plan, chain, ctx)
    if (is.null(capture) || !isTRUE(capture$valid)) next
    valid_count <- valid_count + 1L
    alpha_new <- matrix(vapply(seq_len(n), function(s)
      as.numeric(capture$rows[[s]][1L, chain$par_names]), numeric(p)), p)
    lhat_new <- EMC2:::.emc_collapsed_log_lhat(alpha_new, ctx$fits)
    lhat_old <- EMC2:::.emc_collapsed_log_lhat(x0$alpha, ctx$fits)
    delta <- (capture$ll - lhat_new) - (x0$ll - lhat_old)
    log_correction <- sum(delta)
    records[[valid_count]] <- data.frame(target = label, proposal = valid_count,
      attempts_to_find = attempts, subjects = n,
      log_correction = log_correction,
      expected_acceptance = if (is.finite(log_correction))
        min(1, exp(log_correction)) else 0,
      mean_subject_delta = mean(delta), sd_subject_delta = stats::sd(delta),
      median_abs_subject_delta = median(abs(delta)),
      max_abs_subject_delta = max(abs(delta)),
      lower_bound_violations = capture$bounds["lower"],
      rt_bound_violations = capture$bounds["upper"])
    subject_records[[valid_count]] <- data.frame(target = label,
      proposal = valid_count, subject_index = seq_len(n),
      subject_id = subject_ids, correction = delta)
  }
  elapsed <- proc.time()[["elapsed"]] - started
  if (valid_count < n_valid) stop(label, " obtained only ", valid_count,
    " t0-valid proposals in ", attempts, " attempts")
  out <- do.call(rbind, records)
  attr(out, "subject_corrections") <- do.call(rbind, subject_records)
  attr(out, "setup") <- data.frame(target = label, subjects = n,
    group_parameters = p, chart_dim = plan$dim, scale = 1, eps = plan$eps,
    valid_proposals = valid_count, attempts = attempts, elapsed_seconds = elapsed)
  out
}

fullc <- EMC2:::restore_custom_kernel_pointers(
  EMC2:::restore_duplicates(readRDS(archive)), quiet = TRUE)[[1L]]
lba_prepared <- readRDS(lba_start)
lba <- EMC2:::restore_custom_kernel_pointers(
  EMC2:::restore_duplicates(lba_prepared$start), quiet = TRUE)[[2L]]

fullc_state <- make_state(fullc)
lba_state <- make_state(lba)
fullc_plan <- fullc$collapsed$plan
lba_plan <- lba$collapsed$plan
fullc_rows <- run_target("T30", fullc_state, fullc_plan, 20261021L)
lba_rows <- run_target("LBA", lba_state, lba_plan, 20261022L)
rows <- rbind(fullc_rows, lba_rows)
setup <- rbind(attr(fullc_rows, "setup"), attr(lba_rows, "setup"))
subject_corrections <- rbind(attr(fullc_rows, "subject_corrections"),
                             attr(lba_rows, "subject_corrections"))
write.csv(rows, file.path(artifact_dir,
  "fullc_lba_correction_decomposition.csv"), row.names = FALSE)
write.csv(subject_corrections, file.path(artifact_dir,
  "fullc_lba_subject_corrections.csv"), row.names = FALSE)
write.csv(setup, file.path(artifact_dir,
  "fullc_lba_correction_decomposition_setup.csv"), row.names = FALSE)
summary <- do.call(rbind, lapply(split(rows, rows$target), function(x)
  data.frame(target = x$target[1L], proposals = nrow(x),
    median_expected_acceptance = median(x$expected_acceptance),
    mean_expected_acceptance = mean(x$expected_acceptance),
    median_log_correction = median(x$log_correction),
    median_abs_subject_delta = median(x$median_abs_subject_delta),
    median_sd_subject_delta = median(x$sd_subject_delta))))
print(setup, row.names = FALSE)
print(summary, row.names = FALSE)
