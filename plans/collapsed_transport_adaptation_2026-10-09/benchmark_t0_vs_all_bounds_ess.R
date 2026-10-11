pkgload::load_all(".", quiet = TRUE)
production_mode <- identical(Sys.getenv("EMC2_ESS_PRODUCTION", "FALSE"), "TRUE")
# Controlled mode fixes one collapsed round and disables subject HMC. Production
# mode rebuilds the allocator from the archived sample history with HMC enabled.
options(emc2.subject_hmc = production_mode, emc2.collapsed_move = TRUE,
  emc2.sampler_profile = FALSE)

archive_root <- Sys.getenv(
  "EMC2_FULLC_ROOT",
  ".Rtmp/claude-1000/-data-work-EMC2-dev-oo/3baace33-74ca-4333-ba0b-800d175b2295/scratchpad/runs/final")
case <- Sys.getenv("EMC2_ESS_CASE", "T30")
case_file <- switch(case, T30 = "T30", corr15 = "corr15",
  Forstmann = "forstmann", stop("unknown case"))
default_replicate <- if (case == "Forstmann") 2L else 1L
replicate <- as.integer(Sys.getenv("EMC2_FULLC_REPLICATE", default_replicate))
file <- sprintf("%s_fullc_%d_emc.rds", case_file, replicate)
n_draws <- as.integer(Sys.getenv("EMC2_ESS_DRAWS", "200"))
start <- EMC2:::restore_custom_kernel_pointers(
  EMC2:::restore_duplicates(readRDS(file.path(archive_root, file))), quiet = TRUE)
stopifnot(length(start) >= 2L)
pn <- start[[1L]]$par_names
ptype <- EMC2:::get_p_types(pn)
t0_idx <- which(ptype == "t0")
stopifnot(length(t0_idx) > 0L)

# Reconstruct the historical production selector for a fair comparison after
# the package default changes. This belongs to the audit benchmark only.
legacy_all_hard_fixed_idx <- function(sampler) {
  model <- EMC2:::.oo_model_list(sampler$model)
  bounds <- model$bound$minmax
  transform <- model$transform
  constrained <- character()
  if (is.matrix(bounds) && nrow(bounds) == 2L && !is.null(colnames(bounds))) {
    pars <- colnames(bounds)
    func <- stats::setNames(rep("identity", length(pars)), pars)
    if (length(transform$func)) {
      use <- intersect(pars, names(transform$func))
      func[use] <- transform$func[use]
    }
    lower <- stats::setNames(rep(-Inf, length(pars)), pars)
    upper <- stats::setNames(rep(Inf, length(pars)), pars)
    lower[func %in% c("exp", "pnorm")] <- 0
    upper[func == "pnorm"] <- 1
    if (length(transform$lower)) {
      use <- intersect(pars, names(transform$lower))
      lower[use] <- transform$lower[use]
    }
    if (length(transform$upper)) {
      use <- intersect(pars, names(transform$upper))
      upper[use] <- transform$upper[use]
    }
    constrained <- pars[bounds[1L, pars] > lower | bounds[2L, pars] < upper]
  }
  has_t0 <- "t0" %in% names(model$p_types) ||
    "t0" %in% names(transform$func)
  fixed_pars <- unique(c(constrained, if (has_t0) "t0"))
  if (!length(fixed_pars)) return(integer())
  which(EMC2:::get_p_types(sampler$par_names) %in% fixed_pars)
}

prepare <- function(x) {
  out <- unserialize(serialize(x, NULL))
  for (j in seq_along(out)) {
    ch <- out[[j]]
    if (!production_mode) ch$collapsed$plan$rounds <- 1
    ch$collapsed$plan$block <- EMC2:::.emc_collapsed_record(NULL, ch$n_pars)
    ch$collapsed$plan$last <- EMC2:::.emc_collapsed_record(NULL, ch$n_pars)
    ch$collapsed$plan$pool <- EMC2:::.emc_collapsed_record(NULL, ch$n_pars)
    ch$sample_proposals_next <- if (production_mode) 0L else Inf
    ch$production_options$emc2.subject_hmc <- production_mode
    ch$production_options$emc2.collapsed_move <- TRUE
    out[[j]] <- ch
  }
  out
}

variable_outputs <- list()

run_arm <- function(label) {
  arm <- prepare(start)
  ns <- asNamespace("EMC2")
  binding <- ".emc_collapsed_fixed_idx"
  old <- get(binding, envir = ns, inherits = FALSE)
  unlockBinding(binding, ns)
  selector <- switch(label,
    all_hard_bounds = legacy_all_hard_fixed_idx,
    t0_only = function(sampler) which(
      EMC2:::get_p_types(sampler$par_names) == "t0"),
    unconditioned = function(sampler) integer(),
    stop("unknown arm"))
  assign(binding, selector, envir = ns)
  lockBinding(binding, ns)
  on.exit({
    unlockBinding(binding, ns)
    assign(binding, old, envir = ns)
    lockBinding(binding, ns)
  }, add = TRUE)

  set.seed(20261010L)
  started <- proc.time()[["elapsed"]]
  fit <- EMC2::fit(arm, stage = "sample",
    stop_criteria = list(sample = list(iter = n_draws)), step_size = 100L,
    cores_for_chains = 1L, cores_per_chain = 4L, particle_factor = 10,
    max_tries = 1L, trim = FALSE, verbose = FALSE, verboseProgress = FALSE)
  elapsed <- proc.time()[["elapsed"]] - started

  nch <- length(fit); p <- fit[[1L]]$n_pars; nsub <- fit[[1L]]$n_subjects
  gd <- array(NA_real_, c(n_draws, nch, 2L * p),
    dimnames = list(NULL, NULL,
      c(paste0("mu:", pn), paste0("logsd:", pn))))
  subject_names <- names(fit[[1L]]$data)
  if (is.null(subject_names) || length(subject_names) != nsub)
    subject_names <- as.character(seq_len(nsub))
  alpha_names <- paste0("alpha:", rep(pn, times = nsub), "@",
                        rep(subject_names, each = p))
  ad <- array(NA_real_, c(n_draws, nch, p * nsub),
    dimnames = list(NULL, NULL, alpha_names))
  for (j in seq_len(nch)) {
    ch <- fit[[j]]
    idx <- tail(which(ch$samples$stage[seq_len(ch$samples$idx)] == "sample"), n_draws)
    stopifnot(length(idx) == n_draws)
    gd[, j, seq_len(p)] <- t(ch$samples$theta_mu[, idx, drop = FALSE])
    gd[, j, p + seq_len(p)] <- t(vapply(idx, function(k)
      0.5 * log(diag(ch$samples$theta_var[, , k])), numeric(p)))
    ad[, j, ] <- t(matrix(ch$samples$alpha[, , idx], nrow = p * nsub))
  }
  ge <- posterior::summarise_draws(posterior::as_draws_array(gd),
    posterior::ess_bulk, posterior::rhat)
  ae <- posterior::summarise_draws(posterior::as_draws_array(ad),
    posterior::ess_bulk)
  names(ge)[names(ge) == "posterior::ess_bulk"] <- "ess_bulk"
  names(ge)[names(ge) == "posterior::rhat"] <- "rhat"
  names(ae)[names(ae) == "posterior::ess_bulk"] <- "ess_bulk"
  mu <- ge[startsWith(ge$variable, "mu:"), , drop = FALSE]
  logsd <- ge[startsWith(ge$variable, "logsd:"), , drop = FALSE]
  t0 <- ge[startsWith(ge$variable, "mu:t0"), , drop = FALSE]
  all_ess <- c(ge$ess_bulk, ae$ess_bulk)
  block <- lapply(fit, function(ch) ch$collapsed$plan$block)
  allocated_rounds <- vapply(fit, function(ch) {
    x <- ch$collapsed$plan$rounds
    if (length(x)) sum(x) else 0
  }, numeric(1))
  hmc_rates <- unlist(lapply(fit, function(ch) {
    if (is.null(ch$subject_hmc$alloc$rate)) numeric() else
      ch$subject_hmc$alloc$rate
  }), use.names = FALSE)
  out <- data.frame(case = case, arm = label, chains = nch,
    replicate = replicate, production_mode = production_mode,
    draws_per_chain = n_draws, seconds = elapsed,
    group_mu_median_bulk_ess = median(mu$ess_bulk),
    group_mu_median_ess_per_second = median(mu$ess_bulk) / elapsed,
    group_logsd_median_bulk_ess = median(logsd$ess_bulk),
    group_logsd_median_ess_per_second = median(logsd$ess_bulk) / elapsed,
    t0_group_median_bulk_ess = if (nrow(t0)) median(t0$ess_bulk) else NA_real_,
    t0_group_median_ess_per_second = if (nrow(t0)) median(t0$ess_bulk) / elapsed else NA_real_,
    group_mu_median_rhat = median(mu$rhat, na.rm = TRUE),
    group_mu_max_rhat = max(mu$rhat, na.rm = TRUE),
    group_hyper_median_ess_per_second = median(ge$ess_bulk) / elapsed,
    subject_median_bulk_ess = median(ae$ess_bulk),
    subject_median_ess_per_second = median(ae$ess_bulk) / elapsed,
    all_parameter_median_ess_per_second = median(all_ess, na.rm = TRUE) / elapsed,
    all_parameter_q10_ess_per_second = unname(stats::quantile(
      all_ess / elapsed, 0.1, na.rm = TRUE)),
    mean_allocated_collapsed_rounds = mean(allocated_rounds),
    median_subject_hmc_rate = if (length(hmc_rates)) median(hmc_rates) else NA_real_,
    collapsed_expected_acceptance = mean(vapply(block, function(b)
      if (b$n) b$acc / b$n else NA_real_, 0.0), na.rm = TRUE),
    collapsed_moved_rate = mean(vapply(block, function(b)
      if (b$n) b$moved / b$n else NA_real_, 0.0), na.rm = TRUE))
  hyper_rows <- data.frame(family = ifelse(startsWith(ge$variable, "mu:"),
      "group_mean", "group_logsd"), variable = ge$variable,
    ess_bulk = ge$ess_bulk, rhat = ge$rhat, ess_per_second = ge$ess_bulk / elapsed)
  subject_rows <- data.frame(family = "subject", variable = ae$variable,
    ess_bulk = ae$ess_bulk, rhat = NA_real_, ess_per_second = ae$ess_bulk / elapsed)
  variable_outputs[[paste(case, label, sep = "_")]] <<- transform(
    rbind(hyper_rows, subject_rows), case = case, arm = label)
  out
}

mode <- Sys.getenv("EMC2_ESS_MODE", "t0_vs_all_bounds")
arms <- switch(mode,
  t0_vs_all_bounds = c("all_hard_bounds", "t0_only"),
  t0_vs_unconditioned = c("unconditioned", "t0_only"),
  stop("unknown ESS mode"))
results <- do.call(rbind, lapply(arms, run_arm))
variable_file <- file.path("plans/collapsed_transport_adaptation_2026-10-09",
  paste0(mode, "_all_parameter_ess_", tolower(case), "_run", replicate,
    if (production_mode) "_production.csv" else ".csv"))
write.csv(results, variable_file, row.names = FALSE)
write.csv(do.call(rbind, variable_outputs), sub("[.]csv$", "_variables.csv", variable_file),
  row.names = FALSE)
print(results, row.names = FALSE)
