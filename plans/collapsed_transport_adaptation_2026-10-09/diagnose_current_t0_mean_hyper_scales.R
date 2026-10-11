pkgload::load_all(".", quiet = TRUE, compile = FALSE)
artifact_dir <- "plans/collapsed_transport_adaptation_2026-10-09"
prepared <- readRDS("/tmp/collapsed_transport_followup/length_adaptive_allocation_start.rds")
fit <- EMC2:::restore_custom_kernel_pointers(
  EMC2:::restore_duplicates(prepared$start), quiet = TRUE)
chain <- fit[[2L]]
p <- chain$n_pars
idx <- chain$samples$idx
ctx <- EMC2:::.emc_collapsed_context(chain)
chart <- ctx$chart
mu0 <- chain$samples$theta_mu[, idx]
V0 <- as.matrix(chain$samples$theta_var[, , idx])
alpha0 <- matrix(chain$samples$alpha[, , idx], p)
ll0 <- as.numeric(chain$samples$subj_ll[, idx])
cur <- EMC2:::.emc_shear_state(mu0, V0, ctx$fits, chain$prior, ctx$spec)
phi0 <- EMC2:::.emc_collapsed_to_phi(mu0, V0, chart)
residual <- vapply(seq_len(chain$n_subjects), function(s) {
  f <- cur$frames[[s]]
  drop(f$R %*% (alpha0[, s] - f$m))
}, numeric(p))
if (is.null(dim(residual))) residual <- matrix(residual, p)

t0_mean <- which(startsWith(chain$par_names, "t0"))
t0_all <- EMC2:::.emc_collapsed_t0_chart_idx(chain$par_names, chart)
t0_hyper <- setdiff(t0_all, t0_mean)
C0 <- chain$collapsed$plan$C
scale <- chain$collapsed$plan$scale
path_length <- chain$collapsed$plan$length
base_sd <- sqrt(pmax(diag(C0 %*% t(C0)), 1e-12))
mean_ratios <- c(0.25, 0.5, 1, 2)
hyper_ratios <- c(0.25, 0.5, 1)
n_probe <- as.integer(Sys.getenv("EMC2_T0_DIAG_PROBES", "100"))
records <- list()
g <- 0L

for (mean_ratio in mean_ratios) for (hyper_ratio in hyper_ratios) {
  g <- g + 1L
  rows <- rep(scale, chart$dim)
  rows[t0_mean] <- rows[t0_mean] * mean_ratio
  rows[t0_hyper] <- rows[t0_hyper] * hyper_ratio
  C <- sweep(C0, 1L, rows, `*`)
  set.seed(20261010L + g)
  eps <- EMC2:::.emc_collapsed_tune_eps(phi0, C,
    chain$collapsed$plan$eps, ctx, path_length)
  set.seed(30261010L + g)
  acc <- total <- mean_esjd <- hyper_esjd <- other_esjd <- numeric(n_probe)
  elapsed <- numeric(n_probe)
  lower_violation <- logical(n_probe)
  rt_violation <- logical(n_probe)
  for (k in seq_len(n_probe)) {
    started <- proc.time()[["elapsed"]]
    tr <- EMC2:::collapsed_trajectories_cpp(phi0, C, eps, 3L,
      path_length, ctx$cpp)
    phi <- as.numeric(tr$phi)
    st <- EMC2:::.emc_collapsed_from_phi(phi, chart)
    new <- EMC2:::.emc_shear_state(st$mu, st$V, ctx$fits,
      chain$prior, ctx$spec)
    if (is.null(new) || !is.finite(new$lp)) {
      elapsed[k] <- proc.time()[["elapsed"]] - started
      next
    }
    alpha_new <- vapply(seq_len(chain$n_subjects), function(s) {
      f <- new$frames[[s]]
      f$m + drop(backsolve(f$R, residual[, s]))
    }, numeric(p))
    for (s in seq_len(chain$n_subjects)) {
      pars <- EMC2:::get_pars_oo(
        matrix(alpha_new[, s], 1L,
          dimnames = list(NULL, chain$par_names)), chain$data[[s]], chain$model)
      t0 <- pars[, "t0"]
      rt <- chain$data[[s]]$rt
      stopifnot(length(t0) == length(rt))
      lower_violation[k] <- lower_violation[k] ||
        any(is.finite(t0) & t0 < 0.05)
      rt_violation[k] <- rt_violation[k] || any(is.finite(t0) &
        is.finite(rt) & t0 >= rt)
    }
    # In this LBA fit pGuess and pContaminant are fixed at zero, so either
    # t0 support crossing gives an exact zero. Skip the expensive likelihood
    # evaluation for these proposals; their accepted jump is exactly zero.
    if (lower_violation[k] || rt_violation[k]) {
      elapsed[k] <- proc.time()[["elapsed"]] - started
      next
    }
    ll_new <- vapply(seq_len(chain$n_subjects), function(s) {
      a <- matrix(alpha_new[, s], 1L,
        dimnames = list(NULL, chain$par_names))
      as.numeric(EMC2:::calc_ll_manager(a, dadm = chain$data[[s]],
        model = chain$model))
    }, numeric(1))
    correction <- ll_new - ll0 - vapply(seq_len(chain$n_subjects), function(s) {
      f <- ctx$fits[[s]]
      if (is.null(f)) return(0)
      d <- alpha_new[, s] - alpha0[, s]
      sum(d * (f$lin - 0.5 * drop(f$prec %*%
        (alpha_new[, s] + alpha0[, s]))))
    }, numeric(1))
    log_accept <- if (all(is.finite(ll_new)) && is.finite(sum(correction)))
      min(0, sum(correction)) else -Inf
    acc[k] <- exp(log_accept)
    dphi <- phi - phi0
    z <- dphi / base_sd
    total[k] <- acc[k] * sum(forwardsolve(C0, dphi)^2)
    mean_esjd[k] <- acc[k] * sum(z[t0_mean]^2)
    hyper_esjd[k] <- acc[k] * sum(z[t0_hyper]^2)
    other_esjd[k] <- acc[k] * sum(z[-t0_all]^2)
    elapsed[k] <- proc.time()[["elapsed"]] - started
  }
  records[[g]] <- data.frame(chain = 2L, scale = scale,
    t0_mean_ratio = mean_ratio, t0_hyper_ratio = hyper_ratio,
    effective_t0_mean_scale = scale * mean_ratio,
    effective_t0_hyper_scale = scale * hyper_ratio,
    path_length = path_length, eps = eps,
    exact_acceptance = mean(acc), t0_mean_esjd = mean(mean_esjd),
    t0_hyper_esjd = mean(hyper_esjd), other_esjd = mean(other_esjd),
    total_esjd = mean(total), seconds = mean(elapsed),
    total_esjd_per_second = mean(total) / mean(elapsed),
    t0_mean_esjd_per_second = mean(mean_esjd) / mean(elapsed),
    proposals_with_t0_below_005 = sum(lower_violation),
    proposals_with_t0_at_or_above_rt = sum(rt_violation))
  cat("mean", mean_ratio, "hyper", hyper_ratio,
    "accept", round(mean(acc), 3), "t0 mean ESJD/s",
    signif(mean(mean_esjd) / mean(elapsed), 3), "\n")
  flush.console()
}

out <- do.call(rbind, records)
write.csv(out, file.path(artifact_dir,
  "current_t0_mean_hyper_scale_diagnostic.csv"), row.names = FALSE)
print(out[order(-out$t0_mean_esjd_per_second), ], row.names = FALSE)
