pkgload::load_all(".", quiet = TRUE, compile = FALSE)

fit <- EMC2:::restore_custom_kernel_pointers(
  EMC2:::restore_duplicates(readRDS("/tmp/lba_sampling_audit/snapshot.rds")),
  quiet = TRUE)
lengths <- c(1.5, 3, 6, 12, 24)
n_probe <- 24L
records <- list()

tune_eps_for_length <- function(phi, C, eps0, ctx, len) {
  mu_da <- log(10 * eps0)
  h_bar <- 0
  log_eps <- log(eps0)
  log_bar <- log_eps
  for (k in seq_len(EMC2:::.EMC_COLLAPSED_PILOT)) {
    tr <- EMC2:::collapsed_trajectories_cpp(
      phi, C, exp(log_eps), 1L, len, ctx$cpp)
    phi <- tr$phi
    h_bar <- (1 - 1 / (k + 10)) * h_bar +
      (EMC2:::.EMC_COLLAPSED_TARGET - tr$acc) / (k + 10)
    log_eps <- mu_da - sqrt(k) / 0.05 * h_bar
    log_eps <- min(log(EMC2:::.EMC_COLLAPSED_EPS_RANGE[2L]),
                   max(log(EMC2:::.EMC_COLLAPSED_EPS_RANGE[1L]), log_eps))
    wk <- k^(-0.75)
    log_bar <- wk * log_eps + (1 - wk) * log_bar
  }
  min(EMC2:::.EMC_COLLAPSED_EPS_RANGE[2L],
      max(EMC2:::.EMC_COLLAPSED_EPS_RANGE[1L], exp(log_bar)))
}

for (j in seq_along(fit)) {
  ch <- fit[[j]]
  plan <- ch$collapsed$plan
  ctx <- EMC2:::.emc_collapsed_context(ch)
  p <- ch$n_pars
  idx <- ch$samples$idx
  mu0 <- ch$samples$theta_mu[, idx]
  V0 <- as.matrix(ch$samples$theta_var[, , idx])
  alpha0 <- matrix(ch$samples$alpha[, , idx], p)
  ll0 <- as.numeric(ch$samples$subj_ll[, idx])
  cur <- EMC2:::.emc_shear_state(mu0, V0, ctx$fits, ch$prior, ctx$spec)
  phi0 <- EMC2:::.emc_collapsed_to_phi(mu0, V0, ctx$chart)
  residual <- vapply(seq_len(ch$n_subjects), function(s) {
    f <- cur$frames[[s]]
    drop(f$R %*% (alpha0[, s] - f$m))
  }, numeric(p))
  if (is.null(dim(residual))) residual <- matrix(residual, p)
  subj_names <- dimnames(ch$samples$alpha)[[2L]][seq_len(ch$n_subjects)]
  t0_idx <- EMC2:::.emc_collapsed_t0_chart_idx(ch$par_names, ctx$chart)
  C0 <- plan$C
  scale <- c(0.03125, 0.03125, 0.04419417)[j]
  t0_ratio <- c(0.25, 0.25, 1)[j]
  C <- EMC2:::.emc_collapsed_preconditioner(C0, scale, t0_ratio, t0_idx)
  marginal_sd <- sqrt(pmax(diag(C0 %*% t(C0)), 1e-12))

  for (len in lengths) {
    set.seed(20261010L + 100L * j + match(len, lengths))
    eps <- tune_eps_for_length(phi0, C, 1.5, ctx, len)
    set.seed(30261010L + 100L * j + match(len, lengths))
    elapsed <- exact_acc <- inner_acc <- total_esjd <-
      t0_esjd <- other_esjd <- numeric(n_probe)
    for (k in seq_len(n_probe)) {
      started <- proc.time()[["elapsed"]]
      tr <- EMC2:::collapsed_trajectories_cpp(phi0, C, eps, 3L, len, ctx$cpp)
      phi <- as.numeric(tr$phi)
      st <- EMC2:::.emc_collapsed_from_phi(phi, ctx$chart)
      new <- EMC2:::.emc_shear_state(st$mu, st$V, ctx$fits, ch$prior, ctx$spec)
      if (!is.null(new) && is.finite(new$lp)) {
        alpha_new <- vapply(seq_len(ch$n_subjects), function(s) {
          f <- new$frames[[s]]
          f$m + drop(backsolve(f$R, residual[, s]))
        }, numeric(p))
        ll_new <- vapply(seq_len(ch$n_subjects), function(s) {
          a <- matrix(alpha_new[, s], 1L,
                      dimnames = list(NULL, ch$par_names))
          as.numeric(EMC2:::calc_ll_manager(a, dadm = ch$data[[s]],
                                             model = ch$model))
        }, numeric(1))
        correction <- ll_new - ll0 - vapply(seq_len(ch$n_subjects), function(s) {
          f <- ctx$fits[[s]]
          if (is.null(f)) return(0)
          d <- alpha_new[, s] - alpha0[, s]
          sum(d * (f$lin - 0.5 * drop(f$prec %*%
                                        (alpha_new[, s] + alpha0[, s]))))
        }, numeric(1))
        la <- sum(correction)
        exact_acc[k] <- if (is.finite(la)) exp(min(0, la)) else 0
        dphi <- phi - phi0
        standardized <- dphi / marginal_sd
        total_esjd[k] <- exact_acc[k] * sum(forwardsolve(C0, dphi)^2)
        t0_esjd[k] <- exact_acc[k] * sum(standardized[t0_idx]^2)
        other_esjd[k] <- exact_acc[k] * sum(standardized[-t0_idx]^2)
      }
      inner_acc[k] <- mean(tr$acc)
      elapsed[k] <- proc.time()[["elapsed"]] - started
    }
    records[[length(records) + 1L]] <- data.frame(
      chain = j, scale = scale, t0_ratio = t0_ratio,
      length = len, eps = eps,
      exact_acceptance = mean(exact_acc), inner_acceptance = mean(inner_acc),
      total_esjd = mean(total_esjd), t0_esjd = mean(t0_esjd),
      other_esjd = mean(other_esjd), seconds = mean(elapsed),
      total_esjd_per_second = mean(total_esjd) / mean(elapsed),
      conservative_log_score = EMC2:::.emc_collapsed_log_efficiency(
        ifelse(total_esjd > 0, log(total_esjd), -Inf), mean(elapsed)))
    cat("chain", j, "length", len, "eps", eps,
        "exact acceptance", mean(exact_acc), "\n")
  }
}

out <- do.call(rbind, records)
write.csv(out,
  "plans/collapsed_transport_adaptation_2026-10-09/trajectory_length_diagnostic.csv",
  row.names = FALSE)
print(out, row.names = FALSE)
