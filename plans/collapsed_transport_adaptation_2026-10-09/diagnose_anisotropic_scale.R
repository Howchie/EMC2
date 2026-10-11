library(EMC2)

snapshot <- "/tmp/lba_sampling_audit/snapshot.rds"
fit <- EMC2:::restore_custom_kernel_pointers(
  EMC2:::restore_duplicates(readRDS(snapshot)), quiet = TRUE)

# Compare isotropic scale with the same overall scale applied more strongly
# to t0's group mean and variance coordinates. Every candidate is a valid
# fixed preconditioner for the existing reversible surrogate HMC kernel; the
# exact correction remains the production correction.
grid <- expand.grid(scale = c(.03, .06, .12, .24, .48, 1),
                    t0_ratio = c(.25, .5, 1))
grid <- grid[order(grid$scale, grid$t0_ratio), , drop = FALSE]
records <- list()
subjects <- list()

for (j in seq_along(fit)) {
  ch <- fit[[j]]
  p <- ch$n_pars
  idx <- ch$samples$idx
  ctx <- EMC2:::.emc_collapsed_context(ch)
  chart <- ctx$chart
  mu0 <- ch$samples$theta_mu[, idx]
  V0 <- as.matrix(ch$samples$theta_var[, , idx])
  alpha0 <- matrix(ch$samples$alpha[, , idx], p)
  ll0 <- as.numeric(ch$samples$subj_ll[, idx])
  cur <- EMC2:::.emc_shear_state(mu0, V0, ctx$fits, ch$prior, ctx$spec)
  phi0 <- EMC2:::.emc_collapsed_to_phi(mu0, V0, chart)
  residual <- vapply(seq_len(ch$n_subjects), function(s) {
    f <- cur$frames[[s]]
    drop(f$R %*% (alpha0[, s] - f$m))
  }, numeric(p))
  if (is.null(dim(residual))) residual <- matrix(residual, p)

  t0_par <- which(grepl("^t0", ch$par_names))
  t0_coord <- c(t0_par, chart$gamma_at[match(t0_par, chart$gamma)])
  t0_coord <- unique(t0_coord[is.finite(t0_coord)])
  subj_names <- dimnames(ch$samples$alpha)[[2L]][seq_len(ch$n_subjects)]
  marginal_sd <- sqrt(pmax(diag(ch$collapsed$plan$C %*%
                                  t(ch$collapsed$plan$C)), 1e-12))
  bounds <- ch$model()$bound$minmax
  exc <- ch$model()$bound$exception
  C0 <- ch$collapsed$plan$C

  for (g in seq_len(nrow(grid))) {
    overall <- grid$scale[g]
    t0_ratio <- grid$t0_ratio[g]
    row_scale <- rep(overall, chart$dim)
    row_scale[t0_coord] <- overall * t0_ratio
    C <- sweep(C0, 1L, row_scale, `*`)

    # A changed metric gets its own integrator step size. This prevents a
    # poor leapfrog acceptance from being mistaken for a bad outer scale.
    set.seed(20261010L + 10000L * j + g)
    eps <- EMC2:::.emc_collapsed_tune_eps(phi0, C,
      ch$collapsed$plan$eps, ctx)
    set.seed(30261010L + 10000L * j + g)
    n_probe <- 8L
    elapsed <- numeric(n_probe)
    exact_acc <- numeric(n_probe)
    inner_acc <- numeric(n_probe)
    total_esjd <- numeric(n_probe)
    t0_esjd <- numeric(n_probe)
    other_esjd <- numeric(n_probe)
    t0_bad <- integer(n_probe)
    bound_bad <- integer(n_probe)
    worst <- character(n_probe)
    for (k in seq_len(n_probe)) {
      started <- proc.time()[["elapsed"]]
      tr <- EMC2:::collapsed_trajectories_cpp(phi0, C, eps, 3L, 1.5,
                                               ctx$cpp)
      phi <- as.numeric(tr$phi)
      st <- EMC2:::.emc_collapsed_from_phi(phi, chart)
      new <- EMC2:::.emc_shear_state(st$mu, st$V, ctx$fits,
                                     ch$prior, ctx$spec)
      if (is.null(new) || !is.finite(new$lp)) {
        elapsed[k] <- proc.time()[["elapsed"]] - started
        next
      }
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
      correction <- ll_new - ll0 -
        vapply(seq_len(ch$n_subjects), function(s) {
          f <- ctx$fits[[s]]
          if (is.null(f)) return(0)
          d <- alpha_new[, s] - alpha0[, s]
          sum(d * (f$lin - .5 * drop(f$prec %*%
                                      (alpha_new[, s] + alpha0[, s]))))
        }, numeric(1))
      la <- sum(correction)
      exact_acc[k] <- if (is.finite(la)) exp(min(0, la)) else 0
      inner_acc[k] <- mean(tr$acc)
      dphi <- phi - phi0
      standardized <- dphi / marginal_sd
      total_esjd[k] <- exact_acc[k] * sum(forwardsolve(C0, dphi)^2)
      t0_esjd[k] <- exact_acc[k] * sum(standardized[t0_coord]^2)
      other_esjd[k] <- exact_acc[k] * sum(standardized[-t0_coord]^2)
      worst[k] <- subj_names[which.min(correction)]

      for (s in seq_len(ch$n_subjects)) {
        pars <- EMC2:::get_pars_oo(
          matrix(alpha_new[, s], 1L, dimnames = list(NULL, ch$par_names)),
          ch$data[[s]], ch$model)
        nm <- intersect(colnames(bounds), colnames(pars))
        bad <- sweep(pars[, nm, drop = FALSE], 2L, bounds[1L, nm], `<`) |
          sweep(pars[, nm, drop = FALSE], 2L, bounds[2L, nm], `>`)
        for (name in intersect(names(exc), nm))
          bad[, name] <- bad[, name] & pars[, name] != exc[name]
        bound_bad[k] <- bound_bad[k] + as.integer(any(rowSums(bad) > 0))
        rt <- ch$data[[s]]$rt
        t0_bad[k] <- t0_bad[k] + as.integer(any(
          is.finite(rt) & pars[, "t0"] >= rt))
      }
      elapsed[k] <- proc.time()[["elapsed"]] - started
      subjects[[length(subjects) + 1L]] <- data.frame(
        chain = j, candidate = g, scale = overall, t0_ratio = t0_ratio,
        probe = k, subject = subj_names,
        correction = correction)
    }
    records[[length(records) + 1L]] <- data.frame(
      chain = j, candidate = g, scale = overall, t0_ratio = t0_ratio,
      t0_effective = overall * t0_ratio, eps = eps,
      inner_acceptance = mean(inner_acc), exact_acceptance = mean(exact_acc),
      total_esjd = mean(total_esjd), t0_esjd = mean(t0_esjd),
      other_esjd = mean(other_esjd), seconds = mean(elapsed),
      total_esjd_per_second = mean(total_esjd) / mean(elapsed),
      probes_with_t0_support = sum(t0_bad > 0L),
      probes_with_bound_violation = sum(bound_bad > 0L),
      worst_subject = if (all(!nzchar(worst))) NA_character_ else
        names(sort(table(worst), decreasing = TRUE))[1L])
    cat("chain", j, "candidate", g, "scale", overall, "t0 ratio",
        t0_ratio, "accept", signif(mean(exact_acc), 3), "\n")
    flush.console()
  }
}

out <- do.call(rbind, records)
write.csv(out, "plans/collapsed_transport_adaptation_2026-10-09/anisotropic_scale_summary.csv",
          row.names = FALSE)
write.csv(do.call(rbind, subjects),
          "plans/collapsed_transport_adaptation_2026-10-09/anisotropic_scale_subjects.csv",
          row.names = FALSE)
print(aggregate(cbind(exact_acceptance, total_esjd,
                     total_esjd_per_second, t0_esjd, other_esjd) ~
                  scale + t0_ratio, out, median), row.names = FALSE)
