library(EMC2)
fit <- EMC2:::restore_custom_kernel_pointers(
  EMC2:::restore_duplicates(readRDS("/tmp/lba_sampling_audit/snapshot.rds")),
  quiet = TRUE)

ratios <- c(.1, .25, .5, 1)
records <- list()
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
  t0_par <- which(grepl("^t0", ch$par_names))
  t0_mean <- t0_par
  t0_var <- chart$gamma_at[match(t0_par, chart$gamma)]
  sd <- sqrt(pmax(diag(ch$collapsed$plan$C %*%
                        t(ch$collapsed$plan$C)), 1e-12))

  for (mode in c("mean", "variance")) for (ratio in ratios) {
    C <- sweep(ch$collapsed$plan$C, 1L, rep(.06, chart$dim), `*`)
    selected <- if (mode == "mean") t0_mean else t0_var
    C <- sweep(C, 1L, replace(rep(1, chart$dim), selected, ratio), `*`)
    set.seed(41261010L + 1000L * j + match(mode, c("mean", "variance")) * 10L +
               match(ratio, ratios))
    eps <- EMC2:::.emc_collapsed_tune_eps(phi0, C,
      ch$collapsed$plan$eps, ctx)
    set.seed(51261010L + 1000L * j + match(mode, c("mean", "variance")) * 10L +
               match(ratio, ratios))
    n_probe <- 6L
    acc <- inner <- total <- mean_esjd <- variance_esjd <- other_esjd <- numeric(n_probe)
    seconds <- numeric(n_probe)
    for (k in seq_len(n_probe)) {
      started <- proc.time()[["elapsed"]]
      tr <- EMC2:::collapsed_trajectories_cpp(phi0, C, eps, 3L, 1.5, ctx$cpp)
      phi <- as.numeric(tr$phi)
      st <- EMC2:::.emc_collapsed_from_phi(phi, chart)
      new <- EMC2:::.emc_shear_state(st$mu, st$V, ctx$fits, ch$prior, ctx$spec)
      if (is.null(new) || !is.finite(new$lp)) next
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
      corr <- ll_new - ll0 - vapply(seq_len(ch$n_subjects), function(s) {
        f <- ctx$fits[[s]]
        if (is.null(f)) return(0)
        d <- alpha_new[, s] - alpha0[, s]
        sum(d * (f$lin - .5 * drop(f$prec %*%
                                    (alpha_new[, s] + alpha0[, s]))))
      }, numeric(1))
      exact <- if (is.finite(sum(corr))) exp(min(0, sum(corr))) else 0
      dphi <- (phi - phi0) / sd
      acc[k] <- exact
      inner[k] <- mean(tr$acc)
      total[k] <- exact * sum(forwardsolve(ch$collapsed$plan$C, phi - phi0)^2)
      mean_esjd[k] <- exact * sum(dphi[t0_mean]^2)
      variance_esjd[k] <- exact * sum(dphi[t0_var]^2)
      other_esjd[k] <- exact * sum(dphi[-c(t0_mean, t0_var)]^2)
      seconds[k] <- proc.time()[["elapsed"]] - started
    }
    records[[length(records) + 1L]] <- data.frame(
      chain = j, component = mode, ratio = ratio, effective_t0_scale = .06 * ratio,
      eps = eps, inner_acceptance = mean(inner), exact_acceptance = mean(acc),
      total_esjd = mean(total), t0_mean_esjd = mean(mean_esjd),
      t0_variance_esjd = mean(variance_esjd), other_esjd = mean(other_esjd),
      seconds = mean(seconds), esjd_per_second = mean(total) / mean(seconds))
    cat("chain", j, "component", mode, "ratio", ratio,
        "accept", signif(mean(acc), 3), "\n")
    flush.console()
  }
}
out <- do.call(rbind, records)
write.csv(out,
  "plans/collapsed_transport_adaptation_2026-10-09/t0_component_summary.csv",
  row.names = FALSE)
print(aggregate(cbind(exact_acceptance, total_esjd, esjd_per_second,
                      t0_mean_esjd, t0_variance_esjd, other_esjd) ~
                  component + ratio, out, median), row.names = FALSE)
