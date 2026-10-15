# Offline subject coordinates in the conditional frames used by the shear.
conditional_measure <- function(emc, ms, hst, ev) {
  ns <- asNamespace("EMC2")
  p <- emc[[1]]$n_pars; n <- emc[[1]]$n_subjects
  idx <- lapply(emc, function(ch) ev$.emc_alloc_draws(ch, max(2L, ch$alloc$window$n)))
  frames <- lapply(emc, function(ch) ns$.emc_shear_pack(ns$.emc_shear_fits(ch), p))
  Z <- lapply(seq_along(emc), function(j) {
    ch <- emc[[j]]; fits <- frames[[j]]
    k <- idx[[j]]
    Vbar <- apply(ch$samples$theta_var[, , k, drop = FALSE], 1:2, mean)
    Vi <- solve(Vbar)
    ref <- ns$emc_shear_frames_cpp(fits$cube, fits$lin, Vi, drop(Vi %*% rowMeans(ch$samples$theta_mu[, k, drop = FALSE])))
    rotation <- lapply(seq_len(n), function(s) ms[[j]]$W[[s]] %*% backsolve(ref$R[, , s], diag(p)))
    z <- array(0, c(p, n, length(k)))
    for (t in seq_along(k)) {
      i <- k[t]; Vi <- solve(ch$samples$theta_var[, , i])
      fr <- ns$emc_shear_frames_cpp(fits$cube, fits$lin, Vi, drop(Vi %*% ch$samples$theta_mu[, i]))
      for (s in seq_len(n)) z[, s, t] <- rotation[[s]] %*% fr$R[, , s] %*% (ch$samples$alpha[, s, i] - fr$m[, s])
    }
    z
  })
  vars <- t(vapply(seq_len(n), function(s) ev$.emc_alloc_row_var(do.call(cbind, lapply(Z, function(z) matrix(z[, s, ], p)))), numeric(p)))
  adjusted <- lapply(seq_along(ms), function(j) {
    m <- ms[[j]]; w <- emc[[j]]$alloc$window
    m$s_var <- vars
    for (s in seq_len(n)) {
      z <- matrix(Z[[j]][, s, ], p)
      rho <- ev$.emc_alloc_rho(rowMeans((z[, -1L, drop = FALSE] - z[, -ncol(z), drop = FALSE])^2), vars[s, ])
      P <- m$p_fac[s, ]
      if (w$nP[s] >= ev$.EMC_ALLOC_MIN_PAIRS) P <- ev$.emc_alloc_rho(rowSums((m$W[[s]] %*% w$P[[s]]) * m$W[[s]]) / w$nP[s], vars[s, ])
      H <- ev$.emc_alloc_hmc_factor(hst[[j]]$evidence[[s]], m$W[[s]], vars[s, ], ns$.EMC_HMC_MIN_EVIDENCE)
      if (is.null(H)) H <- rep(1, p)
      factor <- ev$.emc_alloc_dev_factor(m$q_old[s], m$r_old[s], P, H)
      m$p_fac[s, ] <- P
      m$t_rho[s, ] <- pmin(1, rho / pmax(factor, 1e-6))
      m$s_rho[s, ] <- m$t_rho[s, ] * (1 - m$q_old[s] * (1 - P))
    }
    m
  })
  attr(adjusted, "residual_variance") <- vars
  adjusted
}
