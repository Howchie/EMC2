# ============================================================================
# BAwD: the ballistic accumulator with power-law drive and clearance decay
#
#   h_rho(u) = (1 + k u / rho)^(-rho), with h_Inf(u) = exp(-k u)
#   Xdot(u) = V h_rho(u) - ell h_rho(u)^gamma, z ~ U(0, A), static
#
# The state does not leak (that is BAwL).  The fixed options are rho in
# {1, 2, 4, Inf} and gamma in {0, 1/2, 2/3, 3/4, 1}; neither is estimated.
#
# Numerical kernels in src/model_BAwD.h are shared by the wrappers and
# sampled likelihood.  Launch codes must match the C++ constants and adapter
# suffixes.
# ============================================================================

# 0 = truncated normal launch (v, sv); 1 = lognormal launch (mu, sigma);
# 2 = median-parameterised continuous split-lognormal launch (mu, sigma,
# delta); 3 = Weibull launch (shape, mean).  Must match BAWD_LAUNCH_* in src/model_BAwD.h and the values the
# adapter puts in ctx->bawd_launch from the c_name suffix.

# Split-normal geometry in log V.  The two half-normal widths have geometric
# mean sigma and log width ratio delta.  The split point c is chosen so that
# mu, rather than c, is the exact median.  Keeping the delta == 0 branch
# explicit preserves the ordinary lognormal path bit-for-bit.
.bawd_split_params <- function(mu, sigma, delta) {
  sigma_L <- sigma * exp(delta / 2)
  sigma_R <- sigma * exp(-delta / 2)
  a <- plogis(delta)
  c <- mu
  left <- delta > 0
  right <- delta < 0
  if (any(left)) c[left] <- mu[left] - sigma_L[left] *
    qnorm(1 / (4 * a[left]))
  if (any(right)) c[right] <- mu[right] - sigma_R[right] *
    qnorm(0.75 - 0.25 * exp(delta[right]))
  list(c = c, sigma_L = sigma_L, sigma_R = sigma_R, a = a)
}

# Draw V from the continuous split-lognormal.  At delta == 0 this deliberately
# calls rlnorm(), so the split variant has exactly the ordinary lognormal
# simulator when its extra parameter is zero.
.bawd_split_rlnorm <- function(mu, sigma, delta) {
  n <- max(length(mu), length(sigma), length(delta))
  mu <- rep_len(mu, n); sigma <- rep_len(sigma, n); delta <- rep_len(delta, n)
  out <- numeric(n)
  zero <- delta == 0
  if (any(zero)) out[zero] <- rlnorm(sum(zero), mu[zero], sigma[zero])
  nz <- !zero
  if (any(nz)) {
    sp <- .bawd_split_params(mu[nz], sigma[nz], delta[nz])
    left <- runif(sum(nz)) < sp$a
    z <- abs(rnorm(sum(nz)))
    y <- sp$c + z * sp$sigma_R
    y[left] <- sp$c[left] - z[left] * sp$sigma_L[left]
    out[nz] <- exp(y)
  }
  out
}

# Fixed clearance exponents gamma used by the kernels and constructor. Each is
# a distinct model. For finite rho and gamma < 1, the frozen-mass exponent is
# alpha = (rho - 1) / (rho * (1 - gamma)); rho = 1 is the logarithmic limit.
# rho = Inf retains the exponential member with alpha = 1 / (1 - gamma).
.bawd_gamma_values <- c(0, 0.5, 2 / 3, 0.75, 1)
.bawd_check_gamma <- function(gamma) {
  if (length(gamma) != 1L || !is.finite(gamma))
    stop("BAwD gamma must be one of 0, 1/2, 2/3, 3/4, 1; got ",
         paste(gamma, collapse = ", "))
  hit <- which(abs(gamma - .bawd_gamma_values) < 1e-12)
  if (length(hit) != 1L)
    stop("BAwD gamma must be one of 0, 1/2, 2/3, 3/4, 1; got ",
         paste(gamma, collapse = ", "))
  .bawd_gamma_values[hit]
}
.bawd_gamma_suffix <- function(gamma) {
  gamma <- .bawd_check_gamma(gamma)
  c("", "_GAM12", "_GAM23", "_GAM34", "_GAM100")[
    which(abs(gamma - .bawd_gamma_values) < 1e-12)
  ]
}

.bawd_rho_values <- c(1, 2, 4, Inf)
.bawd_check_rho <- function(rho) {
  if (length(rho) != 1L || is.na(rho) ||
      (!is.infinite(rho) && !is.finite(rho)) ||
      (is.infinite(rho) && rho < 0)) {
    stop("BAwD rho must be one of 1, 2, 4, Inf; got ",
         paste(rho, collapse = ", "))
  }
  if (is.infinite(rho)) return(Inf)
  hit <- which(abs(rho - c(1, 2, 4)) < 1e-12)
  if (length(hit) != 1L)
    stop("BAwD rho must be one of 1, 2, 4, Inf; got ",
         paste(rho, collapse = ", "))
  c(1, 2, 4)[hit]
}
.bawd_rho_suffix <- function(rho) {
  rho <- .bawd_check_rho(rho)
  if (!is.finite(rho)) "" else paste0("_RHO", rho)
}



.bawd_check_cols <- function(pars, launch) {
  need <- c(.ba_par_names(launch), "b", "A", "t0", "k", "ell")
  missing <- setdiff(need, colnames(pars))
  if (length(missing))
    stop("BAwD requires parameter columns ", paste(missing, collapse = ", "))
  need
}

dBAwD <- function(rt, pars, launch = 1L, posdrift = TRUE, gamma = 0,
                  rho = Inf) {
  gamma <- .bawd_check_gamma(gamma)
  rho <- .bawd_check_rho(rho)
  nm <- .bawd_check_cols(pars, launch)
  dt <- rt - pars[, "t0"]
  eta <- .tw_eta(pars)
  ok <- (rt > 0) & (dt > 0) & is.finite(dt) & (pars[, "b"] >= pars[, "A"])
  ok[is.na(ok)] <- FALSE
  out <- numeric(length(dt))
  if (any(ok)) {
    s <- .tw_fwd(dt[ok], eta[ok])
    out[ok] <- dbawd(t = s, A = pars[ok, "A"], b = pars[ok, "b"],
                     p1 = pars[ok, nm[1]], p2 = pars[ok, nm[2]],
                     delta = if (launch == 2L) pars[ok, "delta"] else 0,
                     k = pars[ok, "k"], ell = pars[ok, "ell"],
                     launch = as.integer(launch), posdrift = posdrift,
                     gamma = gamma, rho = rho) * .tw_jac(dt[ok], eta[ok])
  }
  out
}

pBAwD <- function(rt, pars, launch = 1L, posdrift = TRUE, gamma = 0,
                  rho = Inf) {
  gamma <- .bawd_check_gamma(gamma)
  rho <- .bawd_check_rho(rho)
  nm <- .bawd_check_cols(pars, launch)
  dt <- rt - pars[, "t0"]
  eta <- .tw_eta(pars)
  # rt = Inf is deliberately kept: the CDF there is F_max (the complement of the
  # never-finish mass), not one.  The compiled kernel returns the frozen branch.
  ok <- (rt > 0) & (dt > 0) & (pars[, "b"] >= pars[, "A"])
  ok[is.na(ok)] <- FALSE
  out <- numeric(length(dt))
  if (any(ok)) {
    s <- .tw_fwd(dt[ok], eta[ok])
    out[ok] <- pbawd(t = s, A = pars[ok, "A"], b = pars[ok, "b"],
                     p1 = pars[ok, nm[1]], p2 = pars[ok, nm[2]],
                     delta = if (launch == 2L) pars[ok, "delta"] else 0,
                     k = pars[ok, "k"], ell = pars[ok, "ell"],
                     launch = as.integer(launch), posdrift = posdrift,
                     gamma = gamma, rho = rho)
  }
  out
}


.bawdp_check_cols <- function(pars, launch) {
  need <- c(.ba_par_names(launch), "b", "A", "t0", "k", "lambda")
  missing <- setdiff(need, colnames(pars))
  if (length(missing))
    stop("BAwDp requires parameter columns ", paste(missing, collapse = ", "))
  need
}

dBAwDp <- function(rt, pars, launch = 1L, posdrift = TRUE) {
  nm <- .bawdp_check_cols(pars, launch)
  dt <- rt - pars[, "t0"]
  eta <- .tw_eta(pars)
  ok <- (rt > 0) & (dt > 0) & is.finite(dt) & (pars[, "b"] >= pars[, "A"])
  ok[is.na(ok)] <- FALSE
  out <- numeric(length(dt))
  if (any(ok)) {
    s <- .tw_fwd(dt[ok], eta[ok])
    out[ok] <- dbawdp(t = s, A = pars[ok, "A"], b = pars[ok, "b"],
                      p1 = pars[ok, nm[1]], p2 = pars[ok, nm[2]],
                      delta = if (launch == 2L) pars[ok, "delta"] else 0,
                      k = pars[ok, "k"], lambda = pars[ok, "lambda"],
                      launch = as.integer(launch), posdrift = posdrift) *
      .tw_jac(dt[ok], eta[ok])
  }
  out
}

pBAwDp <- function(rt, pars, launch = 1L, posdrift = TRUE) {
  nm <- .bawdp_check_cols(pars, launch)
  dt <- rt - pars[, "t0"]
  eta <- .tw_eta(pars)
  ok <- (rt > 0) & (dt > 0) & (pars[, "b"] >= pars[, "A"])
  ok[is.na(ok)] <- FALSE
  out <- numeric(length(dt))
  if (any(ok)) {
    s <- .tw_fwd(dt[ok], eta[ok])
    out[ok] <- pbawdp(t = s, A = pars[ok, "A"], b = pars[ok, "b"],
                      p1 = pars[ok, nm[1]], p2 = pars[ok, nm[2]],
                      delta = if (launch == 2L) pars[ok, "delta"] else 0,
                      k = pars[ok, "k"], lambda = pars[ok, "lambda"],
                      launch = as.integer(launch), posdrift = posdrift)
  }
  out
}

# First (rising-limb) crossing of distance d, or Inf if the peak falls short.
# Mirrors bawd_hit_time_r() in src/model_rng.h.
.bawd_hit_time <- function(V, d, k, ell, gamma = 0, rho = Inf) {
  gamma <- .bawd_check_gamma(gamma)
  rho <- .bawd_check_rho(rho)
  if (!isTRUE(d > 0)) return(0)
  if (is.na(V) || !isTRUE(V > 0)) return(Inf)
  if (is.finite(rho)) {
    if (k <= 1e-10) return(if (V > ell) d / (V - ell) else Inf)
    if (ell <= 1e-12) {
      if (abs(rho - 1) <= 1e-12)
        return(expm1(k * d / V) / k)
      S <- k * d * (rho - 1) / (V * rho)
      if (!isTRUE(S < 1)) return(Inf)
      return(rho * ((1 - S)^(-1 / (rho - 1)) - 1) / k)
    }
    if (gamma >= 1 - 1e-12) {
      if (!isTRUE(V > ell)) return(Inf)
      if (abs(rho - 1) <= 1e-12)
        return(expm1(k * d / (V - ell)) / k)
      S <- k * d * (rho - 1) / ((V - ell) * rho)
      if (!isTRUE(S < 1)) return(Inf)
      return(rho * ((1 - S)^(-1 / (rho - 1)) - 1) / k)
    }
    if (!isTRUE(V > ell)) return(Inf)
    pdiv <- function(c, L) {
      y <- c * L
      L * if (abs(y) < 1e-8) 1 + y * (0.5 + y / 6) else expm1(y) / y
    }
    kq <- function(x) {
      L <- log1p(x / rho)
      rho * pdiv(1 - rho, L)
    }
    kr <- function(x) {
      if (gamma <= 1e-12) return(x)
      L <- log1p(x / rho)
      rho * pdiv(1 - rho * gamma, L)
    }
    x_p <- rho * expm1(log(V / ell) / (rho * (1 - gamma)))
    if (V * kq(x_p) - ell * kr(x_p) < k * d) return(Inf)
    x <- k * d / (V - ell)
    if (!isTRUE(x > 0) || !is.finite(x)) x <- 1e-12
    lo <- 0
    hi <- x_p
    for (it in seq_len(100)) {
      f <- V * kq(x) - ell * kr(x) - k * d
      if (f > 0) hi <- x else lo <- x
      L <- log1p(x / rho)
      h <- exp(-rho * L)
      fp <- V * h - ell * h^gamma
      xn <- if (fp > 0) x - f / fp else (lo + hi) / 2
      if (!isTRUE(xn > lo) || !isTRUE(xn < hi) || !is.finite(xn))
        xn <- (lo + hi) / 2
      if (abs(xn - x) <= 1e-13 * max(1, xn)) return(xn / k)
      x <- xn
    }
    return(x / k)
  }
  qf <- function(u) if (k <= 1e-10) u else -expm1(-k * u) / k
  cf <- function(u) {
    if (k <= 1e-10 || gamma <= 1e-12) u
    else if (gamma >= 1 - 1e-12) qf(u)
    else -expm1(-gamma * k * u) / (gamma * k)
  }
  if (k <= 1e-10) return(if (V > ell) d / (V - ell) else Inf)
  if (ell <= 1e-12) {
    x <- 1 - k * d / V
    return(if (x > 0) -log(x) / k else Inf)
  }
  if (gamma >= 1 - 1e-12) {
    if (!isTRUE(V > ell)) return(Inf)
    x <- 1 - k * d / (V - ell)
    return(if (x > 0) -log(x) / k else Inf)
  }
  if (!isTRUE(V > ell)) return(Inf)
  u_p <- log(V / ell) / ((1 - gamma) * k)
  if (V * qf(u_p) - ell * cf(u_p) < d) return(Inf)
  u <- d / (V - ell)
  if (!isTRUE(u > 0) || !is.finite(u)) u <- 1e-12
  for (it in seq_len(100)) {
    f <- V * qf(u) - ell * cf(u) - d
    fp <- V * exp(-k * u) - ell * exp(-gamma * k * u)
    if (!isTRUE(fp > 0)) break
    un <- u - f / fp
    if (!isTRUE(un > 0)) un <- 0.5 * u
    if (un > u_p) un <- u_p
    done <- abs(un - u) <= 1e-13 * max(1, un)
    u <- un
    if (done) break
  }
  u
}

# Pure-R reference simulator; .rfun_BAwD uses the C++ kernel unless
# options(emc2.cpp_rfun = FALSE).  An all-Inf trial column is the package's
# omission convention (R = NA, rt = Inf), which make_data() already handles.
rBAwD <- function(lR, pars, ok = rep(TRUE, length(lR)), launch = 1L,
                  posdrift = TRUE, gamma = 0, rho = Inf) {
  gamma <- .bawd_check_gamma(gamma)
  rho <- .bawd_check_rho(rho)
  nm <- .bawd_check_cols(pars, launch)
  nr <- length(levels(lR))
  bad <- rep(NA, length(lR) / nr)
  out <- data.frame(R = bad, rt = bad)
  n_trials <- nrow(pars) / nr
  dt <- matrix(Inf, nrow = nr, ncol = n_trials)
  t0 <- pars[, "t0"]

  idx <- which(ok)
  if (length(idx)) {
    p <- pars[idx, , drop = FALSE]
    V <- if (launch == 3L) {
      .rweibull_mean(nrow(p), p[, "shape"], p[, "mean"])
    } else if (launch == 1L) {
      rlnorm(nrow(p), p[, "mu"], p[, "sigma"])
    } else if (launch == 2L) {
      .bawd_split_rlnorm(p[, "mu"], p[, "sigma"], p[, "delta"])
    } else {
      msm::rtnorm(nrow(p), p[, "v"], p[, "sv"],
                  lower = if (posdrift) 0 else -Inf)
    }
    z <- p[, "A"] * runif(nrow(p))
    hit <- mapply(.bawd_hit_time, V, p[, "b"] - z, p[, "k"], p[, "ell"],
                  MoreArgs = list(gamma = gamma, rho = rho))
    hit <- .tw_inv(hit, .tw_eta(p))
    hit[!is.finite(hit) | hit < 0] <- Inf
    dt[idx] <- hit
  }
  dt <- dt + matrix(t0, nrow = nr)

  bad_col <- colSums(!is.infinite(dt)) == 0L
  R <- apply(dt, 2, which.min)
  pick <- cbind(R, seq_len(ncol(dt)))
  rt <- dt[pick]
  R <- factor(levels(lR)[R], levels = levels(lR))
  R[bad_col] <- NA
  rt[bad_col] <- Inf
  ok <- matrix(ok, nrow = nr)[1, ]
  out$R[ok] <- levels(lR)[R][ok]
  out$R <- factor(out$R, levels = levels(lR))
  out$rt[ok] <- rt[ok]
  out <- .apply_timed_guess_winner(out, levels(lR))
  out
}

# Pure-R reference simulator for BAwDp.  The C++ implementation uses the same
# monotone clock inversion; this path is retained for the package's explicit
# emc2.cpp_rfun = FALSE fallback.
.bawdp_hit_time <- function(V, d, k, lambda) {
  if (!isTRUE(d > 0) || is.na(V) || !isTRUE(V > 0) ||
      !isTRUE(k >= 0) || !isTRUE(lambda >= 0) || !isTRUE(lambda < 1)) return(Inf)
  if (k <= 1e-10) return(d / (V * (1 - lambda)))
  target <- d / V
  if (lambda == 0) {
    x <- 1 - k * target
    return(if (x > 0) -log(x) / k else Inf)
  }
  u_star <- -log(lambda) / k
  m_max <- (1 - lambda + lambda * log(lambda)) / k
  if (!isTRUE(target > 0) || !isTRUE(target < m_max)) return(Inf)
  m <- function(u) -expm1(-k * u) / k - lambda * u
  lo <- 0; hi <- u_star
  u <- target / (1 - lambda)
  if (!isTRUE(u > lo) || !isTRUE(u < hi) || !is.finite(u)) u <- (lo + hi) / 2
  for (it in seq_len(100)) {
    f <- m(u) - target
    if (f > 0) hi <- u else lo <- u
    fp <- exp(-k * u) - lambda
    un <- if (fp > 0) u - f / fp else (lo + hi) / 2
    if (!isTRUE(un > lo) || !isTRUE(un < hi) || !is.finite(un)) un <- (lo + hi) / 2
    if (abs(un - u) <= 1e-13 * max(1, un)) return(un)
    u <- un
  }
  u
}

rBAwDp <- function(lR, pars, ok = rep(TRUE, length(lR)), launch = 1L,
                   posdrift = TRUE) {
  nm <- .bawdp_check_cols(pars, launch)
  nr <- length(levels(lR))
  bad <- rep(NA, length(lR) / nr)
  out <- data.frame(R = bad, rt = bad)
  n_trials <- nrow(pars) / nr
  dt <- matrix(Inf, nrow = nr, ncol = n_trials)
  idx <- which(ok)
  if (length(idx)) {
    p <- pars[idx, , drop = FALSE]
    V <- if (launch == 3L) {
      .rweibull_mean(nrow(p), p[, nm[1]], p[, nm[2]])
    } else if (launch == 1L) {
      rlnorm(nrow(p), p[, nm[1]], p[, nm[2]])
    } else if (launch == 2L) {
      .bawd_split_rlnorm(p[, "mu"], p[, "sigma"], p[, "delta"])
    } else {
      msm::rtnorm(nrow(p), p[, nm[1]], p[, nm[2]],
                  lower = if (posdrift) 0 else -Inf)
    }
    z <- p[, "A"] * runif(nrow(p))
    hit <- mapply(.bawdp_hit_time, V, p[, "b"] - z, p[, "k"], p[, "lambda"])
    hit <- .tw_inv(hit, .tw_eta(p))
    hit[!is.finite(hit) | hit < 0] <- Inf
    dt[idx] <- hit + p[, "t0"]
  }
  bad_col <- colSums(!is.infinite(dt)) == 0L
  R <- apply(dt, 2, which.min)
  pick <- cbind(R, seq_len(ncol(dt)))
  rt <- dt[pick]
  R <- factor(levels(lR)[R], levels = levels(lR))
  R[bad_col] <- NA; rt[bad_col] <- Inf
  ok_trial <- matrix(ok, nrow = nr)[1, ]
  out$R[ok_trial] <- levels(lR)[R][ok_trial]
  out$R <- factor(out$R, levels = levels(lR)); out$rt[ok_trial] <- rt[ok_trial]
  out <- .apply_timed_guess_winner(out, levels(lR))
  out
}
#' The Ballistic Accumulator with Drive Decay (BAwD)
#'
#' A race model in which each accumulator is driven by a power-law transient
#' signal and opposed by clearance with the same base kernel:
#' \verb{h_rho(u) = (1 + k u / rho)^(-rho)}, with the exponential limit
#' \verb{h_Inf(u) = exp(-k u)}. The fixed kernel option `rho` is one of
#' `1`, `2`, `4`, or `Inf`.
#'
#' @details
#' The fixed clearance exponent `gamma` can be `0`, `1/2`, `2/3`, `3/4`, or
#' `1`. For `gamma < 1` and finite `rho > 1`, the frozen-mass exponent is
#' \verb{alpha = (rho - 1) / (rho (1 - gamma))}; `rho = 1` is its logarithmic
#' limit. The `rho = Inf` member retains the exponential kernel.
#'
#' The trajectory is
#' \verb{X(u) = z + V Q_rho(u) - ell R_rho(u)}, where
#' \verb{Q_rho(u) = integral_0^u h_rho(s) ds} and
#' \verb{R_rho(u) = integral_0^u h_rho(s)^gamma ds}.
#' For `gamma = 1`, \verb{Xdot = (V - ell) h_rho}; there is no finite endpoint.
#' For `gamma < 1`, weak launch strengths can miss the threshold permanently
#' and a finite `T_max` occurs when `k > 0` and `ell > 0`.
#'
#' The start point is static and weak drive can produce an omission. The
#' fixed `rho` and `gamma` choices are model options, never estimated
#' parameters, so they do not appear in `p_types`.
#'
#' With `drift_distribution = "lognormal"` (the default) `log V ~ N(mu, sigma^2)`.
#' With `drift_distribution = "splitlognormal"`, `Y = log V` has the continuous
#' split-normal density with left and right widths
#' \verb{sigma_L = sigma exp(delta/2)} and
#' \verb{sigma_R = sigma exp(-delta/2)}.  Here `delta` is an unbounded real
#' log-width ratio and `mu` is the median of `Y`; the split point `c` is derived
#' from these parameters so that `P(Y <= mu) = 1/2` exactly.  Setting `delta = 0`
#' reduces exactly to the ordinary lognormal launch.
#' With `drift_distribution = "weibull"`, `V ~ Weibull(shape, scale)` on the
#' positive launch scale, parameterised publicly by its arithmetic `mean` and
#' shape. The conventional scale is `mean / Gamma(1 + 1 / shape)`. The Weibull
#' BAwD likelihood is closed form, including
#' its finite-`Tmax` frozen branch.
#'
#' | **Parameter** | **Transform** | **Natural scale** | **Default** | **Mapping** | **Interpretation** |
#' |---|---|---|---|---|---|
#' | *mu* | identity | \[-Inf, Inf\] | 0 | | Median log launch location. |
#' | *sigma* | log | \[0, Inf\] | log(1) | | Geometric-mean width of log launch. |
#' | *shape* | log | \[0, Inf\] | log(1) | | Weibull shape (Weibull launch only). |
#' | *mean* | log | \[0, Inf\] | log(1) | | Weibull arithmetic mean (Weibull launch only). |
#' | *B* | log | \[0, Inf\] | log(1) | *b* = *B* + *A* | Threshold distance. |
#' | *A* | log | \[0, Inf\] | log(0) | | Start-point range. |
#' | *t0* | log | \[0, Inf\] | log(0) | | Non-decision time. |
#' | *k* | log | \[0, Inf\] | log(0) | | Drive-decay rate. |
#' | *ell* | log | \[0, Inf\] | log(1) | | Clearance rate. |
#' | *eta* | identity | \[-Inf, Inf\] | 0 | | Operational-time warp parameter. |
#'
#' The trial-dependent transform derives `Tmax` from `A`, `b`, `k`, and
#' `ell`, and also reports `rt_max` after adding `t0`.
#' Optional fitting parameters: `pContaminant` is the omission probability and
#' `pGuess` is the uniform-outlier probability.
#' For `drift_distribution = "normal"`, use the `v` and `sv` rows instead of
#' `mu` and `sigma`; the normal launch is truncated at zero by default. The
#' split-lognormal option adds the optional identity-scale `delta` parameter
#' (default `0`).
#' For `drift_distribution = "weibull"`, use positive `shape` and arithmetic
#' `mean` (both use log/exp transforms).
#'
#' **Operational-time warp.** `eta` is a trailing free parameter (identity
#' transform, default `0`, unbounded on the natural scale) and is excluded from
#' `p_types_canonical`. For physical accumulation time `u = rt - t0`, let
#' `omega = exp(eta)` and `s = ((1 + u)^omega - 1) / omega`. The CDF and
#' survivor are the parent CDF/survivor evaluated at `s`; the density is the
#' parent density times the Jacobian `c_eta'(u) = (1 + u)^(omega - 1)`.
#' Simulation maps a parent internal finish `s` back with
#' `u = (1 + omega * s)^(1 / omega) - 1` and returns `rt = t0 + u`, so `t0`
#' remains additive. `eta = 0` is the exact identity warp.
#'
#' For BAwL, this contract applies only to the clock-free, uncorrelated
#' constructor (`erlang_type = "none"`, `correlated = FALSE`). BAwL clock or
#' correlated variants, `LogicalRulesLBA`, and all non-ballistic models reject
#' `eta`.
#'
#' @param drift_distribution Distribution of trialwise launch strength:
#'   `"lognormal"` (default), `"splitlognormal"`, `"weibull"`, or `"normal"`.
#' @param posdrift Logical. If `TRUE` (default), truncate normal launch
#'   strengths to be positive; if `FALSE`, append `IO` to the compiled model
#'   name. Only meaningful for normal launches.
#' @param gamma Fixed clearance exponent: `0`, `1/2`, `2/3`, `3/4`, or `1`.
#' @param rho Fixed base-kernel shape: `1`, `2`, `4`, or `Inf` (default).
#'   `Inf` is the exponential kernel and has no rho suffix; finite values
#'   append `"_RHO1"`, `"_RHO2"`, or `"_RHO4"` after the gamma suffix.
#' @return A model list defining the BAwD race model.
#' @examples
#' ADmat <- matrix(c(-1/2, 1/2), ncol = 1, dimnames = list(NULL, "d"))
#' matchfun <- function(d) d$S == d$lR
#' # ell is left out of the formula, so it stays at its default of 1 and fixes
#' # the evidence scale for the lognormal launch.
#' design_BAwD <- design(data = forstmann, model = BAwD, matchfun = matchfun,
#'                       formula = list(mu ~ lM, sigma ~ 1, B ~ E, A ~ 1,
#'                                      t0 ~ 1, k ~ 1),
#'                       contrasts = list(mu = list(lM = ADmat)))
#' @export
BAwD <- function(drift_distribution = c("lognormal", "normal", "splitlognormal", "weibull"),
                 posdrift = TRUE, gamma = 0, rho = Inf) {
  drift_distribution <- match.arg(drift_distribution)
  gamma <- .bawd_check_gamma(gamma)
  rho <- .bawd_check_rho(rho)

  launch <- .ba_launch_code(drift_distribution, "BAwD")
  lognormal <- launch %in% c(1L, 2L)
  weibull <- launch == 3L
  splitlognormal <- launch == 2L
  if ((lognormal || weibull) && !isTRUE(posdrift)) {
    stop("BAwD: posdrift only applies to drift_distribution = \"normal\"; a ",
         "lognormal and Weibull launch strengths are positive by construction.")
  }

  if (weibull) {
    p_types <- c("shape" = log(1), "mean" = log(1))
    transform <- c(shape = "exp", mean = "exp")
    minmax <- cbind(shape = c(1e-4, Inf), mean = c(1e-4, Inf))
  } else if (lognormal) {
    p_types <- c("mu" = 0, "sigma" = log(1))
    transform <- c(mu = "identity", sigma = "exp")
    minmax <- cbind(mu = c(-Inf, Inf), sigma = c(1e-4, Inf))
    if (splitlognormal) {
      p_types <- c(p_types, "delta" = 0)
      transform <- c(transform, delta = "identity")
      minmax <- cbind(minmax, delta = c(-Inf, Inf))
    }
  } else {
    p_types <- c("v" = 1, "sv" = log(1))
    transform <- c(v = "identity", sv = "exp")
    minmax <- cbind(v = c(-Inf, Inf), sv = c(1e-4, Inf))
  }
  # ell defaults to log(1): leaving it out of the formula fixes the evidence
  # scale, which is the recommended convention for the lognormal launch.
  p_types <- c(p_types, "B" = log(1), "A" = log(0), "t0" = log(0),
               "k" = log(0), "ell" = log(1))
  transform <- c(transform, B = "exp", A = "exp", t0 = "exp", k = "exp",
                 ell = "exp")
  minmax <- cbind(minmax, A = c(1e-4, Inf), B = c(1e-4, Inf),
                  t0 = c(0.05, Inf), k = c(1e-4, Inf), ell = c(1e-4, Inf))
  # ell = 0 (the static-start BAwL limit) and k = 0 (the LBA limit) must stay
  # exactly reachable, so both are bound exceptions rather than clamped.
  exception <- c(A = 0, k = 0, ell = 0)

  # pContaminant (omission) and pGuess (uniform outlier); see add_nuisance_pars().
  # Operational-time warp; eta = 0 is exactly this model.
  .tw <- add_time_warp_par(p_types, transform, minmax, exception)
  p_types <- .tw$p_types; transform <- .tw$transform
  minmax <- .tw$minmax; exception <- .tw$exception
  .nuis <- add_nuisance_pars(p_types, transform, minmax, exception)
  p_types <- .nuis$p_types; transform <- .nuis$transform
  minmax <- .nuis$minmax; exception <- .nuis$exception

  # "_SPLIT" follows "_LOGN"; the IO suffix is only reachable for normal
  # launches, and neither lognormal name contains the substring "IO".
  c_name <- paste0("BAwD", if (weibull) "_WEIB"
                           else if (lognormal) "_LOGN"
                           else if (!posdrift) "IO" else "",
                   if (splitlognormal) "_SPLIT" else "",
                   .bawd_gamma_suffix(gamma), .bawd_rho_suffix(rho))
  list(
    type = "RACE",
    c_name = c_name,
    drift_distribution = drift_distribution,
    gamma = gamma,
    rho = rho,
    p_types = p_types,
    p_types_canonical = setdiff(names(p_types),
                                c(.time_warp_par_name, .nuisance_par_names)),
    transform = list(func = transform),
    bound = list(minmax = minmax, exception = exception),
    Ttransform = function(pars, dadm) {
      b <- pars[, "B"] + pars[, "A"]
      Tmax_op <- bawd_tmax_vec(pars[, "A"], b, pars[, "k"], pars[, "ell"],
                               gamma = gamma, rho = rho)
      Tmax <- .tw_inv(Tmax_op, .tw_eta(pars))
      cbind(pars, b = b, Tmax = Tmax, rt_max = pars[, "t0"] + Tmax)
    },
    rfun = function(data, pars) {
      .rfun_BAwD(data$lR, pars, ok = attr(pars, "ok"), launch = launch,
                 posdrift = posdrift, gamma = gamma, rho = rho)
    },
    dfun = function(rt, pars) dBAwD(rt, pars, launch = launch,
                                    posdrift = posdrift, gamma = gamma,
                                    rho = rho),
    pfun = function(rt, pars) pBAwD(rt, pars, launch = launch,
                                    posdrift = posdrift, gamma = gamma,
                                    rho = rho),
    log_likelihood = function(pars, dadm, model, min_ll = log(1e-10)) {
      stop("BAwD: the likelihood is implemented in the compiled race path; ",
           "the R likelihood route is not supported.")
    }
  )
}

#' The Ballistic Accumulator with a Proportional-Clearance Drive (BAwDp)
#'
#' BAwDp has a separable transient evidence profile
#' \code{Xdot = V(exp(-k u) - lambda), X(0) = z},
#' where \code{z ~ Uniform(0, A)} and \code{b = B + A}.  Its internal evidence
#' clock is \code{m(u) = (1 - exp(-k u))/k - lambda*u}.  For
#' \code{0 < lambda < 1} the clock reaches its maximum at the common freeze
#' time \code{u* = log(1/lambda)/k} for every launch strength and start point.
#' The likelihood is therefore the ordinary closed-form LBA likelihood at
#' \code{m(u)}, multiplied by \code{mprime(u)} before the freeze time, and is
#' exactly zero after it.  The \code{lambda = 0} boundary is the pure drive-decay
#' model with an asymptotic (rather than finite) internal clock.  The \code{k = 0}
#' boundary is the standard LBA limit with effective drift scaled by \code{1 - lambda}.
#'
#' This is not BAwL: BAwL decays the accumulated start point, whereas BAwDp
#' keeps \code{z} static and scales the whole temporal drive profile by \code{V}.
#' A favourable start point consequently changes the omission probability.
#' All PDF/CDF evaluations use closed-form LBA primitives; no quadrature is
#' used.
#'
#' With `drift_distribution = "lognormal"` (the default) `log V ~ N(mu, sigma^2)`.
#' With `drift_distribution = "splitlognormal"`, `Y = log V` is continuous
#' split normal with \verb{sigma_L = sigma exp(delta/2)} and
#' \verb{sigma_R = sigma exp(-delta/2)}; `mu` is its exact median and `delta`
#' is an unbounded log-width ratio.  The derived split point is chosen so that
#' `P(Y <= mu) = 1/2`.  `delta = 0` is exactly the ordinary lognormal launch.
#' With `drift_distribution = "weibull"`, `V ~ Weibull(shape, scale)` with
#' public arithmetic `mean`; the conventional scale is derived internally and
#' the closed-form BAwDp likelihood uses positive `shape` and `mean` directly.
#'
#' | **Parameter** | **Transform** | **Natural scale** | **Default** | **Mapping** | **Interpretation** |
#' |---|---|---|---|---|---|
#' | *mu* | identity | \[-Inf, Inf\] | 0 | | Median log launch location. |
#' | *sigma* | log | \[0, Inf\] | log(1) | | Geometric-mean log width. |
#' | *shape* | log | \[0, Inf\] | log(1) | | Weibull shape (Weibull launch only). |
#' | *mean* | log | \[0, Inf\] | log(1) | | Weibull arithmetic mean (Weibull launch only). |
#' | *B* | log | \[0, Inf\] | log(1) | *b* = *B* + *A* | Threshold distance. |
#' | *A* | log | \[0, Inf\] | log(0) | | Start-point range. |
#' | *t0* | log | \[0, Inf\] | log(0) | | Non-decision time. |
#' | *k* | log | \[0, Inf\] | log(1) | | Drive-decay rate. |
#' | *lambda* | probit | \[0, 1\] | qnorm(.5) | | Proportional clearance fraction. |
#' | *eta* | identity | \[-Inf, Inf\] | 0 | | Operational-time warp parameter. |
#'
#' For `drift_distribution = "normal"`, use `v` and `sv` instead of the
#' lognormal launch rows; the normal launch is truncated at zero by default.
#' The split-lognormal option adds the optional identity-scale `delta` parameter
#' (default `0`).
#' For `drift_distribution = "weibull"`, use `shape` and `mean` instead of
#' the lognormal launch rows.
#' Optional fitting parameters: `pContaminant` is the omission probability and
#' `pGuess` is the uniform-outlier probability.
#'
#' **Operational-time warp.** `eta` is a trailing free parameter (identity
#' transform, default `0`, unbounded on the natural scale) and is excluded from
#' `p_types_canonical`. For physical accumulation time `u = rt - t0`, let
#' `omega = exp(eta)` and `s = ((1 + u)^omega - 1) / omega`. The CDF and
#' survivor are the parent CDF/survivor evaluated at `s`; the density is the
#' parent density times the Jacobian `c_eta'(u) = (1 + u)^(omega - 1)`.
#' Simulation maps a parent internal finish `s` back with
#' `u = (1 + omega * s)^(1 / omega) - 1` and returns `rt = t0 + u`, so `t0`
#' remains additive. `eta = 0` is the exact identity warp.
#'
#' For BAwL, this contract applies only to the clock-free, uncorrelated
#' constructor (`erlang_type = "none"`, `correlated = FALSE`). BAwL clock or
#' correlated variants, `LogicalRulesLBA`, and all non-ballistic models reject
#' `eta`.
#' @param drift_distribution Distribution of the launch strength:
#'   `"lognormal"` (default), `"splitlognormal"`, `"weibull"`, or `"normal"`.
#' @param posdrift Logical; truncate the normal launch at zero when `TRUE`.
#'   It has no effect for the positive lognormal and split-lognormal launches.
#' @return A model list defining the BAwDp race model.
#' @export
BAwDp <- function(drift_distribution = c("lognormal", "normal", "splitlognormal", "weibull"),
                  posdrift = TRUE) {
  drift_distribution <- match.arg(drift_distribution)
  launch <- .ba_launch_code(drift_distribution, "BAwD")
  lognormal <- launch %in% c(1L, 2L)
  weibull <- launch == 3L
  splitlognormal <- launch == 2L
  if ((lognormal || weibull) && !isTRUE(posdrift)) {
    stop("BAwDp: posdrift only applies to drift_distribution = \"normal\"; a ",
         "lognormal and Weibull launch strengths are positive by construction.")
  }

  if (weibull) {
    p_types <- c(shape = log(1), mean = log(1))
    transform <- c(shape = "exp", mean = "exp")
    minmax <- cbind(shape = c(1e-4, Inf), mean = c(1e-4, Inf))
  } else if (lognormal) {
    p_types <- c(mu = 0, sigma = log(1))
    transform <- c(mu = "identity", sigma = "exp")
    minmax <- cbind(mu = c(-Inf, Inf), sigma = c(1e-4, Inf))
    if (splitlognormal) {
      p_types <- c(p_types, delta = 0)
      transform <- c(transform, delta = "identity")
      minmax <- cbind(minmax, delta = c(-Inf, Inf))
    }
  } else {
    p_types <- c(v = 1, sv = log(1))
    transform <- c(v = "identity", sv = "exp")
    minmax <- cbind(v = c(-Inf, Inf), sv = c(1e-4, Inf))
  }
  p_types <- c(p_types, B = log(1), A = log(0), t0 = log(0),
               k = log(1), lambda = qnorm(0.5))
  transform <- c(transform, B = "exp", A = "exp", t0 = "exp",
                 k = "exp", lambda = "pnorm")
  minmax <- cbind(minmax, B = c(1e-4, Inf), A = c(1e-4, Inf),
                  t0 = c(0.05, Inf), k = c(1e-4, Inf),
                  lambda = c(0, 1 - 1e-8))
  # A = 0, k = 0, and lambda = 0 are exact, useful boundaries.  k = 0 gives the
  # standard LBA limit (with effective drift scaled by 1 - lambda), while
  # lambda = 0 gives the pure drive-decay model with unbounded support.
  exception <- c(A = 0, k = 0, lambda = 0)
  .tw <- add_time_warp_par(p_types, transform, minmax, exception)
  p_types <- .tw$p_types; transform <- .tw$transform
  minmax <- .tw$minmax; exception <- .tw$exception
  .nuis <- add_nuisance_pars(p_types, transform, minmax, exception)
  p_types <- .nuis$p_types; transform <- .nuis$transform
  minmax <- .nuis$minmax; exception <- .nuis$exception

  launch_pars <- .ba_par_names(launch)
  list(
    type = "RACE",
    c_name = paste0("BAwDp", if (weibull) "_WEIB"
                              else if (lognormal) "_LOGN"
                              else if (!posdrift) "IO" else "",
                    if (splitlognormal) "_SPLIT" else ""),
    drift_distribution = drift_distribution,
    p_types = p_types,
    p_types_canonical = setdiff(names(p_types),
                                c(.time_warp_par_name, .nuisance_par_names)),
    transform = list(func = transform),
    bound = list(minmax = minmax, exception = exception),
    Ttransform = function(pars, dadm) {
      b <- pars[, "B"] + pars[, "A"]
      k <- pars[, "k"]
      lambda <- pars[, "lambda"]
      finite_mask <- (k > 0) & (lambda > 0) & (lambda < 1)
      finite_mask[is.na(finite_mask)] <- FALSE
      Tmax_op <- rep(Inf, nrow(pars))
      Tmax_op[finite_mask] <- -log(lambda[finite_mask]) / k[finite_mask]
      Tmax <- .tw_inv(Tmax_op, .tw_eta(pars))
      cbind(pars, b = b, Tmax = Tmax, rt_max = pars[, "t0"] + Tmax)
    },
    rfun = function(data, pars) {
      .rfun_BAwDp(data$lR, pars, ok = attr(pars, "ok"), launch = launch,
                  posdrift = posdrift)
    },
    dfun = function(rt, pars) dBAwDp(rt, pars, launch = launch,
                                     posdrift = posdrift),
    pfun = function(rt, pars) pBAwDp(rt, pars, launch = launch,
                                     posdrift = posdrift),
    log_likelihood = function(pars, dadm, model, min_ll = log(1e-10)) {
      log_likelihood_race_missing(pars = pars, dadm = dadm, model = model,
                                   min_ll = min_ll)
    }
  )
}
