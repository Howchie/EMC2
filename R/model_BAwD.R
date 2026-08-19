# ============================================================================
# BAwD: the ballistic accumulator with drive decay
#
#   U(t) = V exp(-k t)                              transient drive
#   Xdot(t) = U(t) - ell exp(-gamma k t)            z ~ U(0, A), static
#
# The state does not leak (that is BAwL). Clearance is constant at gamma = 0,
# fades more slowly than the drive at each interior clearance exponent, and
# co-decays with the drive at gamma = 1. Weak launch strengths can therefore
# miss the threshold permanently. For gamma < 1 this also creates a finite
# right endpoint T_max; co-decay retains omissions but has no finite endpoint.
#
# Numerical kernels in src/model_BAwD.h are shared by the wrappers and
# sampled likelihood.  Launch codes must match the C++ constants and adapter
# suffixes.
# ============================================================================

# 0 = truncated normal launch (v, sv); 1 = lognormal launch (mu, sigma).
# Must match BAWD_LAUNCH_* in src/model_BAwD.h and the value the adapter puts
# in ctx->bawd_launch from the "_LOGN" c_name suffix.
.bawd_launch_code <- function(drift_distribution) {
  switch(drift_distribution, lognormal = 1L, normal = 0L,
         stop("Unknown BAwD drift_distribution: ", drift_distribution))
}
# Fixed clearance exponents gamma used by the kernels and constructor. Each is
# a distinct model: the frozen-mass factor at saturation is 1-(ell/w)^r with
# r = 1/(1-gamma), so gamma = 0,1/2,2/3,3/4,1 give r = 1,2,3,4,Inf. gamma is
# the clearance decay exponent; r is a derived quantity, never a constructor
# option (see BAwD() docs).
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


.bawd_par_names <- function(launch) {
  if (launch == 1L) c("mu", "sigma") else c("v", "sv")
}

.bawd_check_cols <- function(pars, launch) {
  need <- c(.bawd_par_names(launch), "b", "A", "t0", "k", "ell")
  missing <- setdiff(need, colnames(pars))
  if (length(missing))
    stop("BAwD requires parameter columns ", paste(missing, collapse = ", "))
  need
}

dBAwD <- function(rt, pars, launch = 1L, posdrift = TRUE, gamma = 0) {
  gamma <- .bawd_check_gamma(gamma)
  nm <- .bawd_check_cols(pars, launch)
  dt <- rt - pars[, "t0"]
  ok <- (rt > 0) & (dt > 0) & is.finite(dt) & (pars[, "b"] >= pars[, "A"])
  ok[is.na(ok)] <- FALSE
  out <- numeric(length(dt))
  if (any(ok)) {
    out[ok] <- dbawd(t = dt[ok], A = pars[ok, "A"], b = pars[ok, "b"],
                     p1 = pars[ok, nm[1]], p2 = pars[ok, nm[2]],
                     k = pars[ok, "k"], ell = pars[ok, "ell"],
                     launch = as.integer(launch), posdrift = posdrift,
                     gamma = gamma)
  }
  out
}

pBAwD <- function(rt, pars, launch = 1L, posdrift = TRUE, gamma = 0) {
  gamma <- .bawd_check_gamma(gamma)
  nm <- .bawd_check_cols(pars, launch)
  dt <- rt - pars[, "t0"]
  # rt = Inf is deliberately kept: the CDF there is F_max (the complement of the
  # never-finish mass), not one.  The compiled kernel returns the frozen branch.
  ok <- (rt > 0) & (dt > 0) & (pars[, "b"] >= pars[, "A"])
  ok[is.na(ok)] <- FALSE
  out <- numeric(length(dt))
  if (any(ok)) {
    out[ok] <- pbawd(t = dt[ok], A = pars[ok, "A"], b = pars[ok, "b"],
                     p1 = pars[ok, nm[1]], p2 = pars[ok, nm[2]],
                     k = pars[ok, "k"], ell = pars[ok, "ell"],
                     launch = as.integer(launch), posdrift = posdrift,
                     gamma = gamma)
  }
  out
}


.bawdp_check_cols <- function(pars, launch) {
  need <- c(.bawl_par_names(launch), "b", "A", "t0", "k", "lambda")
  missing <- setdiff(need, colnames(pars))
  if (length(missing))
    stop("BAwDp requires parameter columns ", paste(missing, collapse = ", "))
  need
}

dBAwDp <- function(rt, pars, launch = 1L, posdrift = TRUE) {
  nm <- .bawdp_check_cols(pars, launch)
  dt <- rt - pars[, "t0"]
  ok <- (rt > 0) & (dt > 0) & is.finite(dt) & (pars[, "b"] >= pars[, "A"])
  ok[is.na(ok)] <- FALSE
  out <- numeric(length(dt))
  if (any(ok)) {
    out[ok] <- dbawdp(t = dt[ok], A = pars[ok, "A"], b = pars[ok, "b"],
                      p1 = pars[ok, nm[1]], p2 = pars[ok, nm[2]],
                      k = pars[ok, "k"], lambda = pars[ok, "lambda"],
                      launch = as.integer(launch), posdrift = posdrift)
  }
  out
}

pBAwDp <- function(rt, pars, launch = 1L, posdrift = TRUE) {
  nm <- .bawdp_check_cols(pars, launch)
  dt <- rt - pars[, "t0"]
  ok <- (rt > 0) & (dt > 0) & (pars[, "b"] >= pars[, "A"])
  ok[is.na(ok)] <- FALSE
  out <- numeric(length(dt))
  if (any(ok)) {
    out[ok] <- pbawdp(t = dt[ok], A = pars[ok, "A"], b = pars[ok, "b"],
                      p1 = pars[ok, nm[1]], p2 = pars[ok, nm[2]],
                      k = pars[ok, "k"], lambda = pars[ok, "lambda"],
                      launch = as.integer(launch), posdrift = posdrift)
  }
  out
}

# First (rising-limb) crossing of distance d, or Inf if the peak falls short.
# Mirrors bawd_hit_time_r() in src/model_rng.h: Newton from the LBA-limit time,
# which lies below the root because q(u) <= u and the trajectory is concave and
# increasing up to its peak.
# Mirrors bawd_hit_time_r() in src/model_rng.h.
.bawd_hit_time <- function(V, d, k, ell, gamma = 0) {
  gamma <- .bawd_check_gamma(gamma)
  if (!isTRUE(d > 0)) return(0)
  if (is.na(V) || !isTRUE(V > 0)) return(Inf)
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
                  posdrift = TRUE, gamma = 0) {
  gamma <- .bawd_check_gamma(gamma)
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
    V <- if (launch == 1L) {
      rlnorm(nrow(p), p[, "mu"], p[, "sigma"])
    } else {
      msm::rtnorm(nrow(p), p[, "v"], p[, "sv"],
                  lower = if (posdrift) 0 else -Inf)
    }
    z <- p[, "A"] * runif(nrow(p))
    hit <- mapply(.bawd_hit_time, V, p[, "b"] - z, p[, "k"], p[, "ell"],
                  MoreArgs = list(gamma = gamma))
    hit[!is.finite(hit) | hit < 0] <- Inf
    dt[idx] <- hit
  }
  dt <- dt + matrix(t0, nrow = nr)

  bad_col <- apply(dt, 2, function(x) all(is.infinite(x)))
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
    V <- if (launch == 1L) rlnorm(nrow(p), p[, nm[1]], p[, nm[2]]) else
      msm::rtnorm(nrow(p), p[, nm[1]], p[, nm[2]], lower = if (posdrift) 0 else -Inf)
    z <- p[, "A"] * runif(nrow(p))
    hit <- mapply(.bawdp_hit_time, V, p[, "b"] - z, p[, "k"], p[, "lambda"])
    hit[!is.finite(hit) | hit < 0] <- Inf
    dt[idx] <- hit + p[, "t0"]
  }
  bad_col <- apply(dt, 2, function(x) all(is.infinite(x)))
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
#' A race model in which each accumulator is driven by a transient signal
#' opposed by fading clearance \verb{ell exp(-gamma k t)}. The fixed clearance
#' decay exponent `gamma` can be `0`, `1/2`, `2/3`, `3/4`, or `1`.
#'
#' @details
#' `gamma` is the clearance decay exponent, not the frozen-mass factor's
#' exponent. For `0 <= gamma < 1`, the frozen mass is
#' \verb{1 - (ell / w)^r} with the derived exponent
#' \verb{r = 1 / (1 - gamma)}. Thus `gamma = 0, 1/2, 2/3, 3/4` give
#' `r = 1, 2, 3, 4`; `gamma = 1` is the co-decay limit represented by
#' `r = Inf`, but it has no finite saturation wall. `r` is never a
#' constructor option; it is fully determined by `gamma`.
#'
#' For `gamma > 0`, the trajectory is
#' \verb{X(u) = z + V(1 - exp(-k u))/k -
#' ell(1 - exp(-gamma k u))/(gamma k)},
#' with the `gamma = 0` limit
#' \verb{X(u) = z + V(1 - exp(-k u))/k - ell u}.
#' The required launch by time `u` is
#' \verb{V^*(u,z) = [k(b-z) + ell(1-exp(-gamma k u))/gamma] /
#' [1-exp(-k u)]}
#' for `gamma > 0`, with the corresponding `k ell u` term at `gamma = 0`.
#' For `gamma < 1`, saturation uses `w = ell exp((1-gamma) k u)` and
#' \verb{|dz/dw| = k^{-1}[1 - (ell/w)^r]}, where `r = 1/(1-gamma)`.
#' At `gamma = 1`, \verb{X(u) = z + (V-ell)(1-exp(-k u))/k}; the eventual-hit
#' condition is `V > ell + k(b - z)`, so responses can be arbitrarily late.
#'
#' The start point is static and weak drive can produce an omission.
#' For `gamma < 1`, `T_max` is finite when `k > 0` and `ell > 0`; the
#' `gamma = 1` endpoint (co-decay) retains omissions but has no finite
#' endpoint. The five regimes are fixed constructor options, not estimated
#' parameters, so all variants have the same number of free parameters
#' (`gamma` never appears in `p_types`).
#'
#' Default values are used for all parameters that are not explicitly listed in
#' the `formula` argument of `design()`. They can also be accessed with
#' `BAwD()$p_types`.
#'
#' With `drift_distribution = "lognormal"` (the default) `log V ~ N(mu, sigma^2)`.
#'
#' | **Parameter** | **Transform** | **Natural scale** | **Default** | **Mapping** | **Interpretation** |
#' |---|---|---|---|---|---|
#' | *mu* | identity | \[-Inf, Inf\] | 0 | | Mean of the log launch strength. |
#' | *sigma* | log | \[0, Inf\] | log(1) | | SD of log launch strength. |
#' | *B* | log | \[0, Inf\] | log(1) | *b* = *B* + *A* | Threshold distance. |
#' | *A* | log | \[0, Inf\] | log(0) | | Start-point range. |
#' | *t0* | log | \[0, Inf\] | log(0) | | Non-decision time. |
#' | *k* | log | \[0, Inf\] | log(0) | | Drive-decay rate. |
#' | *ell* | log | \[0, Inf\] | log(1) | | Clearance rate. |
#'
#' With `drift_distribution = "normal"`, `mu` and `sigma` are replaced by
#' `v` and `sv`; the launch is truncated positive when `posdrift = TRUE`.
#'
#' @param drift_distribution Distribution of trialwise launch strength:
#'   `"lognormal"` (default) or `"normal"`.
#' @param posdrift Logical. If `TRUE` (default), truncate normal launch
#'   strengths to be positive; if `FALSE`, append `IO` to the compiled model
#'   name. Only meaningful for normal launches.
#' @param gamma Fixed clearance decay exponent: `0` (default), `1/2`, `2/3`,
#'   `3/4`, or `1`. A model option, never an estimated parameter; the
#'   corresponding suffix appended to the compiled model name is `""`,
#'   `"_GAM12"`, `"_GAM23"`, `"_GAM34"`, and `"_GAM100"` respectively.
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
BAwD <- function(drift_distribution = c("lognormal", "normal"),
                 posdrift = TRUE, gamma = 0) {
  drift_distribution <- match.arg(drift_distribution)
  gamma <- .bawd_check_gamma(gamma)

  launch <- .bawd_launch_code(drift_distribution)
  lognormal <- (launch == 1L)
  if (lognormal && !isTRUE(posdrift)) {
    stop("BAwD: posdrift only applies to drift_distribution = \"normal\"; a ",
         "lognormal launch strength is positive by construction.")
  }

  if (lognormal) {
    p_types <- c("mu" = 0, "sigma" = log(1))
    transform <- c(mu = "identity", sigma = "exp")
    minmax <- cbind(mu = c(-Inf, Inf), sigma = c(1e-4, Inf))
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
  .nuis <- add_nuisance_pars(p_types, transform, minmax, exception)
  p_types <- .nuis$p_types; transform <- .nuis$transform
  minmax <- .nuis$minmax; exception <- .nuis$exception

  # "_LOGN" (not "_LN": resolve_race_model_adapter dispatches by substring and
  # "LNR" is an existing key).  The IO suffix is only reachable for the normal
  # launch, and "BAwD_LOGN" deliberately contains no "IO".
  c_name <- paste0("BAwD", if (lognormal) "_LOGN"
                           else if (!posdrift) "IO" else "",
                   .bawd_gamma_suffix(gamma))
  list(
    type = "RACE",
    c_name = c_name,
    drift_distribution = drift_distribution,
    gamma = gamma,
    p_types = p_types,
    p_types_canonical = setdiff(names(p_types), .nuisance_par_names),
    transform = list(func = transform),
    bound = list(minmax = minmax, exception = exception),
    Ttransform = function(pars, dadm) {
      b <- pars[, "B"] + pars[, "A"]
      Tmax <- bawd_tmax_vec(pars[, "A"], b, pars[, "k"], pars[, "ell"],
                            gamma = gamma)
      cbind(pars, b = b, Tmax = Tmax, rt_max = pars[, "t0"] + Tmax)
    },
    rfun = function(data, pars) {
      .rfun_BAwD(data$lR, pars, ok = attr(pars, "ok"), launch = launch,
                 posdrift = posdrift, gamma = gamma)
    },
    dfun = function(rt, pars) dBAwD(rt, pars, launch = launch,
                                    posdrift = posdrift, gamma = gamma),
    pfun = function(rt, pars) pBAwD(rt, pars, launch = launch,
                                    posdrift = posdrift, gamma = gamma),
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
#' @param drift_distribution Distribution of the launch strength: `"lognormal"`
#'   (default), or `"normal"`.
#' @param posdrift Logical; truncate the normal launch at zero when `TRUE`.
#'   It has no effect for the lognormal launch.
#' @return A model list defining the BAwDp race model.
#' @export
BAwDp <- function(drift_distribution = c("lognormal", "normal"),
                  posdrift = TRUE) {
  drift_distribution <- match.arg(drift_distribution)
  launch <- .bawl_launch_code(drift_distribution)
  lognormal <- launch == 1L
  if (lognormal && !isTRUE(posdrift)) {
    stop("BAwDp: posdrift only applies to drift_distribution = \"normal\"; a ",
         "lognormal launch strength is positive by construction.")
  }

  if (lognormal) {
    p_types <- c(mu = 0, sigma = log(1))
    transform <- c(mu = "identity", sigma = "exp")
    minmax <- cbind(mu = c(-Inf, Inf), sigma = c(1e-4, Inf))
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
  .nuis <- add_nuisance_pars(p_types, transform, minmax, exception)
  p_types <- .nuis$p_types; transform <- .nuis$transform
  minmax <- .nuis$minmax; exception <- .nuis$exception

  launch_pars <- .bawl_par_names(launch)
  list(
    type = "RACE",
    c_name = paste0("BAwDp", if (lognormal) "_LOGN"
                              else if (!posdrift) "IO" else ""),
    drift_distribution = drift_distribution,
    p_types = p_types,
    p_types_canonical = setdiff(names(p_types), .nuisance_par_names),
    transform = list(func = transform),
    bound = list(minmax = minmax, exception = exception),
    Ttransform = function(pars, dadm) {
      b <- pars[, "B"] + pars[, "A"]
      k <- pars[, "k"]
      lambda <- pars[, "lambda"]
      finite_mask <- (k > 0) & (lambda > 0) & (lambda < 1)
      finite_mask[is.na(finite_mask)] <- FALSE
      Tmax <- rep(Inf, nrow(pars))
      Tmax[finite_mask] <- -log(lambda[finite_mask]) / k[finite_mask]
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
