# ============================================================================
# BAwD: the ballistic accumulator with drive decay
#
#   U(t) = V exp(-k t)                        transient drive
#   X(t) = z + (V/k)(1 - exp(-k t)) - ell t   z ~ U(0, A) and STATIC
#
# The state does not leak (that is BAwL); the drive decays while a constant
# clearance ell opposes it, so each accumulator rises to a peak and then falls.
# A launch strength that has not crossed b = B + A by its peak never will, which
# gives the model a hard right endpoint T_max and genuine omissions when both
# the drive decays and the clearance is positive (k > 0, ell > 0).
#
# All of the numerics live in src/model_BAwD.h.  The dfun/pfun below call the
# SAME compiled kernels as the sampled likelihood, so make_data()/predict() and
# the fit cannot disagree -- but note that those kernels take the launch
# distribution as an explicit argument, because the exported entry points bypass
# the race context.  .bawd_launch_code() is the single place that value is
# derived, and it also drives the c_name suffix.
# ============================================================================

# 0 = truncated normal launch (v, sv); 1 = lognormal launch (mu, sigma).
# Must match BAWD_LAUNCH_* in src/model_BAwD.h and the value the adapter puts
# in ctx->bawd_launch from the "_LOGN" c_name suffix.
.bawd_launch_code <- function(drift_distribution) {
  switch(drift_distribution, lognormal = 1L, normal = 0L,
         stop("Unknown BAwD drift_distribution: ", drift_distribution))
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

dBAwD <- function(rt, pars, launch = 1L, posdrift = TRUE) {
  nm <- .bawd_check_cols(pars, launch)
  dt <- rt - pars[, "t0"]
  ok <- (rt > 0) & (dt > 0) & is.finite(dt) & (pars[, "b"] >= pars[, "A"])
  ok[is.na(ok)] <- FALSE
  out <- numeric(length(dt))
  if (any(ok)) {
    out[ok] <- dbawd(t = dt[ok], A = pars[ok, "A"], b = pars[ok, "b"],
                     p1 = pars[ok, nm[1]], p2 = pars[ok, nm[2]],
                     k = pars[ok, "k"], ell = pars[ok, "ell"],
                     launch = as.integer(launch), posdrift = posdrift)
  }
  out
}

pBAwD <- function(rt, pars, launch = 1L, posdrift = TRUE) {
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
                     launch = as.integer(launch), posdrift = posdrift)
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

# Single-accumulator diagnostics receive the reporting columns produced by
# BAwD(parameterization = "reduced").  The likelihood adapter itself never
# relies on these derived columns.
dBAwD_reduced <- function(rt, pars) dBAwD(rt, pars, launch = 1L, posdrift = TRUE)
pBAwD_reduced <- function(rt, pars) pBAwD(rt, pars, launch = 1L, posdrift = TRUE)

# First (rising-limb) crossing of distance d, or Inf if the peak falls short.
# Mirrors bawd_hit_time_r() in src/model_rng.h: Newton from the LBA-limit time,
# which lies below the root because q(u) <= u and the trajectory is concave and
# increasing up to its peak.
.bawd_hit_time <- function(V, d, k, ell) {
  if (!isTRUE(d > 0)) return(0)
  if (is.na(V) || !isTRUE(V > 0)) return(Inf)
  qf <- function(u) if (k <= 1e-10) u else -expm1(-k * u) / k
  if (k <= 1e-10) return(if (V > ell) d / (V - ell) else Inf)
  if (ell <= 1e-12) {
    x <- 1 - k * d / V
    return(if (x > 0) -log(x) / k else Inf)
  }
  if (!isTRUE(V > ell)) return(Inf)
  u_p <- log(V / ell) / k
  if (V * qf(u_p) - ell * u_p < d) return(Inf)
  u <- d / (V - ell)
  if (!isTRUE(u > 0) || !is.finite(u)) u <- 1e-12
  for (it in seq_len(100)) {
    f <- V * qf(u) - ell * u - d
    fp <- V * exp(-k * u) - ell
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
                  posdrift = TRUE) {
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
    hit <- mapply(.bawd_hit_time, V, p[, "b"] - z, p[, "k"], p[, "ell"])
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
      !isTRUE(k > 0) || !isTRUE(lambda >= 0) || !isTRUE(lambda < 1)) return(Inf)
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
#' A race model in which each accumulator is driven by a *transient* signal that
#' decays, and is opposed by a constant clearance. For accumulator `i` the start
#' point is `z_i = A * U_i` with `U_i ~ Uniform(0, 1)`, the threshold is
#' `b = B + A`, and with launch strength `V_i` the trajectory is
#' \deqn{X_i(t) = z_i + (V_i/k)(1 - e^{-kt}) - \ell t.}
#' The drive `V_i e^{-kt}` decays, so the trajectory rises to a peak at
#' `log(V_i/\ell)/k` and then falls. An accumulator that has not crossed `b` by
#' its peak never will.
#'
#' @details
#'
#' **This is not BAwL.** [BAwL()] applies leak to the accumulated evidence, so
#' the *state* (including the start point) decays and weak drive produces an
#' arbitrarily slow response. BAwD decays the *drive* while a constant `ell`
#' removes evidence, so the start point is static and weak drive produces an
#' *omission* instead. The two coincide only at `ell = 0` **and** `A = 0`; with
#' `ell = 0` and `A > 0` they are materially different models.
#'
#' At `k = 0`, the trajectory is the ordinary ballistic trajectory with
#' effective drift `D = V - ell`. For the normal launch with `posdrift = FALSE`,
#' this is exactly the unrestricted-normal LBA with mean drift `v - ell` and
#' SD `sv`. With `posdrift = TRUE`, the launch is truncated at `V > 0`, so
#' `D > -ell` rather than `D > 0`; the lognormal launch gives a shifted-
#' lognormal effective drift and is not the standard LBA drift family.
#'
#' When `k > 0` and `ell > 0`, the behavioural signature is a hard right
#' endpoint `T_max` (the peak time of the `z = 0` accumulator at its critical
#' launch strength), with a genuine never-finish mass. For a point start
#' (`A = 0`) the density approaches that endpoint linearly; with a nonzero
#' start-point range (`A > 0`) the collapsing live-start interval adds another
#' factor and the density approaches it quadratically. If `ell = 0`, the model
#' has no finite endpoint (although it can remain defective when the decayed
#' drive asymptote is below threshold); if `k = 0`, it has the ordinary LBA-like
#' infinite support. Thus the finite-window interpretation applies only to
#' `k > 0` and `ell > 0`.
#'
#' Default values are used for all parameters that are not explicitly listed in
#' the `formula` argument of `design()`. They can also be accessed with
#' `BAwD()$p_types`.
#'
#' With `drift_distribution = "lognormal"` (the default) `log V ~ N(mu, sigma^2)`:
#'
#' | **Parameter** | **Transform** | **Natural scale** | **Default** | **Mapping** | **Interpretation** |
#' |---|---|---|---|---|---|
#' | *mu* | identity | \[-Inf, Inf\] | 0 | | Mean of the log launch strength. |
#' | *sigma* | log | \[0, Inf\] | log(1) | | Between-trial SD of the log launch strength. |
#' | *B* | log | \[0, Inf\] | log(1) | *b* = *B* + *A* | Distance from the upper start-point range to the threshold. |
#' | *A* | log | \[0, Inf\] | log(0) | | Start-point range; the start point does **not** decay. |
#' | *t0* | log | \[0, Inf\] | log(0) | | Non-decision time. |
#' | *k* | log | \[0, Inf\] | log(0) | | Drive-decay rate; `k = 0` is the LBA limit. |
#' | *ell* | log | \[0, Inf\] | log(1) | | Tonic clearance / minimum effective sampling rate. |
#' | *pContaminant* | probit | \[0, 1\] | qnorm(0) | | Optional *omission* contaminant probability: mass at `rt = Inf` only, handled by the data pipeline |
#' | *pGuess* | probit | \[0, 1\] | qnorm(0) | | Optional uniform *guess* (outlier) probability, mixed into observed RT densities over the guess window |
#'
#' With `drift_distribution = "normal"`, `mu` and `sigma` are replaced by `v`
#' (identity, default 1) and `sv` (log, default `log(1)`), and
#' `V ~ N(v, sv^2)` truncated to be positive when `posdrift = TRUE`. This is the
#' variant that contains the exact BAwL corner at `ell = 0` and `A = 0` and
#' the exact unrestricted-normal LBA limit at `k = 0` only when `posdrift = FALSE`.
#' With `posdrift = TRUE`, the `k = 0` limit instead has the shifted truncation
#' `V - ell > -ell`. The lognormal variant is the default because its
#' likelihood needs no quadrature.
#'
#' **Fixing the evidence scale.** The evidence axis is defined only up to a
#' scale: `(V, ell, b, A) -> (cV, c*ell, cb, cA)` leaves every crossing time
#' unchanged, so exactly one of those must be fixed. With
#' `drift_distribution = "normal"` the usual `constants = c(sv = log(1))` does
#' it. With `drift_distribution = "lognormal"` it does **not**: `sigma` is
#' dimensionless. Fix `ell` instead --- which is what leaving `ell` out of the
#' `formula` does, since its default is `log(1)` --- and the spec's
#' dimensionless quantities become literal, `r = V` and `c = k b`. Omitting this
#' produces a ridge in the posterior rather than an error.
#'
#' **Identifiability.** At `k = 0` only `v - ell` (or the lognormal analogue) is
#' identified: `ell` is identified purely through the curvature of the
#' trajectory and the width of the response window. Prefer `ell ~ 1`, or a
#' coarse condition factor, and do not cross `ell` with the same factors as the
#' launch-strength mean.
#'
#' **Bounded support.** When `k > 0` and `ell > 0`, an observed response time
#' above `t0 + T_max` has density exactly zero and floors that trial's likelihood.
#' This is a constraint on where the posterior can live, to be handled by
#' initialisation and priors; `pContaminant` is a Bernoulli *omission* rate and
#' does not address late responses. `pGuess` does: it mixes a uniform density
#' over the guess window into observed RTs, so a response past `t0 + T_max` gets
#' a likelihood floor instead of a zero. The `ell = 0` and `k = 0` limits do not
#' have this finite support.
#'
#' @param drift_distribution Distribution of the trialwise launch strength:
#'   `"lognormal"` (the default) for `log V ~ N(mu, sigma^2)`, or `"normal"` for
#'   `V ~ N(v, sv^2)`.
#' @param posdrift Logical. If `TRUE` (default), truncate normal launch
#'   strengths to be positive; if `FALSE`, use the untruncated normal and append
#'   `IO` to the compiled model name. Only meaningful for
#'   `drift_distribution = "normal"`, since a lognormal launch strength is
#'   positive by construction.
#' @param parameterization Character. `"rate"` (default) uses the ordinary
#'   `(mu, sigma, B, A, t0, k, ell)` chart. `"reduced"` uses the identified
#'   lognormal chart `(y0, T_max, A, delta, sigma, t0)`, with `A` the absolute
#'   start-point range and `ell = 1` as the scale convention.
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
                 posdrift = TRUE,
                 parameterization = c("rate", "reduced")) {
  drift_distribution <- match.arg(drift_distribution)
  parameterization <- match.arg(parameterization)

  if (parameterization == "reduced") {
    if (drift_distribution != "lognormal")
      stop("BAwD(parameterization = \"reduced\") currently requires the lognormal launch.")
    if (!isTRUE(posdrift))
      stop("BAwD(parameterization = \"reduced\") has a positive lognormal launch by construction.")

    p_types <- c(y0 = log(1), T_max = log(1), A = log(0),
                 delta = 0, sigma = log(1), t0 = log(0))
    transform <- c(y0 = "exp", T_max = "exp", A = "exp",
                   delta = "identity", sigma = "exp", t0 = "exp")
    minmax <- cbind(y0 = c(1e-4, Inf), T_max = c(0.05, Inf),
                    A = c(1e-4, Inf), delta = c(-Inf, Inf),
                    sigma = c(1e-4, Inf), t0 = c(0.05, Inf))
    exception <- c(A = 0)
    .nuis <- add_nuisance_pars(p_types, transform, minmax, exception)
    p_types <- .nuis$p_types; transform <- .nuis$transform
    minmax <- .nuis$minmax; exception <- .nuis$exception

    # This reporting transform exposes the ordinary BAwD coordinates, but the
    # compiled likelihood maps the six-column prefix itself (see the
    # BAwD_REDUCED adapter).  A is absolute, so threshold changes do not
    # induce a start-point-range change.
    Ttransform_reduced <- function(pars, dadm) {
      y0 <- pars[, "y0"]; Tmax <- pars[, "T_max"]
      h <- expm1(y0) - y0
      k <- y0 / Tmax
      b <- h / k
      A <- pars[, "A"]
      cbind(pars,
            mu = y0 + pars[, "sigma"] * pars[, "delta"],
            B = b - A, k = k, ell = 1, b = b,
            Tmax = Tmax, rt_max = pars[, "t0"] + Tmax)
    }

    return(list(
      type = "RACE",
      c_name = "BAwD_REDUCED",
      drift_distribution = drift_distribution,
      parameterization = parameterization,
      p_types = p_types,
      p_types_canonical = setdiff(names(p_types), .nuisance_par_names),
      transform = list(func = transform),
      bound = list(minmax = minmax, exception = exception),
      Ttransform = Ttransform_reduced,
      rfun = function(data, pars) {
        .rfun_BAwD_reduced(data$lR, pars, ok = attr(pars, "ok"))
      },
      dfun = function(rt, pars) dBAwD_reduced(rt, pars),
      pfun = function(rt, pars) pBAwD_reduced(rt, pars),
      log_likelihood = function(pars, dadm, model, min_ll = log(1e-10)) {
        stop("BAwD: the likelihood is implemented in the compiled race path; ",
             "the R likelihood route is not supported.")
      }
    ))
  }

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
                           else if (!posdrift) "IO" else "")

  list(
    type = "RACE",
    c_name = c_name,
    drift_distribution = drift_distribution,
    p_types = p_types,
    p_types_canonical = setdiff(names(p_types), .nuisance_par_names),
    transform = list(func = transform),
    bound = list(minmax = minmax, exception = exception),
    Ttransform = function(pars, dadm) {
      b <- pars[, "B"] + pars[, "A"]
      # T_max is the estimable window quantity: k and ell are individually
      # near-degenerate along a manifold that holds it fixed, so it is what
      # should be reported and interpreted rather than either rate.  Inf
      # whenever k = 0 or ell = 0, where nothing saturates and support is
      # unbounded.  rt_max is the observable ceiling on this accumulator.
      Tmax <- bawd_tmax_vec(pars[, "A"], b, pars[, "k"], pars[, "ell"])
      cbind(pars, b = b, Tmax = Tmax, rt_max = pars[, "t0"] + Tmax)
    },
    rfun = function(data, pars) {
      .rfun_BAwD(data$lR, pars, ok = attr(pars, "ok"), launch = launch,
                 posdrift = posdrift)
    },
    dfun = function(rt, pars) dBAwD(rt, pars, launch = launch,
                                    posdrift = posdrift),
    pfun = function(rt, pars) pBAwD(rt, pars, launch = launch,
                                    posdrift = posdrift),
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
#' exactly zero after it.  The \code{lambda = 0}
#' boundary is the pure drive-decay model with an asymptotic (rather than finite)
#' internal clock.
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
  # A = 0 and lambda = 0 are exact, useful boundaries.  k remains positive so
  # BAwDp's finite-clock branch is the model used by the race adapter.
  exception <- c(A = 0, lambda = 0)
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
      cbind(pars, b = pars[, "B"] + pars[, "A"])
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
