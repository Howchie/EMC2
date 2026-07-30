# ============================================================================
# BAwD: the ballistic accumulator with drive decay
#
#   U(t) = V exp(-k t)                        transient drive
#   X(t) = z + (V/k)(1 - exp(-k t)) - ell t   z ~ U(0, A) and STATIC
#
# The state does not leak (that is BAwL); the drive decays while a constant
# clearance ell opposes it, so each accumulator rises to a peak and then falls.
# A launch strength that has not crossed b = B + A by its peak never will, which
# gives the model a hard right endpoint T_max and genuine omissions.
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
#' `ell = 0` and `A > 0` they are materially different models. At `k = 0`, BAwD
#' is the LBA with drift `V - ell`.
#'
#' The behavioural signature is the support of the response times. BAwD has a
#' hard right endpoint `T_max` (the peak time of the `z = 0` accumulator at its
#' critical launch strength), the density falls linearly to zero there, and the
#' remaining probability is a genuine never-finish mass. It therefore suits data
#' with a visible response window, substantial omissions and a rapid terminal
#' drop in hit density; it cannot produce a heavy or indefinitely extended tail.
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
#' | *pContaminant* | probit | \[0, 1\] | qnorm(0) | | Optional contamination probability handled by the data pipeline. |
#'
#' With `drift_distribution = "normal"`, `mu` and `sigma` are replaced by `v`
#' (identity, default 1) and `sv` (log, default `log(1)`), and
#' `V ~ N(v, sv^2)` truncated to be positive when `posdrift = TRUE`. This is the
#' variant that nests BAwL and the LBA, so it is the one to use when those
#' comparisons matter; the lognormal variant is the default because its
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
#' **Bounded support.** Because the supported window is finite, an observed
#' response time above `t0 + T_max` has density exactly zero and floors that
#' trial's likelihood. This is a constraint on where the posterior can live, to
#' be handled by initialisation and priors; `pContaminant` is a Bernoulli
#' *omission* rate and does not address late responses.
#'
#' @param drift_distribution Distribution of the trialwise launch strength:
#'   `"lognormal"` (the default) for `log V ~ N(mu, sigma^2)`, or `"normal"` for
#'   `V ~ N(v, sv^2)`.
#' @param posdrift Logical. If `TRUE` (default), truncate normal launch
#'   strengths to be positive; if `FALSE`, use the untruncated normal and append
#'   `IO` to the compiled model name. Only meaningful for
#'   `drift_distribution = "normal"`, since a lognormal launch strength is
#'   positive by construction.
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
                 posdrift = TRUE) {
  drift_distribution <- match.arg(drift_distribution)
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

  p_types <- c(p_types, pContaminant = qnorm(0))
  transform <- c(transform, pContaminant = "pnorm")
  minmax <- cbind(minmax, pContaminant = c(0.001, 0.999))
  exception <- c(exception, pContaminant = 0)

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
    p_types_canonical = names(p_types)[names(p_types) != "pContaminant"],
    transform = list(func = transform),
    bound = list(minmax = minmax, exception = exception),
    Ttransform = function(pars, dadm) {
      cbind(pars, b = pars[, "B"] + pars[, "A"])
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
