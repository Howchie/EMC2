# BAwDD: BAwD at gamma = 0 with a free reciprocal shape parameter alpha = 1/rho.

.bawdd_check_cols <- function(pars, launch) {
  need <- c(.ba_par_names(launch), "b", "A", "t0", "k", "ell", "alpha")
  missing <- setdiff(need, colnames(pars))
  if (length(missing))
    stop("BAwDD requires parameter columns ", paste(missing, collapse = ", "))
  need
}

dBAwDD <- function(rt, pars, launch = 1L, posdrift = TRUE) {
  nm <- .bawdd_check_cols(pars, launch)
  dt <- rt - pars[, "t0"]
  eta <- .tw_eta(pars)
  ok <- (rt > 0) & (dt > 0) & is.finite(dt) & (pars[, "b"] >= pars[, "A"])
  ok[is.na(ok)] <- FALSE
  out <- numeric(length(dt))
  if (any(ok)) {
    s <- .tw_fwd(dt[ok], eta[ok])
    out[ok] <- dbawdd(
      t = s, A = pars[ok, "A"], b = pars[ok, "b"],
      p1 = pars[ok, nm[1]], p2 = pars[ok, nm[2]], k = pars[ok, "k"],
      ell = pars[ok, "ell"], alpha = pars[ok, "alpha"], launch = as.integer(launch),
      posdrift = posdrift,
      delta = if (launch == 2L) pars[ok, "delta"] else 0
    ) * .tw_jac(dt[ok], eta[ok])
  }
  out
}

pBAwDD <- function(rt, pars, launch = 1L, posdrift = TRUE) {
  nm <- .bawdd_check_cols(pars, launch)
  dt <- rt - pars[, "t0"]
  eta <- .tw_eta(pars)
  ok <- (rt > 0) & (dt > 0) & (pars[, "b"] >= pars[, "A"])
  ok[is.na(ok)] <- FALSE
  out <- numeric(length(dt))
  if (any(ok)) {
    s <- .tw_fwd(dt[ok], eta[ok])
    out[ok] <- pbawdd(
      t = s, A = pars[ok, "A"], b = pars[ok, "b"],
      p1 = pars[ok, nm[1]], p2 = pars[ok, nm[2]], k = pars[ok, "k"],
      ell = pars[ok, "ell"], alpha = pars[ok, "alpha"], launch = as.integer(launch),
      posdrift = posdrift,
      delta = if (launch == 2L) pars[ok, "delta"] else 0
    )
  }
  out
}

rBAwDD <- function(lR, pars, ok = rep(TRUE, length(lR)), launch = 1L,
                   posdrift = TRUE) {
  # For BAwDD, the simulator is identical to rbawdd_cpp. We do not provide
  # a pure R fallback for BAwDD now that it includes ell; just use the C++ one.
  res <- rbawdd_cpp(pars, levels(lR), ok, as.integer(launch), posdrift)
  out <- .rfun_cpp_pack(res, levels(lR), length(lR) / length(levels(lR)))
  .apply_timed_guess_winner(out, levels(lR))
}

#' The Ballistic Accumulator with Drive Decay and free Alpha (BAwDD)
#'
#' BAwDD is the `gamma = 0` member of [BAwD] with the power-law shape
#' estimated rather than fixed.  Its trajectory is
#' `X(u) = z + V Q_alpha(u) - ell u`, where `alpha = 1/rho`.  For `alpha > 0`,
#' `h_alpha(u) = (1 + alpha k u)^(-1/alpha)`; `alpha = 0` is the exponential
#' limit `h_0(u) = exp(-k u)`.  The reciprocal power-law shape `alpha` is
#' estimated for each accumulator; all density and distribution evaluations
#' are closed form through the shared BAwD/LBA primitives.  The special value
#' `alpha = Inf` is accepted as the no-decay limit.
#'
#' `alpha` is unbounded above.  `alpha > 1` (equivalently `rho < 1`) is an
#' ordinary member: the drive clock `Q_alpha` is then unbounded, but with
#' `gamma = 0` the constant clearance still outruns it, so the trajectory keeps
#' a finite peak and a finite endpoint `Tmax` whenever `k > 0` and `ell > 0`.
#' The frozen-branch closed forms carry the sign of `rho - 1` explicitly and so
#' hold on both sides of `alpha = 1`.
#'
#' With `drift_distribution = "lognormal"` (the default), `V` is lognormal,
#' sampled on the natural scale as its arithmetic `mean = E[V]` and
#' `cv = SD(V)/E[V]`.  The `splitlognormal`, `weibull`, and `normal` launch
#' choices use the same parameterisations as [BAwD].
#'
#' | **Parameter** | **Transform** | **Natural scale** | **Default** | **Mapping** | **Interpretation** |
#' |---|---|---|---|---|---|
#' | *mean* | log | \[0, Inf\] | log(1) | | Arithmetic mean E[V] of the launch strength (lognormal launch). |
#' | *cv* | log | \[0, Inf\] | log(1) | | Coefficient of variation SD(V)/E[V] (lognormal launch). |
#' | *shape* | log | \[0, Inf\] | log(1) | | Weibull shape (Weibull launch only). |
#' | *mean* | log | \[0, Inf\] | log(1) | | Weibull arithmetic mean (Weibull launch only). |
#' | *v* | identity | \[-Inf, Inf\] | 1 | | Mean normal launch strength (normal launch only). |
#' | *sv* | log | \[0, Inf\] | log(1) | | SD of normal launch strength (normal launch only). |
#' | *B* | log | \[0, Inf\] | log(1) | *b* = *B* + *A* | Threshold distance. |
#' | *A* | log | \[0, Inf\] | log(0) | | Start-point range. |
#' | *t0* | log | \[0, Inf\] | log(0) | | Non-decision time. |
#' | *k* | log | \[0, Inf\] | log(0) | | Drive-decay rate. |
#' | *ell* | log | \[0, Inf\] | log(0) | | Constant clearance rate; `ell = 0` removes it. |
#' | *alpha* | log | \[0, Inf\] | log(1) | | Reciprocal power-law shape `alpha = 1/rho`; `alpha = 0` is exponential. |
#' | *eta* | identity | \[-Inf, Inf\] | 0 | | Operational-time warp parameter. |
#'
#' `alpha` is bounded on the natural scale by `[0, Inf)` -- unbounded above --
#' and has an exact zero exception selecting the exponential member.  `ell`
#' likewise has an exact zero exception removing clearance.  Optional fitting
#' parameters are `pContaminant` (omission probability) and `pGuess`
#' (uniform-outlier probability).
#'
#' @param drift_distribution Distribution of trialwise launch strength:
#'   `"lognormal"` (default), `"splitlognormal"`, `"weibull"`, or `"normal"`.
#' @param posdrift Logical. If `TRUE` (default), truncate normal launch
#'   strengths below zero; if `FALSE`, retain unrestricted normal launches.
#' @return A model list defining the BAwDD race model.
#' @export
BAwDD <- function(drift_distribution = c("lognormal", "normal",
                                         "splitlognormal", "weibull"),
                  posdrift = TRUE) {
  drift_distribution <- match.arg(drift_distribution)
  launch <- .ba_launch_code(drift_distribution, "BAwDD")
  lognormal <- launch %in% c(1L, 2L)
  weibull <- launch == 3L
  splitlognormal <- launch == 2L
  if ((lognormal || weibull) && !isTRUE(posdrift)) {
    stop("BAwDD: posdrift only applies to drift_distribution = \"normal\"; a ",
         "lognormal and Weibull launch strengths are positive by construction.")
  }

  if (weibull) {
    p_types <- c(shape = log(1), mean = log(1))
    transform <- c(shape = "exp", mean = "exp")
    minmax <- cbind(shape = c(1e-4, Inf), mean = c(1e-4, Inf))
  } else if (splitlognormal) {
    # The split launch keeps log-scale (mu, sigma, delta): mu is the exact
    # median and (m, cv, delta) has no closed-form inverse.
    p_types <- c("mu" = 0, "sigma" = log(1), "delta" = 0)
    transform <- c(mu = "identity", sigma = "exp", delta = "identity")
    minmax <- cbind(mu = c(-Inf, Inf), sigma = c(1e-4, Inf),
                    delta = c(-Inf, Inf))
  } else if (lognormal) {
    # Natural-scale launch moments; the kernel converts to (mu, sigma) via
    # launch_lognormal_pair() in src/wald_functions.h.
    p_types <- c("mean" = log(1), "cv" = log(1))
    transform <- c(mean = "exp", cv = "exp")
    minmax <- cbind(mean = c(1e-4, Inf), cv = c(1e-4, Inf))
  } else {
    p_types <- c(v = 1, sv = log(1))
    transform <- c(v = "identity", sv = "exp")
    minmax <- cbind(v = c(-Inf, Inf), sv = c(1e-4, Inf))
  }
  p_types <- c(p_types, B = log(1), A = log(0), t0 = log(0),
               k = log(0), ell = log(0), alpha = log(1))
  transform <- c(transform, B = "exp", A = "exp", t0 = "exp", k = "exp",
                 ell = "exp", alpha = "exp")
  minmax <- cbind(minmax, B = c(1e-4, Inf), A = c(1e-4, Inf),
                  t0 = c(0.05, Inf), k = c(1e-4, Inf), ell = c(0, Inf), alpha = c(0, Inf))
  exception <- c(A = 0, k = 0, ell = 0, alpha = 0)

  .tw <- add_time_warp_par(p_types, transform, minmax, exception)
  p_types <- .tw$p_types; transform <- .tw$transform
  minmax <- .tw$minmax; exception <- .tw$exception
  .nuis <- add_nuisance_pars(p_types, transform, minmax, exception)
  p_types <- .nuis$p_types; transform <- .nuis$transform
  minmax <- .nuis$minmax; exception <- .nuis$exception

  list(
    type = "RACE",
    c_name = paste0("BAwDD", if (weibull) "_WEIB"
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
      # `ell` is sampled here, so the endpoint is only infinite where the
      # geometry says so (ell = 0 or k = 0).  Map alpha onto the shared
      # kernel's rho exactly as the compiled adapter does: alpha = 0 is the
      # exponential member (rho = Inf) and alpha = Inf the no-decay member
      # (k = 0), which bawd_tmax_vec then reports as Tmax = Inf itself.
      alpha <- pars[, "alpha"]
      k <- pars[, "k"]
      alpha_inf <- alpha > 0 & !is.finite(alpha)
      rho <- ifelse(alpha > 0 & is.finite(alpha), 1 / alpha, Inf)
      k[alpha_inf] <- 0
      Tmax_op <- bawd_tmax_vec(pars[, "A"], b, k, pars[, "ell"],
                               gamma = 0, rho = rho)
      Tmax <- .tw_inv(Tmax_op, .tw_eta(pars))
      cbind(pars, b = b, Tmax = Tmax, rt_max = pars[, "t0"] + Tmax)
    },
    rfun = function(data, pars) {
      .rfun_BAwDD(data$lR, pars, ok = attr(pars, "ok"), launch = launch,
                  posdrift = posdrift)
    },
    dfun = function(rt, pars) dBAwDD(rt, pars, launch = launch,
                                     posdrift = posdrift),
    pfun = function(rt, pars) pBAwDD(rt, pars, launch = launch,
                                     posdrift = posdrift),
    log_likelihood = function(pars, dadm, model, min_ll = log(1e-10)) {
      stop("BAwDD: the likelihood is implemented in the compiled race path; ",
           "the R likelihood route is not supported.")
    }
  )
}
