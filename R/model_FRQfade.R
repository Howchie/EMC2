# ============================================================================
# FRQfade: the FRQ counter (at h = 1) on the fading time
#   A(u) = (1 - exp(-kappa u)) / kappa,
# so F(u) = F_1(A(u)) and the completion probability h = F_1(1/kappa) is
# derived.  Kernel: src/model_FRQfade.h over src/model_FRQ.h.
# ============================================================================

.frqfade_check_cols <- function(pars) {
  need <- c("K", "v", "sv", "kappa", "t0")
  missing <- setdiff(need, colnames(pars))
  if (length(missing))
    stop("FRQfade requires parameter columns ", paste(missing, collapse = ", "))
  invisible(NULL)
}

.frqfade_call <- function(fun, dt, pars, ok, ...) {
  fun(dt[ok], pars[ok, "K"], pars[ok, "v"], pars[ok, "sv"], pars[ok, "kappa"],
      plo = .frq_opt(pars, ok, "plo", 0),
      clock = .frq_opt(pars, ok, "clock", 1), ...)
}

dFRQfade <- function(rt, pars) {
  .frqfade_check_cols(pars)
  dt <- rt - pars[, "t0"]
  ok <- (rt > 0) & (dt > 0) & is.finite(dt) & is.finite(pars[, "t0"])
  ok[is.na(ok)] <- FALSE
  out <- numeric(length(dt))
  if (any(ok)) out[ok] <- .frqfade_call(dfrqfade, dt, pars, ok)
  out
}

pFRQfade <- function(rt, pars) {
  .frqfade_check_cols(pars)
  dt <- rt - pars[, "t0"]
  # rt = Inf is kept: the CDF there is the derived h.
  ok <- (rt > 0) & (dt > 0) & is.finite(pars[, "t0"])
  ok[is.na(ok)] <- FALSE
  out <- numeric(length(dt))
  if (any(ok)) out[ok] <- .frqfade_call(pfrqfade, dt, pars, ok)
  out
}

sFRQfade <- function(rt, pars) {
  .frqfade_check_cols(pars)
  dt <- rt - pars[, "t0"]
  out <- rep(1, length(dt))
  ok <- (rt > 0) & (dt > 0) & is.finite(pars[, "t0"])
  ok[is.na(ok)] <- FALSE
  if (any(ok)) out[ok] <- .frqfade_call(pfrqfade, dt, pars, ok,
                                        lower_tail = FALSE)
  out
}

# Pure-R reference simulator; .rfun_FRQfade uses the C++ kernel unless
# options(emc2.cpp_rfun = FALSE).  Draws the h = 1 FRQ decision time x and
# maps it back through the fading clock, u = -log(1 - kappa x) / kappa;
# x >= 1/kappa is an omission.
rFRQfade <- function(lR, pars, ok = rep(TRUE, length(lR))) {
  .frqfade_check_cols(pars)
  ok <- .frq_check_race(lR, pars, ok, "FRQfade")
  kappa <- pars[, "kappa"]
  dt <- rep(Inf, nrow(pars))
  idx <- which(ok & is.finite(pars[, "t0"]) & is.finite(kappa) & kappa >= 0)
  if (length(idx)) {
    base <- cbind(pars[idx, , drop = FALSE], h = 1)
    x <- .frq_draw_dt(base)
    k <- kappa[idx]
    u <- x
    fade <- k > 0
    u[fade] <- Inf
    live <- fade & k * x < 1
    u[live] <- -log1p(-k[live] * x[live]) / k[live]
    dt[idx] <- u
  }
  .frq_race(lR, dt, pars[, "t0"], ok)
}

#' The FRQ Process with a Fading Registration Opportunity (FRQfade)
#'
#' Model file to estimate [FRQ()] with a fading registration opportunity in
#' EMC2.
#'
#' @details
#'
#' Each accumulator's decision time `u = rt - t0` has distribution function
#'
#' `F(u) = F_1(A(u))`, with `A(u) = (1 - exp(-kappa * u)) / kappa`,
#'
#' where `F_1` is the [FRQ()] distribution function at `h = 1` with the same
#' `K`, `v`, `sv`, `plo` and `clock`, and `A(u) = u` at `kappa = 0`. The
#' density is `f(u) = f_1(A(u)) * exp(-kappa * u)`. An accumulator never
#' finishes with probability `1 - h`, where `h = F_1(1/kappa)` is derived; a
#' race omission has probability `prod(1 - h)` over accumulators.
#'
#' | **Parameter** | **Transform** | **Bounds** | **Default** | **Mapping** | **Interpretation** |
#' |---|---|---|---|---|---|
#' | *K* | 1 + exp | \[1, Inf\] | log(2) (K = 3) | | Event threshold |
#' | *v* | log | \[1e-4, Inf\] | log(6) | | Mean initial input rate |
#' | *sv* | log | \[1e-4, Inf\] or 0 | log(0) | *beta* = (*v*/*sv*)^2 | Between-trial SD of the input rate |
#' | *kappa* | log | \[1e-4, Inf\] or 0 | log(0.5) | *h* = *F_1*(1/*kappa*) | Decay rate of the registration opportunity |
#' | *t0* | log | \[0.05, Inf\] or 0 | log(0) | | Non-decision time |
#' | *plo* | probit | \[1e-6, 1\] or 0 | qnorm(0) | | Probability of the low-threshold state |
#' | *clock* | log | \[1e-4, Inf\] or 0 | log(1) | | Shape of the registration law |
#'
#' The values `K = 1`, `sv = 0`, `kappa = 0`, `t0 = 0`, `plo = 0` and
#' `clock = 0` are exact bound exceptions and can be set as constants. At
#' `kappa = 0` the model is [FRQ()] with `h = 1`.
#'
#' `mapped_pars()` adds the derived columns `h`, `Etot` (expected total
#' evidence, `beta * G / (1 - G)` with `G` the registration law at `1/kappa`,
#' or `v / kappa` at `sv = 0`; `Inf` at `kappa = 0`), `beta`, `lambda` and
#' `cv` = `sv / v`.
#'
#' Optional nuisance parameters are `pContaminant`, an omission probability
#' outside the accumulators, and `pGuess`, a uniform-outlier probability.
#'
#' @return A model list defining the FRQfade race model.
#' @examples
#' design_FRQfade <- design(data = forstmann, model = FRQfade,
#'   matchfun = function(d) d$S == d$lR,
#'   formula = list(K ~ 1, v ~ lM + E, sv ~ lM, kappa ~ 1, t0 ~ 1))
#' @export
FRQfade <- function() {
  p_types <- c("K" = log(2), "v" = log(6), "sv" = log(0),
               "kappa" = log(0.5), "t0" = log(0), "plo" = qnorm(0),
               "clock" = log(1))
  transform <- c(K = "exp", v = "exp", sv = "exp", kappa = "exp",
                 t0 = "exp", plo = "pnorm", clock = "exp")
  transform_lower <- c(K = 1)
  # Bounds as in FRQ(); kappa's lower guard keeps a sampled kappa off the
  # exact kappa = 0 (h = 1) boundary, which is reachable as a constant.
  minmax <- cbind(K = c(1, Inf), v = c(1e-4, Inf), sv = c(1e-4, Inf),
                  kappa = c(1e-4, Inf), t0 = c(0.05, Inf),
                  plo = c(1e-6, 1), clock = c(1e-4, Inf))
  exception <- c(K = 1, sv = 0, kappa = 0, t0 = 0, plo = 0, clock = 0)
  nuis <- add_nuisance_pars(p_types, transform, minmax, exception)
  p_types <- nuis$p_types
  transform <- nuis$transform
  minmax <- nuis$minmax
  exception <- nuis$exception

  list(
    type = "RACE",
    c_name = "FRQ_FADE",
    p_types = p_types,
    p_types_canonical = c("K", "v", "sv", "kappa", "t0"),
    transform = list(func = transform, lower = transform_lower),
    bound = list(minmax = minmax, exception = exception),
    Ttransform = function(pars, dadm) {
      # Reporting only: the compiled likelihood derives h itself.
      pl <- if ("plo" %in% colnames(pars)) pars[, "plo"] else rep(0, nrow(pars))
      cu <- if ("clock" %in% colnames(pars)) pars[, "clock"] else rep(1, nrow(pars))
      fs <- frq_fade_summary(pars[, "K"], pars[, "v"], pars[, "sv"],
                             pars[, "kappa"], pl, cu)
      g <- frq_rate(pars[, "K"], pars[, "v"], pars[, "sv"], 1, pl, cu)
      add <- cbind(h = as.numeric(fs[, "h"]), Etot = as.numeric(fs[, "Etot"]),
                   beta = as.numeric(g[, "beta"]),
                   lambda = as.numeric(g[, "lambda"]),
                   cv = as.numeric(pars[, "sv"] / pars[, "v"]))
      rownames(add) <- NULL
      cbind(pars, add)
    },
    rfun = function(data, pars) {
      ok <- attr(pars, "ok")
      if (is.null(ok)) ok <- rep(TRUE, nrow(pars))
      .rfun_FRQfade(data$lR, pars, ok = ok)
    },
    dfun = function(rt, pars) dFRQfade(rt, pars),
    pfun = function(rt, pars) pFRQfade(rt, pars),
    log_likelihood = function(pars, dadm, model, min_ll = log(1e-10)) {
      stop("FRQfade: the likelihood is implemented in the compiled race path; ",
           "the R likelihood route is not supported.")
    }
  )
}
