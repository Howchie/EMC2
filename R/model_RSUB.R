# ============================================================================
# RSUB: racing subordinators (Math/RDMSWTN_Derivation.tex)
#
# Each accumulator's evidence is a nondecreasing pure-jump process with mean
# v u and variance s^2 v u at decision time u, and the accumulator finishes
# when it first reaches the distance d ~ U[B, B + A].  evidence = "dense" is
# the inverse-Gaussian subordinator; "packets" is compound Poisson evidence
# with exponential packets.  The kernels (src/model_RSUB.h) take the
# decision time and (v, B, A, s).
# ============================================================================

.rsub_check_cols <- function(pars) {
  missing <- setdiff(c("v", "B", "A", "t0", "s"), colnames(pars))
  if (length(missing))
    stop("RSUB requires parameter columns ", paste(missing, collapse = ", "))
  invisible(NULL)
}

dRSUB <- function(rt, pars, packets = FALSE) {
  .rsub_check_cols(pars)
  dt <- rt - pars[, "t0"]
  ok <- is.finite(dt) & dt > 0 & is.finite(pars[, "t0"])
  ok[is.na(ok)] <- FALSE
  out <- numeric(length(dt))
  if (any(ok)) {
    out[ok] <- drsub(dt[ok], pars[ok, "v"], pars[ok, "B"], pars[ok, "A"],
                     pars[ok, "s"], packets = packets)
    out[ok][is.na(out[ok])] <- 0
  }
  out
}

pRSUB <- function(rt, pars, packets = FALSE) {
  .rsub_check_cols(pars)
  dt <- rt - pars[, "t0"]
  ok <- !is.na(dt) & dt > 0 & is.finite(pars[, "t0"])
  out <- numeric(length(dt))
  if (any(ok)) {
    out[ok] <- prsub(dt[ok], pars[ok, "v"], pars[ok, "B"], pars[ok, "A"],
                     pars[ok, "s"], packets = packets)
    out[ok][is.na(out[ok])] <- 0
  }
  out
}

# Reference simulator.  Dense: the decision time is (s / v) max_{0 <= w <= d}
# (W(w) + w / s), drawn from the endpoint of the drifted Brownian motion and
# the maximum of the Brownian bridge through it.  Packets: 1 + Pois(2 d / s^2)
# packets are needed, so the decision time is their gamma arrival time at
# rate 2 v / s^2.
rRSUB <- function(lR, pars, packets = FALSE, ok = rep(TRUE, nrow(pars))) {
  .rsub_check_cols(pars)
  nr <- length(levels(lR))
  if (nr < 1L || nrow(pars) %% nr != 0L)
    stop("RSUB: pars rows must be a multiple of accumulators")
  ok <- as.logical(ok)
  ok[is.na(ok)] <- FALSE
  v <- pars[, "v"]; s <- pars[, "s"]; t0 <- pars[, "t0"]
  ok <- ok & is.finite(t0) & is.finite(v) & v > 0 & is.finite(s) & s > 0
  n <- sum(ok)
  dt <- rep(Inf, nrow(pars))
  d <- pars[ok, "B"] + pars[ok, "A"] * runif(n)
  if (packets) {
    c2 <- s[ok]^2
    dt[ok] <- rgamma(n, 1 + rpois(n, 2 * d / c2), rate = 2 * v[ok] / c2)
  } else {
    x <- d / s[ok] + sqrt(d) * rnorm(n)
    mx <- 0.5 * (x + sqrt(x^2 - 2 * d * log(runif(n))))
    dt[ok] <- s[ok] * mx / v[ok]
  }
  n_tri <- nrow(pars) / nr
  finish <- matrix(dt + ifelse(is.finite(t0), t0, 0), nrow = nr, ncol = n_tri)
  R <- max.col(-t(finish), ties.method = "first")
  rt <- finish[cbind(R, seq_len(n_tri))]
  R <- factor(levels(lR)[R], levels = levels(lR))
  R[!is.finite(rt)] <- NA
  out <- cbind.data.frame(R = R, rt = rt)
  .apply_timed_guess_winner(out, levels(lR))
}

#' The Racing Subordinator Model (RSUB)
#'
#' Race model in which each accumulator's evidence `G(u)` is a nondecreasing
#' pure-jump Levy process (a subordinator) with `E[G(u)] = v u` and
#' `Var[G(u)] = s^2 v u`. The accumulator finishes when `G` first reaches the
#' distance `d`, which is uniform on `[B, B + A]`. The finish time is
#' `t0 + u`, and the observed response is the accumulator with the smallest
#' finish time. `evidence` selects the jump law; both laws share the
#' parameters and the first two moments of `G(u)`.
#'
#' `evidence = "dense"` is the inverse-Gaussian subordinator: `G(u)` is the
#' first time a Brownian motion with unit diffusion coefficient and drift
#' `1 / s` reaches the level `(v / s) u`. With `m = v u`, `z = s sqrt(d)`,
#' `alpha = (m - d) / z` and `beta = (m + d) / z`, the point-start
#' first-passage CDF and density are
#'
#' `F(u) = Phi(alpha) - exp(2 m / s^2) Phi(-beta)`,
#'
#' `f(u) = (2 v / z) phi(alpha) - (2 v / s^2) exp(2 m / s^2) Phi(-beta)`.
#'
#' Equivalently the decision time is `(s / v) max_{0 <= w <= d} (W(w) + w / s)`
#' for a standard Brownian motion `W`.
#'
#' `evidence = "packets"` is compound Poisson evidence: events arrive at rate
#' `2 v / s^2` and each adds an exponential packet of mean `s^2 / 2`. With
#' `a = 2 v u / s^2`, `b = 2 d / s^2` and independent `N ~ Poisson(a)`,
#' `K ~ Poisson(b)`, the point-start first-passage CDF and density are
#'
#' `F(u) = P(N - K >= 1)`,
#'
#' `f(u) = (2 v / s^2) exp(-a - b) I_0(2 sqrt(a b))`,
#'
#' with `I_0` the modified Bessel function of the first kind. Equivalently the
#' decision time is the arrival time of `1 + K` events.
#'
#' Under both laws, start-point variability is integrated analytically through
#' the stop-loss transform of the evidence state.
#'
#' | **Parameter** | **Transform** | **Natural scale** | **Default** | **Mapping** | **Interpretation** |
#' |---|---|---|---|---|---|
#' | *v* | log | \[1e-3, Inf\] | log(1) | | Mean evidence rate. |
#' | *B* | log | \[0, Inf\] | log(1) | *b* = *B* + *A* | Distance from *A* to *b* (response threshold). |
#' | *A* | log | \[1e-4, Inf\] | log(0) | | Range of the uniform start point. |
#' | *t0* | log | \[0.05, Inf\] | log(0) | | Non-decision time. |
#' | *s* | log | \[0, Inf\] | log(1) | | Evidence SD per square root of mean evidence (`s^2` is the variance-to-mean ratio). |
#'
#' `mapped_pars()` adds the threshold `b`.
#'
#' Optional nuisance parameters are `pContaminant`, an omission probability
#' outside the accumulators, and `pGuess`, a uniform-outlier probability.
#'
#' @param evidence `"dense"` (the default) for inverse-Gaussian subordinator
#'   evidence, `"packets"` for compound Poisson evidence with exponential
#'   packets.
#' @return A model list defining the RSUB race model.
#' @examples
#' matchfun <- function(d) d$S == d$lR
#' design_RSUB <- design(data = forstmann, model = RSUB, matchfun = matchfun,
#'                       formula = list(v ~ lM, B ~ E, A ~ 1, t0 ~ 1, s ~ 1),
#'                       constants = c(s = log(1)))
#' design_RSUB_pk <- design(data = forstmann,
#'                          model = function() RSUB(evidence = "packets"),
#'                          matchfun = matchfun,
#'                          formula = list(v ~ lM, B ~ E, A ~ 1, t0 ~ 1, s ~ 1),
#'                          constants = c(s = log(1)))
#' @export
RSUB <- function(evidence = c("dense", "packets")) {
  evidence <- match.arg(evidence)
  packets <- evidence == "packets"
  p_types <- c(v = log(1), B = log(1), A = log(0), t0 = log(0), s = log(1))
  transform <- c(v = "exp", B = "exp", A = "exp", t0 = "exp", s = "exp")
  minmax <- cbind(v = c(1e-3, Inf), B = c(0, Inf), A = c(1e-4, Inf),
                  t0 = c(0.05, Inf), s = c(0, Inf))
  exception <- c(A = 0, v = 0, t0 = 0)
  .nuis <- add_nuisance_pars(p_types, transform, minmax, exception)
  list(
    type = "RACE",
    c_name = if (packets) "RSUB_PK" else "RSUB",
    p_types = .nuis$p_types,
    p_types_canonical = c("v", "B", "A", "t0", "s"),
    transform = list(func = .nuis$transform),
    bound = list(minmax = .nuis$minmax, exception = .nuis$exception),
    Ttransform = function(pars, dadm) {
      cbind(pars, b = pars[, "B"] + pars[, "A"])
    },
    rfun = function(data = NULL, pars) .rfun_RSUB(data$lR, pars, packets,
                                                  ok = attr(pars, "ok")),
    dfun = function(rt, pars) dRSUB(rt, pars, packets),
    pfun = function(rt, pars) pRSUB(rt, pars, packets),
    log_likelihood = function(pars, dadm, model, min_ll = log(1e-10)) {
      stop("RSUB: the likelihood is implemented in the compiled race path; ",
           "the R likelihood route is not supported.")
    }
  )
}
