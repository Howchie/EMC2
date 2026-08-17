# ============================================================================
# FRQ: the Finite Reservoir Quorum process
#
#   N potential evidence units, each independently available with probability p
#   and, when available, registering at an Exponential(lambda) latency.  The
#   accumulator responds once K units have registered, so its decision time is
#   the K-th order statistic -- and is +Inf whenever fewer than K units happen
#   to be available at all.
#
# With alpha = K and beta = N - K + 1 the whole thing is closed form,
#
#   F(x) = I_{q(x)}(alpha, beta),   q(x) = p(1 - e^{-lambda x}),
#
# and the continuous relaxation is just alpha, beta > 0 (bounded at 1 in the
# fitted model; see FRQ()).  The shapes are named `alpha`/`beta` rather than
# K and N - K + 1 because U ~ Beta(alpha, beta) is the exact latent
# representation of the quorum, and because the natural name for the second one
# (`R`) is a reserved data column in EMC2 (R/design.R:251).
#
# The estimated coordinates are (alpha, beta, h, tau, t0) rather than
# (alpha, beta, p, lambda, t0): h = I_p(alpha, beta) is the eventual completion
# probability and tau the conditional median decision time.  All of the
# numerics, including that inversion, live in src/model_FRQ.h -- the dfun/pfun
# below call the SAME compiled kernels as the sampled likelihood, so
# make_data()/predict() and the fit cannot disagree.
# ============================================================================

.frq_check_cols <- function(pars) {
  need <- c("alpha", "beta", "h", "tau", "t0", "delta")
  missing <- setdiff(need, colnames(pars))
  if (length(missing))
    stop("FRQ requires parameter columns ", paste(missing, collapse = ", "))
  need
}

dFRQ <- function(rt, pars) {
  .frq_check_cols(pars)
  dt <- rt - pars[, "t0"]
  ok <- (rt > 0) & (dt > 0) & is.finite(dt)
  ok[is.na(ok)] <- FALSE
  out <- numeric(length(dt))
  if (any(ok)) {
    out[ok] <- dfrq(t = dt[ok], alpha = pars[ok, "alpha"],
                    beta = pars[ok, "beta"], h = pars[ok, "h"],
                    tau = pars[ok, "tau"], delta = pars[ok, "delta"])
  }
  out
}

pFRQ <- function(rt, pars) {
  .frq_check_cols(pars)
  dt <- rt - pars[, "t0"]
  # rt = Inf is deliberately kept: the CDF there is h (the complement of the
  # never-finish mass), not one.  The compiled kernel returns the saturated
  # branch, which is what the defective-tail bookkeeping needs.
  ok <- (rt > 0) & (dt > 0)
  ok[is.na(ok)] <- FALSE
  out <- numeric(length(dt))
  if (any(ok)) {
    out[ok] <- pfrq(t = dt[ok], alpha = pars[ok, "alpha"],
                    beta = pars[ok, "beta"], h = pars[ok, "h"],
                    tau = pars[ok, "tau"], delta = pars[ok, "delta"])
  }
  out
}

# Survivor 1 - F, saturating at the defective mass 1 - h.  Computed through the
# reflection I_{1-q}(beta, alpha) rather than by subtraction, so it stays
# accurate in the deep tail where 1 - F is the interesting quantity.
sFRQ <- function(rt, pars) {
  .frq_check_cols(pars)
  dt <- rt - pars[, "t0"]
  out <- rep(1, length(dt))
  ok <- (rt > 0) & (dt > 0)
  ok[is.na(ok)] <- FALSE
  if (any(ok)) {
    out[ok] <- pfrq(t = dt[ok], alpha = pars[ok, "alpha"],
                    beta = pars[ok, "beta"], h = pars[ok, "h"],
                    tau = pars[ok, "tau"], delta = pars[ok, "delta"],
                    lower_tail = FALSE)
  }
  out
}

# Pure-R reference simulator; .rfun_FRQ uses the C++ kernel unless
# options(emc2.cpp_rfun = FALSE).  This is the exact order-statistic
# representation, not an approximation: drawing the latent quorum
# U ~ Beta(alpha, beta) and inverting the registration CDF at U/p is
# distributionally identical to building N cues and waiting for the K-th,
# because F(x) = I_{q(x)}(alpha, beta) makes U <= q(x) equivalent to T <= x.
# U > p means the reservoir saturates below the quorum, so the accumulator
# never terminates; an all-Inf trial column is the package's omission
# convention (R = NA, rt = Inf), which make_data() already handles.
rFRQ <- function(lR, pars, ok = rep(TRUE, length(lR))) {
  .frq_check_cols(pars)
  nr <- length(levels(lR))
  bad <- rep(NA, length(lR) / nr)
  out <- data.frame(R = bad, rt = bad)
  n_trials <- nrow(pars) / nr
  dt <- matrix(Inf, nrow = nr, ncol = n_trials)
  t0 <- pars[, "t0"]

  idx <- which(ok)
  if (length(idx)) {
    p <- pars[idx, , drop = FALSE]
    # Same compiled inversion the likelihood uses, so the simulator cannot
    # drift away from the density it is meant to be sampling from.
    pl <- frq_rate(p[, "alpha"], p[, "beta"], p[, "h"], p[, "tau"],
                   p[, "delta"])
    # With threshold variability the latent quorum percentile is no longer
    # Beta: draw uniformly on the CDF scale and pull it back through H before
    # inverting the incomplete beta.  delta == 0 keeps the plain rbeta draw so
    # that seeded simulations predating delta still reproduce exactly.
    dl <- p[, "delta"]
    U <- rbeta(nrow(p), p[, "alpha"], p[, "beta"])
    vary <- is.finite(dl) & dl > 0
    if (any(vary)) {
      z <- frq_h_inv_r(runif(sum(vary)), dl[vary])
      U[vary] <- qbeta(z, p[vary, "alpha"], p[vary, "beta"])
    }
    hit <- rep(Inf, nrow(p))
    live <- !is.na(pl[, "p"]) & (U <= pl[, "p"])
    live[is.na(live)] <- FALSE
    if (any(live)) {
      hit[live] <- -log1p(-U[live] / pl[live, "p"]) / pl[live, "lambda"]
    }
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

#' The Finite Reservoir Quorum Process (FRQ)
#'
#' Model file to estimate the Finite Reservoir Quorum (FRQ) process in EMC2.
#'
#' @details
#' Model files are almost exclusively used in `design()`.
#'
#' Each accumulator draws on a finite pool of `N` potential evidence units and
#' responds once a quorum of `K` of them has registered. A unit is usable on a
#' given trial with probability `p`, and a usable unit registers after an
#' exponentially distributed latency with rate `lambda`. If fewer than `K` of
#' the `N` units happen to be usable, the quorum is never reached and the
#' accumulator never responds, so the model produces omissions without needing
#' a contaminant parameter.
#'
#' Default values are used for all parameters that are not explicitly listed in
#' the `formula` argument of `design()`. They can also be accessed with
#' `FRQ()$p_types`.
#'
#' | **Parameter** | **Transform** | **Natural scale** | **Default** | **Mapping** | **Interpretation** |
#' |---|---|---|---|---|---|
#' | *alpha* | log | \[1, Inf\] | log(1) | | Quorum size (`K`): how many evidence units are needed to respond |
#' | *beta* | log | \[1, Inf\] | log(1) | | Spare capacity (`N - K + 1`): how many units beyond the quorum the pool holds |
#' | *h* | probit | \[0, 1\] | qnorm(0.95) | | Probability that the accumulator ever responds |
#' | *tau* | log | \[0, Inf\] | log(0.5) | | Median decision time on the trials where it does respond |
#' | *t0* | log | \[0, Inf\] | log(0) | | Non-decision time |
#' | *delta* | log | \[0, 6\] | log(0) | | Between-trial variability in threshold; 0 (the default) turns it off |
#' | *pContaminant* | probit | \[0, 1\] | qnorm(0) | | Optional *omission* contaminant probability: mass at `rt = Inf` only, handled by the data pipeline |
#' | *pGuess* | probit | \[0, 1\] | qnorm(0) | | Optional uniform *guess* (outlier) probability, mixed into observed RT densities over the guess window |
#'
#' The sampled parameters are not the generative ones. `alpha` and `beta` are
#' the quorum size `K` and the spare capacity `N - K + 1`, treated as
#' continuous and bounded below at 1, so `alpha = beta = 1` is a single-unit
#' pool. In place of the availability probability `p` and the registration
#' rate `lambda`, which are difficult to estimate, the model is parameterised
#' by two directly interpretable quantities: `h`, the probability that the
#' accumulator ever responds, and `tau`, its median decision time among the
#' trials on which it does. `tau` is a decision time, so the median observed
#' RT of a lone accumulator is `t0 + tau`. `mapped_pars()` reports `p` and
#' `lambda`, along with the pool size `N` = `alpha` + `beta` - 1 and the
#' quorum fraction `d` = `alpha`/`N`, next to the sampled parameters.
#'
#' `h` and `tau` are the two dimensions of evidence strength and are where the
#' design normally goes: `h` is whether enough evidence is ultimately
#' available, `tau` is how quickly the decision is reached when it is. A
#' design on `tau` alone is a pure speed effect, on `h` alone a pure
#' availability effect, and on both is the model that contains the other two.
#' Note that holding `tau` fixed while `h` varies fixes the completion time
#' *given* completion, which is not the same as holding `lambda` fixed.
#'
#' `delta` adds between-trial variability in how much evidence the accumulator
#' demands before responding, playing the role that start-point or threshold
#' variability plays in the ballistic models. On a trial where `delta` is
#' larger, the requirement is more often unusually lenient *or* unusually
#' strict, and less often typical; the average requirement is unchanged, so
#' `h` and `tau` keep their meanings exactly. At the default of zero there is
#' no such variability and the model reduces to the description above.
#' `mapped_pars()` reports `sQ` = `delta`/sqrt(3), the standard deviation of
#' the trial-to-trial shift on a log-odds scale. Note that `delta` is a
#' property of a single accumulator rather than a state shared by the whole
#' race, and that it is only weakly identified: values below about 2 are
#' difficult to distinguish from a small change in `alpha` and `beta`, so it
#' is worth fitting only with a lot of data or a strong prior.
#'
#' `alpha` and `beta` control the shape of the distribution: `alpha` how
#' sharply the density rises at the leading edge, `beta` how heavy the late
#' tail is. Both are much more weakly informed than `h` and `tau`, so the
#' recommended default is `alpha ~ 1` and `beta ~ 1`, shared across
#' accumulators along with `t0`. Crossing the shapes with the same factors as
#' `tau` gives a poorly identified model rather than an error. Unlike the
#' ballistic models, FRQ has no evidence-scale parameter that must be fixed
#' for identifiability: `h` is a probability and `tau` is a time.
#'
#' A race produces an omission only when every accumulator in it fails, with
#' probability `prod(1 - h)` over accumulators, so an appreciable failure rate
#' per accumulator is compatible with very few observed omissions: two
#' accumulators at `h = 0.98` give a race-level omission rate of 0.0004.
#' `pContaminant` is offered for consistency with the other models but is
#' redundant here unless an omission mechanism *outside* the decision process
#' is wanted. With a single accumulator and `beta` fixed at 1 the two are not
#' identified at all.
#'
#' Because the FRQ is a race model, it has one accumulator per response
#' option. EMC2 automatically constructs a factor representing the
#' accumulators `lR` (i.e., the latent response) with level names taken from
#' the `R` column in the data. For race models, the `design()` argument
#' `matchfun` can be provided, a function that takes the `lR` factor (defined
#' in the augmented data (d) in the following function) and returns a logical
#' defining the correct response. In the example below, the match is simply
#' such that the `S` factor equals the latent response factor:
#' `matchfun=function(d)d$S==d$lR`. Then `matchfun` is used to automatically
#' create a latent match (`lM`) factor with levels `FALSE` (i.e., the stimulus
#' does not match the accumulator) and `TRUE` (i.e., the stimulus does match
#' the accumulator). This is added internally and can also be used in model
#' formula, typically for parameters related to the rate of accumulation.
#'
#' @return A model list defining the FRQ race model.
#' @examples
#' ADmat <- matrix(c(-1/2, 1/2), ncol = 1, dimnames = list(NULL, "d"))
#' matchfun <- function(d) d$S == d$lR
#' # The shapes are held constant; the design acts on the two evidence
#' # dimensions, tau (speed given a response) and h (whether one occurs).
#' design_FRQ <- design(data = forstmann, model = FRQ, matchfun = matchfun,
#'                      formula = list(alpha ~ 1, beta ~ 1, h ~ lM,
#'                                     tau ~ lM + E, t0 ~ 1),
#'                      contrasts = list(h = list(lM = ADmat),
#'                                       tau = list(lM = ADmat)))
#' # For all parameters that are not defined in the formula, default values are
#' # assumed (see Table above).
#' @export
FRQ <- function() {
  # Shapes on the log scale: alpha = 1 (K = 1) and beta = 1 (K = N) are both
  # exactly reachable defaults, and together they are the one-unit reservoir.
  p_types <- c("alpha" = log(1), "beta" = log(1),
               # A completion probability near but not at one: the default must
               # not silently declare a 50% omission rate for any parameter the
               # user leaves out of the formula.
               "h" = qnorm(0.95),
               # A sub-second conditional median: tau is a decision time in
               # seconds, so log(1) would be an implausible starting scale.
               "tau" = log(0.5), "t0" = log(0),
               # Threshold variability OFF by default: the base FRQ is the
               # delta = 0 member of the family, and log(0) = -Inf reaches it
               # exactly (the same idiom t0 uses to default to zero).
               "delta" = log(0))
  transform <- c(alpha = "exp", beta = "exp", h = "pnorm", tau = "exp",
                 t0 = "exp", delta = "exp")
  # alpha, beta >= 1 is the conservative continuous FRQ relaxation: it is the
  # region a literal finite reservoir can reach (K >= 1 and N - K + 1 >= 1),
  # and it keeps the density bounded at the leading edge, where
  # f(x) ~ x^(alpha - 1).  The unrestricted transformed-Beta family
  # (shapes down to ~0) is mathematically valid but is NOT offered: the
  # (h, tau) coordinates are not representable there in double precision.
  # qbeta(h, alpha, beta) underflows to 0 for both p and u at small h (the
  # kernel then rejects the point, giving the sampler an artificial cliff)
  # and rounds to exactly 1 at large h -- at alpha = beta = 0.05, h = 0.99
  # already gives p = 1, i.e. a silently PROPER distribution with none of the
  # requested defect.  The dead zone reaches alpha = 0.5 at h = 1 - 1e-9.
  # h's upper bound stops just short of one so that log(1 - h), the score of an
  # intrinsic no-response trial, stays finite.  h = 1 is the non-defective
  # limit and is approached, not attained.
  # delta is capped at 6 rather than left unbounded.  It is a log-odds
  # half-width, so H'(0)/H'(1/2) = cosh(delta/2)^2 -- already 101 at delta = 6,
  # i.e. extreme quorum percentiles a hundred times more likely than middling
  # ones.  The cap also protects the inversion: H^{-1} maps 1 - h to roughly
  # (1 - h) * delta/sinh(delta), so an unbounded delta would drive p's
  # complement below double precision at the top of h's range.
  minmax <- cbind(alpha = c(1, Inf), beta = c(1, Inf),
                  h = c(1e-6, 1 - 1e-9), tau = c(1e-4, Inf),
                  t0 = c(0.05, Inf), delta = c(1e-4, 6))
  # delta = 0 is the nested base model and must stay reachable even though it
  # is below the lower bound, exactly as t0 = 0 is.
  exception <- c(t0 = 0, delta = 0)

  # pContaminant (omission) and pGuess (uniform outlier); see add_nuisance_pars().
  # FRQ already produces omissions intrinsically through 1 - h, so
  # pContaminant is genuinely redundant here unless a mechanism *outside* the
  # decision process is wanted; it is offered for consistency, not because the
  # model needs it.
  .nuis <- add_nuisance_pars(p_types, transform, minmax, exception)
  p_types <- .nuis$p_types; transform <- .nuis$transform
  minmax <- .nuis$minmax; exception <- .nuis$exception

  list(
    type = "RACE",
    c_name = "FRQ",
    p_types = p_types,
    p_types_canonical = setdiff(names(p_types), .nuisance_par_names),
    transform = list(func = transform),
    bound = list(minmax = minmax, exception = exception),
    Ttransform = function(pars, dadm) {
      # Reporting only: Ttransform does not run on the compiled likelihood
      # path (see src/utils.h), so this costs nothing during sampling and the
      # kernel derives p/lambda itself.  It routes through the same exported
      # inversion, so the reported generative parameters are by construction
      # the ones the likelihood used.
      pl <- frq_rate(pars[, "alpha"], pars[, "beta"], pars[, "h"],
                     pars[, "tau"], pars[, "delta"])
      # N = alpha + beta - 1 is the reservoir size of the literal finite-cue
      # process, and d = alpha/N the fraction of it required to reach quorum:
      # small d is Poisson-counter-like, appreciable d is where the finite
      # reservoir actually bites.  Both are only literal at integer shapes.
      # The alpha, beta >= 1 bounds keep N >= 1 and d in (0, 1].
      N <- pars[, "alpha"] + pars[, "beta"] - 1
      # Built as a bare matrix first: cbind()ing named scalars onto a one-row
      # `pars` makes R invent row names from the argument names, which would
      # then travel with the parameter matrix.
      # sQ = delta/sqrt(3) is SD(eps) for eps ~ U(-delta, delta), i.e. the
      # trialwise SD of the log-odds shift in caution -- the scale a reader
      # actually has intuitions about, unlike the half-width itself.
      add <- cbind(p = as.numeric(pl[, "p"]),
                   lambda = as.numeric(pl[, "lambda"]),
                   N = as.numeric(N),
                   d = as.numeric(pars[, "alpha"] / N),
                   sQ = as.numeric(pars[, "delta"]) / sqrt(3))
      rownames(add) <- NULL
      cbind(pars, add)
    },
    rfun = function(data, pars) {
      .rfun_FRQ(data$lR, pars, ok = attr(pars, "ok"))
    },
    dfun = function(rt, pars) dFRQ(rt, pars),
    pfun = function(rt, pars) pFRQ(rt, pars),
    log_likelihood = function(pars, dadm, model, min_ll = log(1e-10)) {
      stop("FRQ: the likelihood is implemented in the compiled race path; ",
           "the R likelihood route is not supported.")
    }
  )
}
