# ============================================================================
# FRQ: finite-reservoir quorum process
#
# N units are independently available with probability p and, when available,
# register after delays whose mean rate is lambda.  With cv_u = 0 the delays are
# exponential; cv_u > 0 integrates an independent Gamma-distributed unit rate.
# The accumulator responds at the K-th registration and otherwise has an
# infinite finishing time.
#
# With alpha = K, beta = N - K + 1:
#   F(x) = I_{q(x)}(alpha, beta),  q(x) = p G_cv(x),
#   G_0(x) = 1 - exp(-lambda x),
#   G_cv(x) = 1 - (1 + cv_u^2 lambda x)^(-1/cv_u^2).
#
# The fitted coordinates are (alpha, beta, p, lambda, t0, delta, cv_u), where
# p is the per-unit availability probability and lambda is the registration
# rate. The eventual completion probability h is a derived report only. The R
# wrappers use the kernels in src/model_FRQ.h.
# ============================================================================

.frq_check_cols <- function(pars) {
  need <- c("alpha", "beta", "p", "lambda", "t0")
  missing <- setdiff(need, colnames(pars))
  if (length(missing))
    stop("FRQ requires parameter columns ", paste(missing, collapse = ", "))
  need
}

dFRQ <- function(rt, pars) {
  .frq_check_cols(pars)
  dt <- rt - pars[, "t0"]
  ok <- (rt > 0) & (dt > 0) & is.finite(dt) & is.finite(pars[, "t0"])
  ok[is.na(ok)] <- FALSE
  out <- numeric(length(dt))
  if (any(ok)) {
    delta <- if ("delta" %in% colnames(pars)) pars[ok, "delta"] else 0
    cv_u <- if ("cv_u" %in% colnames(pars)) pars[ok, "cv_u"] else 0
    out[ok] <- dfrq(t = dt[ok], alpha = pars[ok, "alpha"],
                    beta = pars[ok, "beta"], p = pars[ok, "p"],
                    lambda = pars[ok, "lambda"], delta = delta, cv_u = cv_u)
  }
  out
}

pFRQ <- function(rt, pars) {
  .frq_check_cols(pars)
  dt <- rt - pars[, "t0"]
  # rt = Inf is deliberately kept: the CDF there is h (the complement of the
  # never-finish mass), not one.  The compiled kernel returns the saturated
  # branch, which is what the defective-tail bookkeeping needs.
  ok <- (rt > 0) & (dt > 0) & is.finite(pars[, "t0"])
  ok[is.na(ok)] <- FALSE
  out <- numeric(length(dt))
  if (any(ok)) {
    delta <- if ("delta" %in% colnames(pars)) pars[ok, "delta"] else 0
    cv_u <- if ("cv_u" %in% colnames(pars)) pars[ok, "cv_u"] else 0
    out[ok] <- pfrq(t = dt[ok], alpha = pars[ok, "alpha"],
                    beta = pars[ok, "beta"], p = pars[ok, "p"],
                    lambda = pars[ok, "lambda"], delta = delta, cv_u = cv_u)
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
  ok <- (rt > 0) & (dt > 0) & is.finite(pars[, "t0"])
  ok[is.na(ok)] <- FALSE
  if (any(ok)) {
    delta <- if ("delta" %in% colnames(pars)) pars[ok, "delta"] else 0
    cv_u <- if ("cv_u" %in% colnames(pars)) pars[ok, "cv_u"] else 0
    out[ok] <- pfrq(t = dt[ok], alpha = pars[ok, "alpha"],
                    beta = pars[ok, "beta"], p = pars[ok, "p"],
                    lambda = pars[ok, "lambda"], delta = delta, cv_u = cv_u,
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
  if (nr < 1L || nrow(pars) %% nr != 0L)
    stop("FRQ: pars rows must be a multiple of accumulators")
  if (length(ok) != nrow(pars))
    stop("FRQ: ok must have one value per parameter row")
  ok <- as.logical(ok)
  ok[is.na(ok)] <- FALSE
  n_trials <- nrow(pars) / nr
  dt <- matrix(Inf, nrow = nr, ncol = n_trials)
  t0 <- pars[, "t0"]
  idx <- which(ok & is.finite(t0))
  if (length(idx)) {
    p <- pars[idx, , drop = FALSE]
    dl <- if ("delta" %in% colnames(p)) p[, "delta"] else rep(0, nrow(p))
    cu <- if ("cv_u" %in% colnames(p)) p[, "cv_u"] else rep(0, nrow(p))
    # With threshold variability the latent quorum percentile is no longer
    # Beta: draw uniformly on the CDF scale and pull it back through H before
    # inverting the incomplete beta.
    U <- rbeta(nrow(p), p[, "alpha"], p[, "beta"])
    vary <- is.finite(dl) & dl > 0
    if (any(vary)) {
      z <- frq_h_inv_r(runif(sum(vary)), dl[vary])
      U[vary] <- qbeta(z, p[vary, "alpha"], p[vary, "beta"])
    }
    hit <- rep(Inf, nrow(p))
    live <- is.finite(p[, "p"]) & is.finite(p[, "lambda"]) &
      is.finite(U) & (U <= p[, "p"])
    live[is.na(live)] <- FALSE
    if (any(live)) {
      ratio <- U[live] / p[live, "p"]
      # A proper model has p = 1.  An exact U = 1 is a zero-probability
      # endpoint, but a finite RNG can return it; keep it finite rather than
      # creating an artificial omission through an infinite inverse.
      proper <- p[live, "p"] == 1
      ratio[proper & ratio >= 1] <- 1 - .Machine$double.eps / 2
      c2 <- cu[live] * cu[live]
      z <- numeric(length(ratio))
      exp_branch <- c2 == 0
      if (any(exp_branch))
        z[exp_branch] <- -log1p(-ratio[exp_branch]) /
          p[live, "lambda"][exp_branch]
      if (any(!exp_branch))
        z[!exp_branch] <- expm1(-c2[!exp_branch] *
                                  log1p(-ratio[!exp_branch])) /
          (c2[!exp_branch] * p[live, "lambda"][!exp_branch])
      hit[live] <- z
    }
    hit[!is.finite(hit) | hit < 0] <- Inf
    dt[idx] <- hit
  }
  # Invalid/non-finite t0 rows stay at +Inf rather than poisoning the race
  # column with NA when the additive shift is applied.
  t0_safe <- t0
  t0_safe[!is.finite(t0_safe)] <- 0
  dt <- dt + matrix(t0_safe, nrow = nr)

  bad_col <- colSums(!is.infinite(dt)) == 0L
  R <- max.col(-t(dt), ties.method = "first")
  pick <- cbind(R, seq_len(ncol(dt)))
  rt <- dt[pick]
  R <- factor(levels(lR)[R], levels = levels(lR))
  R[bad_col] <- NA
  rt[bad_col] <- Inf
  # A trial remains simulable when at least one accumulator is valid; using
  # only the first row here made omission depend on response-level ordering.
  ok_col <- colSums(matrix(ok, nrow = nr)) > 0L
  out <- data.frame(R = rep(NA_character_, n_trials), rt = rep(NA_real_, n_trials))
  out$R[ok_col] <- levels(lR)[R][ok_col]
  out$R <- factor(out$R, levels = levels(lR))
  out$rt[ok_col] <- rt[ok_col]
  .apply_timed_guess_winner(out, levels(lR))
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
#' given trial with probability `p`, and a usable unit registers after a delay
#' with mean rate `lambda`. At `cv_u = 0` the delay is exponential; at positive
#' `cv_u` its latent unit rate is Gamma-distributed. If fewer than `K` of
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
#' | *p* | probit | \[0, 1\] | qnorm(0.95) | | Probability that an evidence unit is available |
#' | *lambda* | log | \[0, Inf\] | log(log(2) / 0.5) | | Registration rate of an available evidence unit |
#' | *t0* | log | \[0, Inf\] | log(0) | | Non-decision time |
#' | *delta* | log | \[0, Inf\] | log(0) | | Half-width of continuous threshold variability |
#' | *cv_u* | log | \[0, Inf\] | log(0) | | CV of an individual unit's latent registration rate |
#'
#' The fitted parameters are the generative ones. `alpha` and `beta` are the
#' quorum size `K` and the spare capacity `N - K + 1`, treated as continuous
#' and bounded below at 1, so `alpha = beta = 1` is a single-unit pool. `p` is
#' the probability that each evidence unit is available and `lambda` is the
#' registration rate of an available unit. `mapped_pars()` reports the
#' eventual completion probability `h` in addition to the pool size `N` =
#' `alpha` + `beta` - 1, the quorum fraction `d` = `alpha`/`N`, and the
#' unit-rate Gamma shape `a_u = 1/cv_u^2` (`Inf` at zero).
#'
#' `p` and `lambda` are the two mechanistic evidence dimensions: `p` changes
#' whether the finite pool can supply a quorum at all, while `lambda` changes
#' the speed of registration conditional on availability. This makes direct
#' rate versus quorum/threshold tests possible. The derived `h` is useful for
#' reporting omission risk, but it is not used as a substitute for either
#' mechanism in a design formula.
#'
#' `delta` adds between-trial variability in how much evidence the accumulator
#' demands before responding, playing the role that start-point or threshold
#' variability plays in the ballistic models. On a trial where `delta` is
#' larger, the requirement is more often unusually lenient *or* unusually
#' strict, and less often typical; the average requirement is unchanged, so
#' `p` and `lambda` keep their meanings exactly. At the default of zero there is
#' no such variability and the model reduces to the description above.
#' `mapped_pars()` reports `sQ` = `delta`/sqrt(3), the standard deviation of
#' the trial-to-trial shift on a log-odds scale. Note that `delta` is a
#' property of a single accumulator rather than a state shared by the whole
#' race, and that it is only weakly identified: values below about 2 are
#' difficult to distinguish from a small change in `alpha` and `beta`, so it
#' is worth fitting only with a lot of data or a strong prior.
#'
#' `cv_u` is independent unit-level rate heterogeneity, not a shared trialwise
#' drift state. It is best held fixed or shared across racers initially. The
#' proper/non-defective boundary is obtained by fixing `p = 1`, which makes
#' every unit available. This boundary is exact and leaves `lambda` as the
#' rate parameter rather than changing coordinates.
#'
#' `alpha` and `beta` control the shape of the distribution: `alpha` how
#' sharply the density rises at the leading edge, `beta` how heavy the late
#' tail is. Both are more weakly informed than `p` and `lambda`, so the
#' recommended default is `alpha ~ 1` and `beta ~ 1`, shared across
#' accumulators along with `t0`.
#'
#' A race produces an omission only when every accumulator in it fails, with
#' probability `prod(1 - h)` over accumulators, so an appreciable failure rate
#' per accumulator is compatible with very few observed omissions: two
#' accumulators at `h = 0.98` give a race-level omission rate of 0.0004.
#' Optional fitting parameters are `pContaminant`, the omission probability,
#' and `pGuess`, the uniform-outlier probability. The optional
#' threshold-variability `delta` and unit-rate `cv_u` parameters are log/exp
#' transformed with default `log(0)`.
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
#' # The shapes are held constant; the design acts on the two mechanistic
#' # evidence dimensions, lambda (registration rate) and p (availability).
#' design_FRQ <- design(data = forstmann, model = FRQ, matchfun = matchfun,
#'                      formula = list(alpha ~ 1, beta ~ 1, p ~ lM,
#'                                     lambda ~ lM + E, t0 ~ 1),
#'                      contrasts = list(p = list(lM = ADmat),
#'                                       lambda = list(lM = ADmat)))
#' # For all parameters that are not defined in the formula, default values are
#' # assumed (see Table above).
#' @export
FRQ <- function() {
  # Shapes on the log scale: alpha = 1 (K = 1) and beta = 1 (K = N) are both
  # exactly reachable defaults, and together they are the one-unit reservoir.
  p_types <- c("alpha" = log(1), "beta" = log(1),
               # A high, but not proper, availability probability. The
               # default remains defective so omissions are available when a
               # design leaves p out of its formula.
               "p" = qnorm(0.95),
               # For alpha = beta = 1, lambda = log(2) / 0.5 gives a
               # one-unit median registration time of about 0.5 seconds.
               "lambda" = log(log(2) / 0.5), "t0" = log(0),
               # Threshold variability OFF by default: the base FRQ is the
               # delta = 0 member of the family, and log(0) = -Inf reaches it
               # exactly (the same idiom t0 uses to default to zero).
               "delta" = log(0), "cv_u" = log(0))
  transform <- c(alpha = "exp", beta = "exp", p = "pnorm", lambda = "exp",
                 t0 = "exp", delta = "exp", cv_u = "exp")
  # alpha, beta >= 1 is the conservative continuous FRQ relaxation: it is the
  # region a literal finite reservoir can reach (K >= 1 and N - K + 1 >= 1),
  # and it keeps the density bounded at the leading edge, where
  # f(x) ~ x^(alpha - 1).  The unrestricted transformed-Beta family
  # (shapes down to ~0) is mathematically valid but is NOT offered because
  # it leaves the finite-reservoir interpretation and leading-edge behavior.
  # The direct (p, lambda) coordinates have no quantile-inversion dead zone.
  # The free p upper bound stops just short of one so that a sampler can keep
  # a finite defective mass; p = 1 is the explicit non-defective boundary and
  # is reachable through the bound exception below.
  # delta is capped at 6 rather than left unbounded. It is a log-odds
  # half-width, so H'(0)/H'(1/2) = cosh(delta/2)^2 -- already 101 at delta = 6,
  # i.e. extreme quorum percentiles a hundred times more likely than middling
  # ones.
  # Unit-rate CV is positive on the free scale but zero is an exact nested
  # boundary.  Keeping a small positive lower bound avoids a near-zero
  # transformed region with essentially no shape information.
  minmax <- cbind(alpha = c(1, Inf), beta = c(1, Inf),
                  p = c(1e-6, 1 - 1e-9), lambda = c(1e-4, Inf),
                  t0 = c(0.05, Inf), delta = c(1e-4, 6),
                  cv_u = c(1e-4, Inf))
  exception <- c(t0 = 0, delta = 0, cv_u = 0, p = 1)

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
    p_types_canonical = c("alpha", "beta", "p", "lambda", "t0"),
    transform = list(func = transform),
    bound = list(minmax = minmax, exception = exception),
    Ttransform = function(pars, dadm) {
      # Reporting only: Ttransform does not run on the compiled likelihood
      # path (see src/utils.h), so this costs nothing during sampling. The
      # fitted p/lambda columns are already the generative coordinates; only
      # h and the structural summaries are derived here.
      dl <- if ("delta" %in% colnames(pars)) pars[, "delta"] else rep(0, nrow(pars))
      cu <- if ("cv_u" %in% colnames(pars)) pars[, "cv_u"] else rep(0, nrow(pars))
      h <- pFRQ(rep(Inf, nrow(pars)), pars)
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
      add <- cbind(h = as.numeric(h),
                   N = as.numeric(N),
                   d = as.numeric(pars[, "alpha"] / N),
                   sQ = as.numeric(dl) / sqrt(3),
                   a_u = ifelse(cu == 0, Inf, 1 / cu^2))
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
