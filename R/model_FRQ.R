# ============================================================================
# FRQ: the Finite Reservoir Quorum process (Math/FRQ.tex)
#
#   N potential evidence units, each independently available with probability p
#   and, when available, registering at an Exponential(lambda) latency.  The
#   accumulator responds once K units have registered, so its decision time is
#   the K-th order statistic -- and is +Inf whenever fewer than K units happen
#   to be available at all.
#
# With A = K and R = N - K + 1 the whole thing is closed form,
#
#   F(x) = I_{q(x)}(A, R),   q(x) = p(1 - e^{-lambda x}),
#
# and the continuous relaxation is just A, R > 0.  Those two shapes are called
# `alpha` and `beta` here: `R` is a reserved data column name in EMC2
# (R/design.R:251), and U ~ Beta(alpha, beta) is the exact latent
# representation of the quorum, so the Beta shape names are the honest ones.
#
# The estimated coordinates are (alpha, beta, h, tau, t0) rather than
# (alpha, beta, p, lambda, t0): h = I_p(alpha, beta) is the eventual completion
# probability and tau the conditional median decision time.  All of the
# numerics, including that inversion, live in src/model_FRQ.h -- the dfun/pfun
# below call the SAME compiled kernels as the sampled likelihood, so
# make_data()/predict() and the fit cannot disagree.
# ============================================================================

.frq_check_cols <- function(pars) {
  need <- c("alpha", "beta", "h", "tau", "t0")
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
                    tau = pars[ok, "tau"])
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
                    tau = pars[ok, "tau"])
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
                    tau = pars[ok, "tau"], lower_tail = FALSE)
  }
  out
}

# Pure-R reference simulator; .rfun_FRQ uses the C++ kernel unless
# options(emc2.cpp_rfun = FALSE).  This is the exact order-statistic
# representation of Math/FRQ.tex Sec. 6, not an approximation: drawing the
# latent quorum U ~ Beta(alpha, beta) and inverting the registration CDF at U/p
# is distributionally identical to building N cues and waiting for the K-th.
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
    pl <- frq_rate(p[, "alpha"], p[, "beta"], p[, "h"], p[, "tau"])
    U <- rbeta(nrow(p), p[, "alpha"], p[, "beta"])
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
#' A race model in which each accumulator draws on a *finite reservoir* of
#' potential evidence units and responds once a quorum of them has registered.
#' There are `N` potential units; each is independently available with
#' probability `p`, and an available unit registers at an `Exponential(lambda)`
#' latency. The accumulator responds when `K` distinct units have registered,
#' so its decision time is the `K`-th order statistic of the effective arrival
#' times --- and is infinite whenever fewer than `K` units happen to be
#' available at all.
#'
#' @details
#'
#' **Closed form.** By time `x` a fraction `q(x) = p G(x)` of the reservoir has
#' registered, so the number registered is `Binomial(N, q(x))` and a response
#' has occurred exactly when that count reaches `K`. Writing the shapes as
#' `alpha = K` and `beta = N - K + 1`,
#' \deqn{F(x) = I_{q(x)}(\alpha, \beta), \qquad
#'       f(x) = \frac{p\,g(x)}{B(\alpha,\beta)}
#'              q(x)^{\alpha-1}[1 - q(x)]^{\beta-1},}
#' with `q(x) = p(1 - exp(-lambda x))` for exponential registration. There is
#' no first-passage solve and no quadrature anywhere in the likelihood. The
#' model is a *continuous relaxation*: `alpha` and `beta` are real-valued, and
#' the integer case is the literal finite-cue process.
#'
#' **Parameters are `(alpha, beta, h, tau, t0)`, not `(N, K, p, lambda, t0)`.**
#' The generative coordinates are badly conditioned to estimate, so the exposed
#' ones are the two Beta shapes plus two interpretable summaries: `h`, the
#' probability that the accumulator *ever* completes, and `tau`, its
#' **conditional median decision time**, `P(T <= tau | T < Inf) = 0.5`. The
#' inversion back to the generative pair is exact,
#' \deqn{p = I^{-1}_{h}(\alpha, \beta), \qquad
#'       \lambda = -\log(1 - u/p)/\tau, \qquad
#'       u = I^{-1}_{h/2}(\alpha, \beta),}
#' and is performed inside the compiled kernel on every evaluation. `Ttransform`
#' reports `p` and `lambda` (and `N = alpha + beta - 1`) through that same
#' inversion, so `mapped_pars()` shows the generative description without the
#' sampler ever having to live in it. Note `tau` is a decision-time median; the
#' observed-RT conditional median is `t0 + tau`.
#'
#' Default values are used for all parameters that are not explicitly listed in
#' the `formula` argument of `design()`. They can also be accessed with
#' `FRQ()$p_types`.
#'
#' | **Parameter** | **Transform** | **Natural scale** | **Default** | **Mapping** | **Interpretation** |
#' |---|---|---|---|---|---|
#' | *alpha* | log | \[1, Inf\] | log(1) | | Quorum depth (`K`). Controls the order of the leading edge: `f(x) ~ x^(alpha-1)`. |
#' | *beta* | log | \[1, Inf\] | log(1) | | Residual redundancy (`N - K + 1`). For `p = 1` the late hazard approaches `beta * lambda`. |
#' | *h* | probit | \[0, 1\] | qnorm(0.95) | | Probability the accumulator *ever* completes; `1 - h` is its never-finish mass. |
#' | *tau* | log | \[0, Inf\] | log(0.5) | | Conditional median decision time, in seconds. |
#' | *t0* | log | \[0, Inf\] | log(0) | | Non-decision time. |
#' | *pContaminant* | probit | \[0, 1\] | qnorm(0) | | Optional *omission* contaminant probability: mass at `rt = Inf` only, handled by the data pipeline |
#' | *pGuess* | probit | \[0, 1\] | qnorm(0) | | Optional uniform *guess* (outlier) probability, mixed into observed RT densities over the guess window |
#'
#' **Two dimensions of evidence strength.** `h` and `tau` separate mechanisms
#' that most race models confound: `h` is whether sufficient evidence is
#' ultimately *attainable*, and `tau` is how fast the process completes *when
#' it succeeds*. Crossing only `tau` with a condition gives a timing-only
#' model, only `h` an availability-only model, and both a dual model; these are
#' nested in one likelihood. The recommended default is to share `alpha`,
#' `beta` and `t0` structurally across accumulators and let `h` and `tau` carry
#' the design.
#'
#' **Defectiveness is intrinsic.** Every accumulator fails permanently with
#' probability `1 - h`, so an FRQ race produces omissions with no contaminant
#' parameter at all. This is *not* an additive omission mixture: `p` also
#' reshapes the finite-RT distribution, and only in the corner `alpha + beta -
#' 1 = alpha` (that is, `beta = 1`, `K = N`) does conditioning on a response
#' remove `p` from the RT distribution entirely. A `J`-accumulator race fails
#' only if every accumulator does, with probability `prod(1 - h_i)`, so
#' appreciable per-accumulator defect is compatible with essentially no
#' observable omissions --- two accumulators at `h = 0.98` give a race-level
#' omission rate of 0.0004. Under a deadline `D` the omission probability is
#' `prod(S_i(D))` instead, which the package's censoring machinery handles.
#'
#' **Identifiability.** `alpha` and `beta` are shape parameters read off the
#' leading edge and the late hazard respectively, and are far more weakly
#' informed than `h` and `tau`. Prefer `alpha ~ 1` and `beta ~ 1` and let the
#' design act on `h` and `tau`; crossing the shapes with the same factors as
#' `tau` will produce a ridge rather than an error. There is no separate
#' evidence-scale parameter to fix --- `h` and `tau` are already on
#' probability and time scales --- so, unlike the ballistic models, FRQ needs
#' no scaling constant.
#'
#' @param relax Logical. If `FALSE` (default), the shapes are bounded below at
#'   1, the *conservative continuous relaxation*: this keeps the model close to
#'   the geometry of a literal finite reservoir, where `alpha = K >= 1` and
#'   `beta = N - K + 1 >= 1`. If `TRUE`, the bound drops to a small epsilon,
#'   giving the unrestricted transformed-Beta family. `alpha < 1` is
#'   mathematically valid but produces a density singularity at the leading
#'   edge, since `f(x) ~ x^(alpha-1)`.
#' @return A model list defining the FRQ race model.
#' @examples
#' ADmat <- matrix(c(-1/2, 1/2), ncol = 1, dimnames = list(NULL, "d"))
#' matchfun <- function(d) d$S == d$lR
#' # The shapes are shared structurally; the design acts on the two evidence
#' # dimensions, tau (speed given completion) and h (attainability).
#' design_FRQ <- design(data = forstmann, model = FRQ, matchfun = matchfun,
#'                      formula = list(alpha ~ 1, beta ~ 1, h ~ lM,
#'                                     tau ~ lM + E, t0 ~ 1),
#'                      contrasts = list(tau = list(lM = ADmat)))
#' @export
FRQ <- function(relax = FALSE) {
  # Shapes on the log scale: alpha = 1 (K = 1) and beta = 1 (K = N) are both
  # exactly reachable defaults, and together they are the one-unit reservoir.
  p_types <- c("alpha" = log(1), "beta" = log(1),
               # A completion probability near but not at one: the default must
               # not silently declare a 50% omission rate for any parameter the
               # user leaves out of the formula.
               "h" = qnorm(0.95),
               # A sub-second conditional median: tau is a decision time in
               # seconds, so log(1) would be an implausible starting scale.
               "tau" = log(0.5), "t0" = log(0))
  transform <- c(alpha = "exp", beta = "exp", h = "pnorm", tau = "exp",
                 t0 = "exp")
  # relax = FALSE keeps alpha, beta >= 1 (the conservative continuous FRQ
  # relaxation, Math/FRQ.tex Sec. 13); relax = TRUE admits the unrestricted
  # transformed-Beta family, whose alpha < 1 corner has a leading-edge
  # singularity.  The kernel is identical either way -- this is only a bound.
  shape_lo <- if (isTRUE(relax)) 1e-4 else 1
  # h's upper bound stops just short of one so that log(1 - h), the score of an
  # intrinsic no-response trial, stays finite.  h = 1 is the non-defective
  # limit and is approached, not attained.
  minmax <- cbind(alpha = c(shape_lo, Inf), beta = c(shape_lo, Inf),
                  h = c(1e-6, 1 - 1e-9), tau = c(1e-4, Inf),
                  t0 = c(0.05, Inf))
  exception <- c(t0 = 0)

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
    # `relax` changes only the R-side bounds, so both variants share one
    # compiled kernel and one c_name.
    c_name = "FRQ",
    relax = isTRUE(relax),
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
                     pars[, "tau"])
      # N = alpha + beta - 1 is the reservoir size of the literal finite-cue
      # process, and d = alpha/N the fraction of it required to reach quorum:
      # small d is Poisson-counter-like, appreciable d is where the finite
      # reservoir actually bites.  Both are only literal at integer shapes.
      N <- pars[, "alpha"] + pars[, "beta"] - 1
      # Built as a bare matrix first: cbind()ing named scalars onto a one-row
      # `pars` makes R invent row names from the argument names, which would
      # then travel with the parameter matrix.
      add <- cbind(p = as.numeric(pl[, "p"]),
                   lambda = as.numeric(pl[, "lambda"]),
                   N = as.numeric(N),
                   d = as.numeric(pars[, "alpha"] / N))
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
